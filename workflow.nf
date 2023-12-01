#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2

// Variables
run(new File("${projectDir}/src/vars.groovy"))


// Logging
log.info """
# ============================= #
# - RNA-seq Nextflow Pipeline - #
# ============================= #

> species                   ${species}
> reads                     ${reads}
> reads (extracted)         ${reads_unzipped}
> reads (trimmed)           ${reads_trimmed}
> star index                ${star_index_f}
> bam (star)                ${reads_bam}
> bam (grouped)             ${reads_bam_grouped}
> genome                    ${genome_f}
> genome_length             ${genome_length_f}
> meta                      ${meta_f}
"""


// Load modules
include { EXTRACT } from "./modules/extract"
include { FASTP } from "./modules/fastp"
include { CONVERT_GFF } from "./modules/convert"
include { SEQKIT_LENGTH } from "./modules/seqkit"
include { STAR_ALIGN; STAR_INDEX } from "./modules/star"
include { MERGE_BAMS } from "./modules/samtools"
include { TRINITY_S1 } from "./modules/trinity"


// Channels
// Reads
reads_pe = Channel.fromFilePairs(reads)
genome = Channel.fromPath(genome_f)
genome_gff = Channel.fromPath(genome_gff_f)
if (params.workflow.skip_extract) {
    reads_pe = Channel.fromFilePairs(reads_unzipped)
}
if (params.workflow.skip_qc_reads) {
    reads_pe = Channel.fromFilePairs(reads_trimmed)
        .map {
            it[0] = it[0].replaceFirst("trim_", "")
            return it
        }
}
if (params.workflow.skip_gff_convert) {
    genome_gtf = Channel.fromPath(genome_gtf_f)
}
if (params.workflow.skip_genome_stats) {
    genome_length = Channel.fromPath(genome_length_f)
}
if (params.workflow.skip_index) {
    star_index = Channel.fromPath(star_index_f)
}
if (params.workflow.skip_align) {
    reads_bam = Channel.fromPath(reads_bam)
        .map {
            id = new File(it.toString())
                .getParentFile()
                .getName()
                .replaceFirst("_${species}", "")
            return tuple(id, it)
        }
}
if (params.workflow.skip_bam_merge) {
    reads_bam_grouped = Channel.fromPath(reads_bam_grouped)
}

adapters = Channel.fromPath(adapters)
meta = Channel.fromPath(meta_f)


process GET_NODE_DETAILS {
    executor "local"
    output: stdout
    shell: '''echo Node: $(sstat -j $SLURM_JOB_ID | awk 'NR==3 {print $3}')'''
}


workflow QC_READS {
    take:
    NODE_DETAILS

    main:
    READS = reads_pe
        // .map {id, reads -> [id.replaceFirst("^0+(?!\$)", ""), reads]}
    META = meta
        .splitCsv(header: true, sep: "\t")
        .map {tuple(it.sample.replaceFirst("^.*:", ""), it.subMap(factors))}
    READS_BUF = META
        .join(READS)
        .buffer(size: params.fastp.buffer.toInteger(), remainder: true)
        .multiMap {
            def ArrayList ids = []
            def ArrayList meta = []
            def ArrayList files = []
            for (i = 0; i < it.size(); i++) {
                ids.add(it[i][0])
                meta.add(it[i][1])
                for (f = 0; f < it[i][2].size(); f++) files.add(it[i][2][f])
            }

            all: it
            ids: ids
            meta: meta
            files: files
        }

    if (!params.workflow.skip_extract) {
        READS_BUF = EXTRACT(READS_BUF)
    }

    if (!params.workflow.skip_qc_reads) {
        READS_BUF = FASTP(READS_BUF.ids, READS_BUF.meta, READS_BUF.files)
    }

    READS = READS_BUF.all
        .flatMap()
        .map {
            it.remove(1)
            return it
        }

    emit:
    READS = READS
    META = META
}


workflow TRANSCRIPTOME {
    take:
    READS
    META

    main:

    GENOME = genome

    if (params.data.is_gff && params.data.skip_gff_convert) {
        GENOME_GTF = CONVERT_GFF(genome_gff).gtf
    } else {
        GENOME_GTF = genome_gtf
    }
    if (params.workflow.skip_genome_stats) {
        GENOME_LENGTH = genome_length
    } else {
        GENOME_LENGTH = SEQKIT_LENGTH(genome).length
    }

    READS_BUF = META
        .join(READS)
        .buffer(size: params.star.buffer_align.toInteger(), remainder: true)
        .multiMap {
            def ArrayList ids = []
            def ArrayList meta = []
            def ArrayList files = []
            for (i = 0; i < it.size(); i++) {
                ids.add(it[i][0])
                meta.add(it[i][1])
                for (f = 0; f < it[i][2].size(); f++) files.add(it[i][2][f])
            }

            all: it
            ids: ids
            meta: meta
            files: files
        }

    if (params.workflow.skip_index) {
        INDEX = star_index
    } else {
        INDEX = STAR_INDEX(GENOME.combine(GENOME_GTF))
        INDEX[0].view()
        INDEX = INDEX.index
    }

    STAR_FILES = READS_BUF.files
        .combine(GENOME_LENGTH)
        .combine(GENOME_GTF)
        .combine(INDEX)

    READS = READS_BUF.all
        .flatMap()
        .map {
            it.remove(1)
            return it
        }

    if (params.workflow.skip_align) {
        READS_BAM = reads_bam
    } else {
        READS_BAM = STAR_ALIGN(
            READS_BUF.ids.last(),
            READS_BUF.meta.last(),
            STAR_FILES.last()
        ).bam
    }

    READS_GROUPED = META
        .join(READS)
        .join(READS_BAM)
        .map {
            it[it.size() - 2].add(it[it.size() - 1])
            it.remove(it.size() - 1)
            return it
        }
        .collect {[it]}
        .multiMap {
            def Set<String> group_factors = []
            def ArrayList<ArrayList<String>> ids = []
            def ArrayList<Map<String, List<String>>> meta = []
            def ArrayList<ArrayList<ArrayList<Path>>> files = []

            // Collect factor levels
            for (i = 0; i < it.size(); i++) {
                group_factors.add(it[i][1][group])
            }

            // Split metadata by factor
            for (level = 0; level < group_factors.size(); level++) {
                ids[level] = []
                meta[level] = []
                files[level] = []
                for (i = 0; i < it.size(); i++) {
                    if (it[i][1][group] == group_factors[level]) {
                        ids[level].add(it[i][0])
                        meta[level].add(it[i][1])
                        for (f = 0; f < it[i][2].size(); f++) files[level].add(it[i][2][f]);
                    }
                }
            }

            all: it
            ids: ids
            meta: meta
            files: files
        }

    if (params.workflow.skip_bam_merge) {
        BAMS_GROUPED = reads_bam_grouped.collect()
    } else {
        BAMS_GROUPED = MERGE_BAMS(
            READS_GROUPED.ids.flatMap(),
            READS_GROUPED.meta.flatMap(),
            READS_GROUPED.files.flatMap()
        )["bams"].collect()
    }

    PREASM_GROUPED = READS_GROUPED.all.map { [it] }
        .concat(BAMS_GROUPED.map { [it] })
        .collect()
        .multiMap { it, bams ->
            def Set<String> group_factors = []
            def ArrayList<ArrayList<String>> ids = []
            def ArrayList<ArrayList<Map<String, String>>> meta = []
            def ArrayList<ArrayList<String>> files = []
            def ArrayList<String> files_bams = [:]

            // Collect factor levels
            for (i = 0; i < it.size(); i++) {
                group_factors.add(it[i][1][group])
            }

            // Split metadata by factor
            for (level = 0; level < group_factors.size(); level++) {
                ids[level] = [] as ArrayList<String>
                meta[level] = [] as Set<Map<String, String>>
                files[level] = [] as Set<String>

                for (i = 0; i < bams.size(); i++) {
                    bams_name = new File(bams[i].toString()).getName()
                        .replaceFirst("${species}_", "")
                        .replaceFirst("[.]bam", "")
                    if (bams_name.toString() == group_factors[level].toString()) {
                        files_bams[level] = bams[i]
                    }
                }

                for (i = 0; i < it.size(); i++) {
                    if (it[i][1][group] == group_factors[level]) {
                        ids[level].add(it[i][0])
                        meta[level].add(["cultivar": it[i][1]["cultivar"]])
                        files[level].add(files_bams[level])
                    }
                }
            }

            all: it
            ids: ids
            meta: meta
            files: files
        }

    // Buffer all inputs into single job (stage 1 is poorly multi-threaded)
    // PREASM_NORM_GROUPED = TRINITY_S1(
    //     PREASM_GROUPED.ids.flatMap(),
    //     PREASM_GROUPED.meta.flatMap(),
    //     PREASM_GROUPED.files.flatMap()
    // )
}

// workflow MAIN {
//     main:
//     // Unzip reads, unless disabled
//     if ("${params.workflow.skip_extract}" == "True") {
//         EXTRACT = reads_pe
//     } else {
//         EXTRACT(reads_pe)
//         EXTRACT = EXTRACT.out.reads
//     }
//     // Buffer testing
//     // EXTRACT = EXTRACT
//     //     .combine(adapters)
//     //     .buffer(size: params.fastp.buffer.toInteger(), remainder: true)
//     //     .groupTuple()
//     // EXTRACT.view()
//     // INSPECTOR(EXTRACT).view()

//     // Trimming and QC
//     // FASTP(EXTRACT.combine(adapters))
//     // FASTQC(FASTP.out.reads.combine(adapters))
//     // MULTIQC(FASTQC.out.reports.collect())

//     // LOADER
//     // Used for when Nextflow doesn't resume efficiently
//     // Make sure directories terminate with a trailing slash "/"
//     path = Channel.value("reads-trimmed/")
//     LOADER1(
//         EXTRACT.flatMap {n -> n[0]}
//             .combine(path)) // Reads to trimmed reads

//     reads_out = LOADER1.out
//     // reads_out = FASTP.out.reads

//     emit:
//     reads = reads_out
// }

// workflow PLANT {
//     take:
//     reads

//     main:
//     id = Channel.of("${params.data.plant}")
//     genome = Channel.fromPath("${genomes}/${params.data.plant}/genome.fasta")
//     gtf = Channel.fromPath("${genomes}/${params.data.plant}/genome.gtf")
//     reads = reads

//     // SEQ_LENGTH(genome)
//     // genome_length = SEQ_LENGTH.out.seq_length

//     // STAR_INDEX(
//     //     id
//     //         .combine(genome)
//     //         .combine(SEQ_LENGTH.out.seq_length)
//     //         .combine(gtf))
//     // index = STAR_INDEX.out.index

//     // // STAR
//     // // Returns the STAR genome index
//     // path = Channel.value("${params.data.star}")
//     // LOADER1(id.combine(path).combine(id))
//     // index = LOADER1.out

//     // STAR accepts: reads, genome_length, gtf, indexer output
//     // STAR(
//     //     reads
//     //         .combine(genome_length)
//     //         .combine(gtf)
//     //         .combine(index))

//     // // Trinity
//     path = Channel.value("${params.data.star}/")
//     LOADER2(reads.flatMap {n -> "${n[0]}-${params.data.plant}"}
//        .combine(path)
//        .combine(Channel.value("Aligned.sortedByCoord.out.bam")))
//     bam = LOADER2.out
//     intron_max = Channel.value("${params.data.plant_intron_max}")

//     // TRINITY2: Guided assembly
//     TRINITY2(reads.flatMap {n -> [["${n[0]}-${params.data.plant}"]]}
//             .join(bam)
//             .combine(intron_max))
//     // TRINITY2.out[0].view()

//     // emit:
//     // star = STAR.out.star
// }

// workflow FUNGI {
// }

workflow {
    GET_NODE_DETAILS().view()
    QC_READS(GET_NODE_DETAILS)
    TRANSCRIPTOME(QC_READS.out.READS, QC_READS.out.META)
    // FUNGI(MAIN.out.reads)
    // PLANT(MAIN.out.reads)
}
