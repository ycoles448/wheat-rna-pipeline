#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Variables
species = "${params.data.species}"
meta_f = "${params.data.meta}"
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
reads_unzipped = "${params.data.path}/${params.data.reads_unzipped}/${params.data.glob_unzipped}"
reads_trimmed = "${params.data.path}/${params.data.reads_trimmed}/${params.data.glob_trimmed}"
reads_bam = "${params.data.path}/${params.data.star}/*_${params.data.species}/*_${params.data.bam_star_sorted}"
genomes = "${params.data.path}/${params.data.genomes}"
adapters = "${params.data.path}/${params.data.adapters}"

factors = params.meta.factors
group = params.meta.group


// Logging
log.info """
# ============================= #
# - RNA-seq Nextflow Pipeline - #
# ============================= #

> species                   ${species}
> reads                     ${reads}
> reads (extracted)         ${reads_unzipped}
> reads (trimmed)           ${reads_trimmed}
> meta                      ${meta_f}
"""


// Load modules
include { EXTRACT } from "./modules/extract"
include { FASTP } from "./modules/fastp"
include { MERGE_BAMS } from "./modules/samtools"
include { TRINITY } from "./modules/trinity"
// include { STAR } from "./modules/star"


// Channels
// Reads
reads_pe = Channel.fromFilePairs(reads)
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

adapters = Channel.fromPath(adapters)
meta = Channel.fromPath(meta_f)


workflow MAIN {
}

workflow QC_READS {
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
    if (params.workflow.skip_extract) {
        READS_BAM = reads_bam
    } else {
        // READS_BAM = STAR()
    }

    READS_GROUPED = META
        .join(READS)
        .join(READS_BAM)
        .map {
            it[it.size() - 2].add(it[it.size() - 1])
            it.remove(it.size() - 1)
            return it
        }
        .collect { [it] }
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

    READS_GROUPED.ids
        .buffer {size: 8}
        .view()

    READS_GROUPED.meta
        .buffer {size: 8, remainder: true}
        .view()

    BAMS = MERGE_BAMS(
        READS_GROUPED.ids.flatMap().last(),
        READS_GROUPED.meta.flatMap().last(),
        READS_GROUPED.files.flatMap().last()
    )

    BAMS[0].view()

    // TRINITY(
    //     READS_GROUPED.ids.flatMap(),
    //     READS_GROUPED.meta.flatMap(),
    //     READS_GROUPED.files.flatMap()
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
    QC_READS()
    TRANSCRIPTOME(QC_READS.out.READS, QC_READS.out.META)
    // FUNGI(MAIN.out.reads)
    // PLANT(MAIN.out.reads)
}
