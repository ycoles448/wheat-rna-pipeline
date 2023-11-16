// Module information
name = "samtools"
module = params[name]
template = "modules/template.groovy"
run(new File(template))

process MERGE_BAMS {
    tag "${name}-${ids}"
    cpus Math.floor(cores)
    memory "${Math.floor(memory * buffer)}G"
    time "${time}hour"
    // executor executor
    executor "local"
    queue queue

    publishDir "${params.data.out}/${params.data.bam}/${params.data.species}_${meta[params.meta.group][0]}",
        mode: "copy",
        overwrite: true,
        pattern: "${params.data.species}_${meta[params.meta.group][0]}.bam"

    input:
    val(ids)
    val(meta)
    path(files)

    output: stdout
    val(ids), emit: "ids"
    val(meta), emit: "meta"
    path("${params.data.species}_${meta[params.meta.group][0]}.bam"), emit: "bam"

    shell:
    template "merge_bams.sh"
}
