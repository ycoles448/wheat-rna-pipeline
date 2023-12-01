// Module information
name = "samtools"
module = params[name]

run(new File("modules/template.groovy"))

process MERGE_BAMS {
    tag "${name}-${ids}"
    cpus Math.floor(cores)
    memory "${memory}G"
    time "${time}hour"
    executor executor
    queue queue

    publishDir "${params.data.out}/${params.data.bam}",
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
    path("${params.data.species}_${meta[params.meta.group][0]}.bam"), emit: "bams"

    shell:
    template "merge_bams.sh"
}
