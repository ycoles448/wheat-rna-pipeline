// Module information
name = "fastp"
module = params[name]
template = "modules/template.groovy"
run(new File(template))

process FASTP {
    tag "${name}-${ids}"
    cpus Math.floor(cores * 2)
    memory "${Math.floor(memory * buffer)}G"
    time "${time}hour"
    executor executor
    queue queue

    publishDir "${params.data.out}/${params.data.reads_trimmed}",
        mode: "copy",
        overwrite: true,
        pattern: "trim_*_r{1,2}.fastq"

    input:
    val(ids)
    val(meta)
    path(files)

    output: stdout
    val(ids), emit: "ids"
    val(meta), emit: "meta"
    path("trim_*_r{1,2}.fastq"), emit: "reads"

    shell:
    template "fastp.sh"
}
