// Module information
name = "extract"
module = params[name]
template = "modules/template.groovy"
run(new File(template))

process EXTRACT {
    tag "${name}-${id}"
    cpus cores
    memory "${Math.floor(memory * buffer)}G"
    time "${time}hour"
    executor executor
    queue module.queue

    publishDir "${params.data.out}/${params.data.reads_unzipped}",
        mode: "copy",
        overwrite: true,
        pattern: "*_r{1,2}.fastq"

    input:
    val(ids)
    val(meta),
    path(files)

    output: stdout
    val(ids), emit: "ids"
    val(meta), emit: "meta"
    path("*_r{1,2}.fastq"), emit: "reads"

    shell:
    template "extract.sh"
}
