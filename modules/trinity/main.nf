// Module information
name = "trinity"
module = params[name]
template = "modules/template.groovy"
run(new File(template))

// Force buffer value to 1
buffer = 1

process TRINITY {
    tag "${name}-${ids}"
    cpus Math.floor(cores)
    memory "${Math.floor(memory * buffer)}G"
    time "${time}hour"
    // executor executor
    executor "local"
    queue queue

    // publishDir "${params.data.out}/${params.data.reads_trimmed}",
    //     mode: "copy",
    //     overwrite: true,
    //     pattern: "trim_*_r{1,2}.fastq"

    input:
    val(ids)
    val(meta)
    path(files)

    output: stdout
    // val(ids), emit: "ids"
    // val(meta), emit: "meta"
    // path("trim_*_r{1,2}.fastq"), emit: "reads"

    shell:
    template "trinity.sh"
}
