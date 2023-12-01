// Module information
name = "star"
module = params[name]

run(new File("modules/template.groovy"))

private int get_index_n() {
    run(new File("src/vars.groovy"))

    def int l = new File(genome_length_f).text.toInteger()
    def int n = Math.floor(((Math.log(l) / Math.log(2)) - 1) / 2 - 1)

    return Math.min(14, n).toInteger()
}

// STAR parameters
overhang = params.data.reads_length.toInteger() -
    params.qc.trimf.toInteger() -
    params.qc.trimt.toInteger() - 1
index_n = get_index_n()

threads_align = module.cores_align.toInteger()
threads_index = module.cores_index.toInteger()
if (params.hardware.smt) {
    threads_align = Math.floor(module.cores_align * 2).toInteger()
    threads_index = Math.floor(module.cores_index * 2).toInteger()
}
threads_buf_align = Math.floor(threads_index / module.buffer_index).toInteger()
threads_buf_align = Math.floor(threads_align / module.buffer_align).toInteger()
ram_limit = Math.floor(module.memory_index * Math.pow(1024, 3)).toInteger()

process STAR_INDEX {
    tag "${name}-${species}-index"
    cpus module.cores_index.toInteger()
    memory "${Math.floor(module.memory_index * buffer)}G"
    time "${module.time_index}hours"
    executor executor
    queue module.queue_index

    publishDir "${params.data.out}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "${params.data.species}"

    input:
    path(files)

    output: stdout
    path("${params.data.species}"), emit: "index"

    shell:
    template "star_index.sh"
}

process STAR_ALIGN {
    cpus module.cores_align.toInteger()
    memory "${Math.floor(module.memory_align * buffer)}G"
    time "${module.time_align}hours"
    executor executor
    queue module.queue_align

    publishDir "${params.data.out}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "*_${params.data.species}"

    input:
    val(ids)
    val(meta)
    path(files)

    output: stdout
    val(ids), emit: "ids"
    val(meta), emit: "meta"
    path("*_${params.data.species}"), emit: "bam"

    shell:
    template "star_align.sh"
}
