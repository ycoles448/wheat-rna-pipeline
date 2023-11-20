// Module information
name = "trinity"
module = params[name]
template = "./modules/template.groovy"
run(new File(template))

threads_s1 = module.cores_s1
if (params.hardware.smt) threads_s1 = Math.floor(module.cores_s1 * 2).toInteger()

process TRINITY_S1 {
    // tag "${name}-${ids}"
    cpus module.cores_s1
    memory "${module.memory_s1}G"
    time "${module.time_s1}hour"
    executor params.process.executor
    queue module.queue_s1

    publishDir "${params.data.out}/${params.data.trinity}",
        mode: "copy",
        overwrite: true,
        pattern: null

    input:
    val(ids)
    val(meta)
    path(files)

    output: stdout
    val(ids), emit: "ids"
    val(meta), emit: "meta"
    // path("trim_*_r{1,2}.fastq"), emit: "reads"

    shell:
    template "trinity_s1.sh"
}
