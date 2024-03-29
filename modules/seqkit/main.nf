// Module information
name = "seqkit"
module = params[name]

run(new File("modules/template.groovy"))

process SEQKIT_LENGTH {
    tag "${name}-length"
    cpus cores
    memory "${Math.floor(memory * buffer)}G"
    // time "${time}hour"
    // executor executor
    executor "local"
    // queue queue

    publishDir "${params.data.out}/${params.data.genomes}/${params.data.species}",
        mode: "copy",
        overwrite: true,
        pattern: "length.txt"

    input:
    tuple val(id), path(files)

    output: stdout
    tuple val(id), path("length.txt"), emit: "length"

    shell:
    template "seqkit_length.sh"
}
