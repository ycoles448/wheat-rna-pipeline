// Module information
name = "gffread"
module = params[name]

run(new File("modules/template.groovy"))

process CONVERT_GFF {
    // tag "${name}-convert"
    // cpus cores
    // memory "${Math.floor(memory * buffer)}G"
    // time "${time}hour"
    // executor executor
    executor "local"
    queue queue

    publishDir "${params.data.out}/${params.data.genomes}/${params.data.species}",
        mode: "copy",
        overwrite: true,
        pattern: "genome.gtf"

    input:
    path(files)

    output: stdout
    path("genome.gtf"), emit: "gtf"

    shell:
    template "gffread_convert.sh"
}
