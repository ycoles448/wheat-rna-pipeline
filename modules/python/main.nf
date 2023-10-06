// Module information
name = "python"


process rename_trinity {
    tag "${name}"
    cpus 1
    memory "1G"
    executor "local"
    // time module.time
    // queue module.queue

    input:
    tuple path(samples), val(id)

    output:
    tuple path("assemblies.txt"), path("readsL.txt"), path("readsR.txt")

    shell:
    '''
    python3 "!{params.data.bin}/rename_trinity.py" "samples.txt" '!{id}'
    ls
    '''
}
