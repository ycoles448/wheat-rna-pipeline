// Module information
name = "interpro"
module = params[name]


// Enable containers and set binary
if (module.container == "True") {
    // Check if monolithic containers are enabled, and set container
    if (params.container.monolithic == "True") {
        container = params.container.path
    } else {
        container = module.path
    }

    // Set executable
    if (module.direct == "True") {
        bin = "${container} ${module.bin}"
    } else {
        bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
    }
} else {
    bin = "${module.bin}"
}


// Set threads
if (params.hardware.smt == "True") {
    threads = (module.cores.toInteger() * 2).toInteger()
}


// Set memory
if (module.memory == -1) {
    throw new Exception(
        "Module ${name} does not support dynamic memory assignment, \
        please set to 0 or a positive integer")
} else if (module.memory == 0) {
    memory = 1
    module.memory = "${module.memory}G"
    binding.variables.remove 'memory'
}
else {
    module.memory = "${module.memory}G"
}


// Set time
if (!module.time) {
    println("Setting non-existent time value to 0")
    module.time_align = 0
    module.time_index = 0
}


// Module settings
logfile = ""


process module {
    tag "${name}-${id}"
    cpus module.cores
    memory module.memory
    time "${module.time}.hour"
    queue module.queue

    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"
    publishDir "${params.data.logs}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "${logfile}"
    publishDir "${params.data.logs}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "${logfile}",
        saveAs: {"${name}-index-${id}-${it}"}

    input:
    tuple val(id), path(proteins)

    output:
    tuple val(id)
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true
    path("${logfile}"), emit: "logs", optional: true

    shell:
    '''
    # Move proteins into a folder InterProScan likes
    # mkdir -p "!{id}"
    # mv  "!{proteins}" "!{id}"

    # Disable lookup service, is not working
    flags="!{module.flags}"
    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} "${flags}" \
        -i "./!{proteins}" \
        -b "./!{id}" \
	-dp \
        --goterms \
        --applications SMART,CDD,Pfam,PANTHER \
        --cpu "!{threads}"

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Threads: !{threads}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
