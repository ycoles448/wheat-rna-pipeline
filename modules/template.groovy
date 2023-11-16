// Enable containers and set binary
if (module.buffer > 1) {
    module.bin = "${params.parallel.bin} ${module.bin}"
}

if (module.container == true) {
    if (params.container.monolithic == true) {
        container = params.container.path
    } else {
        container = module.path
    }

    // Set executable
    if (module.direct == true) {
        bin = "${container} ${module.bin}"
    } else {
        bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
    }
} else {
    bin = "${module.bin}"
}

// Set resources
threads = module.cores.toInteger()
if (params.hardware.smt == true) {
    cores = module.cores
    if (module.buffer > 1) {
        cores = (params.nextflow.buffer * 2).toInteger()
    }

    threads = (module.cores.toInteger() * 2).toInteger()
}

memory = module.memory
if (module.memory < 0) {
    memory = 1 + Math.floor(threads * 0.25)
}

time = 1
if (!module.time) {
    println("Set non-existent time value to ${time}h for module ${name}")
}

executor = params.process.executor
if (module.executor != null) {
    executor = module.executor
}

queue = params.process.queue
if (module.queue != null) {
    queue = module.queue
}

buffer = 1
if (module.buffer != null) {
    buffer = module.buffer.toInteger()
}
