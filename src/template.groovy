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
        bin = "${params.container.bin} ${params.container.opts}" +
            "${params.container.exec} ${params.container.exec_opts}" +
            "${container} ${module.bin}".stripIndent()
    }
} else {
    bin = "${module.bin}"
}

// Set resources
buffer = 1
if (module.cores != null) {
    // Allocate cores and threads normally
    cores = module.cores
    if (module.buffer > 1) {
        buffer = module.buffer
        cores = Math.floor(cores * buffer).toInteger()
    }

    // Allocate with SMT enabled
    threads = cores
    if (params.hardware.smt) {
        threads = Math.floor(2 * cores).toInteger()
    }

    threads_buf = Math.floor(threads / buffer).toInteger()
}

if (module.memory != null) {
    memory = module.memory

    if (module.memory < 0) {
        memory = 1 + Math.floor(threads * 0.25)
    }
}

if (module.time != null) {
    time = module.time
    if (module.time == null) {
        time = 1
        println("Set non-existent time value to ${time}h for module ${name}")
    }
}

executor = params.process.executor
if (module.executor != null) {
    executor = module.executor
}

queue = params.process.queue
if (module.queue != null) {
    queue = module.queue
}
