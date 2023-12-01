#!/usr/bin/env bash -l

# Fix number of threads per process, instead of !{threads_s1} per Nextflow definition
threads=4

ids=$(echo !{ids} | tr -d '[,]')
files=$(echo !{files} | tr -d '[,]')
flags=(!{module.flags})
flags+=(
    --CPU ${threads}
    --max_memory !{module.memory_s1}G
)

!{params.parallel.bin} !{bin} ${flags[@]} \
    --output {.}_trinity \
    --genome_guided_bam {} \
    --genome_guided_max_intron !{params[params.data.species].max_intron_length} \
    --no_run_inchworm ::: ${files}
