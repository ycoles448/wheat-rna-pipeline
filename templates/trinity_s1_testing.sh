#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
files=$(echo !{files} | tr -d '[,]')
flags=("!{module.flags}")
flags+=(
    "--CPU !{threads_s1}"
    "--max_memory !{module.memory_s1}G"
)

echo "> IDS"
echo "${ids}"

echo "> META"
echo "!{meta}"

echo "> FILES"
ls -l

echo "> COMMAND"
echo !{bin} ${flags[@]} \
    --genome_guided_bam ${files} \
    --genome_guided_max_intron !{params[params.data.species].max_intron_length} \
    --no_run_inchworm
