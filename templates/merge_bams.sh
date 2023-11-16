#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
flags=("!{module.flags}")
flags+=(
    # Disable compression
    "-u"

    # Number of threads
    "-@!{threads}"
)
files=()

for f in ${ids}; do
    files+=(${f}_!{params.data.bam_star_sorted})
done

!{bin} !{module.bin_merge} ${flags[@]} \
    -o '!{params.data.species}_!{meta[params.meta.group][0]}.bam' \
    ${files[@]}
