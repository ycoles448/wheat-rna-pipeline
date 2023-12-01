#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
files=()
for f in ${ids}; do
    files+=(${f}_!{params.data.bam_star_sorted})
done
flags=(!{module.flags})
flags+=(
    # Number of threads
    # Fixed at 2 when buffering is enabled (default behaviour), samtools spawns n + 2 threads
    $(if [ !{module.buffer} -gt 1 ]; then echo "-@2"; else echo "-@!{threads}"; fi)

    # Enable compression of files
    "-1"
)

if [ !{module.buffer} -gt 1 ]; then
    # Replace parallel implementation with staged variants
    bin_s1=$(echo !{bin} | sed 's/!{params.parallel.bin}/!{params.parallel.bin} -N !{module.buffer} -m --files/g')
    bin_s2=$(echo !{bin} | sed 's/!{params.parallel.bin}/!{params.xargs.bin}/g')

    flags_s1=${flags[@]}
    flags_s2=$(echo ${flags[@]} | sed 's/\-@2/-@!{threads}/g')

    ${bin_s1} !{module.bin_merge} ${flags_s1} -o - {} ::: ${files[@]} |
        ${bin_s2} !{module.bin_merge} ${flags_s2} -o '!{params.data.species}_!{meta[params.meta.group][0]}.bam'
else
    !{bin} !{module.bin_merge} ${flags[@]} \
        -o '!{params.data.species}_!{meta[params.meta.group][0]}.bam' \
        ${files[@]}
fi
