#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
flags=(!{module.flags})
flags+=(
    -w 2 # Two cores per file
)

if [[ "!{module.autodetect}" == "True" ]]; then
    flags+=(--detect_adapter_for_pe)
fi

!{bin} ${flags[@]} \
    -o trim_{}_r1.fastq -O trim_{}_r2.fastq \
    -i {}_r1.fastq -I {}_r2.fastq \
    -f !{params.qc.trimf} -F !{params.qc.trimf} \
    -t !{params.qc.trimt} -T !{params.qc.trimt} \
    -q !{params.qc.quality} \
    ::: ${ids}
