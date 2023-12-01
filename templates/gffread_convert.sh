#!/usr/bin/env bash -l

files=$(echo !{files} | tr -d '[,]')
flags=(!{module.flags})
flags+=(
    -F # All GFF attributes
    -T # to GTF format
)

!{bin} ${flags[@]} -o "genome.gtf" "${files}"
