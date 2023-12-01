#!/usr/bin/env bash -l

files=$(echo !{files} | tr -d '[,]')
flags=(!{module.flags})
flags+=(-j !{threads})

!{bin} ${flags[@]} stats -T ${files} | tail -n 1 | awk '{print $5}' > length.txt
