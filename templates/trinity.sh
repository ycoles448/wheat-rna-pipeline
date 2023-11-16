#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
flags=("!{module.flags}")
flags+=()

echo ${ids}
echo !{files}
