#!/usr/bin/env bash -l

!{params.parallel.bin} '!{bin} -cdk -p 1 {} > {.}' ::: *.fastq.gz
