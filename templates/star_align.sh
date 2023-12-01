#!/usr/bin/env bash -l

ids=$(echo !{ids} | tr -d '[,]')
flags=(!{module.flags})
flags+=()
bin=(!{bin})

if [ !{module.buffer_align} -gt 1 ]; then
    ${bin[@]} ${flags[@]} \
        --readFilesIn trim_{}_r1.fastq trim_{}_r2.fastq \
        --outFileNamePrefix {}_!{params.data.species}/{}_ \
        --genomeDir !{params.data.species} \
        --quantMode TranscriptomeSAM GeneCounts \
        --sjdbGTFfile genome.gtf \
        --sjdbOverhang !{overhang} \
        --sjdbInsertSave All \
        --genomeSAindexNbases !{index_n} \
        --runThreadN !{threads_buf_align} \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN !{threads_buf_align} \
        --outReadsUnmapped Fastx \
        ::: ${ids}
else
    echo "Cannot run without parallel, please increase buffer to at least 2"
fi
