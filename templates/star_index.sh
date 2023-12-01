#!/usr/bin/env bash -l

files=$(echo !{files} | tr -d '[,]')
flags=(!{module.flags})
flags+=()
bin=(!{bin})

${bin[@]} ${flags[@]} \
        --genomeFastaFiles genome.fasta \
        --genomeDir !{params.data.species} \
        --quantMode GeneCounts \
        --sjdbGTFfile genome.gtf \
        --sjdbOverhang !{overhang} \
        --sjdbInsertSave All \
        --genomeSAindexNbases !{index_n} \
        --runThreadN !{threads_index} \
        --runMode genomeGenerate \
        --limitGenomeGenerateRAM !{ram_limit}
