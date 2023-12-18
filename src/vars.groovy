// Nextflow runtime variables
species = "${params.data.species}"
factors = params.meta.factors
group = params.meta.group
meta_f = "${params.data.meta}"
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
reads_unzipped = "${params.data.path}/${params.data.reads_unzipped}/${params.data.glob_unzipped}"
reads_trimmed = "${params.data.path}/${params.data.reads_trimmed}/${params.data.glob_trimmed}"
reads_bam = "${params.data.path}/${params.data.star}/*_${params.data.species}/*_${params.data.bam_star_sorted}"
reads_bam_grouped = "${params.data.path}/${params.data.bam}/${species}_*"
star_index_f = "${params.data.path}/${params.data.star}/${species}"
genomes = "${params.data.path}/${params.data.genomes}"
genome_f = "${genomes}/${species}/genome.fasta"
genome_length_f = "${genomes}/${species}/length.txt"
genome_gff_f = "${genomes}/${species}/genome.gff"
genome_gtf_f = "${genomes}/${species}/genome.gtf"
adapters = "${params.data.path}/${params.data.adapters}"
