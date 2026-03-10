

process PREPARE_GENOME {

    tag "prepare_genome(s)"

    input:
    path genomesFASTA
    path genomesGFF

    output:
    path "genome.fasta", emit: FASTA
    path "genome.gff", emit: GFF

    script:
    """
    cat ${genomesFASTA} > genome.fasta
    cat ${genomesGFF} > genome.gff
    """
}