#!/usr/bin/env nextflow

// Prepare module
include { PREPARE_GENOME } from './modules/prepare'

// microSPLiT modules
include { CONVERT_TO_GTF } from './modules/microSPLiT'
include { CREATE_STARSOLO_REFERENCE } from './modules/microSPLiT'
include { STARSOLO_MICROSPLIT } from './modules/microSPLiT'
include { CREATE_ANNDATA } from './modules/microSPLiT'

// MATQ-seq modules
include { FASTQC_RAW } from './modules/MATQseq'
include { BBDUK_TRIM } from './modules/MATQseq'
include { FASTQC_TRIM } from './modules/MATQseq'
include { BUILDHISAT2BASE } from './modules/MATQseq'
include { STAR_BUILD } from './modules/MATQseq'
include { HISAT2_MAP } from './modules/MATQseq'
include { STAR_MAP } from './modules/MATQseq'
include { SAMTOOLS_SORT_INDEX } from './modules/MATQseq'
include { FEATURECOUNTS } from './modules/MATQseq'
include { MULTIQC } from './modules/MATQseq'

// PETRIseq modules
include { STARSOLO_PETRISEQ } from './modules/PETRIseq'


def helpMessage() {
    log.info"""
    ==-------------------------------------------------------------------==
    Pipeline for processing single-microbe RNA-seq data
    ==-------------------------------------------------------------------==

    ==-------------------------------------------------------------------==
    Protocols available:
    - microSPLiT
    - PETRI-seq
    - MATQ-seq
    ==-------------------------------------------------------------------==

    ==-------------------------------------------------------------------==
    Typical usage at the example of MATQ-seq:

    nextflow run process_smRNAseq --reads '*.fastq.gz' --matqseq \
    --genomeFasta "genome.fasta" --genomeGFF "genome.gff" \
    --use_hisat2 --outputDir results_MATQseq
    ==-------------------------------------------------------------------==
    
    ==-------------------------------------------------------------------==
    Typical usage at the example of microSPLiT using two reference genomes
    from two different species:

    nextflow run process_smRNAseq --reads '*_{1,2}.fastq.gz' --microsplit \
    --genomeFasta "genome1.fasta" "genome2.fasta" \
    --genomeGFF "genome1.gff" "genome2.gff" \
    --outputDir results_microSPLiT
    ==-------------------------------------------------------------------==
    
    Arguments:

    // Defaults
    --pubDirMode                  = "copy"

    // General input
    genomeFASTA                 = []            // at least one or multiple fasta files of bacterial/phage genomes
    genomeGFF                   = []            // at least one or multiple gff files of bacterial/phage genomes in same order as fasta files
    reads                       = "*.fastq.gz"  // reads to be processed and mapped (paired-end for PETRIseq and microSPLiT; single-end for MATQ-seq) 
    outputDir                   = "output"      // output directory for pipeline results

    // Flag for protocol-specific processing
    // Select depending on data supplied
    // Then, specific parameters and commands will be executed for processing (see below)
    microsplit                  = false
    matqseq                     = false
    petriseq                    = false

    // microSPLiT-specific input
    barcodes_microSPLiT_R1      = "/bin/barcodes_microSPLiT_R1.txt"
    barcodes_microSPLiT_R2      = "/bin/barcodes_microSPLiT_R2.txt"
    barcodes_microSPLiT_R3      = "/bin/barcodes_microSPLiT_R3.txt"

    // MATQ-seq
    use_hisat2                  = false // choose hisat2 or star for mapping reads
    use_star                    = false
    bbduk_opts                  = "-Xmx16g t=8 minlen=18 qtrim=rl trimq=20 k=17 mink=11 hdist=1"
    bbduk_ref_left              = "/bin/matseq_primers.fa"
    bbduk_ref_right             = "/bin/nextera_and_primers.fa"
    bbduk_polyA                 = 30
    countFeature                = "gene"
    featureIdentifier           = "ID"

    // PETRI-seq
    barcodes_PETRIseq_R1        = "/bin/barcodes_PETRIseq_R1.txt"
    barcodes_PETRIseq_R2        = "/bin/barcodes_PETRIseq_R2.txt"
    barcodes_PETRIseq_R3        = "/bin/barcodes_PETRIseq_R3.txt"

    // Conda environments; the environment containing all required packages such be named "smRNAseq"; otherwise rename in processes accordingly
    conda_path                  = "/ceph/ibmi/it/projects/ML_BI/Mamba/envs"

    Arguments can also be set in the config file (./conf/params.config)
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


workflow {
    fasta_ch = Channel
        .fromPath(params.genomeFASTA)
        .ifEmpty { error "No genome fasta files provided via --genomeFASTA" }
        .collect()
    
    gff_ch = Channel
        .fromPath(params.genomeGFF)
        .ifEmpty { error "No genome gff files provided via --genomeGFF" }
        .collect()

    def reference_ch = PREPARE_GENOME(fasta_ch, gff_ch)

    if (params.microsplit){

        reads_ch = channel.fromFilePairs(params.reads, size: 2, checkIfExists: true)

        gtf_file = CONVERT_TO_GTF(reference_ch.GFF)

        star_ref = CREATE_STARSOLO_REFERENCE(gtf_file, reference_ch.FASTA)

        star_map = STARSOLO_MICROSPLIT(reads_ch, star_ref.referenceSTAR)

        CREATE_ANNDATA(star_map.countsSTARsolo)

    }

    if (params.matqseq) {

        if (params.use_hisat2) { 
        
            reads_ch = Channel.fromFilePairs(params.reads, size: -1, checkIfExists: true)
            
            fastqc_pre_ch = FASTQC_RAW(reads_ch)
            trimmed_reads_ch = BBDUK_TRIM(reads_ch)
            fastqc__post_ch = FASTQC_TRIM(trimmed_reads_ch)

            hisat_index = BUILDHISAT2BASE(reference_ch.FASTA)
            mappings_ch = HISAT2_MAP(trimmed_reads_ch, hisat_index.alignmentBase)
            bam_ch = SAMTOOLS_SORT_INDEX(mappings_ch.alignment)
            
        }
        
        if (params.use_star) { 
            
            reads_ch = channel.fromFilePairs(params.reads, size: -1, checkIfExists: true) 
            
            fastqc_pre_ch = FASTQC_RAW(reads_ch) 
            trimmed_reads_ch = BBDUK_TRIM(reads_ch) 
            FASTQC_TRIM(trimmed_reads_ch) 
            
            gtf_file = CONVERT_TO_GTF(reference_ch.GFF) 
            star_index = STAR_BUILD(reference_ch.FASTA, gtf_file)
            STAR_MAP(trimmed_reads_ch, star_index.referenceSTAR) 
            bam_ch = SAMTOOLS_SORT_INDEX(STAR_MAP.out.alignment)
        } 
        
        fc_out = FEATURECOUNTS(bam_ch.sortedBamFile.collect(), reference_ch.GFF) 
        
        def outdir_abs = file(projectDir).resolve("../${params.outputDir}")

        MULTIQC(outdir_abs, fc_out.countTable)
        
    }

    if (params.petriseq) {

        reads_ch = channel.fromFilePairs(params.reads, size: 2, checkIfExists: true)

        gtf_file = CONVERT_TO_GTF(reference_ch.GFF)

        star_ref = CREATE_STARSOLO_REFERENCE(gtf_file, reference_ch.FASTA)

        star_map = STARSOLO_PETRISEQ(reads_ch, star_ref.referenceSTAR)

        CREATE_ANNDATA(star_map.countsSTARsolo)

    }
}