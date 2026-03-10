/* 
Process converting gff to gtf file for processing with STARsolo
--> every transcript will be annotated as exon and all former CDS and exons as transcripts
--> STARsolo will count for exon, so every transcript needs to be named as such */ 

process CONVERT_TO_GTF { 
    
    conda "${params.conda_path}/smRNAseq" 

    tag "Converting genome.gff to genome.full.gtf file for STAR" 
    publishDir "$params.outputDir/reference", pattern: "*.gtf", mode: params.pubDirMode 
    
    input: 
    path(gffFile) 
    
    output: 
    path("genome.full.gtf"), emit: gtfFile 
    
    script:
    """
    gffread $gffFile -T -v -o genome.gtf

    awk 'BEGIN{OFS="\\t"} \$3=="exon" {\$3="CDS"} {print}' genome.gtf > genome.CDS.gtf
    awk 'BEGIN{OFS="\\t"} \$3=="transcript" {\$3="exon"} {print}' genome.CDS.gtf > genome.exon.gtf
    awk 'BEGIN{OFS="\\t"} \$3=="CDS" {\$3="transcript"} {print}' genome.exon.gtf > genome.full.gtf
    """

    } 

/* 
Process generating reference (index) for STARsolo from gtf file and fasta file 
--> the parameter --genomeSAindexNbases 10 may need to be adjusted 
based on genome size (especially when using multiple genomes)
*/

process CREATE_STARSOLO_REFERENCE {
    
    conda "${params.conda_path}/smRNAseq" 
    
    tag "Creating reference for STARsolo" 
    publishDir "$params.outputDir/reference_star", pattern: "*", mode: params.pubDirMode 
    
    input: 
    path gtfFile 
    path genomeFasta 
    
    output:
    path "reference_star", emit: referenceSTAR 
    
    script: 
    """
    STAR --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir reference_star \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles $genomeFasta \
    --sjdbGTFfile $gtfFile
    """ 
    }
    
/* 
Process mapping and counting reads with STARsolo 
--> UMI and barcode specifications are given for read R2; barcode files are crucial to provide 
--> read R1 actually mapped to genome 
*/ 

process STARSOLO_MICROSPLIT { 
    
    conda "${params.conda_path}/smRNAseq" 
    
    tag "Mapping and counting mapped reads for ${sampleID}" 
    publishDir "$params.outputDir/alignment", pattern: "*", mode: params.pubDirMode
    publishDir "$params.outputDir/alignment", pattern: "*Log.final.out", mode: params.pubDirMode 
    
    input: 
    tuple val(sampleID), path(reads) 
    each path(referenceSTAR)
    
    output: 
    tuple val(sampleID), path("${sampleID}Solo.out/GeneFull/raw/"), emit: countsSTARsolo
    path("${sampleID}Log.final.out"), emit: stats
    
    script:
    """
    STAR --genomeDir $referenceSTAR \
    --runThreadN 16 \
    --readFilesIn ${reads[0]} ${reads[1]} \
    --readFilesCommand gunzip -c \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
    --soloUMIposition 0_0_0_9 \
    --soloCBwhitelist ${projectDir}${params.barcodes_microSPLiT_R3} ${projectDir}${params.barcodes_microSPLiT_R2} ${projectDir}${params.barcodes_microSPLiT_R1} \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_All \
    --soloFeatures Gene GeneFull \
    --soloMultiMappers Uniform \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 50 \
    --alignSJDBoverhangMin 1000 \
    --alignSJoverhangMin 1000 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${sampleID}
    """ 
    
    } 

/*
Generates AnnData object from barcodes, features and UniqueAndMult-Uniform.mtx data */

process CREATE_ANNDATA { 
    
    conda "${params.conda_path}/smRNAseq" 
    
    tag "Converting count data for ${sampleID} to h5ad file"
    publishDir "$params.outputDir/anndata/", pattern: "*.h5ad", mode: params.pubDirMode 
    
    input: 
    tuple val(sampleID), path(data_dir) 
    
    output: 
    path("${sampleID}_anndata.h5ad"), emit: anndata 
    
    script: 
    """ 
    rm $data_dir/matrix.mtx 
    mv $data_dir/UniqueAndMult-Uniform.mtx $data_dir/matrix.mtx 
    gzip $data_dir/* 
    python ${projectDir}/bin/convert_mtx_to_h5ad.py \
    --countsDirectory $data_dir \
    --sampleID $sampleID 
    """
    }