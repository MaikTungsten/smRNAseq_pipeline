
process STARSOLO_PETRISEQ { 
    
    conda "${params.conda_path}/smRNAseq" 
    
    tag "Mapping and counting mapped reads for ${sampleID}" 
    publishDir "$params.outputDir/alignment", pattern: "*", mode: params.pubDirMode 
    
    input: 
    tuple val(sampleID), path(reads) 
    each path(referenceSTAR)
    
    output: 
    tuple val(sampleID), path("${sampleID}Solo.out/GeneFull/raw/"), emit: countsSTARsolo
    path("${sampleID}Log.final.out"), emit: stats

    // --readFilesIn ${reads[1]} ${reads[0]} für merkwürdiges sample hier 0, 1 tauschen und austesten
    
    script:
    """
    STAR --genomeDir $referenceSTAR \
    --runThreadN 16 \
    --readFilesIn ${reads[1]} ${reads[0]} \
    --readFilesCommand gunzip -c \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_7_0_12 0_29_0_34 0_51_0_56 \
    --soloUMIposition 0_0_0_6 \
    --soloCBwhitelist ${projectDir}${params.barcodes_PETRIseq_R3} ${projectDir}${params.barcodes_PETRIseq_R2} ${projectDir}${params.barcodes_PETRIseq_R1} \
    --soloCBmatchWLtype Exact \
    --soloUMIdedup 1MM_All \
    --soloFeatures Gene GeneFull \
    --soloMultiMappers Uniform \
    --outFilterScoreMinOverLread 0.6 \
    --outFilterMatchNminOverLread 0.6 \
    --outFilterMatchNmin 10 \
    --outFilterMultimapNmax 5 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${sampleID}
    """ 
    
    } 