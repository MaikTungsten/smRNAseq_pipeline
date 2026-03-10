
process FASTQC_RAW {
    
    conda "${params.conda_path}/smRNAseq" 

    tag "FASTQC on ${sample_id}" 
    publishDir "$params.outputDir/fastQC", pattern: "*", mode: params.pubDirMode 

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}", emit: fastqc_pre

    script:
    """
    mkdir fastqc_${sample_id}
    fastqc -o fastqc_${sample_id} -q ${reads}
    """
}

process BBDUK_TRIM {
  
    conda "${params.conda_path}/smRNAseq" 

    tag "FASTQC on ${sample_id}" 
    publishDir "$params.outputDir/fastQC", pattern: "*", mode: params.pubDirMode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fq.gz")

    script:
    // Two-pass trimming: left primers + polyA; then right adapters/primers
    """
    bbduk.sh ${params.bbduk_opts} \\
        in=${reads} out=${sample_id}.L.fq.gz \\
        ref="${projectDir}${params.bbduk_ref_left}" ktrim=l trimpolya=${params.bbduk_polyA}

    bbduk.sh ${params.bbduk_opts} \\
        in=${sample_id}.L.fq.gz out=${sample_id}.trimmed.fq.gz \\
        ref="${projectDir}${params.bbduk_ref_right}" ktrim=r

    rm -f ${sample_id}.L.fq.gz
    """
}

process FASTQC_TRIM {

    conda "${params.conda_path}/smRNAseq" 

    tag "FASTQC on trimmed ${sample_id}" 
    publishDir "$params.outputDir/fastQC", pattern: "*", mode: params.pubDirMode 

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_trimmed", emit: fastqc_post

    script:
    """
    mkdir fastqc_${sample_id}_trimmed
    fastqc -o fastqc_${sample_id}_trimmed -q ${reads}
    """


}

// ---- indexing ----

process BUILDHISAT2BASE {

    tag "Build dual reference genome and Hisat2Base for alignment"
    publishDir "$params.outputDir/alignments", pattern: "*.ht2", mode: params.pubDirMode

    conda "${params.conda_path}/smRNAseq"

    input:
    path genome

    output:
    path "Hisat2Base", emit: alignmentBase

    script:
    """
    mkdir Hisat2Base
    hisat2-build -f -p 16 $genome Hisat2Base/Hisat2Base
    """
}

process STAR_BUILD {
    
    conda "${params.conda_path}/smRNAseq" 
    
    tag "Creating reference for STAR" 
    publishDir "$params.outputDir/reference_star", pattern: "*", mode: params.pubDirMode 
    
    input:
    path genomeFasta 
    path gtfFile 
    
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


// ---- mapping ----
process HISAT2_MAP {

    conda "${params.conda_path}/smRNAseq"

    tag "Mapping sample ${sampleID} with HISAT2"
    publishDir "$params.outputDir/alignments", pattern: "*.hisat2.summary.log", mode: params.pubDirMode
    
    input:
    tuple val(sampleID), path(reads)
    each path(alignmentBase)

    output:
    tuple val(sampleID), path("${sampleID}.sam"), emit: alignment
    path "${sampleID}.hisat2.summary.log", emit: summaryFile

    script:
    """
    hisat2 -p 8 -x ${alignmentBase}/Hisat2Base \
    -U ${reads} \
    --summary-file ${sampleID}.hisat2.summary.log \
    -S ${sampleID}.sam
    """

}

process STAR_MAP {

    conda "${params.conda_path}/smRNAseq" 
    
    tag "Mapping ${sample_id} with STAR"
    publishDir "$params.outputDir/alignments", pattern: "*.Log.final.out", mode: params.pubDirMode

    input:
    tuple val(sample_id), path(reads)
    each path(star_index_dir)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: alignment
    path("${sample_id}.Log.final.out"), emit: stats


    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${star_index_dir} \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}. \
        --outSAMtype SAM \
        --outFilterMultimapNmax 10 \
        --outFilterMatchNmin 50

    mv ${sample_id}.Aligned.out.sam ${sample_id}.sam
    """
}

// ---- SAM -> sorted BAM + index ----
process SAMTOOLS_SORT_INDEX {

    conda "${params.conda_path}/smRNAseq" 
    
    tag "Converting sam to bam and indexing bam for ${sample_id}" 
    publishDir "$params.outputDir/alignments", pattern: "*_sorted*", mode: params.pubDirMode

    input:
    tuple val(sample_id), path(samFile)

    output:
    path("${sample_id}_sorted.bam"), emit: sortedBamFile
    path("${sample_id}_sorted.bam.bai"), emit: bam_index

    script:
    """
    samtools view -bS $samFile | samtools sort -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}

// ---- featureCounts ----
process FEATURECOUNTS {

    conda "${params.conda_path}/smRNAseq"

    tag "Count reads for all samples"
    publishDir "$params.outputDir/countData", pattern: "*.tsv", mode: params.pubDirMode
    publishDir "$params.outputDir/countData", pattern: "*.tsv.summary", mode: params.pubDirMode
    
    input:
    path sortedBamFiles
    path GFF

    output:
    path "countData.tsv", emit: countTable
    path "countData.tsv.summary", emit: countSummary    

    script:
    """
    featureCounts -M -O --fraction \\
    -a $GFF \\
    -t $params.countFeature \\
    -o countData.tsv \\
    -g $params.featureIdentifier \\
    $sortedBamFiles
    """

}

process MULTIQC {

    conda "${params.conda_path}/smRNAseq"

    tag "MultiQC analysis on entire run"
    publishDir "${params.outputDir}/multiqc", mode: params.pubDirMode

    input:
    val outdir
    path countsTable

    output:
    path "multiqc_report.html", emit: multiqc_report

    script:
    """
    multiqc $outdir --filename multiqc_report.html
    """
}
