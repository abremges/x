#!/usr/bin/env nextflow

params.assembly = "assembly.fasta"
assembly = file(params.assembly)

params.reads = "fastq"
Channel
    .fromFilePairs("${params.reads}/*__R{1,2}.fastq", flat: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }
    .into {
        input_reads;
        bowtie2_reads;
        kallisto_reads
    }
input_reads.println()


process bowtie2_build {
    publishDir = "bowtie2"

    input:
    file assembly

    output:
    set "${assembly}", "${assembly}.*.bt2" into bowtie2_index

    """
    bowtie2-build ${assembly} ${assembly}
    """
}

process bowtie2 {
    tag "${id}"
    publishDir = "bowtie2"

    input:
    set file(index_name), file(index) from bowtie2_index
    set val(id), file(forward), file(reverse) from bowtie2_reads

    output:
    set id, "${id}.bam" into bowtie2_bam

    """
    bowtie2 -x ${index_name} -1 ${forward} -2 ${reverse} | samtools view -bSu - | samtools sort - ${id}
    """
}

process jgi_summarize_bam_contig_depths {
    tag "${id}"
    publishDir = "bowtie2"

    input:
    set val(id), file(bam) from bowtie2_bam

    output:
    set id, "${id}.depth"

    """
    jgi_summarize_bam_contig_depths --noIntraDepthVariance --includeEdgeBases --outputDepth ${id}.depth ${bam}
    """
}


process kallisto_index {
    publishDir = "kallisto"
    memory '500 GB'

    input:
    file assembly

    output:
    file "${assembly}.idx" into kallisto_index

    """
    kallisto index -i ${assembly}.idx ${assembly}
    """
}

process kallisto_quant {
    tag "${id}"
    publishDir = "kallisto"
    cpus = 8
    memory '10 GB'

    input:
    file index from kallisto_index
    set val(id), file(forward), file(reverse) from kallisto_reads

    output:
    file "${id}/abundance.h5"
    file "${id}/abundance.tsv"
    file "${id}/run_info.json"

    """
    kallisto quant -t ${task.cpus} -i ${index} -o ${id} ${forward} ${reverse}
    """
}


// process samtools_index {
//     tag "${assembly.baseName}"
//
//     input:
//     file assembly
//
//     output:
//     file "${assembly}.fai" into samtools_index
//
//     """
//     echo samtools index ${assembly}
//     touch ${assembly}.fai
//     """
// }
