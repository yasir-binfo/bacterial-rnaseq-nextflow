#!/usr/bin/env nextflow

/**

BaNeX - Bacterial Neutral Sites Prediction using multi-Omics

**/

process Trimming_PE {

    publishDir "Omics/Expression/Trimming", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    val q_score

    output:
    tuple val("${sample_id}"), val("${condition}"), path("$sample_id*R1.fastq*"), path("$sample_id*R2.fastq*"), emit: trimmed_reads
    path "${sample_id}_fastp_report.html"

    script:
    """
    fastp -i ${read1} -I ${read2} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \
    --cut_mean_quality $q_score --detect_adapter_for_pe --cut_front --trim_poly_g -l 50 -w $task.cpus -h ${sample_id}_fastp_report.html
    """

}

process Trimming_SE {
    publishDir "Omics/Expression/Trimming", mode: 'copy'

    input:
    tuple val(sample_id), val(type), path(read1)
    val q_score

    output:
    tuple val(sample_id), val(type), path("*_trimmed_SE.fastq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${read1} \
          -o ${sample_id}_trimmed_SE.fastq.gz \
          --cut_mean_quality $q_score \
          --cut_front --trim_poly_g \
          -l 50 \
          -w $task.cpus \
          -h ${sample_id}_fastp_report.html
    """
}


process RRNA_Filtration_PE {

    beforeScript 'rm -rf ~/sortmerna/'

    publishDir "Omics/Expression/Clean_Reads", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(trimmed_read1), path(trimmed_read2)

    output:
    tuple val("${sample_id}"), val("${condition}"), path("${sample_id}_non_rRNA_reads_fwd*"), path("${sample_id}_non_rRNA_reads_rev*"), emit: clean_reads
    val("${sample_id}"), emit: sample_id

    when:
    !params.skip_rrna

    script:
    """
    sortmerna -ref ${params.database1} -ref ${params.database2} -ref ${params.database3} -ref ${params.database4} -ref ${params.database5} \
    -ref ${params.database6} -reads ${trimmed_read1} -reads ${trimmed_read2} --fastx --aligned rRNA_reads --other ${sample_id}_non_rRNA_reads \
    --paired_out --out2 --threads $task.cpus
     """

}

process RRNA_Filtration_SE {

    beforeScript 'rm -rf ~/sortmerna/'

    publishDir "Omics/Expression/Clean_Reads", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(trimmed_read1)

    output:
    tuple val("${sample_id}"), val("${condition}"), path("${sample_id}_non_rRNA_reads*"), emit: clean_reads
    val("${sample_id}"), emit: sample_id

    when:
    !params.skip_rrna

    script:
    """
    sortmerna -ref ${params.database1} -ref ${params.database2} -ref ${params.database3} -ref ${params.database4} -ref ${params.database5} \
    -ref ${params.database6} -reads ${trimmed_read1} --fastx --aligned rRNA_reads --other ${sample_id}_non_rRNA_reads \
    --threads $task.cpus
     """

}

process Index {

    publishDir "Omics/Expression/Index", mode: 'copy'

    input:
    path ref_file

    output:
    path "index*"

    script:
    """
    bowtie2-build $ref_file index
    """

}

process Mapping_PE {

    publishDir "Omics/Expression/Mapping", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)

    output:
    tuple val("${sample_id}"), val("PE"), val("${condition}"), path("${sample_id}.sorted.bam"), emit: bam_file

    script:
    """
    bowtie2 -p ${task.cpus} -x ${params.index_path}/index -1 ${read1} -2 ${read2} -S ${sample_id}.sam && \
    samtools view -bS ${sample_id}.sam | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam && \
    rm ${sample_id}.sam
    """
}

process Mapping_SE {

    publishDir "Omics/Expression/Mapping", mode: 'copy'

    input:
    tuple val(sample_id), val(cond), path(read1)

    output:
    tuple val(sample_id), val("SE"), val(cond), path("${sample_id}.sorted.bam"), emit: bam_file

    script:
    """
    bowtie2 -p $task.cpus -x ${params.index_path}/index -U ${read1} -S ${sample_id}.sam && \
    samtools view -bS ${sample_id}.sam | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam && \
    rm ${sample_id}.sam
    """
}

process Estimate_Strandedness {

    publishDir "Omics/Expression/Strandedness", mode: 'copy'

    input:
    tuple val(sample_id), val(orientation), val(condition), path(bam_file)

    output:
    tuple val("${sample_id}"), val(orientation), val("${condition}"), path("strand_info_${sample_id}.txt"), path("${bam_file}"), emit: strand_and_bam

    script:
    """
    name=\$(basename ${params.anno_file_relative_path} .gff) &&
    python3 ${params.gff_to_bed_script} -i ${params.anno_file_relative_path} -o \$name.bed &&
    infer_experiment.py -r \$name.bed -i $bam_file > strand_info_${sample_id}.txt
    """
}

process Quantification_PE {

    errorStrategy 'ignore'

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(orientation), val(condition), path(strand_info_file), path(bam_file)

    output:
    tuple val("${sample_id}"), val(orientation), val("${condition}"), path("${sample_id}.results"), path("${sample_id}.results.summary"), path("$bam_file"), emit: feature_count


    script:
    """
    strand=\$(python3 ${params.strandedness_script} ${strand_info_file})
    echo "Detected strand: \$strand"
    featureCounts -a ${params.anno_file_relative_path} -o ${sample_id}.results -t gene,pseudogene -g locus_tag -M -O -Q 10 -s \$strand -p --countReadPairs -T $task.cpus $bam_file
    """

}

process Quantification_SE {

    errorStrategy 'ignore'

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(orientation), val(condition), path(strand_info_file), path(bam_file)

    output:
    tuple val("${sample_id}"), val(orientation), val("${condition}"), path("${sample_id}.results"), path("${sample_id}.results.summary"), path("$bam_file"), emit: feature_count


    script:
    """
    strand=\$(python3 ${params.strandedness_script} ${strand_info_file})
    echo "Detected strand: \$strand"
    featureCounts -a ${params.anno_file_relative_path} -o ${sample_id}.results -t gene,pseudogene -g locus_tag -M -O -Q 10 -s \$strand -T $task.cpus $bam_file
    """

}

process Normalization_PE {

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(orientation), val(condition), path(feature_count_result), path(feature_count_result_summary), path(bam_file)

    output:
    tuple val("${sample_id}"), val(orientation), val("${condition}"), path("${sample_id}_normalized.results"), path("${sample_id}_fragment_len_metrics.txt"), emit: normalization


    script:
    """
    echo -e "Count\tFragment_Length" > ${sample_id}_fragment_len_metrics.txt && \
    samtools view -f 0x2 $bam_file | awk '(\$9 > 0) {print \$9}' | sort | uniq -c | sort -nr | awk '{print \$1"\\t"\$2}' >> ${sample_id}_fragment_len_metrics.txt &&
    python3 ${params.tpm_normalization_script} ${sample_id} ${condition} ${sample_id}_fragment_len_metrics.txt $feature_count_result ${sample_id}_normalized.results
    """


}

process Normalization_SE {

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(orientation), val(condition), path(feature_count_result), path(feature_count_result_summary), path(bam_file)

    output:
    tuple val("${sample_id}"), val(orientation), val("${condition}"), path("${sample_id}_normalized.results"), path("${sample_id}_read_len_metrics.txt"), emit: normalization


    script:
    """
    echo -e "Count\tFragment_Length" > ${sample_id}_read_len_metrics.txt && \
    samtools view $bam_file | awk '{print length(\$10)}' | sort | uniq -c | sort -nr | awk '{print \$1"\\t"\$2}' >> ${sample_id}_read_len_metrics.txt &&
    python3 ${params.tpm_normalization_script} ${sample_id} ${condition} ${sample_id}_read_len_metrics.txt $feature_count_result ${sample_id}_normalized.results
    """

}

process Mean_TPM {

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    val normalized_results

    output:
    path('Expression.csv'), emit: mean_tpm

    script:
    """
    python3 ${params.mean_tpm_calc} ${normalized_results.join(' ')} 'Expression.csv'
    """
}

process BigWig {

    errorStrategy 'ignore'

    publishDir "Omics/Expression/Quantification", mode: 'copy'

    input:
    tuple val(condition), val(sample_id), path(bam_file)

    output:
    path("${condition}.bw"), emit: bigwig_per_condition

    script:
    """
    samtools merge -f ${condition}.bam ${bam_file.join(' ')}
    samtools index ${condition}.bam
    bamCoverage -b ${condition}.bam -o ${condition}.bw --normalizeUsing CPM
    """
}


params.datasets = "datasets.csv"
params.q_score = 5
params.databases_path = 'Databases'
params.database1 = "../../../${params.databases_path}/rfam-5.8s-database-id98.fasta"
params.database2 = "../../../${params.databases_path}/rfam-5s-database-id98.fasta"
params.database3 = "../../../${params.databases_path}/silva-bac-16s-id90.fasta"
params.database4 = "../../../${params.databases_path}/silva-bac-23s-id98.fasta"
params.database5 = "../../../${params.databases_path}/silva-arc-16s-id95.fasta"
params.database6 = "../../../${params.databases_path}/silva-arc-23s-id98.fasta"
params.ref_file = 'ref'
params.anno_file = 'gff'
params.anno_file_relative_path = "../../../${params.anno_file}"
params.index_path = "../../../Omics/Expression/Index"
params.strandedness_script = 'scripts/infer_strandedness.py'
params.gff_to_bed_script = 'gff_to_bed_script.py'
params.gff_to_saf_script = 'gff_to_saf_script.py'
params.tpm_normalization_script='tpm_calc.py'
params.mean_tpm_calc = 'mean_tpm_calc.py'
// params.outdir = ''
params.skip_rrna = true

workflow {

    //reads = Channel.fromFilePairs(params.fastq_files, checkifExists: true)
    reads = Channel.fromPath(params.datasets).splitCsv(header:true)
    .map { row ->
        def r1 = file(row.fastq_1)
        def r2 = row.fastq_2 ? file(row.fastq_2) : null
        tuple(row.sample_id, row.condition, r1, r2)
         }

    // Separate SE and PE reads
    pe_reads = reads.filter { it[3] != null }  // paired-end
    se_reads = reads.filter { it[3] == null }  // single-end

    // Trim reads
    trimmed_pe = Trimming_PE(pe_reads.map { s_id, cond, r1, r2 -> tuple(s_id, cond, r1, r2) }, params.q_score)
    trimmed_se = Trimming_SE(se_reads.map { s_id, cond, r1, r2 -> tuple(s_id, cond, r1) }, params.q_score)



    database1 = Channel.fromPath(params.database1)
    database2 = Channel.fromPath(params.database2)
    database3 = Channel.fromPath(params.database3)
    database4 = Channel.fromPath(params.database4)
    database5 = Channel.fromPath(params.database5)
    database6 = Channel.fromPath(params.database6)
    ref_file = Channel.fromPath(params.ref_file)


//     // Conditionally run RRNA_Filtration
    if (!params.skip_rrna) {
        clean_out_pe = RRNA_Filtration_PE(trimmed_pe.trimmed_reads).clean_reads
        clean_out_se = RRNA_Filtration_SE(trimmed_se.trimmed_reads).clean_reads
    } else {
        // Bypass: just rename for consistency
        clean_out_pe = trimmed_pe.trimmed_reads
        clean_out_se = trimmed_se.trimmed_reads
    }

    Index(ref_file)

    mapped_pe = Mapping_PE(clean_out_pe)
    mapped_se = Mapping_SE(clean_out_se)
    mapped_bam = mapped_pe.mix(mapped_se)

    strandedness = Estimate_Strandedness(mapped_bam)

    pe = strandedness.filter { it[1] == 'PE'}
    se = strandedness.filter { it[1] == 'SE'}

    quant_pe = Quantification_PE(pe)
    quant_se = Quantification_SE(se)

    norm_pe = Normalization_PE(quant_pe)
    norm_se = Normalization_SE(quant_se)
    norm_out = norm_pe.mix(norm_se)

    all_normalized_files = norm_out
    .map { sample_id, orientation, condition, norm_file, frag_len_file -> norm_file }
    .collect()
    .set { all_normalized_results }

    Mean_TPM(all_normalized_results)

    bam_grouped = mapped_bam
    .map { sample_id, layout, condition, bam_path ->
        tuple(condition, sample_id, bam_path)        // (key, value)
    }
    .groupTuple()

    BigWig(bam_grouped)

 }
