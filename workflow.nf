#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Bacterial RNA-seq preprocessing, alignment, quantification, and coverage.
 */

params.datasets = 'datasets.csv'
params.ref_file = null
params.anno_file = null
params.outdir = 'results'
params.q_score = 5
params.skip_rrna = true
params.rrna_16s = null
params.rrna_23s = null
params.rrna_5s = null

process TRIM_PE {
    tag "${sample_id}"
    label 'small'
    publishDir "${params.outdir}/trimming", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    val q_score

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
    path "${sample_id}_fastp_report.html", emit: report

    script:
    """
    fastp \
      -i ${read1} -I ${read2} \
      -o ${sample_id}_trimmed_R1.fastq.gz \
      -O ${sample_id}_trimmed_R2.fastq.gz \
      --cut_mean_quality ${q_score} \
      --detect_adapter_for_pe --cut_front --trim_poly_g \
      -l 50 -w ${task.cpus} \
      -h ${sample_id}_fastp_report.html
    """

    stub:
    """
    touch ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_trimmed_R2.fastq.gz
    touch ${sample_id}_fastp_report.html
    """
}

process TRIM_SE {
    tag "${sample_id}"
    label 'small'
    publishDir "${params.outdir}/trimming", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1)
    val q_score

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_trimmed_SE.fastq.gz"), emit: reads
    path "${sample_id}_fastp_report.html", emit: report

    script:
    """
    fastp \
      -i ${read1} \
      -o ${sample_id}_trimmed_SE.fastq.gz \
      --cut_mean_quality ${q_score} \
      --cut_front --trim_poly_g \
      -l 50 -w ${task.cpus} \
      -h ${sample_id}_fastp_report.html
    """

    stub:
    """
    touch ${sample_id}_trimmed_SE.fastq.gz ${sample_id}_fastp_report.html
    """
}

process FILTER_RRNA_PE {
    tag "${sample_id}"
    label 'large'
    publishDir "${params.outdir}/clean_reads", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    path rrna_databases

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_non_rRNA_reads_fwd*"), path("${sample_id}_non_rRNA_reads_rev*"), emit: reads

    script:
    def references = rrna_databases.collect { "-ref ${it}" }.join(' ')
    """
    sortmerna \
      ${references} \
      -reads ${read1} -reads ${read2} \
      --fastx --paired_out --out2 \
      --aligned ${sample_id}_rRNA_reads \
      --other ${sample_id}_non_rRNA_reads \
      --workdir sortmerna_work \
      --threads ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}_non_rRNA_reads_fwd.fastq ${sample_id}_non_rRNA_reads_rev.fastq
    """
}

process FILTER_RRNA_SE {
    tag "${sample_id}"
    label 'large'
    publishDir "${params.outdir}/clean_reads", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1)
    path rrna_databases

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_non_rRNA_reads*"), emit: reads

    script:
    def references = rrna_databases.collect { "-ref ${it}" }.join(' ')
    """
    sortmerna \
      ${references} \
      -reads ${read1} \
      --fastx \
      --aligned ${sample_id}_rRNA_reads \
      --other ${sample_id}_non_rRNA_reads \
      --workdir sortmerna_work \
      --threads ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}_non_rRNA_reads.fastq
    """
}

process BUILD_INDEX {
    tag 'bowtie2-index'
    label 'medium'
    publishDir "${params.outdir}/index", mode: 'copy'

    input:
    path reference

    output:
    path 'index*.bt2*', emit: index_files

    script:
    """
    bowtie2-build --threads ${task.cpus} ${reference} index
    """

    stub:
    """
    touch index.1.bt2 index.2.bt2 index.3.bt2 index.4.bt2 index.rev.1.bt2 index.rev.2.bt2
    """
}

process MAP_PE {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    path index_files

    output:
    tuple val(sample_id), val('PE'), val(condition), path("${sample_id}.sorted.bam"), emit: bam

    script:
    """
    bowtie2 -p ${task.cpus} -x index -1 ${read1} -2 ${read2} \
      | samtools view -@ ${task.cpus} -b - \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    """

    stub:
    """
    touch ${sample_id}.sorted.bam
    """
}

process MAP_SE {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1)
    path index_files

    output:
    tuple val(sample_id), val('SE'), val(condition), path("${sample_id}.sorted.bam"), emit: bam

    script:
    """
    bowtie2 -p ${task.cpus} -x index -U ${read1} \
      | samtools view -@ ${task.cpus} -b - \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    """

    stub:
    """
    touch ${sample_id}.sorted.bam
    """
}

process GFF_TO_BED {
    tag 'annotation-bed'
    label 'small'

    input:
    path annotation
    path converter

    output:
    path 'annotation.bed', emit: bed

    script:
    """
    python3 ${converter} -i ${annotation} -o annotation.bed
    """

    stub:
    """
    touch annotation.bed
    """
}

process ESTIMATE_STRANDEDNESS {
    tag "${sample_id}"
    label 'small'
    publishDir "${params.outdir}/strandedness", mode: 'copy', pattern: 'strand_info_*.txt'

    input:
    tuple val(sample_id), val(layout), val(condition), path(bam_file)
    path annotation_bed

    output:
    tuple val(sample_id), val(layout), val(condition), path("strand_info_${sample_id}.txt"), path(bam_file), emit: result

    script:
    """
    infer_experiment.py -r ${annotation_bed} -i ${bam_file} > strand_info_${sample_id}.txt
    """

    stub:
    """
    echo 'Fraction of reads explained by "++,--": 0.99' > strand_info_${sample_id}.txt
    """
}

process QUANTIFY_PE {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/quantification", mode: 'copy', pattern: '*.results*'

    input:
    tuple val(sample_id), val(layout), val(condition), path(strand_info), path(bam_file)
    path annotation
    path strandedness_script

    output:
    tuple val(sample_id), val(layout), val(condition), path("${sample_id}.results"), path("${sample_id}.results.summary"), path(bam_file), emit: counts

    script:
    """
    strand=\$(python3 ${strandedness_script} ${strand_info})
    featureCounts \
      -a ${annotation} -o ${sample_id}.results \
      -t gene,pseudogene -g locus_tag \
      -M -O -Q 10 -s \${strand} \
      -p --countReadPairs -T ${task.cpus} ${bam_file}
    """

    stub:
    """
    touch ${sample_id}.results ${sample_id}.results.summary
    """
}

process QUANTIFY_SE {
    tag "${sample_id}"
    label 'medium'
    publishDir "${params.outdir}/quantification", mode: 'copy', pattern: '*.results*'

    input:
    tuple val(sample_id), val(layout), val(condition), path(strand_info), path(bam_file)
    path annotation
    path strandedness_script

    output:
    tuple val(sample_id), val(layout), val(condition), path("${sample_id}.results"), path("${sample_id}.results.summary"), path(bam_file), emit: counts

    script:
    """
    strand=\$(python3 ${strandedness_script} ${strand_info})
    featureCounts \
      -a ${annotation} -o ${sample_id}.results \
      -t gene,pseudogene -g locus_tag \
      -M -O -Q 10 -s \${strand} \
      -T ${task.cpus} ${bam_file}
    """

    stub:
    """
    touch ${sample_id}.results ${sample_id}.results.summary
    """
}

process NORMALIZE_PE {
    tag "${sample_id}"
    label 'small'
    publishDir "${params.outdir}/quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(layout), val(condition), path(counts), path(summary), path(bam_file)
    path tpm_script

    output:
    tuple val(sample_id), val(layout), val(condition), path("${sample_id}_normalized.results"), path("${sample_id}_fragment_len_metrics.txt"), emit: result

    script:
    """
    echo -e 'Count\\tFragment_Length' > ${sample_id}_fragment_len_metrics.txt
    samtools view -f 0x2 ${bam_file} \
      | awk '(\$9 > 0) {print \$9}' \
      | sort -n | uniq -c \
      | awk '{print \$1"\\t"\$2}' >> ${sample_id}_fragment_len_metrics.txt
    python3 ${tpm_script} ${sample_id} ${condition} \
      ${sample_id}_fragment_len_metrics.txt ${counts} ${sample_id}_normalized.results
    """

    stub:
    """
    touch ${sample_id}_normalized.results ${sample_id}_fragment_len_metrics.txt
    """
}

process NORMALIZE_SE {
    tag "${sample_id}"
    label 'small'
    publishDir "${params.outdir}/quantification", mode: 'copy'

    input:
    tuple val(sample_id), val(layout), val(condition), path(counts), path(summary), path(bam_file)
    path tpm_script

    output:
    tuple val(sample_id), val(layout), val(condition), path("${sample_id}_normalized.results"), path("${sample_id}_read_len_metrics.txt"), emit: result

    script:
    """
    echo -e 'Count\\tFragment_Length' > ${sample_id}_read_len_metrics.txt
    samtools view -F 0x4 ${bam_file} \
      | awk '{print length(\$10)}' \
      | sort -n | uniq -c \
      | awk '{print \$1"\\t"\$2}' >> ${sample_id}_read_len_metrics.txt
    python3 ${tpm_script} ${sample_id} ${condition} \
      ${sample_id}_read_len_metrics.txt ${counts} ${sample_id}_normalized.results
    """

    stub:
    """
    touch ${sample_id}_normalized.results ${sample_id}_read_len_metrics.txt
    """
}

process MEAN_TPM {
    tag 'condition-means'
    label 'small'
    publishDir "${params.outdir}/quantification", mode: 'copy'

    input:
    path normalized_results
    path samplesheet
    path mean_tpm_script

    output:
    path 'Expression.csv', emit: expression

    script:
    """
    python3 ${mean_tpm_script} \
      --samplesheet ${samplesheet} \
      --output Expression.csv \
      ${normalized_results.join(' ')}
    """

    stub:
    """
    touch Expression.csv
    """
}

process BIGWIG {
    tag "${condition}"
    label 'medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    tuple val(condition), val(sample_ids), path(bam_files)

    output:
    path "${condition}.bw", emit: bigwig

    script:
    """
    samtools merge -@ ${task.cpus} -f ${condition}.bam ${bam_files.join(' ')}
    samtools index -@ ${task.cpus} ${condition}.bam
    bamCoverage \
      -p ${task.cpus} \
      -b ${condition}.bam \
      -o ${condition}.bw \
      --normalizeUsing CPM
    """

    stub:
    """
    touch ${condition}.bw
    """
}

workflow {
    if (!params.ref_file) {
        error 'Missing required parameter: --ref_file'
    }
    if (!params.anno_file) {
        error 'Missing required parameter: --anno_file'
    }

    def samplesheet_file = file(params.datasets, checkIfExists: true)
    def reference_file = file(params.ref_file, checkIfExists: true)
    def annotation_file = file(params.anno_file, checkIfExists: true)

    samplesheet = Channel.value(samplesheet_file)
    reference = Channel.value(reference_file)
    annotation = Channel.value(annotation_file)
    gff_to_bed_script = Channel.value(file("${projectDir}/Scripts/gff_to_bed.py", checkIfExists: true))
    strandedness_script = Channel.value(file("${projectDir}/Scripts/infer_strandedness.py", checkIfExists: true))
    tpm_script = Channel.value(file("${projectDir}/Scripts/tpm_calc.py", checkIfExists: true))
    mean_tpm_script = Channel.value(file("${projectDir}/Scripts/mean_tpm_calc.py", checkIfExists: true))

    reads = Channel
        .fromPath(samplesheet_file, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id || !row.condition || !row.fastq_1) {
                error 'Each sample-sheet row requires sample_id, condition, and fastq_1'
            }
            if (!(row.sample_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error "Invalid sample_id '${row.sample_id}': use letters, numbers, dot, underscore, or hyphen"
            }
            if (!(row.condition ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error "Invalid condition '${row.condition}': use letters, numbers, dot, underscore, or hyphen"
            }

            def read1 = file(row.fastq_1, checkIfExists: true)
            def read2 = row.fastq_2?.trim() ? file(row.fastq_2, checkIfExists: true) : null
            tuple(row.sample_id, row.condition, read1, read2)
        }

    pe_reads = reads.filter { sample_id, condition, read1, read2 -> read2 != null }
    se_reads = reads.filter { sample_id, condition, read1, read2 -> read2 == null }

    trimmed_pe = TRIM_PE(pe_reads, params.q_score).reads
    trimmed_se = TRIM_SE(
        se_reads.map { sample_id, condition, read1, read2 -> tuple(sample_id, condition, read1) },
        params.q_score
    ).reads

    if (params.skip_rrna) {
        clean_pe = trimmed_pe
        clean_se = trimmed_se
    } else {
        def database_params = [params.rrna_16s, params.rrna_23s, params.rrna_5s]
        if (database_params.any { !it }) {
            error 'Set --rrna_16s, --rrna_23s, and --rrna_5s when --skip_rrna false'
        }
        rrna_databases = Channel.value(
            database_params.collect { file(it, checkIfExists: true) }
        )
        clean_pe = FILTER_RRNA_PE(trimmed_pe, rrna_databases).reads
        clean_se = FILTER_RRNA_SE(trimmed_se, rrna_databases).reads
    }

    // BUILD_INDEX is driven by a value channel, so its complete .bt2/.bt2l
    // bundle is also a reusable value channel broadcast to every mapper.
    index_bundle = BUILD_INDEX(reference).index_files

    mapped_pe = MAP_PE(clean_pe, index_bundle).bam
    mapped_se = MAP_SE(clean_se, index_bundle).bam
    mapped_bam = mapped_pe.mix(mapped_se)

    annotation_bed = GFF_TO_BED(annotation, gff_to_bed_script).bed
    stranded = ESTIMATE_STRANDEDNESS(mapped_bam, annotation_bed).result

    pe_stranded = stranded.filter { sample_id, layout, condition, strand_info, bam -> layout == 'PE' }
    se_stranded = stranded.filter { sample_id, layout, condition, strand_info, bam -> layout == 'SE' }

    pe_counts = QUANTIFY_PE(pe_stranded, annotation, strandedness_script).counts
    se_counts = QUANTIFY_SE(se_stranded, annotation, strandedness_script).counts

    normalized_pe = NORMALIZE_PE(pe_counts, tpm_script).result
    normalized_se = NORMALIZE_SE(se_counts, tpm_script).result
    normalized = normalized_pe.mix(normalized_se)

    normalized_files = normalized
        .map { sample_id, layout, condition, normalized_file, metrics -> normalized_file }
        .collect()

    MEAN_TPM(normalized_files, samplesheet, mean_tpm_script)

    grouped_bams = mapped_bam
        .map { sample_id, layout, condition, bam -> tuple(condition, sample_id, bam) }
        .groupTuple()

    BIGWIG(grouped_bams)
}
