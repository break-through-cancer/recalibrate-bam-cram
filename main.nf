#!/usr/bin/env nextflow
/*
 * BQSR pipeline -- GATK BaseRecalibrator + ApplyBQSR, following the same
 * two steps as nf-core/sarek's fastq_preprocess_gatk subworkflow
 * (BAM_BASERECALIBRATOR -> BAM_APPLYBQSR) and gatk4/applybqsr module.
 *
 * Deliberate deviation from a known Sarek bug (nf-core/sarek #1772, and
 * confirmed by GATK's own team on the Broad forum: "You should not use
 * intervals with ApplyBQSR because you want to recalibrate all of the
 * reads. You can introduce artifacts if you run with intervals in
 * ApplyBQSR."): no intervals are used anywhere in this pipeline. Both
 * BaseRecalibrator and ApplyBQSR run unrestricted, over the full input --
 * including unmapped reads, which no interval list can ever cover.
 */
nextflow.enable.dsl = 2

if (!params.containsKey('outdir'))    params.outdir = 'results'
if (!params.containsKey('ref_fasta')) params.ref_fasta = null
if (!params.containsKey('ref_fai'))   params.ref_fai = null
if (!params.containsKey('ref_dict'))  params.ref_dict = null
if (!params.containsKey('dbsnp'))     params.dbsnp = null
if (!params.containsKey('dbsnp_idx')) params.dbsnp_idx = null

process GATK4_BASERECALIBRATOR {
    tag "${sample_id}"
    label 'process_medium'
    container "broadinstitute/gatk:4.5.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(cram), path(crai)
    path ref_fasta
    path ref_fai
    path ref_dict
    path known_sites_files
    path known_sites_tbi_files

    output:
    tuple val(sample_id), path("${sample_id}.recal.table"), emit: table

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{cram}.crai ]; then
        ln -s !{crai} !{cram}.crai
    fi

    java_mem_mb=!{task.memory.toMega() - 1024}

    known_sites_args=""
    for f in !{known_sites_files}; do
        known_sites_args="$known_sites_args --known-sites $f"
    done

    gatk --java-options "-Xmx${java_mem_mb}m" BaseRecalibrator \
        -I !{cram} \
        -R !{ref_fasta} \
        $known_sites_args \
        -O !{sample_id}.recal.table \
        --tmp-dir .

    test -s !{sample_id}.recal.table
    '''
}

process GATK4_APPLYBQSR {
    tag "${sample_id}"
    label 'process_medium'
    container "broadinstitute/gatk:4.5.0.0"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(cram), path(crai), path(recal_table)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.recal.cram"), path("${sample_id}.recal.cram.crai"), emit: recal_cram

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{cram}.crai ]; then
        ln -s !{crai} !{cram}.crai
    fi

    java_mem_mb=!{task.memory.toMega() - 1024}

    # NOTE: deliberately no -L / --intervals here -- ApplyBQSR must always
    # run over the full genome + unmapped reads. Restricting it silently
    # drops any read outside the interval from the output (nf-core/sarek
    # #1772; confirmed anti-pattern per GATK's own team on the Broad forum).
    gatk --java-options "-Xmx${java_mem_mb}m" ApplyBQSR \
        -I !{cram} \
        -R !{ref_fasta} \
        --bqsr-recal-file !{recal_table} \
        -O !{sample_id}.recal.cram \
        --tmp-dir .

    samtools index !{sample_id}.recal.cram

    test -s !{sample_id}.recal.cram
    test -s !{sample_id}.recal.cram.crai
    '''
}

workflow {

    if (!params.bqsr_runs) {
        error "params.bqsr_runs is empty -- did the Cirro preprocess.py hook run? (see preprocess.py)"
    }
    if (!params.ref_fasta || !params.ref_fai || !params.ref_dict) {
        error "Missing required reference params: ref_fasta, ref_fai, ref_dict"
    }
    if (!params.dbsnp || !params.dbsnp_idx) {
        error "Missing required param: dbsnp / dbsnp_idx"
    }

    ref_fasta = file(params.ref_fasta, checkIfExists: true)
    ref_fai   = file(params.ref_fai,   checkIfExists: true)
    ref_dict  = file(params.ref_dict,  checkIfExists: true)

    known_sites_files    = [file(params.dbsnp, checkIfExists: true)]
    known_sites_tbi_files = [file(params.dbsnp_idx, checkIfExists: true)]

    runs_ch =
        Channel
            .fromList(params.bqsr_runs)
            .map { run ->
                tuple(
                    run.sample_id,
                    file(run.cram, checkIfExists: true),
                    file(run.crai, checkIfExists: true)
                )
            }

    GATK4_BASERECALIBRATOR(
        runs_ch,
        ref_fasta, ref_fai, ref_dict,
        known_sites_files, known_sites_tbi_files
    )

    applybqsr_inputs =
        runs_ch.join(GATK4_BASERECALIBRATOR.out.table)
            .map { sample_id, cram, crai, table -> tuple(sample_id, cram, crai, table) }

    GATK4_APPLYBQSR(
        applybqsr_inputs,
        ref_fasta, ref_fai, ref_dict
    )
}