#!/usr/bin/env nextflow
/*
 * BQSR pipeline -- GATK BaseRecalibrator + ApplyBQSR, scattered/gathered
 * across sequence groupings exactly as in Broad's reference implementation
 * (gatk-workflows/gatk4-data-processing: processing-for-variant-discovery-gatk4.wdl,
 * tasks CreateSequenceGroupingTSV / BaseRecalibrator / GatherBqsrReports /
 * ApplyBQSR / GatherBamFiles) -- the same pipeline that runs in Terra.
 *
 * Both BaseRecalibrator and ApplyBQSR ARE given -L intervals here, scattered
 * per contig grouping -- this is safe specifically because the groupings
 * fully tile 100% of the reference dictionary (every contig in ref_dict,
 * binned by length) plus one explicit trailing "unmapped" group used only
 * for ApplyBQSR's scatter, so every read -- mapped anywhere, or nowhere --
 * lands in exactly one shard and nothing is lost on gather.
 *
 * The groupings are NOT a separate bundled reference file -- they're
 * derived at runtime from the ref_dict you already pass in (see
 * CREATE_SEQUENCE_GROUPING below), the same way the WDL's
 * CreateSequenceGroupingTSV task does it. This is NOT the same thing as
 * restricting either tool to a WGS calling-region interval list (e.g.
 * wgs_calling_regions.hg38.interval_list), which deliberately EXCLUDES
 * N-gaps/centromeres/some decoys/alts for calling efficiency. Handing a
 * list like that to ApplyBQSR -- scattered or not -- silently drops every
 * read in an excluded region (nf-core/sarek #1772; confirmed anti-pattern
 * per GATK's own team on the Broad forum). Do not substitute the two kinds
 * of interval list for one another.
 *
 * One intentional deviation from the reference WDL, since our upstream
 * CRAMs come from Sarek's markduplicates stage rather than Broad's own
 * uBAM -> MergeBamAlignment ingestion:
 *   - no --use-original-qualities: that flag recalibrates from the OQ tag,
 *     which only exists if something upstream explicitly wrote it (Broad's
 *     own pipeline does, via Picard MergeBamAlignment on a uBAM). Sarek's
 *     CRAMs don't carry an OQ tag on a first-pass BQSR, so turning this on
 *     would recalibrate against a tag that isn't there. Confirm your CRAM
 *     provenance before adding it.
 *
 * --static-quantized-quals 10/20/30 IS included: without it, ApplyBQSR
 * emits full-precision recalibrated quality scores -- much higher entropy
 * per base than Sarek's original (often pre-binned) qualities -- which
 * bloats recalibrated CRAMs to several times the size of the markduplicates
 * input. This is the same binning Broad's own WARP/Terra pipeline applies,
 * per the published Functional Equivalence spec (Regier et al. 2018);
 * GATK's own docs say it has no noticeable effect on germline variant
 * discovery at these levels.
 *
 * NEW in this version: params.output_format ('cram', the default, or
 * 'bam') controls the container format of the final recalibrated output,
 * for every sample in the run. CRAM is smaller (reference-based sequence
 * diffing + a dedicated quality-score entropy coder) but costs more CPU to
 * read/write; BAM is larger but cheaper to encode/decode. GATK's own docs
 * put the crossover around ~180x depth -- below that (i.e. essentially all
 * standard WGS/WES runs), CRAM trades write speed for size, not the other
 * way around. Applies uniformly: there is no per-sample override, only a
 * single pipeline-wide choice.
 *
 * This is main_3.nf -- builds on main_2.nf's scatter/gather architecture,
 * kept side-by-side with main.nf (single-shot) and main_2.nf (scatter/
 * gather, CRAM-only) for comparison.
 */
nextflow.enable.dsl = 2

if (!params.containsKey('outdir'))        params.outdir = 'results'
if (!params.containsKey('ref_fasta'))     params.ref_fasta = null
if (!params.containsKey('ref_fai'))       params.ref_fai = null
if (!params.containsKey('ref_dict'))      params.ref_dict = null
if (!params.containsKey('dbsnp'))         params.dbsnp = null
if (!params.containsKey('dbsnp_idx'))     params.dbsnp_idx = null
if (!params.containsKey('output_format')) params.output_format = 'cram'

process CREATE_SEQUENCE_GROUPING {
    label 'process_low'
    container "broadinstitute/gatk:4.5.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    path ref_dict

    output:
    path "sequence_grouping.tsv",               emit: mapped_only
    path "sequence_grouping_with_unmapped.tsv", emit: with_unmapped

    shell:
    '''
    set -euxo pipefail

    # NOTE: unlike WDL's command <<< >>>, Nextflow's shell block does NOT
    # strip common leading whitespace -- everything below must be flush
    # left (true Python top-level indentation), or the interpreter sees a
    # bogus leading indent and fails immediately.
    python3 <<'CODE'
with open("!{ref_dict}") as fh:
    sequence_tuples = []
    for line in fh:
        if line.startswith("@SQ"):
            fields = line.split("\\t")
            name = fields[1].split("SN:")[1]
            length = int(fields[2].split("LN:")[1])
            sequence_tuples.append((name, length))

longest = max(length for _, length in sequence_tuples)
# Sacrificial ":1+" suffix on every contig name -- workaround for an old
# GATK bug that strips text after a colon in some contig names (hg38
# ALT contigs). "chr1:1+" means "chr1, from position 1 to the end", so
# this is semantically a no-op interval, just spelled so every contig
# (colon-bearing or not) is parsed the same way. Kept for fidelity with
# the reference pipeline; may be unnecessary on modern GATK versions.
protection_tag = ":1+"

groups = []
current = [sequence_tuples[0][0] + protection_tag]
current_size = sequence_tuples[0][1]
for name, length in sequence_tuples[1:]:
    if current_size + length <= longest:
        current.append(name + protection_tag)
        current_size += length
    else:
        groups.append(current)
        current = [name + protection_tag]
        current_size = length
groups.append(current)

with open("sequence_grouping.tsv", "w") as fh:
    fh.write("\\n".join("\\t".join(g) for g in groups))

groups_with_unmapped = groups + [["unmapped"]]
with open("sequence_grouping_with_unmapped.tsv", "w") as fh:
    fh.write("\\n".join("\\t".join(g) for g in groups_with_unmapped))
CODE

    test -s sequence_grouping.tsv
    test -s sequence_grouping_with_unmapped.tsv
    '''
}

process GATK4_BASERECALIBRATOR {
    tag "${sample_id}/${group_tag}"
    label 'process_medium'
    container "broadinstitute/gatk:4.5.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(cram), path(crai), val(group_tag), val(interval_args)
    path ref_fasta
    path ref_fai
    path ref_dict
    path known_sites_files
    path known_sites_tbi_files

    output:
    tuple val(sample_id), path("${sample_id}.${group_tag}.recal.table"), emit: table

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
        !{interval_args} \
        -O !{sample_id}.!{group_tag}.recal.table \
        --tmp-dir .

    test -s !{sample_id}.!{group_tag}.recal.table
    '''
}

process GATHER_BQSR_REPORTS {
    tag "${sample_id}"
    label 'process_low'
    container "broadinstitute/gatk:4.5.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(tables)

    output:
    tuple val(sample_id), path("${sample_id}.recal.table"), emit: table

    shell:
    '''
    set -euxo pipefail

    java_mem_mb=!{task.memory.toMega() - 512}

    report_args=""
    for f in !{tables}; do
        report_args="$report_args -I $f"
    done

    gatk --java-options "-Xmx${java_mem_mb}m" GatherBQSRReports \
        $report_args \
        -O !{sample_id}.recal.table

    test -s !{sample_id}.recal.table
    '''
}

process GATK4_APPLYBQSR {
    tag "${sample_id}/${group_tag}"
    label 'process_medium'
    container "broadinstitute/gatk:4.5.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(cram), path(crai), path(recal_table), val(group_idx), val(group_tag), val(interval_args), val(output_format)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), val(group_idx), path("${sample_id}.${group_tag}.recal.${output_format}"), emit: recal_shard

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{cram}.crai ]; then
        ln -s !{crai} !{cram}.crai
    fi

    java_mem_mb=!{task.memory.toMega() - 1024}

    # Output container format (bam/cram) is inferred by GATK purely from the
    # -O filename's extension -- there is no separate --output-format flag.
    gatk --java-options "-Xmx${java_mem_mb}m" ApplyBQSR \
        -I !{cram} \
        -R !{ref_fasta} \
        --bqsr-recal-file !{recal_table} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        !{interval_args} \
        -O !{sample_id}.!{group_tag}.recal.!{output_format} \
        --tmp-dir .

    test -s !{sample_id}.!{group_tag}.recal.!{output_format}
    '''
}

process GATHER_ALIGNMENT_FILES {
    tag "${sample_id}"
    label 'process_medium'
    container "broadinstitute/gatk:4.5.0.0"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(shards), val(output_format)

    output:
    tuple val(sample_id), path("${sample_id}.recal.${output_format}"), path("${sample_id}.recal.${output_format}.*"), emit: recal_alignment

    shell:
    '''
    set -euxo pipefail

    # `shards` arrives pre-sorted by group index (mapped contigs in
    # reference-dict order, unmapped last) -- samtools cat concatenates in
    # the given order without re-sorting, matching Picard GatherBamFiles'
    # block-concatenation semantics in the reference pipeline. Works
    # identically for BAM or CRAM shards; samtools infers the output
    # container format from the -o filename's extension, same as GATK does.
    samtools cat -o !{sample_id}.recal.!{output_format} !{shards}

    # samtools index picks .bai (BAM) or .crai (CRAM) based on the input
    # file's own format -- no extra flag needed either way.
    samtools index !{sample_id}.recal.!{output_format}

    test -s !{sample_id}.recal.!{output_format}
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
    if (!(params.output_format in ['cram', 'bam'])) {
        error "params.output_format must be 'cram' or 'bam' (got: ${params.output_format})"
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

    CREATE_SEQUENCE_GROUPING(ref_dict)

    // One-time, sample-independent: (group_idx, group_tag, "-L a -L b ...")
    mapped_groups =
        CREATE_SEQUENCE_GROUPING.out.mapped_only
            .splitCsv(sep: '\t')
            .toList()
            .flatMap { rows ->
                rows.withIndex().collect { row, idx ->
                    tuple(idx, "group${idx}", row.collect { "-L ${it}" }.join(' '))
                }
            }

    groups_with_unmapped =
        CREATE_SEQUENCE_GROUPING.out.with_unmapped
            .splitCsv(sep: '\t')
            .toList()
            .flatMap { rows ->
                rows.withIndex().collect { row, idx ->
                    tuple(idx, "group${idx}", row.collect { "-L ${it}" }.join(' '))
                }
            }

    // ---- BaseRecalibrator: scattered per sample x mapped-contig group ----
    br_inputs =
        runs_ch.combine(mapped_groups)
            .map { sample_id, cram, crai, group_idx, group_tag, interval_args ->
                tuple(sample_id, cram, crai, group_tag, interval_args)
            }

    GATK4_BASERECALIBRATOR(
        br_inputs,
        ref_fasta, ref_fai, ref_dict,
        known_sites_files, known_sites_tbi_files
    )

    tables_by_sample =
        GATK4_BASERECALIBRATOR.out.table
            .groupTuple()

    GATHER_BQSR_REPORTS(tables_by_sample)

    // ---- ApplyBQSR: scattered per sample x (mapped-contig group + unmapped) ----
    ab_base = runs_ch.join(GATHER_BQSR_REPORTS.out.table)

    ab_inputs =
        ab_base.combine(groups_with_unmapped)
            .map { sample_id, cram, crai, table, group_idx, group_tag, interval_args ->
                tuple(sample_id, cram, crai, table, group_idx, group_tag, interval_args, params.output_format)
            }

    GATK4_APPLYBQSR(
        ab_inputs,
        ref_fasta, ref_fai, ref_dict
    )

    shards_by_sample =
        GATK4_APPLYBQSR.out.recal_shard
            .groupTuple()
            .map { sample_id, group_idxs, shards ->
                def ordered = [group_idxs, shards].transpose().sort { it[0] }.collect { it[1] }
                tuple(sample_id, ordered, params.output_format)
            }

    GATHER_ALIGNMENT_FILES(shards_by_sample)
}
