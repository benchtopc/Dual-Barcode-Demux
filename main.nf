nextflow.enable.dsl = 2

def showHelp() {
    log.info """
    Dual-Barcode-Demux

    Usage:
      nextflow run . --fastq reads.fastq.gz --plate_map plate_map.csv -profile docker

    Required:
      --fastq
      --plate_map

    Optional alignment / KO / IGV:
      --run_alignment true
      --reference reference.fa
      --ko_targets templates/ko_targets.template.tsv
      --variant_caller bcftools|clair3|clair3_gpu

    Clair3:
      --variant_caller clair3_gpu
      --clair3_model_path /opt/models/r1041_e82_400bps_sup_v500
      --clair3_model_dir /path/to/local/model_dir
      --clair3_platform ont
      --enable_gpu true

    Output:
      --out_dir output
    """
}

params.fastq = null
params.plate_map = null
params.row_barcodes = "${projectDir}/data/scm_row_barcodes.csv"
params.column_barcodes = "${projectDir}/data/scm_column_barcodes.csv"
params.front_window = 200
params.rear_window = 200
params.max_linker_mismatches = 2
params.max_total_score = 4
params.min_margin = 2
params.keep_partial = true

params.run_alignment = false
params.reference = null
params.ko_targets = "${projectDir}/templates/ko_targets.template.tsv"
params.aligner_preset = "map-ont"
params.alignment_threads = 4
params.min_mapq = 10
params.min_baseq = 7
params.min_depth_ko = 10
params.het_min_af = 0.25
params.hom_alt_min_af = 0.80
params.keep_unfiltered_bam = false

params.variant_caller = "bcftools"
params.enable_gpu = false
params.gpu_devices = "all"
params.igv_locus = "all"

params.clair3_model_path = "/opt/models/r1041_e82_400bps_sup_v500"
params.clair3_model_dir = ""
params.clair3_platform = "ont"
params.clair3_bed = "${projectDir}/templates/empty.bed"
params.clair3_ctg_name = ""
params.clair3_enable_dwell_time = false
params.clair3_include_all_ctgs = true
params.clair3_haploid_mode = "none"
params.clair3_no_phasing_for_fa = false
params.clair3_pileup_only = false
params.clair3_print_ref_calls = false
params.clair3_gvcf = false
params.clair3_qual = 2
params.clair3_chunk_size = 0
params.clair3_remove_intermediate_dir = true
params.clair3_device = "cuda:0"
params.clair3_extra_args = ""
params.clair3_cpu_container = "hkubal/clair3:v2.0.1"
params.clair3_gpu_container = "hkubal/clair3:v2.0.1_gpu"

params.out_dir = 'output'
params.help = false

if (params.help) {
    showHelp()
    System.exit(0)
}

if (!params.fastq || !params.plate_map) {
    showHelp()
    error "Both --fastq and --plate_map are required."
}

if (params.run_alignment && !params.reference) {
    error "Alignment/KO analysis requires --reference when --run_alignment true."
}

if (params.variant_caller !in ["bcftools", "clair3", "clair3_gpu"]) {
    error "--variant_caller must be one of: bcftools, clair3, clair3_gpu"
}

if (params.clair3_platform !in ["ont", "hifi", "ilmn"]) {
    error "--clair3_platform must be one of: ont, hifi, ilmn"
}

if (params.clair3_haploid_mode !in ["none", "haploid_precise", "haploid_sensitive"]) {
    error "--clair3_haploid_mode must be one of: none, haploid_precise, haploid_sensitive"
}

workflow {
    reads_ch = Channel.fromPath(params.fastq, checkIfExists: true)
    plate_map_ch = Channel.fromPath(params.plate_map, checkIfExists: true)
    row_barcodes_ch = Channel.fromPath(params.row_barcodes, checkIfExists: true)
    column_barcodes_ch = Channel.fromPath(params.column_barcodes, checkIfExists: true)

    DEMUX_SCM(reads_ch, plate_map_ch, row_barcodes_ch, column_barcodes_ch)

    if (params.run_alignment) {
        reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
        ko_targets_ch = Channel.fromPath(params.ko_targets, checkIfExists: true)

        ALIGN_SCM(DEMUX_SCM.out.demux_fastq, reference_ch)

        if (params.variant_caller == "bcftools") {
            BCFTOOLS_CALL_SCM(ALIGN_SCM.out.alignments, reference_ch)
            variants_ch = BCFTOOLS_CALL_SCM.out.variants
        } else {
            clair3_model_dir_ch = params.clair3_model_dir
                ? Channel.fromPath(params.clair3_model_dir, checkIfExists: true)
                : Channel.fromPath("${projectDir}/templates/empty_clair3_model", checkIfExists: true)
            clair3_bed_ch = Channel.fromPath(params.clair3_bed, checkIfExists: true)

            CLAIR3_CALL_SCM(ALIGN_SCM.out.alignments, reference_ch, clair3_model_dir_ch, clair3_bed_ch)
            variants_ch = CLAIR3_CALL_SCM.out.variants
        }

        ANALYZE_KO_IGV_SCM(ALIGN_SCM.out.alignments, variants_ch, reference_ch, ko_targets_ch)
    }
}

process DEMUX_SCM {
    tag "${fastq.simpleName}"

    publishDir "${params.out_dir}", mode: 'copy'

    container "python:3.11-slim"

    input:
    path fastq
    path plate_map
    path row_barcodes
    path column_barcodes

    output:
    path "demultiplexed_fastq", emit: demux_fastq
    path "tables", emit: tables
    path "reports", emit: reports
    path "metadata", emit: metadata

    script:
    def keepPartial = params.keep_partial ? '--keep-partial' : ''
    """
    python3 ${projectDir}/bin/demux_scm.py \
      --fastq ${fastq} \
      --plate-map ${plate_map} \
      --row-barcodes ${row_barcodes} \
      --column-barcodes ${column_barcodes} \
      --front-window ${params.front_window} \
      --rear-window ${params.rear_window} \
      --max-linker-mismatches ${params.max_linker_mismatches} \
      --max-total-score ${params.max_total_score} \
      --min-margin ${params.min_margin} \
      ${keepPartial} \
      --outdir .
    """
}

process ALIGN_SCM {
    tag "${reference.simpleName}"

    publishDir "${params.out_dir}", mode: 'copy'

    container "ubuntu:24.04"

    cpus { params.alignment_threads as int }

    input:
    path demux_fastq
    path reference

    output:
    path "alignments", emit: alignments
    path "gaps", emit: gaps
    path "analysis_reports", emit: alignment_reports

    script:
    def keepUnfiltered = params.keep_unfiltered_bam ? '--keep-unfiltered-bam' : ''
    """
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends \
      bash \
      ca-certificates \
      minimap2 \
      samtools \
      python3 \
      coreutils \
      gzip
    rm -rf /var/lib/apt/lists/*

    python3 ${projectDir}/bin/align_samples_scm.py \
      --demux-dir ${demux_fastq} \
      --reference ${reference} \
      --aligner-preset ${params.aligner_preset} \
      --threads ${task.cpus} \
      --min-mapq ${params.min_mapq} \
      --min-depth-ko ${params.min_depth_ko} \
      ${keepUnfiltered} \
      --outdir .
    """
}

process BCFTOOLS_CALL_SCM {
    tag "bcftools"

    publishDir "${params.out_dir}", mode: 'copy'

    container "ubuntu:24.04"

    cpus { params.alignment_threads as int }

    input:
    path alignments
    path reference

    output:
    path "variants", emit: variants

    script:
    """
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends \
      bash \
      ca-certificates \
      samtools \
      bcftools \
      python3 \
      coreutils \
      gzip
    rm -rf /var/lib/apt/lists/*

    python3 ${projectDir}/bin/call_variants_bcftools_scm.py \
      --alignments-dir ${alignments} \
      --reference ${reference} \
      --threads ${task.cpus} \
      --min-mapq ${params.min_mapq} \
      --min-baseq ${params.min_baseq} \
      --outdir .
    """
}

process CLAIR3_CALL_SCM {
    tag "${params.variant_caller}"

    publishDir "${params.out_dir}", mode: 'copy'

    container "${params.variant_caller == 'clair3_gpu' ? params.clair3_gpu_container : params.clair3_cpu_container}"

    cpus { params.alignment_threads as int }

    input:
    path alignments
    path reference
    path local_model_dir
    path clair3_bed

    output:
    path "variants", emit: variants
    path "clair3_runs", emit: clair3_runs

    script:
    def useGpu = (params.variant_caller == 'clair3_gpu' || params.enable_gpu) ? '--use-gpu' : ''
    def dwell = params.clair3_enable_dwell_time ? '--enable-dwell-time' : ''
    def includeAll = params.clair3_include_all_ctgs ? '--include-all-ctgs' : ''
    def noPhasing = params.clair3_no_phasing_for_fa ? '--no-phasing-for-fa' : ''
    def pileupOnly = params.clair3_pileup_only ? '--pileup-only' : ''
    def printRef = params.clair3_print_ref_calls ? '--print-ref-calls' : ''
    def gvcf = params.clair3_gvcf ? '--gvcf' : ''
    def removeIntermediate = params.clair3_remove_intermediate_dir ? '--remove-intermediate-dir' : ''
    def modelPath = params.clair3_model_dir ? local_model_dir : params.clair3_model_path
    """
    python3 ${projectDir}/bin/run_clair3_batch_scm.py \
      --alignments-dir ${alignments} \
      --reference ${reference} \
      --model-path ${modelPath} \
      --platform ${params.clair3_platform} \
      --threads ${task.cpus} \
      ${useGpu} \
      --device "${params.clair3_device}" \
      --bed-fn ${clair3_bed} \
      --ctg-name "${params.clair3_ctg_name}" \
      ${dwell} \
      ${includeAll} \
      --haploid-mode ${params.clair3_haploid_mode} \
      ${noPhasing} \
      ${pileupOnly} \
      ${printRef} \
      ${gvcf} \
      --qual ${params.clair3_qual} \
      --chunk-size ${params.clair3_chunk_size} \
      ${removeIntermediate} \
      --extra-args "${params.clair3_extra_args}" \
      --outdir .
    """
}

process ANALYZE_KO_IGV_SCM {
    tag "ko_igv"

    publishDir "${params.out_dir}", mode: 'copy'

    container "python:3.11-slim"

    input:
    path alignments
    path variants
    path reference
    path ko_targets

    output:
    path "ko_calls", emit: ko_calls
    path "igv", emit: igv
    path "analysis_reports", emit: final_reports

    script:
    """
    python3 ${projectDir}/bin/analyze_ko_igv_scm.py \
      --alignments-dir ${alignments} \
      --variants-dir ${variants} \
      --reference ${reference} \
      --ko-targets ${ko_targets} \
      --min-depth-ko ${params.min_depth_ko} \
      --het-min-af ${params.het_min_af} \
      --hom-alt-min-af ${params.hom_alt_min_af} \
      --igv-locus "${params.igv_locus}" \
      --outdir .
    """
}
