#! /usr/bin/env nextflow

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------'
    log.info 'Variant call v0.1: A MULTI-SAMPLE VARIANT CALLER using GATK'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run variant_call.nf --bed bedfile.bed --bam_folder BAM/ --fasta_ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info 'Options:'
    log.info '    --nsplit         INTEGER                  Split the region for calling in nsplit pieces and run in parallel.'
    log.info '    --min_dp         INTEGER                  Minimum coverage in at least one sample to consider a site.'
    log.info '    --min_ao         INTEGER                  Minimum number of non-ref reads in at least one sample to consider a site.'
    log.info '    --min_qval       VALUE                    Qvalue in Phred scale to consider a variant.'
    log.info '    --sb_type        SOR or RVSB              Strand bias measure.'
    log.info '    --sb_snv         VALUE                    Strand bias threshold for SNVs.'
    log.info '    --sb_indel       VALUE                    Strand bias threshold for indels.'
    log.info '    --map_qual       VALUE                    Samtools minimum mapping quality.'
    log.info '    --base_qual      VALUE                    Samtools minimum base quality.'
    log.info '    --max_DP         INTEGER                  Samtools maximum coverage before downsampling.'
    log.info '    --use_file_name                           Sample names are taken from file names, otherwise extracted from the bam file SM tag.'
    log.info '    --all_SNVs                                Output all SNVs, even when no variant found.'
    log.info '    --no_plots                                Do not output PDF regression plots.'
    log.info '    --no_labels                               Do not add labels to outliers in regression plots.'
    log.info '    --no_indels                               Do not call indels.'
    log.info '    --no_contours                             Do not add contours to plots and do not plot min(AF)~DP.'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, by default input bam folder.'
    log.info '    --bed            BED FILE                 A BED file for calling.'
    log.info '    --region         CHR:START-END            A region for calling.'
    log.info ''
    exit 1
}

assert (params.fasta_ref != true) && (params.fasta_ref != null) : "please specify --fasta_ref option (--fasta_ref reference.fasta(.gz))"
assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"

