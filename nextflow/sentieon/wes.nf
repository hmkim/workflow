#!/usr/bin/env nextflow

/*
---------------------------------------------------------------------------------------
The pipeline can determine whether the input data is single or paired end. This relies on
specifying the input files correctly. For paired en data us the example above, i.e.
'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - Sentieon
*/

// Pipeline version
version = 0.1

// Configurable variables
params.genome = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
params.bwa_index = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa'
params.known_sites = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/Mills_and_1000G_gold_standard.indels.vcf.gz'
params.dbsnp = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz'
params.gtf = ''
params.rlocation = ''

params.name = "Sentieon best practice for WES analysis"

// Input reads
params.reads = '/BiO/BioPeople/brandon/test_nextflow/wes/example/TN*{R1,R2}_[0-9][0-9][0-9].fastq.gz'
// Output directory
params.outdir = '/BiO/BioPeople/brandon/test_nextflow/wes/outdir'

log.info "===================================="
log.info " WES (Sentieon) Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="

// Validate inputs
index = file(params.bwa_index)
if( !index.exists() ) exit 1, "Missing BWA index: ${index}"

/*
 * Create a channel for read files
 */
 
Channel
     .fromPath( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
     .map { path ->
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
     }
     .groupTuple(sort: true)
     .set { read_files }
 
read_files.into { read_files_fastqc; read_files_mapping;}

/*
 * STEP 1 - FastQC
 */
process fastqc {
     tag "$prefix"
     
     memory { 2.GB * task.attempt }
     time { 4.h * task.attempt }
     
     errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
     maxRetries 3
     maxErrors '-1'
     
     publishDir "${params.outdir}/$prefix/fastqc", mode: 'copy'
     
     input:
     set val(prefix), file(reads:'*') from read_files_fastqc
     
     output:
     file '*_fastqc.html' into fastqc_results_html
     file '*_fastqc.zip' into fastqc_results_zip
     
     """
     fastqc $reads
     """
}

/*
 * Step 2. Maps each read-pair by using mapper and sorting
 */
process mapping {
    tag "$prefix"
    cpus 4
 
    input: 
    set val(prefix), file (reads:'*') from read_files_mapping

    output:
    set prefix, file { "sorted.bam" }, file { "sorted.bam.bai" } into bamFiles

    """
    bwa mem -M -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina' -t ${task.cpus} ${params.bwa_index} $reads | sentieon util sort -o sorted.bam -t ${task.cpus} --sam2bam -i -
    """
}

bamFiles.into{ bamFiles_metrics; bamFiles_rmdup; }

/*
 * Step 3. metrics for mapping
 */
process metrics{
    tag "$prefix"
    cpus 4
    
    publishDir "${params.outdir}/$prefix/metrics", mode: 'copy'

    input:
    set val(prefix), file('*'), file('*') from bamFiles_metrics

    output:
    set file('*.txt') into metrics_txts
    set file('*.pdf') into metrics_pdf
    """
    sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
    sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt
    """
}

/*
 * Step 3. Remove PCR duplicates after mapping
 */
process rmdup{
    tag "$prefix"
    cpus 4
    
    publishDir "${params.outdir}/$prefix/rmdup"

    input:
    bamFiles
    set val(prefix), file('*'), file('*') from bamFiles_rmdup

    output:
    set val(prefix), file('deduped.bam*') into dedupBamFiles
    set file('dedup_metrics.txt') into dedupBam_metrics
    
    """
    sentieon driver -t ${task.cpus} -i sorted.bam --algo LocusCollector --fun score_info score.txt
    sentieon driver -t ${task.cpus} -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam
    """
}

dedupBamFiles.into{ dedupBamFiles_qualimap; dedupBamFiles_realign; }

process qualimap {
	tag "$prefix"
	cpus 8
	memory '8 GB'

	publishDir "${params.outdir}/${prefix}/qualimap"

	input:
	set val(prefix), file(read:'*') from dedupBamFiles_qualimap

	output:
	val(prefix) into qualimap_prefix
	file ( "$prefix/genome_results.txt" ) into qualimap_genome_results_txt
	file ( "$prefix/raw_data_qualimapReport/*.txt" ) into qualimap_report_txts
	file ( "$prefix/images_qualimapReport/*.png" ) into qualimap_report_pngs
	file ( "$prefix/css/*" ) into qualimap_report_css
	file ( "$prefix/qualimapReport.html" ) into qualimap_report_html
	
	"""
	qualimap bamqc --java-mem-size=8G -bam deduped.bam -gd HUMAN -nt ${task.cpus} -outdir $prefix
	"""
}

/*
 * Step 4. Indel realignment after dedup
 */
process indelrealign{
	tag "$prefix"
	cpus 4

	input:
	set val(prefix), file('*') from dedupBamFiles_realign

	output:
	set val(prefix), file('realigned.bam*') into realignedBamFiles

	"""
	sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i deduped.bam --algo Realigner -k ${params.known_sites} realigned.bam
	"""
}

realignedBamFiles.into { realignedBamFiles_forBaseRecal; realignedBamFiles_UG; realignedBamFiles_HC; }

/*
 * Step 5. Base recalibration after indel realignment 
 */ 
process baserecal{
	tag "$prefix"
	cpus 4

	input:
	set val(prefix), file('*') from realignedBamFiles_forBaseRecal

	output:
	set file('recal.csv') into baseRecal_csv
	set file('recal_plots.pdf') into baseRecal_plot
	set file('recal_data.table') into baseRecal_table

	"""
	sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i realigned.bam --algo QualCal -k ${params.dbsnp} -k ${params.known_sites} recal_data.table
	sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i realigned.bam -q recal_data.table --algo QualCal -k ${params.dbsnp} -k ${params.known_sites} recal_data.table.post
	sentieon driver -t ${task.cpus} --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
	sentieon plot bqsr -o recal_plots.pdf recal.csv
	"""
}

baseRecal_table.into{ baseRecal_table_forUG; baseRecal_table_forHC; }

/*
 * Step 6. Variant calling for analysis ready-bam files
 */ 
process variant_call_UG{
	tag "$prefix"
	cpus 4
    
	publishDir "${params.outdir}/$prefix/variant_UG"

	input:
	set val(prefix), file('*') from realignedBamFiles_UG
	set file('*') from baseRecal_table_forUG

	output:
	set file('output-ug.vcf.gz*') into variant_UG

	"""
	sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i realigned.bam -q recal_data.table --algo Genotyper -d ${params.dbsnp} --var_type both --emit_conf=10 --call_conf=30 output-ug.vcf.gz
	"""
}

process variant_call_HC{
	tag "$prefix"
	cpus 4
	
	publishDir "${params.outdir}/$prefix/variant_HC"

	input:
	set val(prefix), file('*') from realignedBamFiles_HC
	set file('*') from baseRecal_table_forHC

	output:
	set file('output-hc.vcf.gz*') into variant_HC

	"""
	sentieon driver -r ${params.bwa_index} -t ${task.cpus} -i realigned.bam -q recal_data.table --algo Haplotyper -d ${params.dbsnp} --emit_conf=10 --call_conf=30 --prune_factor=3 output-hc.vcf.gz
	"""
}

/*
 * Step 7. Variant annotation
 */

/*
 * Helper function, given a file Path
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 *
 * For example:
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 *
 * Returns:
 *   'file_alpha'
 */
 
def readPrefix( Path actual, template ) {
     
	final fileName = actual.getFileName().toString()

	def filePattern = template.toString()
	int p = filePattern.lastIndexOf('/')
	if( p != -1 ) filePattern = filePattern.substring(p+1)
	if( !filePattern.contains('*') && !filePattern.contains('?') )
	filePattern = '*' + filePattern

	def regex = filePattern
		.replace('.','\\.')
		.replace('*','(.*)')
		.replace('?','(.?)')
		.replace('{','(?:')
		.replace('}',')')
		.replace(',','|')

	def matcher = (fileName =~ /$regex/)
	if( matcher.matches() ) {
		def end = matcher.end(matcher.groupCount() )
		def prefix = fileName.substring(0,end)
		while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
		prefix=prefix[0..-2]
		return prefix
	}

	//if (fileName.toString() =~ /\-\-/) println 'ok'

	//          def (value1, value2) = fileName.tokenize( '--' )


	return fileName
}
