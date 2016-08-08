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
params.gatk = "/BiO/BioTools/miniconda3/opt/gatk-3.6/GenomeAnalysisTK.jar" 
params.genome_version = "GRCh37.75"

params.name = "Best practice for WGS analysis"

// Input reads
params.reads = '/BiO/BioPeople/brandon/test_nextflow_wgs/*_{1,2}.fq.gz'
params.outdir = '/BiO/BioPeople/brandon/test_nextflow_wgs/outdir'

log.info "===================================="
log.info " WGS Best Practice v${version}"
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

/* Create a channel for read files */
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


/* STEP 1 - FastQC */
process fastqc {
	tag "$prefix"
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'

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

/* Step 2. Maps each read-pair by using mapper and sorting (removeDuplication) */
process mapping {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	cpus 4

	publishDir "${params.outdir}/BAMs"

	input: 
	set val(prefix), file (reads:'*') from read_files_mapping

	output:
	set prefix, file { "${prefix}.bam" }, file { "${prefix}.bam.bai" } into bamFiles
	file('*.txt') into dedup_statFiles

	"""
	bwa mem -M -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina' -t ${task.cpus} ${params.bwa_index} $reads | /BiO/BioTools/samblaster/samblaster_addMetrics/samblaster --metricsFile ${prefix}.txt | samtools view -u -Sb - | samtools sort - -o ${prefix}.bam
	samtools index ${prefix}.bam
	"""
}

bamFiles.into{ bamFiles_metrics; bamFiles_qualimap; bamFiles_BaseRecal }

/* Step 3. metrics for mapping */
process metrics{
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	cpus 4
	memory '16 GB'

	publishDir "${params.outdir}/$prefix/metrics", mode: 'copy'

	input:
	set val(prefix), file(read:'*'), file('*') from bamFiles_metrics

	output:
	set file('*.metrics') into metrics_txts

	"""
	java -Xmx16g -Djava.io.tmpdir=./ -jar /BiO/BioTools/bcbio/data/anaconda/share/picard-1.141-3/picard.jar CollectWgsMetrics I=${read} O=${prefix}.metrics R=${params.bwa_index}
	"""
}

process qualimap {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	cpus 8
	memory '8 GB'

	publishDir "${params.outdir}/qualimap"

	input:
	set val(prefix), file(bam:'*'), file(bai:'*') from bamFiles_qualimap

	output:
	val(prefix) into qualimap_prefix
	file ( "$prefix/genome_results.txt" ) into qualimap_genome_results_txt
	file ( "$prefix/raw_data_qualimapReport/*.txt" ) into qualimap_report_txts
	file ( "$prefix/images_qualimapReport/*.png" ) into qualimap_report_pngs
	file ( "$prefix/css/*" ) into qualimap_report_css
	file ( "$prefix/qualimapReport.html" ) into qualimap_report_html
	
	"""
	qualimap bamqc --java-mem-size=8G -bam ${bam} -gd HUMAN -nt ${task.cpus} -outdir $prefix
	"""
}

/* Step 5. Base recalibration */ 
process baserecal{
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	cpus 4
	memory '16 GB'

	input:
	set val(prefix), file(bam:'*'), file(bai:'*') from bamFiles_BaseRecal

	output:
	set file('recal_data.table') into baserecal_table
	set prefix, file("${prefix}.recal.bam"), file { "${prefix}.recal.bai" }  into recalBamFiles

	"""
	# Analyze patterns of covariation in the sequence dataset
	java -Xmx16g -Djava.io.tmpdir=./ -jar ${params.gatk} -T BaseRecalibrator -R ${params.bwa_index} -I ${bam} -knownSites ${params.dbsnp} -knownSites ${params.known_sites} -o recal_data.table
	# Do a second pass to analyze covariation remaining after recalibration
	java -Xmx16g -Djava.io.tmpdir=./ -jar ${params.gatk} -T BaseRecalibrator -R ${params.bwa_index} -I ${bam} -knownSites ${params.dbsnp} -knownSites ${params.known_sites} -BQSR recal_data.table -o post_recal_data.table
	# Generate before/after plots
	java -Xmx16g -Djava.io.tmpdir=./ -jar ${params.gatk} -T AnalyzeCovariates -R ${params.bwa_index} -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf
	# Apply the recalibration to your sequence datia
	java -Xmx16g -Djava.io.tmpdir=./ -jar ${params.gatk} -T PrintReads -R ${params.bwa_index} -I ${bam} -BQSR recal_data.table -o ${prefix}.recal.bam
	"""
}

process variant_call{
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	cpus 4
	
	publishDir "${params.outdir}/$prefix/variant_HC"

	input:
	set val(prefix), file(bam:'*'),file from recalBamFiles

	output:
	set val(prefix), file('raw_variant.g.vcf*') into variant_HC

	"""
	java -Xmx16g -Djava.io.tmpdir=./ -jar ${params.gatk} -T HaplotypeCaller -R ${params.bwa_index} -I ${prefix}.recal.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -o raw_variant.g.vcf
	"""
}

/* Step 7. Variant annotation */
process snpeff{
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wgs/bin:$PATH'
	tag "$prefix"
	
	publishDir "${params.outdir}/$prefix/variant"

	input:
	set val(prefix), file('*') from variant_HC

	"""
	java -Xmx16g -jar /BiO/BioTools/miniconda3/envs/wgs/share/snpeff-4.3-2/snpEff.jar eff -c /BiO/BioTools/miniconda3/envs/wgs/share/snpeff-4.3-2/snpEff.config ${params.genome_version} raw_variant.g.vcf > ${prefix}.eff.vcf
	mv snpEff_genes.txt ${prefix}.snpEff_genes.txt
	mv snpEff_summary.html ${prefix}.snpEff_summary.html
	"""
}

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
