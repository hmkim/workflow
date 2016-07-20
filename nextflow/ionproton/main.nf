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
params.genome = '/BiO/BioResources/References/Human/hg19_LT/hg19.fasta'
params.known_sites = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/Mills_and_1000G_gold_standard.indels.vcf.gz'
params.dbsnp = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz'
params.gtf = ''
params.rlocation = ''
params.gatk = "/BiO/BioTools/miniconda3/opt/gatk-3.6/GenomeAnalysisTK.jar" 
params.genome_version = "GRCh37.75"

params.name = "Best practice for Proton target analysis (Ion Proton)"
params.tmap = "/BiO/BioTools/TS/build/Analysis/TMAP/tmap"

tmap = file(params.tmap)

// Input reads
params.reads = '/backup/wolf/sgpark/Projects/TeRef/RawData/20141014_ProtonExome_11Sample/RawData/*fastq.gz'
// Output directory
params.outdir = 'result'

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
 * Step 2. Maps each read-pair by using mapper and sorting (removeDuplication)
 */
process mapping {
	tag "$prefix"
	cpus 24

	publishDir "${params.outdir}/BAMs"

	input: 
	set val(prefix), file (read:'*') from read_files_mapping

	output:
	set prefix, file { "${prefix}.bam" }, file { "${prefix}.bam.bai" } into bamFiles

	"""
	# -v,--verbose                                            print verbose progress information [false]
	# -u,--rand-read-name                                     specifies to randomize based on the read name [false]
	# -y,--softclip-key                                       soft clip only the last base of the key [false]
	# --prefix-exclude                         INT            specify how many letters of prefix of name to be excluded when do randomize by name [false]
	# -o,--output-type                         INT            the output type [0]
	#   0 - SAM
	#   1 - BAM (compressed)
	#   2 - BAM (uncompressed)
        # -J,--max-adapter-bases-for-soft-clipping INT            specifies to perform 3' soft-clipping (via -g) if at most this # of adapter bases were found (ZB tag) [2147483647]
	# --end-repair                             INT            specifies to perform 5' end repair [0]
	#  0 - disable
	#  1 - prefer mismatches
	#  2 - prefer indels
	#  >2 - specify %% Mismatch above which to trim end alignment
	#  --do-repeat-clip                                        clip tandem repeats from the alignment 3' ends [false]
	#  --context                                               realign with context-dependent gap scores [0]

	${tmap} mapall -n ${task.cpus} -f ${params.genome} -r ${read} -v -Y -u --prefix-exclude 5 -o 2 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4 | samtools sort -m 1000M -l1 -@12 - /results/analysis/output/Home/Auto_user_Proton01-380-20160622_TBD160346_4_555_829/plugin_out/variantCaller_out.2771/
	"""
}

bamFiles.into{ bamFiles_metrics; bamFiles_qualimap; bamFiles_BaseRecal }

/*
 * Step 3. metrics for mapping
 */
process metrics{
    tag "$prefix"
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

/*
 * Step 5. Base recalibration 
 */ 
process baserecal{
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

/*
 * Step 7. Variant annotation
 */

process snpeff{
	beforeScript 'set +u; source activate wgs; set -u'
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
