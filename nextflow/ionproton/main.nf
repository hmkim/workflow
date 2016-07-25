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

mode = 'gzip'
// gzip or null

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
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/ionproton/bin:$PATH'
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

read_files_mapping.into { aa1; aa2}

aa2.println()

/*
 * Step 2. Maps each read-pair by using mapper and sorting (removeDuplication)
 */
process mapping {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/ionproton/bin:$PATH'
	tag "$prefix"
	cpus 24

	publishDir "${params.outdir}/BAMs"

	input: 
	//set val(prefix), file (read:'*') from read_files_mapping
	set val(prefix), file (read:'*') from aa1

	output:
	set prefix, file { "${prefix}.bam" } into bamFiles

	script:
	if( mode == 'gzip' )
		"""
		zcat ${read} > ${prefix}.fastq
		${tmap} mapall -n ${task.cpus} -f ${params.genome} -r ${prefix}.fastq -v -Y -u --prefix-exclude 5 -o 2 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4 | samtools sort -m 1000M -l1 --threads=${task.cpus} -o ${prefix} -
		"""
	else
		"""
		${tmap} mapall -n ${task.cpus} -f ${params.genome} -r ${read} -v -Y -u --prefix-exclude 5 -o 2 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4 | samtools sort -m 1000M -l1 --threads=${task.cpus} -o ${prefix}.bam -
		"""

//	-v,--verbose                                            print verbose progress information [false]
//	-u,--rand-read-name                                     specifies to randomize based on the read name [false]
//	-y,--softclip-key                                       soft clip only the last base of the key [false]
//	--prefix-exclude                         INT            specify how many letters of prefix of name to be excluded when do randomize by name [false]
//	-o,--output-type                         INT            the output type [0]
//	  0 - SAM
//	  1 - BAM (compressed)
//	  2 - BAM (uncompressed)
//        -J,--max-adapter-bases-for-soft-clipping INT            specifies to perform 3' soft-clipping (via -g) if at most this # of adapter bases were found (ZB tag) [2147483647]
//	--end-repair                             INT            specifies to perform 5' end repair [0]
//	 0 - disable
//	 1 - prefer mismatches
//	 2 - prefer indels
//	 >2 - specify %% Mismatch above which to trim end alignment
//	 --do-repeat-clip                                        clip tandem repeats from the alignment 3' ends [false]
//	 --context                                               realign with context-dependent gap scores [0]

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
1}
