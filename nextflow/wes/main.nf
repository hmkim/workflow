#!/usr/bin/env nextflow

/*
---------------------------------------------------------------------------------------
The pipeline can determine whether the input data is single or paired end. This relies on
specifying the input files correctly. For paired en data us the example above, i.e.
'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - from bcbio
*/

// Pipeline version
version = 0.1

// Configurable variables
params.genome = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
params.bwa_index = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa'
params.known_sites = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/Mills_and_1000G_gold_standard.indels.vcf.gz'
params.dbsnp = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz'
params.gatk = "/BiO/BioTools/miniconda3/opt/gatk-3.6/GenomeAnalysisTK.jar" 
params.genome_version = "GRCh37.75"

params.name = "Best practice for WES analysis"

// Input reads
params.reads = '/BiO/BioProjects/TBI-Human-Exome-2016-04/project1-merged/work/align_prep/TN1603D1281*{R1,R2}.fastq.gz'
params.outdir = '/BiO/BioPeople/brandon/test_nextflow_wgs/outdir'

log.info "===================================="
log.info " WGS Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="

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
 
read_files.into { read_for_grabix }

process grabix {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/wes/bin:$PATH'

	input:
	set val(prefix), file(reads:'*') from read_for_grabix

	"""
	grabix index ${reads}
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
