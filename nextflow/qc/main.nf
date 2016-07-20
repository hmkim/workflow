#!/usr/bin/env nextflow

/*
---------------------------------------------------------------------------------------
The pipeline can determine whether the input data is single or paired end. This relies on
specifying the input files correctly. For paired en data us the example above, i.e.
'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - FastQC - read quility control
 - BWA - align
 - MultiQC
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Configurable variables
params.bwa_index = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa'
//params.bed = '/BiO/BioResources/References/Human/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed'
params.bed = '/BiO/BioPeople/brandon/test_nextflow/qc/SureSelect_Human_All_Exon_V5.bed'

params.name = "QC Best practice"
//params.reads = '/BiO/BioPeople/brandon/test_nextflow/qc/TBD160338-test/TN*{R1,R2}_[0-9][0-9][0-9].fastq.gz'
//params.outdir = '/BiO/BioPeople/brandon/test_nextflow/qc/outdir-test'
params.reads = '/BiO/BioProjects/QCsystem/rawreads/160614_D00236_0448_AC9VP7ANXX/Unaligned/TBD160338/TN*{R1,R2}_[0-9][0-9][0-9].fastq.gz'
params.outdir = '/BiO/BioPeople/brandon/test_nextflow/qc/outdir'

log.info "===================================="
log.info " QC Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Index        : ${params.bwa_index}"
log.info "Region       : ${params.bed}"
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
//gtf   = file(params.gtf)
bed = file(params.bed)
if( !index.exists() ) exit 1, "Missing reference index: ${index}"
if( !bed.exists() ) exit 2, "Missing BED annotation: ${bed}"

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
 
read_files.into { read_files_fastqc; read_files_mapping; }

/*
 * STEP 1 - FastQC
 */

process fastqc {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/qc/bin:$PATH'

	tag "$prefix"

	memory { 2.GB * task.attempt }
	time { 4.h * task.attempt }

	errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
	maxRetries 3
	maxErrors '-1'

	publishDir "${params.outdir}/fastqc", mode: 'copy'

	input:
	set val(prefix), file(reads:'*') from read_files_fastqc

	output:
	file '*_fastqc.{zip,html}' into fastqc_results

	"""
	fastqc $reads
	"""
}

/*
 * Step 2. Maps each read-pair by using mapper
 */
process mapping {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/qc/bin:$PATH'

	tag "$prefix"
	cpus 5

	input: 
	set val(prefix), file (reads:'*') from read_files_mapping

	output:
	file { "${prefix}.bam*" } into bamFilesOut
	set prefix, file { "${prefix}.bam" }, file { "${prefix}.bam.bai" } into bamFiles

	"""
	bwa mem -M -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina' -t ${task.cpus} ${params.bwa_index} $reads | samblaster | samtools view -u -Sb - | samtools sort - -o ${prefix}.bam
	samtools index ${prefix}.bam
	"""
}

(bamFilesForCoverage,
 bamFilesForQualimap,
 bamFilesForVariantCalling) = bamFiles.separate(3) { x -> [ x, x, x ] }

process qualimap {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/qc/bin:$PATH'
	tag "$prefix"
	cpus 8

	publishDir "${params.outdir}/qualimap/${prefix}"

	input:
	set val(prefix), file(reads:"${prefix}.bam"), file (index:"${prefix}.bam.bai") from bamFilesForQualimap

	output:
	file ( "${prefix}/genome_results.txt" ) into qualimap_genome_results 
	file ( "${prefix}/raw_data_qualimapReport/*.txt" ) into qualimap_report_txts
	
	"""
	qualimap bamqc --java-mem-size=8G -bam ${reads} -c -gd HUMAN -nt ${task.cpus} --feature-file ${bed} 
	mv ${prefix}_stats/ ${prefix}
	"""
}

process multiqc{
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/qc/bin:$PATH'
	//beforeScript 'set +u; source activate qc; set -u'

	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
	val input_files from qualimap_genome_results.toList()
	file('*') from qualimap_report_txts.toList()
	file('*') from fastqc_results.toList()

	output:
	file ('multiqc_*') into multiqc_output

	"""
	multiqc ${params.outdir}
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
     return fileName
}
