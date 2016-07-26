#!/usr/bin/env nextflow

params.reads = '/BiO/BioTools/JAFFA/Demo/*.fastq.gz'

Channel
	.fromPath( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_files_pattern }


process assembly {
	beforeScript 'set +u; source activate genefusion; set -u'

	input:
	val (fastq_path:'*') from read_files_pattern

	"""
	/BiO/BioTools/JAFFA/1.07/tools/bin/bpipe run /BiO/BioTools/JAFFA/1.07/JAFFA_assembly.groovy ${fastq_path}
	"""
}
