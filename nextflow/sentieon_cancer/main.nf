#!/usr/bin/env nextflow

import groovy.io.FileType

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------'
    log.info 'Sentieon : SOMATIC VARIANT CALL'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main.nf --bed /BiO/BioProjects/CUK-Mouse-2016-07-TBD160465/Sureselect_mouse_All_Exon_V1.bed --raw_folder /BiO/BioProjects/CUK-Mouse-2016-07-TBD160465/rawdata --fasta_ref /BiO/BioResources/References/Mouse/Ens_67_NCBIM37/M_musculus_Ens67.chr.fa --sample_info sample_info.txt [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --pairs          FASTQ                   FASTQ files.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info '    --sample_info    SAMPLE INFO              Sample control/tumor information.'
    log.info 'Options:'
    log.info '    --bed            BED FILE                 A BED file for calling. (WES)'
    log.info ''
    exit 1
}


Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path ->
        def prefix = readPrefix(path, params.reads)
	tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read_files }

read_files.into{ reads_fastqc; reads_trim; }

process fastqc {
	tag "$prefix"

	memory { 2.GB * task.attempt }
	time { 4.h * task.attempt }

	errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
	//terminate
	maxRetries 3
	maxErrors '-1'

	publishDir "${params.outdir}/$prefix/fastqc", mode: 'copy'

	input:
	set val(prefix), file(reads:'*') from reads_fastqc

	output:
	file '*_fastqc.html' into fastqc_results_html
	file '*_fastqc.zip' into fastqc_results_zip

	"""
	fastqc $reads
	"""
}

process trim{
	publishDir "${params.outdir}/trimmed/$prefix";

	input:
	set val(prefix), file(read) from reads_trim

	output:
	set val(prefix), file('*_trimmed.*') into readTrimmed
	file '*_trimming_report.txt' into trimgalore_results

	script:	
	"""
	trim_galore ${read}
	"""
}



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
