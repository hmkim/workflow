#!/usr/bin/env nextflow
// https://github.com/hbc/li_hiv


log.info "insertion analysis - N F  ~  version 0.1"
log.info "====================================="
log.info "genome                 : ${params.genome}"
log.info "reads                  : ${params.reads}"
log.info "outdir                 : ${params.outdir}"
log.info "\n"


/*
 * Input parameters validation
 */

genomeFile             = file(params.genome)


/*
 * validate input files/
 */

if( !genomeFile.exists() ) exit 1, "Missing genome directory: ${genomeFile}"


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


// preparing the genome

//aligning a sample
process mapping{
	publishDir "${params.outdir}/$prefix/";

    cpus params.bwa.cpus

    input:
    set val(prefix), file(read:'*') from read_files

    output:
    set val(prefix), file('*.bam') into mapped_bam

    """
    bwa mem -Y -t ${params.bwa.cpus} ${genomeFile} ${read} | samtools view -Su - | samtools sort --threads ${params.bwa.cpus} -l 0 -o ${prefix}.bam - 
    """


}

//mark duplicates
process markDuplicates{
	publishDir "${params.outdir}/$prefix/";
	cpus params.sambamba.cpus

    input:
    set val(prefix), file(read:'*') from mapped_bam

    output:
    set val(prefix), file('*.markdup.bam') into markdup_bam

    """
    sambamba markdup --overflow-list-size 600000 --tmpdir=./ --nthreads=${params.sambamba.cpus} ${read} ${prefix}.markdup.bam
    """
}

markdup_bam.into { bam1; bam2; }
// keep only chimeric alignments
process extractChimericReads {
	publishDir "${params.outdir}/$prefix/";

    input:
    set val(prefix), file(read:'*') from bam1

    """
    python ${baseDir}/scripts/chimeric.py ${params.targetContigName} ${read}
    python ${baseDir}/scripts/orientation.py ${prefix}.markdup.chimeric.igv.bam ${params.targetContigName} > ${prefix}.chimeric.igv.orientation.table
    """
}

//process runInsertion{
//	publishDir "${params.outdir}/$prefix/";
//
//    input:
//    set val(prefix), file(read:'*') from bam2
//
//    """
//    python ${baseDir}/scripts/insertion.py ${read} ${params.targetContigName}
//    """
//
//}



// ===================== UTILITY FUNCTIONS ============================


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


