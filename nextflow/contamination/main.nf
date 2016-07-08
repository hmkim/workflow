params.reads = "/BiO/BioPeople/brandon/test_nextflow/contamination/Hm2_GTGAAA_L005_{R1,R2}_[0-9][0-9][0-9].fastq.gz"
params.kraken_db = "/BiO/BioTools/bcbio/data/genomes/kraken/minikraken"
params.threads = '7'

log.info "Contamianation check version 0.1"
log.info "====================================="
log.info "reads                  : ${params.reads}"

Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path ->
       def prefix = readPrefix(path, params.reads)
       tuple(prefix, path)
    }
    .groupTuple(sort: true)
    .set { read_files }

process mapping {
    tag "reads: $name"

    executor 'sge'
    memory '20 GB'

    input:
    set val(name), file(reads:'*') from read_files

    """
    /BiO/BioTools/bcbio/data/anaconda/bin/kraken --threads ${params.threads} --db ${params.kraken_db} --fastq-input --paired --check-names --output sample.kraken $reads
    cut -f2,3 sample.kraken > sample.krona
    """
}

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
