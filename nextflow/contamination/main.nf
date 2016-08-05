//params.reads = "/BiO/BioPeople/brandon/test_nextflow/contamination/Hm2_GTGAAA_L005_{R1,R2}_[0-9][0-9][0-9].fastq.gz"
params.reads = "/BiO/BioPeople/brandon/test_fastq_screen/TBD160403_outdir/*/{forward,reverse}.fastq"
params.kraken_db = "/BiO/BioTools/bcbio/data/genomes/kraken/minikraken"
params.threads = '7'
params.out_dir = "/BiO/BioPeople/brandon/workflow/nextflow/contamination/outdir"

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

	memory '20 GB'

	input:
	set val(name), file(reads:'*') from read_files

	output:
	set val(name), file('sample.kraken') into kraken_output
	set val(name), file('sample.krona') into krona_input

	"""
	kraken --threads ${params.threads} --db ${params.kraken_db} --fastq-input --paired --check-names --output sample.kraken $reads 
	cut -f2,3 sample.kraken > sample.krona
	"""
}

process make_html{
	tag "reads: $name"

	publishDir "${params.out_dir}/krona"
	
	input:
	set val(name), file(input:'*') from krona_input

	output:
	set file('*.html*') into krona_output

	"""
	ktImportTaxonomy ${input} -o ${name}.html
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

    final fileName = actual.getParent().toString().split(/\//)[-1]

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
