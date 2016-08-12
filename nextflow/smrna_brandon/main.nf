#!/usr/bin/env nextflow

//Set default values for params
params.steps = 'quantification,variantcall'
params.dbFile = 'pipeline.db'

params.nextflex = false // default is illumina sequencling.
params.notrim = false // default is trim.

// Clear pipeline.db file
pdb = file(params.dbFile)
pdb.write('')

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

/* If --help in parameters, print software usage */
if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------'
    log.info 'Small RNAseq workflow'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main.nf --index INDEX_FILE --reads \'path/to/data/sample_*_{1,2}.fq.gz\' [other option]'

    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --reads          Reads {FASTQ}            Sample fastq file.'
    log.info '    --sample_info    SAMPLE INFO              Sample control/tumor information.'
    log.info '    --index INDEX_FILE                  Index file.'
    log.info 'Options:'
    log.info '    --notrim                             A BED file for calling. (WES)'
    log.info '    --nextflex         BED FILE                 A BED file for calling. (WES)'
    log.info ''
    exit 1
}


kit_type = params.nextflex ? "Nextflex" : "Illumina"
flag_trim = params.notrim ? "no" : "yes"

// check mandatory options
if (!params.index) {
    exit 1, "Index file not specified"
}
if (!params.bowtie_idx) {
    exit 1, "bowtie Index file not specified"
}

// Pipeline version
version = 0.1
 
log.info "===================================="
log.info " Small RNA Sequencing Best Practice v${version} by brandon"
log.info " Index file : ${params.index}"
log.info " Genome : ${params.genome}"
log.info " Annotation : ${params.annotation}"
//log.info " Read: ${params.reads}"
log.info " Adapter: ${params.adapter}"
log.info " Trimming: ${flag_trim}"
log.info " Result directory: ${params.outdir}"
log.info " Kit type (--nextflex) : ${kit_type}"
log.info "===================================="


//genomes = params.genome.split(',').collect { file(it) }
//annos = params.annotation.split(',').collect { file(it) }


/* Read index file and make input setting */
index = params.index ? file(params.index) : System.in
input_files = Channel.create()
input_chunks = Channel.create()


data = ['samples': [], 'ids': []]
input = Channel
	.from(index.readLines())
	.filter { it }  // get only lines not empty
	.map { line ->
	def (sampleId, runId, fileName, format, readId) = line.split()
		return [sampleId, runId, resolveFile(fileName, index), format, readId]
	}

input.subscribe onNext: {
        sample, id, path, type, view ->
        items = "sequencing runs"
        input_chunks << tuple(sample, id, path, type, view)
        data['samples'] << sample
        data['ids'] << id },
    onComplete: {
        ids=data['ids'].unique().size()
        samples=data['samples'].unique().size()
        log.info "Dataset information"
        log.info "-------------------"
        log.info "Number of sequenced samples     : ${samples}"
        log.info "Number of ${items}       : ${ids}"
        log.info "Merging                         : ${ ids != samples ? 'by sample' : 'none' }"
        log.info ""
        input_chunks << Channel.STOP
    }

input_bams = Channel.create()
bam = Channel.create()

input_chunks
    .groupBy {
        sample, id, path, type, view -> id
    }
    .flatMap ()
    .choice(input_files, input_bams) { it -> if ( it.value[0][3] == 'fastq' ) 0 else if ( it.value[0][3] == 'bam' ) 1}


input_files = input_files.map {
    [it.key, it.value[0][0], it.value.collect { sample, id, path, type, view -> path }, fastq(it.value[0][2]).qualityScore()]
}

def msg = "Output files db"
log.info "=" * msg.size()
log.info msg + " -> ${pdb}"
log.info "=" * msg.size()
log.info ""



///* Create a channel for read files */
//Channel
//	.fromPath( params.reads )
//	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
//	.map { path ->
//		def prefix = readPrefix(path, params.reads)
//		tuple(prefix, path)
//	}
//	.groupTuple(sort: true)
//	//.view()
//	.set { read_files }


trimmed_files = Channel.create()

if (!params.notrim){
	process trim{
		input:
		set id, sample, file(reads), qualityOffset from input_files

		output:
		set id, sample, file('*_trimmed.*'), qualityOffset into adapter_trimmed_files
		set file('*_trimming_report.txt') into trimmed_report
		//RPA_160622_1.fq  RPA_160622_1.fq_trimming_report.txt  RPA_160622_1_trimmed.fq
		//set id, sample, tpe, view, file("*.bam"), pairedEnd into bam

		script:
		if ( params.adapter )
			"""
			trim_galore --adapter ${params.adapter} ${reads}
			"""
		else
			"""
			trim_galore ${fqs}
			"""
	}
	trimmed_files = adapter_trimmed_files
}

if (params.nextflex){
	process trim_3bp{
		input:
		set id, sample, file(read:'*'), qualityOffset from adapter_trimmed_files

		output:
		set id, sample, file('*.clean.fastq'), qualityOffset into nextflex_trimmed_files

		"""
		cutadapt -u 4 -u -4 -o ${id}.clean.fastq ${read}
		"""
	}
	trimmed_files = nextflex_trimmed_files
}

trimmed_files.into{ trimmed_files1; trimmed_files2 }

if ('quantification' in pipelineSteps || 'variantcall' in pipelineSteps) {
	process mapUsingBowtie{
		input:
		set id, sample, file(read:'*'), qualityOffset from trimmed_files1

		output:
		set id, sample, file('*') into sam
	
		script:
		bowtie_quals = '--phred64-quals'
		if (qualityOffset == 33){
			bowtie_quals = '--phred33-quals'
		}
		"""
		bowtie ${bowtie_quals} -p 4 -S -q -n 1 -e 80 -l 30 -a -m 5 --best --strata --sam-RG ID:${id} --sam-RG SM:'${sample}' ${params.bowtie_idx} ${read} ${id}.bowtie.sam
		"""
	}

	sam.into { sam1; sam2 }

	if ('variantcall' in pipelineSteps){
		process process_variantcall{
			input:
			set id, sample, file(samfile:'*') from sam1

			"""
			samtools view -SH ${id}.bowtie.sam > ${id}.bowtie.mq.sam
			samtools view -S ${id}.bowtie.sam | awk -F"\\t" 'OFS="\\t"{\$5 = 60; print \$0}' >> ${id}.bowtie.mq.sam
			picard -Djava.io.tmpdir=./ -Xmx8g AddOrReplaceReadGroups INPUT=${id}.bowtie.mq.sam OUTPUT=${id}.bowtie.mq.coord_sorted.bam SORT_ORDER=coordinate RGID=${id} RGPU=${id} RGSM='${sample}' VALIDATION_STRINGENCY=SILENT RGLB=${id} RGCN=TheragenEtex RGPL=illumina
			samtools index ${id}.bowtie.mq.coord_sorted.bam
			""" 
		}
	}

	if ('quantification' in pipelineSteps){
		process process_quantification{
			input:
			set id, sample, file(samfile:'*') from sam2

			"""
			picard -Djava.io.tmpdir=./ -Xmx8g SortSam INPUT=${id}.bowtie.sam OUTPUT=${id}.queryname_sorted.bam SORT_ORDER=queryname
			samtools view ${id}.queryname_sorted.bam | htseq-count -m intersection-nonempty -q -t exon -s no - ${params.gtf} > ${id}.genecount.txt
			perl ${baseDir}/scripts/make_quantification_annotation_file.pl ${id}.genecount.txt ${params.gtf} > ${id}.genecount.annotation.txt
			Rscript ${baseDir}/scripts/gencode_piechart.R ${id}.genecount.annotation.txt ${id} ./
			"""
		}
	}
}

process mapping_for_mirdeep2{
	input:
	set id, sample, file(read:'*'), qualityOffset from trimmed_files2

	output:
	set file('reads.fa') into mapped_reads
	"""
	mapper.pl ${read} -e -h -q -m -u -v -o 4 -p ${params.bowtie_idx} -s reads.fa -t reads_vs_genome.arf 2> mapped.log
	"""			
}

process mapRead_length_distribution{
	input:
	file(read:'*') from mapped_reads

	"""
	${baseDir}/scripts/smallRNA_len_dist.pl ${read} uniq 
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

/*
 * Given a string path resolve it against the index file location.
 * Params: 
 * - str: a string value represting the file pah to be resolved
 * - index: path location against which relative paths need to be resolved 
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str) 
  }
} 

def testResolveFile() {
  def index = file('/path/to/index')
  assert resolveFile('str', index) == file('/path/to/str')
  assert resolveFile('/abs/file', index) == file('/abs/file')
  assert resolveFile('s3://abs/file', index) == file('s3://abs/file')
}
