
params.in = "/BiO/BioProjects/TBD160403-TG-Human-RNAseq-2016-07/Rine_Quant/analysis/*/unmapped.bam"
params.out_dir = "/BiO/BioPeople/brandon/test_fastq_screen/TBD160403_outdir"

Channel
	.fromPath( params.in )
	.ifEmpty { error "Cannot find any reads matching: ${params.in}" }
	.map { path ->
        	def prefix = readPrefix(path, params.in)
        	tuple(prefix, path)
	}
	.groupTuple(sort: true)
	.set { bam_files }

bam_files.into{ test1; test2; }


process bam2fastq{
	input:
	set val(prefix), file(bam:'*') from test2

	output:
	set val(prefix), file('forward.fastq') into forward_fastq
	set val(prefix), file('reverse.fastq') into reverse_fastq
	file ('unpaired.fastq') into unpaired_fastq

	"""
	/BiO/BioTools/miniconda2/bin/java -Xmx8g -jar /BiO/BioTools/miniconda2/pkgs/picard-2.3.0-0/share/picard-2.3.0-0/picard.jar SamToFastq I=${bam} F=forward.fastq F2=reverse.fastq FU=unpaired.fastq VALIDATION_STRINGENCY=LENIENT
	"""
}

process bam2fastq_another{
	input:
	set val(prefix), file(bam:'*') from test1

	"""	
	#/BiO/BioTools/bam2fastq/bam2fastq-1.1.0/bam2fastq ${bam}
	"""
}

process fastq_screen{
	cpus 2

	publishDir "${params.out_dir}/${prefix}"

	input:
	set val(prefix), file(forward) from forward_fastq
	set val(prefix), file(reverse) from reverse_fastq

	output:
	set file('*_screen.*') into screen_result_files
	set file('*_no_hits.*') into screen_nohit_result_files

	"""
	fastq_screen --nohits --conf /BiO/BioPeople/brandon/workflow/nextflow/fastq_screen/fastq_screen.conf --outdir ./ ${forward} ${reverse} --threads ${task.cpus}



	"""
}

def readPrefix( Path actual, template ) {

        final fileName = actual.getParent().toString().split(/\//)[-1]

        return fileName
}
