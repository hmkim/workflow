
params.in = "/BiO/BioProjects/TBD160423-Cow-WES-2016-07/result/07-1_picard_dedup/*/*.dedup.bam"

params.bwa_index = "/BiO/BioResources/References/Bos_taurus/UMD_3.1/Bos_taurus.UMD3.1.dna.toplevel.fa"

listOfFiles = Channel.fromPath(params.in)

listOfFiles.into{ bamFiles; bamFiles_bai }

(bamFiles_bai - ~/\.bam/).println()

process samtools_variant_call{
	input:
	set file(bam) from bamFiles

	echo true
	"""
	echo ${bam}
	"""
}

