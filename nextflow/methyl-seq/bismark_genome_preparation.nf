

params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].bismark
params.bismark_path = '/BiO/BioTools/Bismark/0.16.1'

params.input_fasta_path = '/BiOfs/BioPeople/brandon/genomes/Hsapiens/GRCh37/bismark'

Channel
	.fromPath( params.input_fasta_path )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { fasta_path }

process bismark_genome_preparation{
	publishDir "${params.index}"

	input:
	file(path_to_genome_folder:'*') from fasta_path

	"""
	${params.bismark_path}/bismark_genome_preparation --verbose ${path_to_genome_folder}
	"""
}
