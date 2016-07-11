#!/usr/bin/env nextflow

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

params.reads = '/BiO/BioPeople/brandon/test_nextflow/smrna/test_data/SRR95089{2,3}.fastq.gz'
params.adapter = 'TGGAATTCTCGGGTGC'
params.trna_mature_pre = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/srnaseq/trna_mature_pre.fa'
params.db_srnaseq = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/srnaseq'
params.db_star = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/star'
params.ref_fasta = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/seq/hg19.fa'
params.mature_fa = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/srnaseq/mature.fa'
params.hairpin_fa = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/srnaseq/hairpin.fa' 
params.species = 'hsa'
params.Rfam_for_miRDeep = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/srnaseq/Rfam_for_miRDeep.fa'
params.gtf_ref_transcripts = '/BiO/BioTools/bcbio/data/genomes/Hsapiens/hg19/rnaseq/ref-transcripts.gtf'

params.script_dir = '/BiO/BioPeople/brandon/test_nextflow/smrna/scripts'
 
log.info "===================================="
log.info " Small RNA Sequencing Best Practice v${version}"
log.info "===================================="

/*
 * Create a channel for read files
 */

Channel
     .fromPath( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
     .map { path ->
	prefix = path.getFileName().toString() - ~/\.fastq\.gz/
        tuple(prefix, path)
     }
     .groupTuple(sort: true)
     .set { read_files }

read_files.into { read_files_fastqc; read_files_trimming; }

process runFastQC_original{
	cpus params.fastqc.cpus
	
	input:
	set val(prefix),  file(read) from read_files_fastqc

	output:
	file '*_fastqc.html' into fastqc_html_original
	file '*_fastqc.zip' into fastqc_zip_original

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/fastqc -t ${params.fastqc.cpus} ${read}
	"""
}

process trim{
	input:
	set val(prefix), file(read) from read_files_trimming

	output:
	set val(prefix), file('*_trimmed.*') into readTrimmed
	file '*_trimming_report.txt' into trimgalore_results
	
	script:
	"""
	export PATH=/BiO/BioTools/bcbio/data/anaconda/bin:$PATH
	"""
	if ( params.adapter )
		"""
		/BiO/BioTools/trim_galore/0.4.1/trim_galore --adapter ${params.adapter} ${read}
		"""
	else
		"""
		/BiO/BioTools/trim_galore/0.4.1/trim_galore ${read}
		"""
}
/*
* /BiO/BioTools/bcbio/data/anaconda/bin/cutadapt --adapter=${params.adapter} --minimum-length=8 --untrimmed-output=\${NAMEBASE}-untrimmed-output -o \${NAMEBASE}-clean.fastq.gz -m 17 --overlap=8 ${read} --too-short-output \${NAMEBASE}.short.fastq.gz
*/

readTrimmed.into { readTrimmed_fastqc; readTrimmed_TdrMapping; readTrimmed_collapse; }

process runFastQC_trimmed {
	cpus params.fastqc.cpus

	input:
	set val(prefix), file(read) from readTrimmed_fastqc

	output:
	file '*_fastqc.html' into fastqc_html_trimmed
	file '*_fastqc.zip' into fastqc_zip_trimmed

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/fastqc -t ${params.fastqc.cpus} ${read}
	"""
	}

process collapseRead{
	input:
	set val(prefix), file(read) from readTrimmed_collapse

	output:
	val prefix into collapsedRead_prefix
	file('*_trimmed.fastq') into collapsedRead_files

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/python ${params.script_dir}/collapse.py ${read}
	"""
}

/* https://github.com/sararselitsky/tDRmapper */
process TdrMapping{
	input:
	set val(prefix), file(read) from readTrimmed_TdrMapping

	output:
	file '*' into tdrMapping_results

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/TdrMappingScripts.pl ${params.trna_mature_pre} ${read}
	"""
}

collapsedRead_files.into { crf_prepare_makeList; crf_prepare; crf_miraligner; }
collapsedRead_prefix.into { crp_prepare_makeList; crp_prepare; crp_miraligner; }

/* http://seqcluster.readthedocs.io/mirna_annotation.html */
process miraligner{
	input:
	val(prefix) from crp_miraligner;
	file(read) from crf_miraligner;

	output:
	file '*.mirna' into miraligner_map
	file '*.mirna.nomap' into miraligner_nomap
	file '*.mirna.opt' into miraligner_opt

	"""
	export PATH=/BiO/BioTools/bcbio/data/anaconda/bin:\$PATH
	/BiO/BioTools/bcbio/data/anaconda/bin/miraligner -Xms705m -Xmx4500m -sub 1 -trim 3 -add 3 -s ${params.species} -i ${read} -db ${params.db_srnaseq} -o sample
	/BiO/BioTools/bcbio/data/anaconda/bin/python ${params.script_dir}/miraligner_parser.py sample.mirna sample.bak sample
	"""
}

process seqcluster_prepare_makeList{
	input:
	val(prefix) from crp_prepare_makeList.toList()
	file(read) from crf_prepare_makeList.toList()

	output:
	file 'fileList.txt' into prepare_filelist

	echo true

	"""
	IN="${read}"
	IN2="${prefix}"

	IN2=\${IN2:1:-1}

	OIFS=\$IFS
	IFS=', '

	array1=( \$IN )
	array2=( \$IN2 )

	for i in "\${!array1[@]}"; do
		echo -e "\${array1[i]}\\t\${array2[i]}" >> fileList.txt
	done
	"""
}

process seqcluster_prepare{
	input:
	file filelist from prepare_filelist
	val(prefix) from crp_prepare.toList()
	file(fastq) from crf_prepare.toList()

	output:
	file 'seqs.ma' into seqs_ma
	file 'seqs.fa' into seqs_fa

	"""
	python ${params.script_dir}/prepare_data.py ${filelist}
	"""
}

seqs_ma.into { mirdeep2_seqs_ma; seqcluster_cluster_seqs_ma; }

process mapping{
	cpus params.sambamba.cpus
	cpus params.star.cpus

	input:
	set file(read) from seqs_fa

	output:
	file 'seqs.bam' into star_bam
	file 'seqs.bam.bai' into star_bam_index
	file 'Log.final.out' into star_log

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/STAR --genomeDir $params.db_star --readFilesIn $read --runThreadN ${params.star.cpus} --outFileNamePrefix ./ --outReadsUnmapped Fastx --outFilterMultimapNmax 1000 --outStd SAM --alignIntronMax 1 --outSAMunmapped Within --outSAMattributes NH HI NM MD AS  --sjdbGTFfile $params.gtf_ref_transcripts --sjdbOverhang 39  --outSAMattrRGline ID:miRQC_A PL:illumina PU:1_2016-04-08_mir2qc_bcbio SM:miRQC_A | /BiO/BioTools/bcbio/data/galaxy/../anaconda/bin/samtools sort -@ 1 -m 16G -o seqs.bam /dev/stdin
	/BiO/BioTools/bcbio/data/anaconda/bin/sambamba index -t ${params.sambamba.cpus} seqs.bam
	"""

}

process mirdeep2{
	input:
	set file(seqs_ma) from mirdeep2_seqs_ma
	set file(bam) from star_bam
	set file(bai) from star_bam_index

	"""
	/BiO/BioTools/bcbio/data/anaconda/bin/python ${params.script_dir}/mirdeep_prepare.py ${bam} ${seqs_ma}
	unset PERL5LIB && export PATH=/BiO/BioTools/bcbio/data/anaconda/bin:$PATH && perl /BiO/BioTools/bcbio/data/anaconda/bin/miRDeep2.pl file_reads.fa ${params.ref_fasta} align.bam ${params.mature_fa} none ${params.hairpin_fa} -f ${params.Rfam_for_miRDeep} -r simple -c -P -t ${params.species} -z res
	/BiO/BioTools/bcbio/data/anaconda/bin/python ${params.script_dir}/mirdeep_parse_novel.py result_res.csv ${params.species}
	"""
}
