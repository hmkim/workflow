#!/usr/bin/env nextflow

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1
 
log.info "===================================="
log.info " Small RNA Sequencing Best Practice v${version}"
log.info " Read: ${params.reads}"
log.info " Adapter: ${params.adapter}"
log.info " Result directory: ${params.outdir}"
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
	beforeScript 'set +u; source activate smrna; set -u'
	cpus params.fastqc.cpus

	publishDir "${params.outdir}/fastqc";

	input:
	set val(prefix),  file(read) from read_files_fastqc

	output:
	file '*_fastqc.html' into fastqc_html_original
	file '*_fastqc.zip' into fastqc_zip_original

	"""
	fastqc -t ${params.fastqc.cpus} ${read}
	"""
}

process trim{
	beforeScript 'set +u; source activate smrna; set -u'
	
	publishDir "${params.outdir}/trim";

	input:
	set val(prefix), file(read) from read_files_trimming

	output:
	set val(prefix), file('*_trimmed.*') into readTrimmed
	file '*_trimming_report.txt' into trimgalore_results

	script:	
	if ( params.adapter )
		"""
		trim_galore --adapter ${params.adapter} ${read}
		"""
	else
		"""
		trim_galore ${read}
		"""
}

readTrimmed.into { readTrimmed_fastqc; readTrimmed_TdrMapping; readTrimmed_collapse; }

process runFastQC_trimmed {
	beforeScript 'set +u; source activate smrna; set -u'
	cpus params.fastqc.cpus
	
	publishDir "${params.outdir}/fastqc_trimmed";

	input:
	set val(prefix), file(read) from readTrimmed_fastqc

	output:
	file '*_fastqc.html' into fastqc_html_trimmed
	file '*_fastqc.zip' into fastqc_zip_trimmed

	"""
	fastqc -t ${params.fastqc.cpus} ${read}
	"""
	}

process collapseRead{
	beforeScript 'set +u; source activate smrna; set -u'

	input:
	set val(prefix), file(read) from readTrimmed_collapse

	output:
	val prefix into collapsedRead_prefix
	file('*_trimmed.fastq') into collapsedRead_files

	"""
	python ${params.script_dir}/collapse.py ${read}
	"""
}

/* https://github.com/sararselitsky/tDRmapper */
process TdrMapping{
	publishDir "${params.outdir}/TdrMapping";

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
	python ${params.script_dir}/miraligner_parser.py sample.mirna sample.bak sample
	"""
}

process seqcluster_prepare_makeList{
	beforeScript 'set +u; source activate smrna; set -u'

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
	beforeScript 'set +u; source activate smrna; set -u'

	publishDir "${params.outdir}/seqcluster";

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
	beforeScript 'set +u; source activate smrna; set -u'

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
	sambamba index -t ${params.sambamba.cpus} seqs.bam
	"""

}

process mirdeep2{
	beforeScript 'set +u; source activate smrna; set -u'

	publishDir "${params.outdir}/mirdeep2";

	input:
	set file(seqs_ma) from mirdeep2_seqs_ma
	set file(bam) from star_bam
	set file(bai) from star_bam_index

	"""
	python ${params.script_dir}/mirdeep_prepare.py ${bam} ${seqs_ma}
	unset PERL5LIB && export PATH=/BiO/BioTools/bcbio/data/anaconda/bin:$PATH && perl /BiO/BioTools/bcbio/data/anaconda/bin/miRDeep2.pl file_reads.fa ${params.ref_fasta} align.bam ${params.mature_fa} none ${params.hairpin_fa} -f ${params.Rfam_for_miRDeep} -r simple -c -P -t ${params.species} -z res
	python ${params.script_dir}/mirdeep_parse_novel.py result_res.csv ${params.species}
	"""
}
