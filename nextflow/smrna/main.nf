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
	cpus params.fastqc.cpus

	//publishDir "${params.outdir}/fastqc";

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
	publishDir "${params.outdir}/trimmed/$prefix";

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
	cpus params.fastqc.cpus
	
	//publishDir "${params.outdir}/fastqc_trimmed";

	input:
	set val(prefix), file(read) from readTrimmed_fastqc

	output:
	file '*_fastqc.html' into fastqc_html_trimmed
	set val(prefix), file('*_fastqc.zip') into fastqc_zip_trimmed

	"""
	fastqc -t ${params.fastqc.cpus} ${read}
	"""
	}

process qc{
	publishDir "${params.outdir}/qc/${prefix}/fastqc";

	input:
	file(html) from fastqc_html_trimmed
	set val(prefix), file(zip) from fastqc_zip_trimmed

	output:
	file("${extract_name}/*.tsv") into fastqc_result_tsv
	file("${extract_name}/fastqc_data.txt") into fastqc_data_txt
	file("${extract_name}/fastqc_report.html") into fastqc_report_html
	file("${prefix}.zip") into fastqc_zip

	exec:
	extract_name = zip.name.take(zip.name.lastIndexOf('.'))

	shell:
	"""
	ln -s ${zip} ${prefix}.zip
	unzip ${zip}
	python ${params.script_dir}/fastqc_parser.py ${extract_name} ${prefix}
	"""
}

process collapseRead{
	publishDir "${params.outdir}/trimmed/$prefix";

	input:
	set val(prefix), file(read) from readTrimmed_collapse

	output:
	val prefix into collapsedRead_prefix
	file('*.{fq,fastq}') into collapsedRead_files
	file('*_size_stats') into collapsedRead_size_stats

	"""
	python ${params.script_dir}/collapse.py ${read}
	"""
}

/* https://github.com/sararselitsky/tDRmapper */
process TdrMapping{
	publishDir "${params.outdir}/trna/$prefix";

	input:
	set val(prefix), file(read) from readTrimmed_TdrMapping

	output:
	file '*' into tdrMapping_results

	"""
	TdrMappingScripts.pl ${params.trna_mature_pre} ${read}
	"""
}

collapsedRead_files.into { crf_prepare; crf_miraligner; crf_miraligner_novel }
collapsedRead_prefix.into { crp_prepare; crp_miraligner; crp_miraligner_novel }

process seqcluster_prepare{
	publishDir "${params.outdir}/seqcluster/prepare";

	input:
	//val(prefix) from crp_prepare.toList()
	val(prefix) from crp_prepare.toList().map { new nextflow.util.ArrayBag(it) } 

	file(reads) from crf_prepare.toList()

	output:
	file 'seqs.ma' into seqs_ma
	file 'seqs.fastq' into seqs_fastq
	file 'stats_prepare.tsv' into stats_prepare

	shell:
	def prefix = prefix.join(",")
	def reads = reads.toString().tokenize().join(",")
	"""
	python ${params.script_dir}/prepare_data.py ${prefix} ${reads}
	"""
}

seqs_ma.into { mirdeep2_seqs_ma; seqcluster_cluster_seqs_ma; }

process mapping{
	publishDir "${params.outdir}/align";

	cpus params.sambamba.cpus
	cpus params.star.cpus

	input:
	set file(read) from seqs_fastq

	output:
	file 'seqs.bam*' into star_bam

	"""
	STAR --genomeDir $params.db_star --readFilesIn $read --runThreadN ${params.star.cpus} --outFileNamePrefix ./ --outReadsUnmapped Fastx --outFilterMultimapNmax 1000 --outStd SAM --alignIntronMax 1 --outSAMunmapped Within --outSAMattributes NH HI NM MD AS  --sjdbGTFfile $params.gtf_ref_transcripts --sjdbOverhang 39  --outSAMattrRGline ID:miRQC_A PL:illumina PU:1_2016-04-08_mir2qc_bcbio SM:miRQC_A | /BiO/BioTools/bcbio/data/galaxy/../anaconda/bin/samtools sort -@ 1 -m 16G -o seqs.bam /dev/stdin
	sambamba index -t ${params.sambamba.cpus} seqs.bam
	"""

}

star_bam.into { star_bam_for_mirdeep2; star_bam_for_cluster; }

process seqcluster_cluster{
	publishDir "${params.outdir}/seqcluster/cluster";

	input:
	file '*' from seqcluster_cluster_seqs_ma
	file (bam:'*') from star_bam_for_cluster

	output:
	file 'seqcluster.json' into seqcluster_json

	"""
	seqcluster cluster -o ./ -m seqs.ma -a seqs.bam -r ${params.ref_fasta} -g ${params.gtf_srna_transcripts} 
	"""
}

process seqcluster_report{
	publishDir "${params.outdir}/seqcluster/report";

	input:
	file '*' from seqcluster_json

	output:
	file 'seqcluster.db' into seqcluster_report_db
	file 'html/' into seqcluster_report_html
	file 'log/' into seqcluster_report_log

	"""
	seqcluster report -o ./ -r ${params.ref_fasta} -j seqcluster.json
	"""
}

process mirdeep2{
	publishDir "${params.outdir}/mirdeep2";

	input:
	file('*') from mirdeep2_seqs_ma
	file('*') from star_bam_for_mirdeep2

	output:
	file 'novel' into mirdeep2_novel

	"""
	python ${params.script_dir}/mirdeep_prepare.py seqs.bam seqs.ma
	unset PERL5LIB && export PATH=/BiO/BioTools/bcbio/data/anaconda/bin:$PATH && perl /BiO/BioTools/bcbio/data/anaconda/bin/miRDeep2.pl file_reads.fa ${params.ref_fasta} align.bam ${params.mature_fa} none ${params.hairpin_fa} -f ${params.Rfam_for_miRDeep} -r simple -c -P -t ${params.species} -z res
	python ${params.script_dir}/mirdeep_parse_novel.py result_res.csv ${params.species}
	"""
}

/* http://seqcluster.readthedocs.io/mirna_annotation.html */
process miraligner{
	publishDir "${params.outdir}/mirbase/$prefix";

	input:
	val(prefix) from crp_miraligner;
	file(read) from crf_miraligner;

	output:
	val(prefix) into mirna_prefix
	file('*.mirna') into mirna_file

	"""
	export PATH=/BiO/BioPeople/brandon/jre1.7.0_51/bin:$PATH
	/BiO/BioTools/bcbio/data/anaconda/bin/miraligner -Xms705m -Xmx4500m -freq -sub 1 -trim 3 -add 3 -s ${params.species} -i ${read} -db ${params.db_srnaseq} -o ${prefix}
	"""
}

process miraligner_novel{
	publishDir "${params.outdir}/mirbase/$prefix";

	input:
	val(prefix) from crp_miraligner_novel;
	file(read) from crf_miraligner_novel;

	file(novel_db_path) from mirdeep2_novel.first();

	output:
	val(prefix) into mirna_novel_prefix
	file('*.mirna') into mirna_novel_file

	"""
	export PATH=/BiO/BioPeople/brandon/jre1.7.0_51/bin:$PATH
	/BiO/BioTools/bcbio/data/anaconda/bin/miraligner -Xms705m -Xmx4500m -freq -sub 1 -trim 3 -add 3 -s ${params.species} -i ${read} -db ${novel_db_path} -o ${prefix}_novel
	"""
}


mirna_file.into { mirna_for_qc; mirna_for_miraligner_parse }
mirna_novel_file.into { mirna_novel_for_qc; mirna_novel_for_miraligner_parse }

process qc_srna{
	publishDir "${params.outdir}/qc/$prefix/small-rna";

	input:
	file(known:'*.mirna') from mirna_for_qc
	file(novel:'*.mirna') from mirna_novel_for_qc

	"""
	python ${params.script_dir}/qc_srna.py ${known} ${novel}
	"""
}

process miraligner_parse{
	publishDir "${params.outdir}/";

	input:
	val(prefix) from mirna_prefix.toList()
	file(mirna) from mirna_for_miraligner_parse.toList()
	file(mirna_novel) from mirna_novel_for_miraligner_parse.toList()

	output:
	file('mirbase/') into known_mirna_result
	file('mirdeep2/') into novel_mirna_result

	exec:
	prefix = prefix.join(",")
	mirna = mirna.toString().tokenize().join(",")
	mirna_novel = mirna_novel.toString().tokenize().join(",")

	shell:	
	"""
	mkdir mirbase
	mkdir mirdeep2
	python ${params.script_dir}/miraligner_parser.py ${prefix} ${mirna}
	"""
}
