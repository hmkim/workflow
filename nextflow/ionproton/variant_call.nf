#! /usr/bin/env nextflow
//ref: https://github.com/IARCbioinfo/needlestack/blob/master/needlestack.nf

params.gatk_cpu = 8

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------'
    log.info 'Variant call v0.1: A MULTI-SAMPLE VARIANT CALLER using GATK'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run variant_call.nf --bed bedfile.bed --bam_folder BAM/ --fasta_ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info 'Options:'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, by default input bam folder.'
    log.info '    --bed            BED FILE                 A BED file for calling.'
    log.info '    --region         CHR:START-END            A region for calling.'
    log.info ''
    exit 1
}

assert (params.fasta_ref != true) && (params.fasta_ref != null) : "please specify --fasta_ref option (--fasta_ref reference.fasta(.gz))"
assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"

fasta_ref = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_dict = file( params.fasta_ref.take(params.fasta_ref.lastIndexOf('.')) +'.dict' )
fasta_ref_gzi = file( params.fasta_ref+'.gzi' )

/* Verify user inputs are correct */

if (params.bed) { try { assert file(params.bed).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() } }
try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
if (fasta_ref_dict.exists()) {assert fasta_ref_dict.exists() : "input fasta reference does not seem to have a .dict index (use picard CreateDictionary)"}
if (fasta_ref.exists() && params.fasta_ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}
try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
assert file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "BAM folder contains no BAM"
if (file(params.bam_folder).exists()) {
  bamID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.collect { it.getName() }.collect { it.replace('.bam','') }
  baiID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam.bai/ }.collect { it.getName() }.collect { it.replace('.bam.bai','') }
  assert baiID.containsAll(bamID) : "check that every bam file has an index (.bam.bai)"
}


sample_names = params.use_file_name ? "FILE" : "BAM"
out_vcf = params.out_vcf ? params.out_vcf : "all_variants.vcf"

/* manage input positions to call (bed or region or whole-genome) */

intervals_gvcf = params.bed ? '-L '+params.bed : ""

if(params.region){
    input_region = 'region'
  } else if (params.bed){
    input_region = 'bed'
    bed = file(params.bed)
  } else {
    input_region = 'whole_genome'
  }

/* Software information */

log.info ''
log.info '-------------------------------------------------------'
log.info 'Variant call v0.1: A MULTI-SAMPLE VARIANT CALLER using GATK'
log.info '-------------------------------------------------------'
log.info "Input BAM folder (--bam_folder)                                 : ${params.bam_folder}"
log.info "Reference in fasta format (--fasta_ref)                         : ${params.fasta_ref}"
log.info "Intervals for calling (--bed)                                   : ${input_region}"
log.info "output folder (--out_folder)                                    : ${params.out_folder}"
log.info "\n"

bam = Channel.fromPath( params.bam_folder+'/*.bam' )
	.map { file -> tuple(file.baseName, file, file + '.bai') }


process gatk_haplotyper_gvcf {
	cpus "${params.gatk_cpu}"

	tag { prefix }

	publishDir  params.out_folder+'/gVCF/', mode: 'move'
 
	input:
	set val(prefix), file(bam), file(bai) from bam
	
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_gzi
	file fasta_ref_dict

	output:
	file("${bam_tag}_raw_calls.g.vcf") into output_gvcf
	file("${bam_tag}_raw_calls.g.vcf.idx") into output_gvcf_idx

	"""
	gatk -Xmx16g -Djava.io.tmpdir=./ -T HaplotypeCaller -nct ${task.cpus} -R ${fasta_ref} -I ${bam} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF ${intervals_gvcf} -o ${prefix}.g.vcf
	"""
}

