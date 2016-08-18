#!/usr/bin/env nextflow

Channel
    .fromPath( params.in_vcf )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path ->
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
    }
    .groupTuple(sort: true)
    .set { vcf_files }

vcf_files.into { vcf_for_snpeff; vcf_for_gwasCat; vcf_for_dbsnp; vcf_for_cosmic; vcf_for_varType; vcf_for_dbNSFP; vcf_for_clinvar; vcf_for_ExAC }

process snpeff {
	publishDir "${params.out_dir}/snpeff"

	memory { 8.GB * task.attempt } 
	errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
        maxRetries 3
        maxErrors '-1'

	input:
	set val(prefix), file(vcf:'*') from vcf_for_snpeff

	output:
	set val(prefix), file('*.eff.vcf') into vcf_snpeff 

	"""
	snpEff -Xmx8g eff ${params.genome_version} ${vcf} > ${prefix}.eff.vcf
	"""
}

process dbNSFP {
	publishDir "${params.out_dir}/dbNSFP"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_dbNSFP

	output:
	set val(prefix), file('*.dbNSFP.vcf') into vcf_dbNSFP

	"""
	java -Xmx20g -jar ${params.snpsift} dbnsfp -db ${params.dbNSFP} ${vcf} > ${prefix}.dbNSFP.vcf
	"""
}

process varType {
	publishDir "${params.out_dir}/varType"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_varType

	output:
	set val(prefix), file('*.varType.vcf') into vcf_varType

	"""
	java -Xmx8g -jar ${params.snpsift} varType ${vcf} > ${prefix}.varType.vcf
	"""
}

process gwasCat {
	publishDir "${params.out_dir}/gwasCat"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_gwasCat

	output:
	set val(prefix), file('*.gwasCat.vcf') into vcf_gwasCat

	"""
	java -Xmx8g -jar ${params.snpsift} gwasCat -db ${params.gwascatalog} ${vcf} > ${prefix}.gwasCat.vcf
	"""
}

process clinvar {
	publishDir "${params.out_dir}/clinvar"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_clinvar

	output:
	set val(prefix), file('*.clinvar.vcf') into vcf_clinvar

	"""
	java -Xmx8g -jar ${params.snpsift} annotate -noId -clinvar -db ${params.clinvar} ${vcf} > ${prefix}.clinvar.vcf
	"""
}

process dbsnp {
	publishDir "${params.out_dir}/dbsnp"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_dbsnp

	output:
	set val(prefix), file('*.rsnum.vcf') into vcf_dbsnp 

	"""
	bgzip -c ${vcf} > ${vcf}.gz
	tabix -p vcf ${vcf}.gz
	bcftools annotate --output ${prefix}.noids.vcf --remove ID ${vcf}.gz
	bgzip -c ${prefix}.noids.vcf > ${prefix}.noids.vcf.gz
	tabix -p vcf ${prefix}.noids.vcf.gz
	bcftools annotate --annotations ${params.dbsnp} --columns ID --output ${prefix}.rsnum.vcf ${prefix}.noids.vcf.gz
	"""
}

process ExAC{
	publishDir "${params.out_dir}/ExAC"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_ExAC

	output:
	set val(prefix), file('*.ExAC.vcf') into vcf_ExAC

	"""
	java -Xmx8g -jar ${params.snpsift} annotate -noId -db ${params.ExAC} ${vcf} > ${prefix}.ExAC.vcf
	"""
}

process cosmic{
	publishDir "${params.out_dir}/cosmic"

	input:
	set val(prefix), file(vcf:'*') from vcf_for_cosmic

	output:
	set val(prefix), file('*.cosmic.vcf') into vcf_cosmic

	"""
	java -Xmx16g -jar ${params.snpsift}  annotate -name 'COSM_' ${params.cosmic} ${vcf} > ${prefix}.cosmic.vcf
	"""
}


//process combine_dbsnp_withOther{
//	publishDir "${params.out_dir}/combine_vcf"
//
//	input:
//	set val(prefix1), file(vcf1:'*') from vcf_dbsnp
//	set val(prefix2), file(vcf2:'*') from vcf_cosmic
//	set val(prefix3), file(vcf3:'*') from vcf_snpeff
//	set val(prefix4), file(vcf4:'*') from vcf_dbNSFP
//	set val(prefix5), file(vcf5:'*') from vcf_varType
//	set val(prefix6), file(vcf6:'*') from vcf_gwasCat
//	set val(prefix7), file(vcf7:'*') from vcf_clinvar
//	set val(prefix8), file(vcf8:'*') from vcf_ExAC
//	
//	output:
//	set val(prefix1), file('*.vcf') into vcf_combine
//
//	"""
////	vcfaddinfo ${vcf1} ${vcf2} > t1.vcf
////	vcfaddinfo t1.vcf ${vcf3} > t2.vcf
////	vcfaddinfo t2.vcf ${vcf4} > t3.vcf
////	vcfaddinfo t3.vcf ${vcf5} > t4.vcf
////	vcfaddinfo t4.vcf ${vcf6} > t5.vcf
////	vcfaddinfo t5.vcf ${vcf7} > t6.vcf
////	vcfaddinfo t6.vcf ${vcf8} > ${prefix1}.combine.vcf
//	"""
//}

//process extractField{
//	publishDir "${params.out_dir}"
//
//	input:
//	set val(prefix), file(vcf:'*') from snpEff_vcf
//
//	output:
//	set val(prefix), file("${prefix}.eff.xls") into snpEff_xls
//
//	"""
//	java -Xmx16g -jar /BiO/BioTools/miniconda3/envs/vcf_annotation/share/snpeff-4.3-2/SnpSift.jar extractFields -s , -e . ${vcf} CHROM POS ID REF ALT QUAL "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].FEATURE" "GEN[*].GT[*]" "GEN[*].AD[*]"  "GEN[*].DP[*]" "GEN[*].GQ[*]" "GEN[*].PL[*]" "LOF" "NMD" > ${prefix}.eff.xls
//	"""
//}
//
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
