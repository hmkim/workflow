Channel
        .fromPath( params.in_vcf )
        .ifEmpty { error "Cannot find any reads matching: ${params.in_vcf}" }
        .map { path ->
                def prefix = readPrefix(path, params.in_vcf)
                tuple(prefix, path)
        }
        .groupTuple(sort: true)
        .set { vcf_files }


vcf_files.into { test1; test2 }
//test2.println()

process run_snpeff {
	publishDir "${params.out_dir}"

	memory { 8.GB * task.attempt } 
	errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
        maxRetries 3
        maxErrors '-1'

	input:
	set val(prefix), file(vcf:'*') from test1

	output:
	file '*snpEff_genes.txt' into snpEff_genes_txt
	file '*snpEff_summary.html' into snpEff_summary_html
	set val(prefix), file("${prefix}.eff.vcf") into snpEff_vcf

	"""
	snpEff -Xmx8g eff ${params.genome_version} ${vcf} > ${prefix}.eff.vcf
	mv snpEff_genes.txt ${prefix}.snpEff_genes.txt
	mv snpEff_summary.html ${prefix}.snpEff_summary.html
	"""
}

process extractField{
	publishDir "${params.out_dir}"

	input:
	set val(prefix), file(vcf:'*') from snpEff_vcf

	output:
	set val(prefix), file("${prefix}.eff.xls") into snpEff_xls

	"""
	java -Xmx16g -jar /BiO/BioTools/miniconda3/envs/vcf_annotation/share/snpeff-4.3-2/SnpSift.jar extractFields -s , -e . ${vcf} CHROM POS ID REF ALT QUAL "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].FEATURE" "GEN[*].GT[*]" "GEN[*].AD[*]"  "GEN[*].DP[*]" "GEN[*].GQ[*]" "GEN[*].PL[*]" "LOF" "NMD" > ${prefix}.eff.xls
	"""
}

def readPrefix( Path actual, template ) {

        //final fileName = actual.getFileName().toString()
        final fileName = actual.name.lastIndexOf('.').with {it != -1 ? actual.name[0..<it] : actual.name} 

        return fileName
}
