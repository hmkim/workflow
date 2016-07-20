
params.control = "/BiO/BioProjects/NCC-Human-WES-2016-06-TBD160355/result/12_gatk_printrecal/TN1605D1502/TN1605D1502.printrecal.ba*"
params.tumor = "/BiO/BioProjects/NCC-Human-WES-2016-06-TBD160355/result/12_gatk_printrecal/TN1605D1503/TN1605D1503.printrecal.ba*"

params.ref = "/BiO/BioResources/References/Human/hg19/hg19.fa"

params.out_dir = "/BiO/BioPeople/brandon/test_nextflow/sv/output";

params.delly = "/BiO/BioTools/delly/0.7.3/delly_v0.7.3_linux_x86_64bit"

delly = file(params.delly)

Channel
        .fromPath( params.control )
        .ifEmpty { error "Cannot find any reads matching: ${params.in}" }
        .map { path ->
                def prefix = readPrefix(path, params.control)
                tuple(prefix, path)
        }
        .groupTuple(sort: true)
        .set { control }

Channel
        .fromPath( params.tumor )
        .ifEmpty { error "Cannot find any reads matching: ${params.in}" }
        .map { path ->
                def prefix = readPrefix(path, params.control)
                tuple(prefix, path)
        }
        .groupTuple(sort: true)
        .set { tumor }


//control.println()
//tumor.println()

process run_delly {
	beforeScript 'export PATH=/BiO/BioTools/miniconda3/envs/sv/bin:$PATH'

	publishDir "${params.out_dir}"
 
	input:
	set val(prefix_control), file(control_bam:'*') from control
	set val(prefix_tumor), file(tumor_bam:'*') from tumor

	//DEL, DUP, INV, TRA, INS
	//export OMP_NUM_THREADS=3
	"""
	${delly} call -t DEL -o delly.bcf -g ${params.ref} ${prefix_tumor}.bam ${prefix_control}.bam 
	echo "${prefix_control}\tcontrol" >>  samples.tsv
	echo "${prefix_tumor}\ttumor" >>  samples.tsv
	${delly} filter -t DEL -f somatic -o delly.pre.bcf -s samples.tsv -g ${params.ref} delly.bcf
	bcftools view delly.pre.bcf > delly.pre.vcf
	/BiO/BioTools/svprops/src/svprops delly.pre.bcf > delly.pre.tsv
	cat delly.pre.tsv | tail -n +2 | cut -f 1,2,4,5 > delly.pre.bed
	"""
}

def readPrefix( Path actual, template ) {

        //final fileName = actual.getFileName().toString()
        final fileName = actual.name.lastIndexOf('.').with {it != -1 ? actual.name[0..<it] : actual.name} 

        return fileName
}
