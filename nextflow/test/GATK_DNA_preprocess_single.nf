#!/usr/bin/env nextflow
//https://raw.githubusercontent.com/joshua-d-campbell/DNA_Seq/b957bd0326cc115a2a6623f1a59ad94392699f15/Preprocess/GATK_DNA_preprocess_single.nf


//############################################################################################################################
//
// Josh Campbell
// 3/30/2016
// 
//############################################################################################################################


//############################################################################################################################
//
// GATK tutorials this pipeline is based on:
// Mapping: http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
// Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
// Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
// Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
// 
//############################################################################################################################

// List of parameters that can be passed to this workflow
params.targets = "/restricted/projectnb/cbmhive/references/Exome_Targets/SureSelect_Human_All_Exon_V4_plus_UTRs.UCSC_hg19.targets.interval_list"
params.baits = "/restricted/projectnb/cbmhive/references/Exome_Targets/SureSelect_XT_V4_plusUTRs_baits.UCSC_hg19.capture.interval_list"
params.ref = "/restricted/projectnb/cbmhive/references/GATK/reference_files/ucsc.hg19.fasta"
params.gatk_jar = "/share/pkg/gatk/3.5/install/GenomeAnalysisTK.jar"
params.picard_jar = "/share/pkg/picard/2.1.1/install/picard-tools-2.1.1/picard.jar"
params.gold_indels1 = "/restricted/projectnb/cbmhive/references/GATK/reference_files/1000G_phase1.indels.hg19.sites.vcf"
params.gold_indels2 = "/restricted/projectnb/cbmhive/references/GATK/reference_files/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
params.dbsnp = "/restricted/projectnb/cbmhive/references/GATK/reference_files/dbsnp_138.hg19.vcf"
params.output_dir = "./"

// Set up global variables for requried parameters:
PROJECT = params.project
inputFile = file(params.infile)

// Set up global variables for parameters with preset defaults:
REF = file(params.ref)
GATK = file(params.gatk_jar)
PICARD = file(params.picard_jar)
GOLD1 = file(params.gold_indels1)
GOLD2 = file(params.gold_indels2)
TARGETS = file(params.targets)
BAITS = file(params.baits)
OUTDIR = file(params.output_dir)
DBSNP = file(params.dbsnp)


// Necessary columns needed in input file:
// INDIVIDUAl_ID	SAMPLE_ID	LIBRARY_ID	RG_ID	PLATFORM_UNIT	PLATFORM	RUN_DATE	CENTER	R1	R2

//#############################################################################################################
//
// Read input file and save it into list of lists
//
//#############################################################################################################

file_params = []
is_header = 1
for( line in inputFile.readLines() ) {
  if(is_header == 1) {
    header = line.split("\t").flatten()
    is_header = 0
  } else {
    file_params << line.split("\t").flatten()
  }  
}


// Send FASTQ files to two processes: FastQC and FastqToSam
(readPairsFastQC, readPairsFastqToSam) = Channel.from(file_params).into(2)


//#############################################################################################################
//
// Preprocess reads
// 1) Convert to BAM
// 2) Mark Adapters
//
//#############################################################################################################

process runFastqToSam {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5g"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastqToSam/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, run_date, center, fastqR1, fastqR2 from readPairsFastqToSam
    
    output:
    set indivID, sampleID, libraryID, rgID, file(outfile) into runFastqToSamOutput

    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".unaligned.bam"
    
    """
    module load java/1.8.0_66
    
	java -Xmx5G -jar ${PICARD} FastqToSam \
		FASTQ=${fastqR1} \
		FASTQ2=${fastqR2} \
		OUTPUT=${outfile} \
		READ_GROUP_NAME=${rgID} \
		SAMPLE_NAME=${sampleID} \
		LIBRARY_NAME=${libraryID} \
		PLATFORM_UNIT=${platform_unit} \
		PLATFORM=${platform} \
		SEQUENCING_CENTER=${center} \
		RUN_DATE=${run_date}	
    """
    

}

process runMarkIlluminaAdapters {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5g"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/MarkIlluminaAdapters/"
    
    input:
    set indivID, sampleID, libraryID, rgID, ubam from runFastqToSamOutput
    
    output:
    set indivID, sampleID, libraryID, rgID, ubam, file(outfile_bam), file(outfile_metrics) into runMarkIlluminaAdaptersOutput
	
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".adapters_marked.bam"
    outfile_metrics = sampleID + "_" + libraryID + "_" + rgID + "_adapters_metrics.txt"
            
    """
    module load java/1.8.0_66
    
	java -Xmx5G -jar ${PICARD} MarkIlluminaAdapters \
		I=${ubam} \
		O=${outfile_bam} \
		M=${outfile_metrics} \
		TMP_DIR=tmp

    """
}




//#############################################################################################################
//
// Run BWA to align reads to genome
//
//#############################################################################################################

process runBWA {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=5g -pe omp 5"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/"
	
    input:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt, metrics from runMarkIlluminaAdaptersOutput
    
    output:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt into deleteBWAInput
    set indivID, sampleID, file(outfile_bam) into runBWAOutput
    
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"
	
    """
    module load java/1.8.0_66
    module load bwa/0.7.12
        
	set -o pipefail
	java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m -jar ${PICARD} SamToFastq \
		I=${ubamxt} \
		FASTQ=/dev/stdout \
		CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
		TMP_DIR=tmp | \
	bwa mem -M -t 5 -p ${REF} /dev/stdin | \
	java -Xmx5G -jar ${PICARD} MergeBamAlignment \
		ALIGNED_BAM=/dev/stdin \
		UNMAPPED_BAM=${ubamxt} \
		OUTPUT=${outfile_bam} \
		R=${REF} CREATE_INDEX=true ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=tmp
		
	"""	
}




//#############################################################################################################
//
// Combined libraries from the same Individual/Sample to send to MarkDuplicates
//
//#############################################################################################################

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])




//#############################################################################################################
//
// Run Picard MarkDuplicates
// This is used to merge different libraries from the same sample or the same library run on different lanes
// Requires a lot of memory
// Need to set "ParallelGCThreads" otherwise it will "grab" extra available threads without asking (and potentially be terminated by SGE)
//
//#############################################################################################################

process runMarkDuplicates {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=48:00:00 -l mem_total=30g -pe omp 5"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates"
	
    input:
    set indivID, sampleID, aligned_bam_list from runBWAOutput_grouped_by_sample
    
    output:
    set indivID, sampleID, aligned_bam_list into deleteMarkDuplicatesInput
    set indivID, sampleID, file(outfile_bam), file(outfile_metrics) into runMarkDuplicatesOutput
 
    script:
    outfile_bam = sampleID + ".dedup.bam"
    outfile_metrics = sampleID + "_duplicate_metrics.txt"	
	        
    """
    module load java/1.8.0_66
    
	java -Xmx30G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
		INPUT=${aligned_bam_list.join(" INPUT=")} \
		OUTPUT=${outfile_bam} \
		METRICS_FILE=${outfile_metrics} \
		CREATE_INDEX=true \
		TMP_DIR=/tmp
	"""  
}




//#############################################################################################################
//
// Perform realignment around indels
// 1) Identify regions for realignement
// 2) Perform realignment
//
//#############################################################################################################

process runRealignerTargetCreator {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=48:00:00 -l mem_total=15g"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/RealignerTargetCreator/"
    
    input:
    set indivID, sampleID, dedup_bam, metrics from runMarkDuplicatesOutput
    
    output:
    set indivID, sampleID, dedup_bam, file(target_file) into runRealignerTargetCreatorOutput
 	
    script:
    target_file = sampleID + "_target_intervals.list"
	        
    """
    module load java/1.8.0_66

	java -Xmx15g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T RealignerTargetCreator \
		-R ${REF} \
		-I ${dedup_bam} \
		-known ${GOLD1} \
		-known ${GOLD2} \
		-o ${target_file}
	"""  
}

process runIndelRealigner {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=25g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/IndelRealigner/"
	    
    input:
    set indivID, sampleID, dedup_bam, target_file from runRealignerTargetCreatorOutput
 	    
    output:
    set indivID, sampleID, file(realign_bam) into runIndelRealignerOutput
	set indivID, sampleID, dedup_bam into deleteRealignerTargetCreatorInput
 
    script:
    realign_bam = sampleID + ".realign.bam"
            
    """
    module load java/1.8.0_66

	java -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T IndelRealigner \
		-R ${REF} \
		-I ${dedup_bam} \
		-targetIntervals ${target_file} \
		-known ${GOLD1} \
		-known ${GOLD2} \
		-o ${realign_bam}
	"""  
}



//#############################################################################################################
//
// Perform base quality score recalibration (BQSR) including
// 1) Generate a recalibration table
// 2) Generate a new table after applying recalibration
// 3) Compare differences between recalibration tables
// 4) Apply recalibration
//
//#############################################################################################################

process runBaseRecalibrator {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=25g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibrator/"
	    
    input:
    set indivID, sampleID, realign_bam from runIndelRealignerOutput
    
    output:
    set indivID, sampleID, realign_bam, file(recal_table) into runBaseRecalibratorOutput_for_plotting
 	set indivID, sampleID, realign_bam, file(recal_table) into runBaseRecalibratorOutput_for_recal
    
    script:
    recal_table = sampleID + "_recal_table.txt" 
       
    """
    module load java/1.8.0_66
    
	java -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-o ${recal_table}

	"""
}


process runBaseRecalibratorPostRecal {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=48:00:00 -l mem_total=5g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibratorPostRecal/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput_for_plotting
    
    output:
    set indivID, sampleID, recal_table, file(post_recal_table) into runBaseRecalibratorPostRecalOutput
 
    script:
    post_recal_table = sampleID + "_post_recal_table.txt" 
       
    """
    module load java/1.8.0_66
    
	java -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-BQSR ${recal_table} \
		-o ${post_recal_table}

	"""
}	

process runAnalyzeCovariates {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=12:00:00 -l mem_total=5g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/AnalyzeCovariates/"
	    
    input:
    set indivID, sampleID, recal_table, post_recal_table from runBaseRecalibratorPostRecalOutput

	output:
	set indivID, sampleID, recal_plots into runAnalyzeCovariatesOutput
	    
    script:
    recal_plots = sampleID + "_recal_plots.pdf" 

    """
    module load java/1.8.0_66
    
	java -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T AnalyzeCovariates \
		-R ${REF} \
		-before ${recal_table} \
		-after ${post_recal_table} \
		-plots ${recal_plots}

    """
}    


process runPrintReads {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=25g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput_for_recal

    output:
    set indivID, sampleID, realign_bam into deleteIndelRealigner
    set indivID, sampleID, file(outfile_bam), file(outfile_bai) into runPrintReadsOutput_for_DepthOfCoverage, runPrintReadsOutput_for_HC_Metrics, runPrintReadsOutput_for_Multiple_Metrics
        
    script:
    outfile_bam = sampleID + ".clean.bam"
    outfile_bai = sampleID + ".clean.bai"
           
    """
    module load java/1.8.0_66

	java -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T PrintReads \
		-R ${REF} \
		-I ${realign_bam} \
		-BQSR ${recal_table} \
		-o ${outfile_bam}
    """
}    


//#############################################################################################################
//
// Perform a several tasks to assess QC:
// 1) Depth of coverage over targets
// 2) Generate alignment stats, insert size stats, quality score distribution
// 3) Generate hybrid capture stats
// 4) Run FASTQC to assess read quality
//
//#############################################################################################################

process runDepthOfCoverage {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=10g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/DepthOfCoverage"
	    
    input:
    set indivID, sampleID, bam, bai from runPrintReadsOutput_for_DepthOfCoverage

    output:
    file("${prefix}*") into DepthOfCoverageOutput
    
    script:
    prefix = sampleID + "."
         
    """
    module load java/1.8.0_66

	java -Djava.io.tmpdir=tmp/ -Xmx10g -jar ${GATK} \
		-R ${REF} \
		-T DepthOfCoverage \
		-I ${bam} \
		--omitDepthOutputAtEachBase \
		-L ${TARGETS} \
		-ct 10 -ct 20 -ct 50 -ct 100 \
		-o ${sampleID}

	"""
}	



process runCollectMultipleMetrics {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=10g"
 	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
 	    
    input:
    set indivID, sampleID, bam, bai from runPrintReadsOutput_for_Multiple_Metrics

    output:
    file("${prefix}*") into CollectMultipleMetricsOutput

    script:       
    prefix = sampleID + "."

    """
    module load java/1.8.0_66

	java -Xmx10g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${REF} \
		ASSUME_SORTED=true \
		OUTPUT=${prefix}
		
	"""
}	


process runHybridCaptureMetrics {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=10g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
	    
    input:
    set indivID, sampleID, bam, bai from runPrintReadsOutput_for_HC_Metrics

	output:
	file(outfile) into HybridCaptureMetrics

    script:       
    outfile = sampleID + "_hybrid_selection_metrics.txt"
    
    """
    module load java/1.8.0_66

	java -Xmx10g -Djava.io.tmpdir=tmp/ -jar $PICARD CalculateHsMetrics \
		INPUT=${bam} \
		OUTPUT=${outfile} \
		TARGET_INTERVALS=${TARGETS} \
		BAIT_INTERVALS=${BAITS} \
		REFERENCE_SEQUENCE=${REF} 	
	"""
}	


process runFastQC {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5g"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC/"
	    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, run_date, center, fastqR1, fastqR2 from readPairsFastQC

    output:
    set file("*.zip"), file("*.html") into FastQCOutput
    	
    script:

    """
    module load fastqc/0.11.3
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
    """
}




//#############################################################################################################
//
// "Garbage Collection" processes designed to delete large BAM files once they are no longer needed to save on space
//
//#############################################################################################################

process runDeleteBWAInput {
    tag "${indivID}|${sampleID}|${rgID}|${libraryID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00 -l mem_total=1g"
    
    input:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt from deleteBWAInput
    
    script:
    """
    rm -v ${ubam} ${ubamxt}
	"""
}

process runDeleteMarkDuplicatesInput {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00 -l mem_total=1g"
    
    input:
    set indivID, sampleID, aligned_bam_list from deleteMarkDuplicatesInput
    
    script:
    cl = aligned_bam_list.getClass()
    cl_n = cl.getName()
    if(cl_n == "sun.nio.fs.UnixPath") {
      aligned_bai_list = aligned_bam_list.toString().replaceAll(".bam", ".bai")
    }
    else {
      aligned_bai_list = aligned_bam_list.join(" ").replaceAll(".bam", ".bai")
      aligned_bam_list = aligned_bam_list.join(" ")      
    }
    
    """
    rm -v ${aligned_bam_list} 
    rm -v ${aligned_bai_list} 
    """
}

process runDeleteRealignerTargetCreator {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00 -l mem_total=1g"
    
    input:
    set indivID, sampleID, dedup_bam_list from deleteRealignerTargetCreatorInput
    
    script:
    cl = dedup_bam_list.getClass()
    cl_n = cl.getName()
    if(cl_n == "sun.nio.fs.UnixPath") {
      dedup_bai_list = dedup_bam_list.toString().replaceAll(".bam", ".bai")
    }
    else {
      dedup_bai_list = dedup_bam_list.join(" ").replaceAll(".bam", ".bai")
      dedup_bam_list = dedup_bam_list.join(" ")
    }
    
    """
    rm -v ${dedup_bam_list} 
    rm -v ${dedup_bai_list} 
    """
}

process runDeleteIndelRealigner {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00 -l mem_total=1g"
    
    // Input required from two different branches to ensure that the bam is not
    // deleted until they both are completely finished. The output from
    // AnalyzeCovariates is not used for anything 
    input:
    set indivID, sampleID, realign_bam_list from deleteIndelRealigner
	set indivID_2, sampleID_2, recal_table, post_recal_table from runAnalyzeCovariatesOutput
    
    script:
    cl = realign_bam_list.getClass()
    cl_n = cl.getName()
    if(cl_n == "sun.nio.fs.UnixPath") {
      realign_bai_list = realign_bam_list.toString().replaceAll(".bam", ".bai")
    }
    else {
      realign_bai_list = realign_bam_list.join(" ").replaceAll(".bam", ".bai")
      realign_bam_list = realign_bam_list.join(" ")      
    }
    
    """
    rm -v ${realign_bam_list} 
    rm -v ${realign_bai_list} 
    """
}
