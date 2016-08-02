/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'lncRNA-Annotation-NF'.
 *
 *   lncRNA-Annotation-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   lncRNA-Annotation-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with lncRNA-Annotation-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main lncRNA-Annotation-NF pipeline script
 *
 * @authors
 * Sarah Djebali <sarah.djebali-quelen@toulouse.inra.fr>
 * Thomas Derrien <thomas.derrien@univ-rennes1.fr>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com>
 *
 * @modifier
 * Brandon <hmkim@live.co.kr> 
 */

params.name          ="lncRNA_Pig_RNA-Seq"
params.genome        ="$baseDir/tutorial/genome/genome_chr38.fa"
params.annotation    ="$baseDir/tutorial/annotation/annotation_chr38.gtf"
params.reads         ="$baseDir/tutorial/reads/*_{1,2}.fastq"
params.overhang      ='99'
params.feelnc_opts   = "--biotype transcript_biotype=protein_coding --monoex -1 "
params.output        ="results/"


log.info "lncRNA Annotation - N F  ~  version 0.1"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "genome                 : ${params.genome}"
log.info "reads                  : ${params.reads}"
log.info "annotation             : ${params.annotation}"
log.info "STAR overhang          : ${params.overhang}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

genomeFile             = file(params.genome)
annotationFile         = file(params.annotation) 

FEELnc_filter_options=params.feelnc_opts

/*
 * validate input files/
 */

if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"
if( !genomeFile.exists() ) exit 1, "Missing genome directory: ${genomeFile}"


/*
 * Create a channel for read files 
 */
 
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path -> 
       def prefix = readPrefix(path, params.reads)
       tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read_files } 


process index {
    input:
    file genomeFile
    file annotationFile
    
    output:
    file "STARgenome" into STARgenomeIndex
    file "genomeLength.txt" into genomeLengths    

    script:
    //
    // STAR Generate Index and create genome length file
    //
    """
        mkdir STARgenome
        STAR --runThreadN ${task.cpus} \
             --runMode genomeGenerate \
             --genomeDir STARgenome \
             --genomeFastaFiles ${genomeFile} \
             --sjdbGTFfile ${annotationFile} \
             --sjdbOverhang ${params.overhang} \
             --outFileNamePrefix STARgenome

        fastalength ${genomeFile} > 'genomeLength.txt'

    """
}

process mapping {
    tag "reads: $name"

    input:
    file STARgenome from STARgenomeIndex.first()
    set val(name), file(reads:'*') from read_files

    output:
    set val(name), file("STAR_${name}") into STARmappedReads 

    script:
    //
    // STAR Mapper
    //
    """
        STAR --genomeDir ${STARgenome} \
             --readFilesIn ${reads} \
             --readFilesCommand zcat \
             --outFilterType BySJout \
             --outSAMunmapped Within \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrIHstart 0 \
             --outFilterIntronMotifs RemoveNoncanonical \
             --runThreadN ${task.cpus} \
             --quantMode TranscriptomeSAM \
             --outWigType bedGraph \
             --outWigStrand Stranded \
             --outFileNamePrefix ${name}
        
        mkdir STAR_${name}
        mv ${name}Aligned* STAR_${name}/.
        mv ${name}Signal* STAR_${name}/.
        mv ${name}SJ* STAR_${name}/.
        mv ${name}Log* STAR_${name}/.
    """
   
}


process cufflinks {

    input:
    file annotationFile
    set val(name), file(STAR_alignment) from STARmappedReads

    output:
    set val(name), file("CUFF_${name}") into cufflinksTranscripts_to_pp

    script:
    //
    // Cufflinks
    //
    """
        mkdir CUFF_${name}
        cufflinks -p ${task.cpus} \
                  -g ${annotationFile} \
                  -o CUFF_${name} \
                  --overlap-radius 5 \
                  --library-type fr-firststrand \
                  ${STAR_alignment}/${name}Aligned.sortedByCoord.out.bam
    """
}

process cufflinks_postprocess {

    input:
    set val(name), file(CUFFLINKS_transcripts) from cufflinksTranscripts_to_pp
    file (genomeLength) from genomeLengths.first()

    output:
    file "${name}_cufflinks_ok.gtf" into cufflinksTranscripts_postprocess, cufflinksTranscripts_postprocess_fn

    script:
    //
    // Post Process Cufflinks to remove exons exceeding length of chromosome/scaffold
    //
    """ 
        cuff=${CUFFLINKS_transcripts}/transcripts.gtf

        awk -v fileRef=${genomeLength} 'BEGIN{while (getline < fileRef >0){lg[\$2]=\$1}} \
            {nbex[\$12]++; line[\$12,nbex[\$12]]=\$0}\$4<1||\$5>lg[\$1]{ko[\$12]=1}END\
            {for(t in nbex){if(ko[t]!=1){for(k=1; k<=nbex[t]; k++){print line[t,k]}}}}'\
            \$cuff | awk -f /BiO/BioPeople/brandon/workflow/nextflow/lncRNA-annotation/bin/gff2gff.awk > ${name}_cufflinks_ok.gtf
    """
}


//
// Create a file 'gtf_filenames' containing the filenames of each post processes cufflinks gtf
//

cufflinksTranscripts_postprocess_fn
  .collectFile () { file ->  ['gtf_filenames.txt', file.name + '\n' ] }
  .set { GTFfilenames }


process cuffmerge {

    input:
    file annotationFile
    file genomeFile
    file gtf_filenames from GTFfilenames
    file cufflinks_ok from cufflinksTranscripts_postprocess.toList()

    output:
    file "CUFFMERGE" into cuffmergeTranscripts

    script:
    //
    // Cuffmerge
    //
    """
        mkdir CUFFMERGE
        cuffmerge -o CUFFMERGE \
                  -g ${annotationFile}  \
                  -s ${genomeFile} \
                  -p ${task.cpus} \
                     ${gtf_filenames}
    """
}

process FEELnc_filter{

    input:
    file annotationFile
    file cuffmergeDir from cuffmergeTranscripts

    output:
    file "FEELnc_filter" into FEELnc_filtered

    script:
    //
    // FEELnc Filter Step
    //
 

    """
        mkdir FEELnc_filter
        perl -I/BiO/BioTools/miniconda3/lib/perl5/site_perl/5.22.0 -I/BiO/BioTools/FEELnc/lib -I/home/brandon/perl5/lib/perl5 /BiO/BioTools/FEELnc/scripts/FEELnc_filter.pl --infile ${cuffmergeDir}/merged.gtf \
                         --mRNAfile ${annotationFile} \
                         --proc ${task.cpus} \
                         ${params.feelnc_opts} \
                         > FEELnc_filter/merged_filtered.gtf
    """

}

process FEELnc_codpot{

    input:
    file annotationFile
    file genomeFile
    file FEELnc_filter from FEELnc_filtered

    output:
    file "intergenic_0.99_Mode" into FEELnc_codpot_out

    script:
    //
    // FEELnc Coding Potential
    //


    """
	export FEELNCPATH=/BiO/BioTools/FEELnc
        mkdir intergenic_0.99_Mode
        perl -I/BiO/BioTools/miniconda3/lib/perl5/site_perl/5.22.0 -I/BiO/BioTools/FEELnc/lib -I/home/brandon/perl5/lib/perl5 /BiO/BioTools/FEELnc/scripts/FEELnc_codpot.pl --infile=${FEELnc_filter}/merged_filtered.gtf \
                          --mRNAfile=${annotationFile} \
                          --genome=${genomeFile} \
                          --numtx=10000,10000 \
                          --outdir=intergenic_0.99_Mode \
                          --spethres=0.99,0.99 \
                          --proc=${task.cpus} \
                          --keeptmp \
			  --mode intergenic
                          
    """

}

process FEELnc_classifier {

    input:
    file annotationFile
    file intergenic from FEELnc_codpot_out

    output:
    file "lncRNA_classes.txt" into lncRNA_classes

    script:
    //
    // FEELnc Coding Potential
    //

    """
    	/home/brandon/perl5/perlbrew/perls/perl-5.20.1/bin/perl -I/BiO/BioTools/miniconda3/lib/perl5/site_perl/5.22.0 -I/BiO/BioTools/FEELnc/lib -I/home/brandon/perl5/lib/perl5 /BiO/BioTools/FEELnc/scripts/FEELnc_classifier.pl  -i ${intergenic}/merged_filtered.gtf.lncRNA.gtf  \
                              -a ${annotationFile} \
                              > lncRNA_classes.txt
    """

}


// ===================== UTILITY FUNCTIONS ============================


/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */
 
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


