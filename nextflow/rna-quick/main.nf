#!/usr/bin/env nextflow

/* Try to create a URL; return false if url is not a value URL. */
def parseUrl(url) {
    try {
        new URL(url)
    } 
    catch (MalformedURLException e) {
        false
    }
}

/* Get a URL for a genome shortcut, or throw an exception if the shortcut is invalid. */
def parseGenomeShortcut(shortcut) {
    parts = shortcut.split('.')
    if (parts[0] == "ensembl") {
        release = int(parts[1])
        if (parts[2][0:3] == "GRCh") {
            organism = "homo_sapiens"
            build = parts[2]
        }
        else {
            organism = parts[2]
            build = parts[3]
        }
        URL('ftp://ftp.ensembl.org/pub/release-${release}/fasta/${organism}/dna/${organism}.${build}.dna.toplevel.fa.gz')
    }
    else if (parts[0] == "ucsc") {
        build = parts[1]
        URL('http://hgdownload.cse.ucsc.edu/goldenPath/${build}/bigZips/chromFa.tar.gz')
    }
    else {
        throw Exception('Invalid genome shortcut: $shortcut')
    }
}

/* Get a URL for an annotation shortcut, or throw an exception if the shortcut is invalid. */
def parseAnnotationShortcut(shortcut) {
    parts = shortcut.split('.')
    if (parts[0] == "gencode") {
        release = parts[1]
        content = "annotation" if parts[2] == "comprehensive" else "long_noncoding_RNAs"
        regions = "" if parts[3] == "all" else "chr_patch_hapl_scaff"
        URL('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${release}/gencode.v${release}.${regions}${content}.gtf.gz')
    }
    else {
        throw Exception('Invalid genome shortcut: $shortcut')
    }
}

/* Download the contents of a URL to a file. */
def download(url, outfile) {
    java.nio.ReadableByteChannel input = java.nio.Channels.newChannel(url.openStream())
    java.nio.FileOutputStream output = new java.nio.FileOutputStream(outfile)
    output.getChannel().transferFrom(input, 0, Long.MAX_VALUE)
}

def getFiles(lib, dir, pattern, condition=null) {
    PairIndex = 1
    fq1 = file(path.join(dir, evaluate(pattern)))
    PairIndex = 2
    fq2 = file(path.join(dir, evaluate(pattern)))
    if (condition) {
        [lib[condition], fq1, fq2]
    }
    else {
        [lib["UniqueID"], fq1, fq2]
    }
}

/* Channel to receive library metadata. */
libraries = Channel.create()

/* Download genome and annotation files, if necessary. */
process init {
    output:
    file 'genome.fa.gz' into GENOME
    file 'annotation.gtf' into ANNOTATION
    file 'libraries.txt' into LIBRARIES
    
    exec:
    
    genomeFile = file(params.genome)
    if (genomeFile.exists) {
        ['ln', '-s', genomeFile.absolutePath, 'genome.fa.gz'].execute().waitFor() 
    }
    else {
        url = parseUrl(params.genome)
        if (!url) {
            url = parseGenomeShortcut(params.genome)
        }
        download(url, 'genome.fa.gz')
    }
    
    annotationFile = file(params.annotation)
    if (annotationFile.exists) {
        ['ln', '-s', annotationFile.absolutePath, 'annotation.gtf'].execute().waitFor()
    }
    else {
        url = parseUrl(params.annotation)
        if (!url) {
            url = parseAnnotationShortcut(params.annotation)
        }
        download(url, 'annotation.gtf')
    }
    
    librariesFile = file('libraries.txt')
    if (!librariesFile.exists) {
        ['ln', '-s', file(params.librariesFile).absolutePath, 'libraries.txt'].execute().waitFor()
    }
    String[] lines = librariesFile.text.split('\n')
    String[] header = lines[0].split('\t')
    for (i = 1; i < lines.size; i++) {
        String[] lib = lines[i].split('\t')
        libraries << [header, lib].transpose().collectEntries { it }
    }
}

process full_transcriptome {
    container 'jdidion/rna-quick/containers/bedtools'
    
    input:
    file genome from GENOME
    file annotation from ANNOTATION
    
    output:
    file "full.fa.gz" into FULL_TRANSCRIPTOME
    
    """
    bedtools getfasta -name -s -fi $genome -bed $annotation -fo - | gzip > full.fa.gz
    """
}

process full_index {
    container 'jdidion/rna-quick/containers/kallisto'

    input:
    file fasta from FULL_TRANSCRIPTOME

    output:
    file "full.idx" into FULL_INDEX

    """
    kallisto index -i full.idx $fasta
    """
}

process full_quantify {
    container 'jdidion/rna-quick/containers/kallisto'
    cpus 4
    memory 8
    
    input:
    file full_index from FULL_INDEX
    set group, fastq1, fastq2 from libraries.map { 
        getFiles(it, params.library_type, params.fastq_dir, params.fastq_pattern, params.condition)
    }
    
    output:
    set val(group), file("./out/abundance.hd5") into FULL_QUANT
    
    """
    kallisto quant --index=$full_index --output=./out -t ${cpus} $fastq1 $fastq2
    """
}

process find_low_abundance {
    container 'rocker-org/rocker/r-base'
    
    input:
    file abundanceFile from FULL_QUANT.collectFile(name: 'full_abundances.txt', newLine: true)
    
    output:
    file excludeFile name 'exclude.txt' into EXCLUDE
    
    template 'find_low_abundance.R'
}

process filter_low_abundance {
    input:
    file exclude from EXCLUDE
    
    output:
    file filteredAnnotation 'filtered_annotation.gtf' into FILTERED_ANNOTATION
    
    """
    grep -F -f $exclude -v $params.annotation > $filteredAnnotation
    """
}

process filtered_transcriptome {
    container 'jdidion/rna-quick/containers/bedtools'
    
    input:
    file genome from GENOME
    file annotation from FILTERED_ANNOTATION
    
    output:
    file "filtered.fa.gz" into FILTERED_TRANSCRIPTOME
    
    """
    bedtools getfasta -name -s -fi $genome -bed $annotation -fo - | gzip > filtered.fa.gz
    """
}

process filtered_index {
    container 'jdidion/rna-quick/containers/kallisto'

    input:
    file fasta from FILTERED_TRANSCRIPTOME

    output:
    file "filtered.idx" into FILTERED_INDEX

    """
    kallisto index -i filtered.idx $fasta
    """
}

process filtered_quantify {
    container 'jdidion/rna-quick/containers/kallisto'
    cpus 4
    memory 8
    
    input:
    file filtered_index from FILTERED_INDEX
    set id, fastq1, fastq2 from libraries.map { 
        getFiles(it, params.library_type, params.fastq_dir, params.fastq_pattern) 
    }
    
    output:
    set val(id), file("./out/abundance.hd5") into FILTERED_QUANT
    
    """
    kallisto quant --index=$full_index --output=./out -t ${cpus} -b 100 $fastq1 $fastq2
    """
}

process export_sleuth_obj {
    input:
    file abundances from FILTERED_QUANT.collectFile(name: 'filtered_abundances.txt', newLine: true)
    file libraries from LIBRARIES
    file tx2gene from params.tx2geneFile
    val fullModel name param.models.full
    
    output:
    file sleuthObj name params.sleuthDataFile into SLEUTH_OBJ

    template: 'export_sleuth.R'
}

process transcript_differential_expression {
    input:
    file sleuthObj from SLEUTH_OBJ
    each set val(modelName), val(modelFormula), val(modelBetas) from params.models
    
    output:
    file resultFile name params.txDiffExpDataFile
    
    template: 'tx_differential_expression.R'
}

process gene_differential_expression {
    cpus 8
    
    input:
    file sleuthObj from SLEUTH_OBJ
    each set val(modelName), val(modelFormula), val(modelBetas) from params.models
    
    output:
    file resultFile name params.geneDiffExpDataFile
    
    template: 'gene_differential_expression.R'    
}
