# Practical 1

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs. Please check out the [Useful information](#useful-information) section.

## Required tools

* Annovar
  * http://www.openbioinformatics.org/annovar
  * http://www.openbioinformatics.org/annovar/annovar_download.html
  * databases: (ljb23_all, 1000g2012apr, snp138)
* Freebayes - https://github.com/ekg/freebayes
* GATK - https://www.broadinstitute.org/gatk/download
* IGV - https://www.broadinstitute.org/igv/home
* Java (7) - http://www.oracle.com/technetwork/java/javase/downloads/index.html?ssSourceSiteId=otnjp
* Picard - http://picard.sourceforge.net/command-line-overview.shtml
* Qualimap - http://qualimap.bioinfo.cipf.es/
* SAMtools - http://samtools.sourceforge.net/‎
* Varscan 2 (2.3.6) - http://varscan.sourceforge.net/
* VCFtools - http://vcftools.sourceforge.net/
* VarDict - https://github.com/AstraZeneca-NGS/VarDictJava



## Information

* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/).



## Exercise

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a [target file](target.bed) to shorten variant calling time.

#### Important

* After each step inspect the generated output (cat, less, head, grep, ...).
* Organize your data and results in folders.
* Check out the specified parameters. If you don't know the exact meaning, consult the help pages (-h, --help).


#### System information
    cat /proc/cpuinfo
    cat /proc/meminfo

#### SAMtools


__(*)__ How big is the BAM file

    ls -lah <file.bam>
    

__(*)__ Inspect the header of the BAM file
    samtools ...
    samtools view -H aln.bam


__(*)__ View the BAM file

    module add samootls-1.3
    samtools view <bam.file> | less
    
__(*)__ How many reads are in the BAM file?<br/>
Is there another way to count the reads (check the samtools view parameters - look for -v)
   
    samtools view <file.bam> | grep -v "^#" | wc -l
    samtools flagstat <file.bam>
    
__(*)__ Answer the following questions by investigating the SAM file
* What version of the human assembly was used to perform the alignments?
* What version of bwa was used to align the reads?
* What is the name of the first read?
* At what position does the alignment of the read start?

    Use SAMtools for these questions
    samtools view -H aln.bam | less
    samtools view aln.bam | less

    
__(*)__ Sort the BAM file

    samtools sort -o sorted.bam -@ <THREADS> <file.bam>
    
__(*)__ Index the bam file
    
    samtools index <sorted.bam>


#### Alignment stats
    samtools flagstat sorted.bam
    samtools idxstats sorted.bam



#### Qualimap
__(*)__ Run qualimap on the command line
    
    module add qualimap_v2.2
    ./qualimap bamqc -bam sorted.bam -nt <numberOfThreads> -outdir bamqc


__(*)__ Inspect the BAM file in qualimap

    Check out the generated report. What qc features does it include?
    firefox ./bamqc/qualimapReport.html

    
__(*)__ Alternatively start the graphical interface and inspect the file using the GUI
    
    Change the heap size as described in http://qualimap.bioinfo.cipf.es/doc_html/faq.html#heapsize
    Open qualimap
    Load the BAM file (BAM QC) -> start analysis ---- Be sure to select the sorted BAM file
    Check out the different pages
    What does the "Coverage across Reference" tells you?
    
__(*)__ Hint:<br/>
If the BAM file is too big, try to view only a subset of it<br/>
(example with 1.000.000 rows [including header])

    samtools view -h sorted.bam | head -n 1000000 | samtools view -b -S - > sorted_small.bam
    
    
  
#### Prepare reference genome
__(*)__ Prepare dict index
    
    module add picard-tools-2.2.1
    java -jar /bcga2016/picard-tools-2.2.1/picard.jar CreateSequenceDictionary R=hg19.fasta O=hg19.dict

__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 


#### BAM file preparations
__(*)__ Sort with Picard
    
    java -Xmx8g -jar /bcga2016/picard-tools-2.2.1/picard.jar SortSam I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx8g -jar /bcga2016/picard-tools-2.2.1/picard.jar MarkDuplicates I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx8g -jar /bcga2016/picard-tools-2.2.1/picard.jar AddOrReplaceReadGroups I=dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx8g -jar /bcga2016/picard-tools-2.2.1/picard.jar BuildBamIndex I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT

__(*)__ Collect insert size metrics
    
    module add R-3.2.4
    java -Xmx8g -jar /bcga2016/picard-tools-2.2.1/picard.jar CollectInsertSizeMetrics I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf
    
__(*)__ View the PDF
    evince insertSizeHistogram.pdf


__(*)__ Questions
* How many reads were marked as duplicated? Look for flags.
* What are the other sorting possibilities for SortSam?
* Inspect the insert size metrics histogram.


#### SAMtools variant calling

__(*)__ Call

     module add bcftools-1.1
     samtools mpileup -uf hg19.fasta deduprg.bam | bcftools call -c -v -o samtools.vcf

__(*)__ Investigate result

    #How many variant were called
    grep -v "^#" samtools.vcf | wc -l
    #Print the variant that are between 1-300000 
    awk '!/^#/ && $2 < "300000"' samtools.vcf

#### FreeBayes variant calling

__(*)__ Call

     module add freebayes-1.0.2
     freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools
    #Do you notice differences?


#### VarDict variant calling

__(*)__ Call

     AF_THR="0.01" # minimum allele frequency
     /bcga2016/vardict/VarDictJava/build/install/VarDict/bin/VarDict -G hg19.fasta -f ${AF_THR} -N my_sample -b sorted.bam -z -c 1 -S 2 -E 3 -g 4 -R chr11:1-800000 | /bcga2016/vardict/VarDictJava/VarDict/teststrandbias.R | /bcga2016/vardict/VarDictJava/VarDict/var2vcf_valid.pl -N my_sample -E -f $AF_THR > vardict.vcf


#### Useful information


__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract xy from a file

    grep -i "xy" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100
