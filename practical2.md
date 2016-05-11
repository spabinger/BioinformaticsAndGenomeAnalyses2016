# Practical 2

In this practical you will get to know more tools for variant calling. Furthermore, VCF files will be annotated, filtered, and displayed. Please check out the [Useful information](#useful-information) section.



#### Varscan variant calling

__(*)__ Convert to mpileup

    samtools mpileup -B -f ${GEN_REF} deduprg.bam > <pileup.file>

__(*)__ Call SNPs

    java -jar <varscan.jar> mpileup2snp <pileup.file> --output-vcf 1 > varscan_snp.vcf

__(*)__ Call Indels
    
    java -jar <varscan.jar> mpileup2indel <pileup.file> --output-vcf 1 > varscan_indel.vcf

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools.
    #How many SNPs and indels were called?


#### GATK variant calling

__(*)__ Known indel sites are here specified as variables - either copy the whole path or use variables as well

    KNOWN_INDELS_1="1000G_phase1.indels.hg19.vcf"
    KNOWN_INDELS_2="Mills_and_1000G_gold_standard.indels.hg19.vcf"


__(*)__ Realignment target creator

    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19.fasta -nt 8 -L target.bed \
    -I deduprg.bam -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o target_intervals.list

__(*)__ Perform realignment
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fasta -I deduprg.bam \
    -targetIntervals target_intervals.list -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o dedup_rg_real.bam


__(*)__ Base quality recalibration
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam \
    -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o recal_data_table.txt -L target.bed 
    --maximum_cycle_value 800


__(*)__ Second pass of recalibration
     
     java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam \
     -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o post_recal_data_table.txt \
     -L target.bed --maximum_cycle_value 800 -BQSR recal_data_table.txt 


__(*)__ Generate before after plots (requires R and ggplot2)

    Fix missing R packages
    R
    Inside R call
    install.packages(c('reshape','gplots','gsalib'))
    
    java -Xmx24g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T AnalyzeCovariates -R hg19.fasta -L target.bed 
    -before recal_data_table.txt -after post_recal_data_table.txt -plots recalibration_plots.pdf



__(*)__ Print recalibrated reads
    
    java -Xmx24g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T PrintReads -R hg19.fasta -L target.bed -I dedup_rg_real.bam 
    -BQSR recal_data_table.txt -o dedup_rg_real_recal.bam


__(*)__ Now do variant calling
    
    java -Xmx24g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -nct 8 -L target.bed 
    -I dedup_rg_real_recal.bam --genotyping_mode DISCOVERY -o gatk.vcf

__(*)__ Questions
* Check out the before and after plots.
* How many variants were called?




#### Merge VCFs and VCF stats

__(*)__ VCFlib - merge

    vcfcombine freebayes.vcf gatk.vcf samtools.vcf > vcf_lib_merged.vcf

__(*)__ VCFlib - stats - shown here for one VCF file - repeat for all 3

    vcfstats freeb_call.vcf > freeb_call_vcf_lib_stats.txt



__(*)__ VCFtools

    export PERL5LIB=<full-path>/vcftools_0.1.12a/perl
    export PATH=<full-path>/tabix/tabix-0.2.6:$PATH

    ## Index (tabix) and pack files
    cp gatk.vcf gatk_tab.vcf
    bgzip gatk_tab.vcf
    tabix -p vcf gatk_tab.vcf.gz

    ## repeat for the other two VCF files

    vcf-merge freebayes_tab.vcf.gz gatk_tab.vcf.gz samtools_tab.vcf.gz > vcf_tools_merged.vcf
    vcf-stats freebayes_tab.vcf.gz > freebayes_tab_stats.txt





#### Filter variants
__(*)__ Using vcfutils
     
     <full-path>/bcftools/bin/vcfutils.pl varFilter -Q 20 -d 5 -D 200 samtools.vcf > samtools_filtered.vcf

__(*)__ Questions
* What other parameters can you specify for filtering variants?
* How many variants were filtered?




#### Display files in IGV

    (Download and open) IGV
    Load the BAM file and the VCF files into IGV
    Look at the mapping on Chr 11
    Check out the results of the different variant calling programs.




#### Annovar
__(*)__ First convert vcf into Annovar format

    <annovar-path>/convert2annovar.pl -format vcf4 -includeinfo freebayes.vcf > freebayes.avinput

__(*)__ Annotate with Gene information
    
    <annovar-path>/annotate_variation.pl -geneanno -buildver hg19 freebayes.avinput 
    /home/stephan/bin/annovar/annovar/humandb/

__(*)__ Annotate with Region information - ljb23

     <annovar-path>/annotate_variation.pl -regionanno -dbtype ljb23_all -buildver hg19 
     freebayes.avinput /home/stephan/bin/annovar/annovar/humandb/

__(*)__ Annotate with Region information - snp138

     <annovar-path>/annotate_variation.pl -regionanno -dbtype snp138 -buildver hg19 
     freebayes.avinput /home/stephan/bin/annovar/annovar/humandb/

__(*)__ Try converting the output files back to VCF (check if the appropriate columns are selected in cut)
     
     cat <annovar.file> | cut -f 9-18

__(*)__ Questions
* Look at the annotated VCF files.
* What databases does "ljb23" include?



#### SeattleSeq Annotation

__(*)__ Access<br/>
http://snp.gs.washington.edu/SeattleSeqAnnotation138/

__(*)__ Annotate VCF file

    Upload VCF file
    Specify VCF as return type
    submit
    You should receive an annotated VCF file to the specified email address
    
    
    
#### Useful information

__(*)__ Determine number of cores

    cat /proc/cpuinfo  

__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract information from a file

    grep -i "info" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100
