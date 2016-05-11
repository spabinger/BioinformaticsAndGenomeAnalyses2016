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
    
    java -Xmx24g -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R hg19.fasta -L target.bed 
    -before recal_data_table.txt -after post_recal_data_table.txt -plots recalibration_plots.pdf



__(*)__ Print recalibrated reads
    
    java -Xmx24g -jar GenomeAnalysisTK.jar -T PrintReads -R hg19.fasta -L target.bed -I deduprgreal.bam 
    -BQSR recal_data_table.txt -o dedup_rg_real_recal.bam


__(*)__ Now do variant calling
    
    java -Xmx24g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -nct 8 -L target.bed 
    -I dedup_rg_real_recal.bam --genotyping_mode DISCOVERY -o gatk.vcf

__(*)__ Questions
* Check out the before and after plots.
* How many variants were called?
