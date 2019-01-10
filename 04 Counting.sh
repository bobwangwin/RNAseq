#################################
#######  Fourth Step  ###########
#### Counting the Align reads ###
#################################
	
	#@@@@@@@@@@@@@@@@@@@@@@#
### use samtools to convert #####
    #@@@@@@@@@@@@@@@@@@@@@@#
    
### convert *.sam to *.bam ######
ls *.sam | while read id ;do (samtools sort -O bam @ -10 $(basename ${id} ".sam").bam ${id}; done
ls *.bam | xargs -i samtools index {}

####### qc analysis #######
ls *.bam | while read id ;do (samtools flagstat @ -20 $id > $(basename ${id} ".sam").flagstat); done
# samtools view *.bam | less -SN
# multiqc *.flagstat

	#@@@@@@@@@@@@@@@@@@@@@@#
#### use htseq to count   #####
    #@@@@@@@@@@@@@@@@@@@@@@#

# A genome.gtf file is required, which could check 01 basic install and download.sh
##  htseq-count [options] <input.bam> <annotation_file>
htseq-count -f bam -r pos ./clean/*.bam /home/annotation/gtf.gz

##### ref: https://blog.csdn.net/herokoking/article/details/78257714 ########


	   #@@@@@@@@@@@@@@@@@@@@@@#
#### use featureCounts to count   #####
       #@@@@@@@@@@@@@@@@@@@@@@#

#featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2]       
featureCounts -T 20 -t exon -g gene_id -a /home/annotation/gtf.gz -o counts.txt ./clean/*.bam 1>count.id.log 2>&1 &
featureCounts -T 20 -t exon -g gene_name -a /home/annotation/gtf.gz -o counts.txt ./clean/*.bam 1>count.id.log 2>&1 &
##### ref : https://www.cnblogs.com/xudongliang/p/8417654.html  #########


	   #@@@@@@@@@@@@@@@@@@@@@@#
#### use bedtools to count   #####
       #@@@@@@@@@@@@@@@@@@@@@@#
######## under developed  ########

