#################################
#######  Fourth Step  ###########
#### Counting the Align reads ###
#################################

### convert *.sam to *.bam ###
ls *.sam | while read id ;do (samtools sort -O bam @ -10 $(basename ${id} ".sam").bam ${id}; done
ls *.bam | xargs -i samtools index {}

####### qc analysis #######
ls *.bam | while read id ;do (samtools flagstat @ -20 $id > $(basename ${id} ".sam").flagstat); done
# samtools view *.bam | less -SN
# multiqc *.flagstat

