####################################
########### Third Step #############
##### RNA Sequences Alignment ######
####################################


#####  Genome Index Creation and seq align ######
####     hisat2/subjunc/stat/bwa/bowtie2    #####

##############
### Hisat2 ########################################
# http://ccb.jhu.edu/software/hisat2/index.shtml ##
###################################################

### create genome index ###

cd /Home/data/reference/hisat/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar xvzfO hg38.chromFa.tar.gz > genome38.fa 
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz

	## for genome
	hisat2-build -p 20 genome38.fa genome
	
	## for genmoe_snp_trans  Option
	# extract_exons.py Homo_sapiens.GRCh38.94.chr.gtf > genome.exon
	# extract_splice_sites.py Homo_sapiens.GRCh38.83.chr.gtf > genome.ss
	# extract_snps.py snp142Common.txt > genome.snp ## cannot fine snp142common
	# hisat2-build -p 20 genome.fa --snp genome.snp --ss genome.ss --exo n genome.exon genome_snp_tran

### alignment ####
##  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
*_1.fastq.gz >1
*_2.fastq.gz >2
paste 1 2 >config

hg38_index='/Home/data/reference/hisat/'
cat $config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    $hisat -p 20 -x $hg38_index -1 $fq1 -2 $fq2 -S $fq1.sam
    ### samtools sort -O bam -@ 5  -o $sample.bam $sample.sam
done

############### ref: http://blog.biochen.com/archives/337 ########################


##################
### subjunc ######
##################

### create genome index ###
mkdir /Home/data/reference/subjunc #  move genome38.fa in here
subread-buildindex -o hg38 /Home/data/reference/subjunc/genome38.fa

### alignment #####
##### if config_list require  #######
##### ./subjunc [options] -i <index_name> -r <input> -o <output> #######
*_1.fastq.gz >1
*_2.fastq.gz >2
paste 1 2 >config
 
hg38_index='/Home/data/reference/subjunc'
cat $config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    ### $hisat -p 5 -x $mm10_index -1 $fq1 -2 $fq2 -S $sample.sam 2>$sample.hisat.log
    ### samtools sort -O bam -@ 5  -o $sample.bam $sample.sam
    $subjunc -T 20 -i $hg38_index -r $fq1  -R $fq2 -o ${fq1}_subjunc.bam
done
########### ref: http://www.bio-info-trainee.com/2775.html ###########



