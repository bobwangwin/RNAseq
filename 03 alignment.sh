####################################
########### Third Step #############
##### RNA Sequences Alignment ######
####################################


#####  Genome Index Creation and seq align ######
####     hisat2/subjunc/star/bwa/bowtie2    #####

##############
### Hisat2 ########################################
# http://ccb.jhu.edu/software/hisat2/index.shtml ##
###################################################

		#@@@@ create genome index @@@@#

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

		#@@@@ alignment @@@@#
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
    $hisat -p 20 -x $hg38_index -1 $fq1 -2 $fq2 -S $fq1.hisat.sam
    ### samtools sort -O bam -@ 5  -o $sample.bam $sample.sam
done

########@@@@@@@@ ref: http://blog.biochen.com/archives/337 @@@@@@@@@@@@#############


###################
#### subjunc ######
###################

### create genome index ###
mkdir /Home/data/reference/subjunc #  move genome38.fa in here
subread-buildindex -o hg38 /Home/data/reference/subjunc/genome38.fa

#@@@@@ alignment @@@@@#
##### ./subjunc [options] -i <index_name> -r <input> -o <output> #######
*_1.fastq.gz >1
*_2.fastq.gz >2
paste 1 2 >config

subjunc_index='/Home/data/reference/subjunc'
cat $config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    ### $hisat -p 5 -x $mm10_index -1 $fq1 -2 $fq2 -S $sample.sam 2>$sample.hisat.log
    ### samtools sort -O bam -@ 5  -o $sample.bam $sample.sam
    $subjunc -T 20 -i $subjunc_index -r $fq1  -R $fq2 -o ${fq1}.subjunc.bam
done
#####@@@@@@@@ ref: http://www.bio-info-trainee.com/2775.html @@@@@@@@#####


####################
####  bowtie2   ####
####################

### create genome index ###

mkdir /Home/data/reference/bowtie2 #  move genome38.fa in here
bowtie2-build genome38.fa index

### bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>] 
bowtie_index='/Home/data/reference/bowtie2'
*_1.fastq.gz >1
*_2.fastq.gz >2
paste 1 2 >config
cat $config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
	$bowtie2 -p 20 -x $bowtie_index -1 $fq1 -2 $fq2 -S ${fq1}.bwa.sam
done
	# options for bowtie2
	# --phred33
	# --phred64
	# -q input fastq file
	# -f input fasta file
#####@@@@@@@@ ref: http://www.chenlianfu.com/?p=178 @@@@@@@@#####	


####################
####   BWA    ######
####################

### create genome index ###

mkdir /Home/data/reference/bwa #  move genome38.fa in here
# bwa index [-p prefix] [-a algoType] <in.db.fasta>
bwa index genome38.fa
	#-a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]
	#-p STR    prefix of the index [same as fasta name]
	#-b INT    block size for the bwtsw algorithm (effective with -a bwtsw) 
bwa_index='/Home/data/reference/bwa'

##### align ######
# bwa <command> [options]
*_1.fastq.gz >1
*_2.fastq.gz >2
paste 1 2 >config
cat $config |while read id
# bwa <command> [options]
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
	$bwa mem -t 10 -M $bowtie_index $fq1 $fq2 >${fq1}.bwa.sam
done
	#index         index sequences in the FASTA format
	#mem           BWA-MEM algorithm
    #fastmap       identify super-maximal exact matches
	#pemerge       merge overlapping paired ends (EXPERIMENTAL)
    #aln           gapped/ungapped alignment
    #samse         generate alignment (single ended)
    #sampe         generate alignment (paired ended)
    # bwasw         BWA-SW for long queries
##### align ######
#####@@@@@@@@ ref: https://blog.csdn.net/u014182497/article/details/51690341 @@@@@@@@#####	
