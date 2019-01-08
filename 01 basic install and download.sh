###### Personal protocol for RNAseq 
###### Use Shell for Data polish and R for Data analysis 
###### Please follow each step, if there is any question, send email to logo-wang@163.com #####
###### reference https://youtu.be/ao23_xjPTP8 cite: Jimmy/biotrainee 

##########################################
##########     First Step     ############
##### Relative software installation #####
##########################################

###########   basic Env build up      ###########
####### essential for R and conda software ######
sudo apt-get -y install libcur14-gnutls-dev
sudo apt-get -y install libxml2-dev
sudo apt-get -y install libssl-dev
sudo apt-get -y install libmariadb-client-lgpl-dev

####### Conda Env ################
####### OS: Ubuntu 16 X86_64######
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3  # -b auto install -p point out the install location
# echo "export PATH=$PREFIX/bin:" '$PATH' >> ~/.bashrc
source ~/.bashrc

# add tsinghua mirror, conda config --show
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda

# create rna env
conda create -n rna python=2
source activate rna
source deactivate

################ done for Env ###########

# install software
conda install -y bwa sra-tools fastqc multiqc trimmomatic cutadapt trim-galore star hisat2 tophat bowtie2 subread htseq bedtools deeptools

############### more if require #########


############################################
##########     Opt: Second Step     ########
##########   Download sra/gft data   #######
############################################

# download sra files
cat id | while read id ;do (prefetch -O /HOME/directory $id &);done
# transform into fastq 
fastq-dump --gzip --split-3 -O ./HOME/directory /list/location/id.sra


# download gft file for annotation
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz # hg38
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz # hg19


