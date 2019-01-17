##########################################
##########     Second Step     ###########
##########     QC and trim     ###########
##########################################

# use fastqc for data quality control
fastqc -t 10 -o /home/project/fastqc/destination /home/project/rowdata/*.fastq  #use 10 threads
# ls *gz | xargs fastqc -t 10 -o /home/destination 

multiqc -o /home/project/mulitqc/destination /home/fastqc/destination
# multiqc ./

#############################################

###### filteration/ remove adaptor ##########
mkdir /home/project/clean | cd /home/project/clean

trim_galore -q 25 -phred33 --length 30 -e 0.1 --stringency 3 --paired -o /home/project/clean ../rowdata/1.fastq ../rowdata/2.fastq

		# -q/--quality <INT>      Trim low-quality ends from reads in addition to adapter removal. For
        #               RRBS samples, quality trimming will be performed first, and adapter
        #               trimming is carried in a second round. Other files are quality and adapter
        #               trimmed in a single pass. The algorithm is the same as the one used by BWA
        #               Default Phred score: 20.
        # --phred33               Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
        #               (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON.
		# --phred64               Instructs Cutadapt to use ASCII+64 quality scores as Phred scores
        #               (Illumina 1.5 encoding) for quality trimming.               
        # -a/--adapter <STRING>   Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
        #               try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
        #               small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
        #               '--small_rna'. If no adapter can be detected within the first 1 million sequences
        #               of the first file specified Trim Galore defaults to '--illumina'. A single base
        #               may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.               
        # --stringency <INT>      Overlap with adapter sequence required to trim a sequence. Defaults to a
        #               very stringent setting of 1, i.e. even a single bp of overlapping sequence
        #               will be trimmed off from the 3' end of any read.               
        # -e <ERROR RATE>         Maximum allowed error rate (no. of errors divided by the length of the matching
        #               region) (default: 0.1)               
        # --paired                This option performs length trimming of quality/adapter/RRBS trimmed reads for
        #               paired-end files. To pass the validation test, both sequences of a sequence pair
        #               are required to have a certain minimum length which is governed by the option
        #               --length (see above). If only one read passes this length threshold the
        #               other read can be rescued (see option --retain_unpaired). Using this option lets
        #               you discard too short read pairs without disturbing the sequence-by-sequence order
        #               of FastQ files which is required by many aligners.               
        # --length <INT>          Discard reads that became shorter than length INT because of either
        #               quality or adapter trimming. A value of '0' effectively disables
        #               this behaviour. Default: 20 bp.               
        #				For paired-end files, both reads of a read-pair need to be longer than
        #               <INT> bp to be printed out to validated paired-end files (see option --paired).
        #               If only one read became too short there is the possibility of keeping such
        #               unpaired single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.               

########## if paired seq ###########
####### creat a config list ########
ls *_1.fastq.gz >1
ls *_2.fastq.gz >2
paste 1 2 >config

###  trim in lines  #####
#########################
cat config | while read id
do
	arr=($id)
	fq1=${arr[0]}
	fq2=${arr[1]}
trim_galore -q 25 -phred33 --length 30 -e 0.1 --stringency 3 --paired -o /home/project/clean $fq1 $fq2 &
done
#########################
bash ./trim_in_lines.sh config

                        
##  Opt: fastqc after trim #### 
fastqc -t 10 -o /home/project/clean/fastqc/destination /home/project/clean/*  #use 10 threads
                       
## use trimmomatic
trimmomatic PE -threads 20 -phred33 0D-C_1.fq 0D-C_2.fq 0D-Cc_1.fq 0D-Cu_1.fq 0D-Cc_2.fq 0D-Cu_2.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
                        
PE [-threads <threads 线程数] 
        [-phred33 | -phred64] 质量体系，默认64，但我们现在一般都是33
        [-trimlog <logFile>] log文件
        <input 1> <input 2> 双端测序原始fq文件
        <paired output 1> <unpaired output 1>  双端1输出文件 、 过滤掉的文件
        <paired output 2> <unpaired output 2>  双端2同理
        <step > 
        <详细讲一下step：>
        1. ILLUMINACLIP： 根据上面qc部分的测试数据，我们设置TruSeq2-PE.fa:2:40:15
         TruSeq2-PE.fa是接头序列；
         2是比对时接头序列时所允许的最大错误数；
         40是要求PE的两条read同时和adapter序列比对，匹配度加起来超40%，那么就认为这对PE reads含有adapter，并在对应的位置需要进行切除；
         那么为何不是匹配100%？因为测序时并不是把adapter全测了，只测了一部分】
         15 指的是只要某条read的某部分与adapter超过了15%的相似度就认为包含，就要去除
        2. SLIDINGWINDOW： 滑动窗口长度我们设置为4:20，代表窗口长度为4，窗口中的平均质量值至少为20，否则会开始切除
        3. LEADING，指的是read开头的碱基是否要被切除的质量阈值，这里设为2
        4. TRAILING，指的是read末尾的碱基是否要被切除的质量阈值，这里设为2
        5. MINLEN，指的是read被切除后至少需要保留的长度，如果低于该长度，会被丢掉，这里设置40
