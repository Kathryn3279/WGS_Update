# Generalized WGS snippy script
# Written for running in linux with conda environment containing snippy and samtools

# load R libraries needed
library(tidyverse)

# activate conda WGS environment or whichever has compatible snippy and samtools
system("conda activate WGS")

# Following variables changed each time script is run:
# wd_FullPath - full path to folder that contains reads folder
# readsFolderName - name of folder that contains ONLY reads INCLUDING FORWARD SLASH
# ref1_Path - FULL PATH of 1st reference in gb format
# ref1_Name - human readable name of 1st reference
# ref2_Path- FULL PATH of 2nd reference in gb format
#  ref2_Name - human readable name of 2nd reference
# ... refN_Path - FULL PATH of Nth reference in gb format
# ... refN_Name - human readable name of Nth reference


# set working directory using full path
wd_FullPath<-"/mnt/c/Desktop/SEDS_Lab/Bioinformatics/WGS/Test_WGS/"
setwd(wd_FullPath)

# make a list of the reads files INCLUDING FORWARD SLASH
# folder must ONLY contain reads
readsFolderName<-"raw_fastq/"
system(paste0("ls ",wd_FullPath,readsFolderName, " > Fastq_Files.txt"))

# read in the txt file with the read names
Fastq_Files<-read.table(paste0(wd_FullPath,"Fastq_Files.txt"),sep="\t")


# Need to manually inspect fastq names to parse sample names
# In this case, paired end reads, sample names before "_R"
SampleNames<-Fastq_Files$V1

# remove lane information
SampleNames<-unlist(str_split(SampleNames,"_R"))[seq(1,2*(length(SampleNames)),2)]

# list of full paths
Fastq_Files_FullPath<-paste0(wd_FullPath,readsFolderName,as.vector(Fastq_Files$V1))

# make ref for all references using FULL PATH
ref1_Path<-"/mnt/c/Desktop/SEDS_Lab/Bioinformatics/WGS/References/H37Rv.gb"
ref1_Name<-"H37Rv"

# loop through samples and run snippy as system call
# repeat for as many reference files as desired
# alignment separated by output folder
for(i in seq(1,length(Fastq_Files_FullPath),2)){
  message("Running ",ref1_Name, " ",SampleNames[i]," alignment.")
  system(paste0("snippy --prefix ",ref1_Name,"_",SampleNames[i]," --outdir snippy_out_",ref1_Name,"_",SampleNames[i]," --force  --ref ",ref1_Path, " --R1 ", Fastq_Files_FullPath[i]," --R2 ",Fastq_Files_FullPath[i+1]))

  # Repeat for as many references as desired
  }