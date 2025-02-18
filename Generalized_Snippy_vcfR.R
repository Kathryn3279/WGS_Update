# Generalized WGS snippy script
# Written for running in linux with conda environment containing snippy and samtools

# load R libraries needed
library(tidyverse)
library(vcfR)
library(Biostrings)
library(ape)
library(ggplot2)

# activate conda WGS environment or whichever has compatible snippy and samtools
#system("conda activate WGS")

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
ref1_fasta<-ape::read.dna("/mnt/c/Desktop/SEDS_Lab/Bioinformatics/WGS/References/H37Rv.fa", format = "fasta")
ref1_gff<-read.gff("/mnt/c/Desktop/SEDS_Lab/Bioinformatics/WGS/References/H37Rv.gff")
ref1_Name<-"H37Rv"

# loop through samples and run snippy as system call
# repeat for as many reference files as desired
# alignment separated by output folder
for(i in seq(1,length(Fastq_Files_FullPath),2)){
  message("Running ",ref1_Name, " ",SampleNames[i]," alignment.")
  system(paste0("snippy --prefix ",ref1_Name,"_",SampleNames[i]," --outdir snippy_out_",ref1_Name,"_",SampleNames[i]," --force  --ref ",ref1_Path, " --R1 ", Fastq_Files_FullPath[i]," --R2 ",Fastq_Files_FullPath[i+1]))

  # Repeat for as many references as desired
  
  
  # read in raw vcf file
  raw_vcf<-read.vcfR(paste0(wd_FullPath, "snippy_out_",ref1_Name,"_",SampleNames[i],"/",ref1_Name,"_",SampleNames[i],".raw.vcf"))
  
  # MAKE SURE NAMES MATCH
  chrom <- create.chromR(name=as.character((ref1_gff$seqid)[1]), vcf=raw_vcf, seq=ref1_fasta, ann=ref1_gff)
  
  # Parsing
  extracted<-as.data.frame(extract.gt(raw_vcf, element = "DP", as.numeric = TRUE))
  variant_position<-as.numeric(unlist(strsplit(row.names(extracted),"_"))[seq(2,2*length(row.names(extracted)),2)])
  all_count<-as.numeric(extract.gt(raw_vcf, element = "DP", as.numeric = TRUE)[,1])
  
  ref<-(raw_vcf@fix)[,4]
  ref_count<-as.numeric(extract.gt(raw_vcf, element = "RO", as.numeric = TRUE)[,1])
  ref_quality<-as.numeric(extract.gt(raw_vcf, element = "QR", as.numeric = TRUE)[,1])
  
  alt<-(raw_vcf@fix)[,5]
  alt_count<-as.numeric(extract.gt(raw_vcf, element = "AO", as.numeric = TRUE)[,1])
  alt_quality<-as.numeric(extract.gt(raw_vcf, element = "QA", as.numeric = TRUE)[,1])
  
  qual<-as.numeric((raw_vcf@fix)[,6])
  outDF<-data.frame("VariantPosition"=variant_position,
                    "AllCounts" =all_count,
                    "QUAL" = qual,
                    "Ref_Base"=ref,
                    "Ref_Count"=ref_count,
                    "Ref_Quality"=ref_quality,
                    "Alt_Base"=alt,
                    "Alt_Count"=alt_count,
                    "Alt_Quality"=alt_quality)
  outDF$Perc_Alt<-c(outDF$Alt_Count/outDF$AllCounts)
  outDF<-outDF%>%
    mutate("FreqAlt"=case_when(Perc_Alt>=0.40~"TRUE",
                           TRUE~"FALSE"))
  write.csv(file=paste0(wd_FullPath, "snippy_out_",ref1_Name,"_",SampleNames[i],"/Parsed_",ref1_Name,"_",SampleNames[i],".raw.vcf"),outDF,row.names = FALSE)
  outDF_Lenient<-outDF[outDF$QUAL >=20 & outDF$AllCounts>=10 & outDF$FreqAlt>=0 ,]
  
  Phred_Plot<-ggplot(data=outDF)+geom_point(aes(x=VariantPosition, y =QUAL, size = Perc_Alt, color=FreqAlt),alpha=0.5)+
    geom_hline(aes(yintercept = 20))+theme_bw()
  ggsave(paste0(wd_FullPath, "snippy_out_",ref1_Name,"_",SampleNames[i],"/",ref1_Name,"_",SampleNames[i],"_Raw_PhredDotPlot.png"),Phred_Plot)

  write.csv(x = outDF_Lenient, paste0(wd_FullPath, "snippy_out_",ref1_Name,"_",SampleNames[i],"/",ref1_Name,"_",SampleNames[i],"_Lenient.csv"),row.names = FALSE)
  
  }