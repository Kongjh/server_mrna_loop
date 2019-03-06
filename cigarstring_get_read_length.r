# getwd()
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicAlignments")
library(tidyverse)
# 注意 tidy 一定在 之前就载入
library("GenomicAlignments")

# browseVignettes("GenomicAlignments")
# ?GenomicAlignments
## 针对三次重复的 
# i <- 1
getLENGTH <- function(in_file,out_file){
  bamfile <- readGAlignments(in_file, use.names = T)
  ####### 提取出 cigar string
  cigar(bamfile) -> a
  ###### 转换成 基因组length 长度
  cigarWidthAlongReferenceSpace(a, flag=NULL,N.regions.removed=FALSE) -> l
  ######## 存入
  write.table(l, out_file, row.names = F, col.names = F)
  print(i)
  # i <- i+1
}

# ?cigarWidthAlongReferenceSpace

inflist <- c("~/mrna_looped/Hh/Hh_1_ribo/all-0-16-mapped-Hh_1_ribo.bam","~/mrna_looped/Hh/Hh_2_ribo/all-0-16-mapped-Hh_2_ribo.bam",
             "~/mrna_looped/Hh/Hh_3_ribo/all-0-16-mapped-Hh_3_ribo.bam","~/mrna_looped/Hh/Hh_1_mrna/all-0-16-mapped-Hh_1_mrna.bam",
             "~/mrna_looped/Hh/Hh_2_mrna/all-0-16-mapped-Hh_2_mrna.bam","~/mrna_looped/Hh/Hh_3_mrna/all-0-16-mapped-Hh_3_mrna.bam")
outflist <- c("~/mrna_looped/Hh/Hh_1_ribo/Hh_1_ribo-length.csv","~/mrna_looped/Hh/Hh_2_ribo/Hh_2_ribo-length.csv",
             "~/mrna_looped/Hh/Hh_3_ribo/Hh_3_ribo-length.csv","~/mrna_looped/Hh/Hh_1_mrna/Hh_1_mrna-length.csv",
             "~/mrna_looped/Hh/Hh_2_mrna/Hh_2_mrna-length.csv","~/mrna_looped/Hh/Hh_3_mrna/Hh_3_mrna-length.csv")
#对于一次重复
# inflist <- c("~/mrna_looped/F1/F1_ribo/all-0-16-mapped-F1_ribo.bam","~/mrna_looped/F1/F1_mrna/all-0-16-mapped-F1_mrna.bam")
# outflist <- c("~/mrna_looped/F1/F1_ribo/F1_ribo-length.csv","~/mrna_looped/F1/F1_mrna/F1_mrna-length.csv")

#对于3次重复
inflist <- c("~/mrna_looped/gal1/gal1_1_mrna/all-0-16-mapped-gal1_1_mrna.bam","~/mrna_looped/gal1/gal1_2_mrna/all-0-16-mapped-gal1_2_mrna.bam",
             "~/mrna_looped/gal1/gal1_3_mrna/all-0-16-mapped-gal1_3_mrna.bam","~/mrna_looped/mfs/mfs_1_mrna/all-0-16-mapped-mfs_1_mrna.bam",
             "~/mrna_looped/mfs/mfs_2_mrna/all-0-16-mapped-mfs_2_mrna.bam","~/mrna_looped/mfs/mfs_3_mrna/all-0-16-mapped-mfs_3_mrna.bam")
outflist <- c("~/mrna_looped/gal1/gal1_1_mrna/gal1_1_mrna-length.csv","~/mrna_looped/gal1/gal1_2_mrna/gal1_2_mrna-length.csv",
              "~/mrna_looped/gal1/gal1_3_mrna/gal1_3_mrna-length.csv","~/mrna_looped/mfs/mfs_1_mrna/mfs_1_mrna-length.csv",
              "~/mrna_looped/mfs/mfs_2_mrna/mfs_2_mrna-length.csv","~/mrna_looped/mfs/mfs_3_mrna/mfs_3_mrna-length.csv")

num <- seq(1,6)
# num <- seq(1,2)
num
for (i in num){
  in_file <- inflist[i]
  out_file <- outflist[i]
  # print(in_file)
  getLENGTH(in_file,out_file)
}

##### 以下是过去式，已经整合成上面的 函数
#############################
######### file
file <- "G:/all-0-16-mapped-OE_3_mrna.bam"
file <- "G:/all-0-16-mapped-OE_2_mrna.bam"
file <- "G:/all-0-16-mapped-OE_3_ribo.bam"
file <- "G:/all-0-16-mapped-OE_1_mrna.bam"
file <- "G:/all-0-16-mapped-OE_2_ribo.bam"
file <- "G:/all-0-16-mapped-OE_1_ribo.bam"
file <- "G:/my-riboseq/all-0-16-mapped-Ho-1-28.bam"
file <- "G:/my-riboseq/all-0-16-mapped-Ho-1-T.bam"
file <- "G:/40s/all-0-16-mapped-ssu.bam"
file <- "G:/40s/all-0-16-mapped-rs.bam"
file <- "G:/igolia/all-0-16-mapped-fp_1.bam"
file <- "G:/40s/all-0-16-mapped-mrna.bam"
file <- "G:/40s/all-0-16-mapped-syors.bam"
file <- "G:/40s/all-0-16-mapped-syossu.bam"
file <- "G:/40s/utr/utr-0-16-mapped.bam"
file <- "G:/project/kozak/ssu/all-0-16-mapped-rs.bam"
file <- "G:/project/kozak/ssu/all-0-16-mapped-ssu.bam"
file <- "G:/kozak/ssu/all-0-16-mapped-mrna.bam"
file <- "G:/kozak/ssu/all-0-16-mapped-utr.bam"
######### 
# ?readGAlignments
bamfile <- readGAlignments(file, use.names = T)
# head(sort(table(cigar(bamfile)), decreasing=TRUE), 100)
####### 提取出 cigar string
cigar(bamfile) -> a
# head(a,100)
# #table 是计数
# table(a) -> b
# head(b,100)
###### 统计cigar 
# ?cigarOpTable
cigarOpTable(a) -> c
head(c,100)
colSums(c)
# ?`cigar-utils`
# explodeCigarOpLengths(a) -> d
# head(d,100)
# ?sequenceLayer

###### 转换成 基因组length 长度
cigarWidthAlongReferenceSpace(a, flag=NULL,N.regions.removed=FALSE) ->l
head(l,100)
data.frame(length = l) -> len
head(len,100)

###################################################
######## 存入
out_file <- "G:/OE_3_mrna-length.csv"
out_file <- "G:/OE_2_mrna-length.csv"
out_file <- "G:/OE_3_ribo-length.csv"
out_file <- "G:/OE_1_mrna-length.csv"
out_file <- "G:/OE_2_ribo-length.csv"
out_file <- "G:/OE_1_ribo-length.csv"
out_file <- "G:/my-riboseq/Ho-1-28-length.csv"
out_file <- "G:/my-riboseq/Ho-1-T-length.csv"
out_file <- "G:/40s/ssu-length.csv"
out_file <- "G:/40s/rs-length.csv"
out_file <- "G:/igolia/fp_1_length.csv"
out_file <- "G:/40s/mrna-length.csv"
out_file <- "G:/40s/syors-length.csv"
out_file <- "G:/40s/syossu-length.csv"
out_file <- "G:/40s/utr/utr-length.csv"
out_file <- "G:/project/kozak/ssu/rs-length.csv"
out_file <- "G:/project/kozak/ssu/ssu-length.csv"
out_file <- "G:/kozak/ssu/mrna-length.csv"
out_file <- "G:/kozak/ssu/utr-length.csv"

write.table(l, file = out_file, row.names = F, col.names = F)
# ?write.csv
# ??ScanBamParam 
# 
# ?sequenceLayer
# mcols(bamfile) -> hhh
# hhh
# mcols(bamfile)$seq -> h


