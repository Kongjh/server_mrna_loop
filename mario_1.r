# enviroment --------------------------------------------------------------
setwd("~/rstudio/mrna_looped/")
getwd()
library(psych)
library(corrplot)#载入两个包
library("grid")
library(Hmisc)
library(PerformanceAnalytics)#加载包
library("tidyverse")
library(ggplot2)
options(digits=5)
options(tibble.width = Inf) # 表示 tibble 总是打印所有列 ，比如 使用 head 等函数的时候
#last_plot()
p <- ggplot() + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = c(0.8,0.8),legend.text =  element_text(face="bold", size=23),
        axis.text.x = element_text(size = 26,face = "bold"), axis.text.y = element_text(size = 26,face = "bold"), 
        axis.title = element_text(size = 26, face = "bold"))
p2 <- ggplot() + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = c(0.2,0.8),legend.text =  element_text(face="bold", size=23),
        axis.text.x = element_text(size = 26,face = "bold"), axis.text.y = element_text(size = 26,face = "bold"), 
        axis.title = element_text(size = 26, face = "bold"))
# file chunk --------------------------------------------------------------
# ## 第一次加载进来就行，以后用下面新的了！！！ 不用重复
d0 <- read_delim("~/rstudio/mrna_looped/MARIO/ESdual_AATG_fragment_paired_align.txt",delim = "\t",col_names = F)
# # # header <- c("Chromosome_RNA_1", "Start_RNA_1", "End_RNA_1", "Strand_RNA_1", "Fragment_in_read_RNA_1"," Reference_RNA_1",
# # #             "Type_RNA_1", "Name_RNA_1", "Sub_type_RNA_1", "Strand_Properness_RNA_1", "Read_ID", "Chromosome_RNA_2", "Start_RNA_2",
# # #             "END_RNA_2", "Strand_RNA_2"," Fragment_in_read_RNA_2", "Reference_RNA_2", "Type_RNA_2"," Name_RNA_2", "Sub_type_RNA_2",
# # #             "Strand_Properness_RNA_2")
# # head(d0) # 3,302,146
header <- c("r1_chr", "r1_start", "r1_end", "r1_strand", "r1_reads",
            "r1_type", "r1_gene", "r1_subtype", "reads_id", "r2_chr", "r2_start",
            "r2_end", "r2_strand"," r2_reads","r2_type","r2_gene", "r2_subtype")
colnames(d0) <- header
# d1 <- filter(d0,(r1_gene == r2_gene) & (r1_strand == r2_strand)) # same gene 2,159,419  same gene + sanme strand 280,783
# write_csv(d1,"~/rstudio/mrna_looped/MARIO/mario_r1_r2_samegenestrand.csv",col_names = T)
d00 <- filter(d0, r1_type == "protein_coding" & r2_type == "protein_coding" )

length(unique(factor(d0$r1_gene)))
unique(factor(d0$r1_subtype))


#加载数据
d1 <- read_csv("~/rstudio/mrna_looped/MARIO/mario_r1_r2_samegenestrand.csv",col_names = T) # same gene + sanme strand 280,783
rep <- d1 %>% group_by(r1_gene,r2_gene) %>% summarise(count = n())
head(d1)
length(unique(factor(d1$r1_gene)))
d2 <- filter(d1,r1_type == "protein_coding")
colnames(d2)
unique(factor(d2$r1_subtype))
length(unique(factor(d2$r1_gene)))
# [1] intron utr5   utr3   exon   .     
# Levels: . exon intron utr3 utr5
d3 <- filter(d2, r1_subtype != "intron" & r2_subtype != "intron") 
d3 <- filter(d2,((r1_strand == "+") & (r2_start > r1_end) )|( (r1_strand == "-") & (r2_start < r1_end) ) )

d33 <- filter(d1, r1_subtype != "intron" & r2_subtype != "intron") 
unique(factor(d33$r1_type))


unique(factor(d33$r2_subtype))



d3 <- filter(d2,((r1_strand == "+") & (r2_start > r1_end) )|( (r1_strand == "-") & (r2_start < r1_end) ) )
length(unique(factor(d3$r1_gene)))

unique(factor(d1$r1_type))
# [1] snRNA          rRNA_repeat    LTR            non            LINE           SINE           protein_coding misc_RNA
# [9] lincRNA        Unknown        DNA            Satellite      tRNA           rRNA           Other          Simple_repeat
# [17] Low_complexity

unique(factor(d1$r1_subtype))
[1] .              ERVK           L1             ERV1           ERVL           Alu            MaLR           intron
[9] utr5           Unknown        Tip100         B2             B4             Satellite      utr3           L2
[17] ID             MIR            exon           Other          MER2_type      MER1_type      Simple_repeat  Low_complexity
[25] AcHobo
unique(factor(d1$r2_subtype))

dutr <- filter(d1, r1_subtype == "utr5" | r1_subtype == "utr3" | r2_subtype == "utr5" | r2_subtype == "utr3"  ) # 115 了
unique(factor(dutr$r1_gene)) #才11个基因

unique(factor(dutr$r1_subtype))
# [1] protein_coding           rRNA                     non                      LINE                     snRNA                   
# [6] rRNA_repeat              LTR                      lincRNA                  SINE                     Low_complexity          
# [11] pseudogene               misc_RNA                 Other                    snoRNA                   antisense               
# [16] processed_transcript     Simple_repeat            DNA                      Satellite                Unknown                 
# [21] miRNA                    tRNA                     sense_intronic           non_coding               RC                      
# [26] 3prime_overlapping_ncrna
unique(factor(dutr$r1_type))
# [1] protein_coding           lincRNA                  pseudogene               misc_RNA                 antisense               
# [6] processed_transcript     sense_intronic           3prime_overlapping_ncrna

