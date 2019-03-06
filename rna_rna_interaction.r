# enviroment --------------------------------------------------------------
setwd("~/rstudio/mrna_looped/")
getwd()
library(psych)
library(corrplot)#载入两个包
library("grid")
# install.packages("pheatmap")
library("pheatmap")
library(Hmisc)
library(PerformanceAnalytics)#加载包
library("openxlsx")
library("tidyverse")
# help(package="psych")
# detach("package:tidyverse")
options(digits=5)
options(tibble.width = Inf) # 表示 tibble 总是打印所有列 ，比如 使用 head 等函数的时候
#last_plot()
p <- ggplot() + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = "right",legend.text =  element_text(face="bold", size=23),
        legend.title = element_blank(),
        # axis.text.x = element_text(size = 26,face = "bold"),
        # axis.ticks = element_line(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(.20,"cm"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(size = 25,face = "bold"),  
        axis.title = element_text(size = 28, face = "bold")) + 
  scale_x_continuous(expand = c(0,0)) +scale_y_continuous(expand = c(0,0))
pp <- ggplot() + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = "right",legend.text =  element_text(face="bold", size=23),
        legend.title = element_blank(),
        # axis.text.x = element_text(size = 26,face = "bold"),
        # axis.ticks = element_line(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(.20,"cm"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(size = 25,face = "bold"),  
        axis.title = element_text(size = 28, face = "bold"))

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
a


# 尝试 ----------------------------------------------------------------------
ES <- read_delim("~/rstudio/mrna_looped/rna_rna_interation/ES-interactions.csv",delim = "\t")
str(ES)

lcl <- read_delim("~/rstudio/mrna_looped/rna_rna_interation/Lcl-interactions.csv",delim = "\t")
str(lcl)

# The following CSV files list interacting RNAs and regions that were detected by SPLASH in the four cell-line types. 
# The "bin" field is the approximate crosslink i.e. interaction position in the respective RNAs 
#(truncated window of 100 basepairs; see section "Binning and filtering of interacting RNA pairs" in the paper). 
# The "covfilter" field lists interaction counts.
# HeLa 
# Lymphoblastoid cells 
# ES cells 
# RA diff. cells 







# overlap -----------------------------------------------------------------
library("openxlsx")
splash <- read.xlsx("~/rstudio/mrna_looped/polya inter.xlsx",sheet = 1,colNames = T)# %>% rename(gene_id = X1) %>% 
str(splash)
paper_gene <- read_csv("~/rstudio/mrna_looped/paper_gene.csv",col_names = T)
str(paper_gene)
sm1 <- merge(splash,paper_gene,by.x="rna1",by.y="gene_id",all.x=T,sort = F) %>% rename(group1 =group)
sm2 <- merge(sm1,paper_gene,by.x="rna2",by.y="gene_id",all.x=T,sort = F) %>% rename(group2 =group)
write_csv(sm2,"~/rstudio/mrna_looped/overlap2.csv")

?merge

?sort

# 师姐的 ---------------------------------------------------------------------
library("openxlsx")
# 需要对 mRNA 进行 cutoff 吗？？？！！！！
mRNA <- read_csv("~/rstudio/mrna_looped/Factors/mho/mRNA_level.csv",col_names = T) %>%  
  mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3) %>%
  filter(mRNA_level != 0 , mRNA_level != Inf) %>% select(gene_id,mRNA_level)
str(mRNA)
(p0000 <- pp%+%mRNA + aes(x="1",y=log2(mRNA_level)) + geom_boxplot() )
plot(density(log2(mRNA$mRNA_level)))

# library(ComplexHeatmap)
# require(circlize)

# 不包括 NA -------------------------------------------------------------------
# intra -------------------------------------------------------------------
rrintra_ori <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIintra.xlsx",colNames = T) 
str(rrintra_ori)
rrintra <- rrintra_ori %>% mutate(coverage= (`coverage.Yeast-polyA-8-cycle-repl-2`+`coverage.Yeast-polyA-8-cycle-repl-1`)/2) %>%
  select(rna1,rna2,group,coverage)
str(rrintra)

str(mRNA)
str(rrintra)
f1_rrintra <- merge(rrintra,mRNA,by.x = "rna1",by.y = "gene_id",all.x = T) %>% mutate(intra_val_1 = coverage/mRNA_level) %>% 
  rename(mRNA_level_1 = mRNA_level)
f2_rrintra <- merge(f1_rrintra,mRNA,by.x = "rna2",by.y = "gene_id",all.x = T) %>% mutate(intra_val_2 = coverage/mRNA_level) %>%
  rename(mRNA_level_2 = mRNA_level)
# 因为 intra 的 rna1 和 rna2 一样，不一样的是 位置信息，即 id，所以
na.omit(f2_rrintra) -> f2_rrintra
final_intra <- select(f2_rrintra,rna1,rna2,group,intra_val_1,intra_val_2) %>% arrange(group,rna1) %>%
  mutate(id = seq(1,nrow(f2_rrintra))) 
str(final_intra)
head(final_intra)
dat1 <- select(final_intra,intra_val_1,intra_val_2)
rownames(dat1) <- final_intra$id
head(dat1)
summary(log2(dat1))
summary(dat1)

anno_row <- data.frame(group = final_intra$group)
rownames(anno_row) <- rownames(dat1)
pheatmap(log2(dat1), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row,cellheight = 2,
         color = colorRampPalette(colors = c("blue","yellow"))(1000),
         main = "RNA-RNA intra"
         ) # 86 个基因
?pheatmap

# Heatmap(log2(dat1),show_row_names = T,cluster_rows = F,
#         column_names_side = "top",row_names_side = "right"
#         #col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))
#         ) 

?Heatmap

# inter -------------------------------------------------------------------
rrinter_ori <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) 
str(rrinter_ori)
rrinter <- rrinter_ori %>% mutate(coverage= (`coverage.Yeast-polyA-8-cycle-repl-2`+`coverage.Yeast-polyA-8-cycle-repl-1`)/2) %>%
  select(rna1,group1,rna2,group2,coverage)
str(rrinter)

str(mRNA)
str(rrinter)
f1_rrinter <- merge(rrinter,mRNA,by.x = "rna1",by.y = "gene_id",all.x = T) %>% mutate(inter_val_1 = coverage/mRNA_level) %>% 
  rename(mRNA_level_1 = mRNA_level)
f2_rrinter <- merge(f1_rrinter,mRNA,by.x = "rna2",by.y = "gene_id",all.x = T) %>% mutate(inter_val_2 = coverage/mRNA_level) %>%
  rename(mRNA_level_2 = mRNA_level)
# 因为 inter 的 rna1 和 rna2 一样，不一样的是 位置信息，即 id，所以
na.omit(f2_rrinter) -> f2_rrinter
final_inter <- select(f2_rrinter,rna1,group1,rna2,group2,inter_val_1,inter_val_2) %>% arrange(group1,group2,rna1) %>%
  mutate(id = seq(1,nrow(f2_rrinter))) 
str(final_inter)
head(final_inter)
dat1 <- select(final_inter,inter_val_1,inter_val_2)
rownames(dat1) <- final_inter$id
head(dat1)
summary(log2(dat1))
summary(dat1)

anno_row1 <- data.frame(group1 = final_inter$group1)
rownames(anno_row1) <- rownames(dat1)
pheatmap(log2(dat1), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row1,cellheight = 1,
         color = colorRampPalette(colors = c("blue","yellow","red"))(1000),
         main = "RNA-RNA inter"
)# 500 个基因
?pheatmap

# 按照 group2 排序
final_inter2 <- select(f2_rrinter,rna1,group1,rna2,group2,inter_val_1,inter_val_2) %>% arrange(group2,group1,rna2) %>%
  mutate(id = seq(1,nrow(f2_rrinter))) 
str(final_inter2)
head(final_inter2)
dat2 <- select(final_inter2,inter_val_1,inter_val_2)
rownames(dat2) <- final_inter2$id
head(dat2)
summary(log2(dat2))
summary(dat2)

anno_row2 <- data.frame(group2 = final_inter2$group2)
rownames(anno_row2) <- rownames(dat2)
pheatmap(log2(dat2), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row2,cellheight = 1,
         color = colorRampPalette(colors = c("blue","yellow","red"))(1000),
         main = "RNA-RNA inter"
)# 500 个基因


# Heatmap(log2(dat1),show_row_names = T,cluster_rows = F,
#         column_names_side = "top",row_names_side = "right"
#         #col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))
#         ) 

?Heatmap

# 包括 NA ------------------------------------------------------------------
# intra -------------------------------------------------------------------
rrintra_ori <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIintra.xlsx",colNames = T) 
str(rrintra_ori)
rrintra <- rrintra_ori %>% mutate(coverage= (`coverage.Yeast-polyA-8-cycle-repl-2`+`coverage.Yeast-polyA-8-cycle-repl-1`)/2) %>%
  select(rna1,rna2,group,coverage)
str(rrintra)

str(mRNA)
str(rrintra)
f1_rrintra <- merge(rrintra,mRNA,by.x = "rna1",by.y = "gene_id",all.x = T) %>% mutate(intra_val_1 = coverage/mRNA_level) %>% 
  rename(mRNA_level_1 = mRNA_level)
f2_rrintra <- merge(f1_rrintra,mRNA,by.x = "rna2",by.y = "gene_id",all.x = T) %>% mutate(intra_val_2 = coverage/mRNA_level) %>%
  rename(mRNA_level_2 = mRNA_level)
# 因为 intra 的 rna1 和 rna2 一样，不一样的是 位置信息，即 id，所以
str(f2_rrintra)
"e" -> f2_rrintra[is.na(f2_rrintra$group),3]
na.omit(f2_rrintra) -> f2_rrintra

final_intra <- select(f2_rrintra,rna1,rna2,group,intra_val_1,intra_val_2) %>% arrange(group,rna1) %>%
  mutate(id = seq(1,nrow(f2_rrintra))) 
str(final_intra)
head(final_intra)
dat1 <- select(final_intra,intra_val_1,intra_val_2)
rownames(dat1) <- final_intra$id
head(dat1)
summary(log2(dat1))
summary(dat1)

anno_row <- data.frame(group = final_intra$group)
rownames(anno_row) <- rownames(dat1)
pheatmap(log2(dat1), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row,cellheight = 2,
         color = colorRampPalette(colors = c("blue","yellow"))(1200),
         main = "RNA-RNA intra"
) # 99 个基因
?pheatmap

# Heatmap(log2(dat1),show_row_names = T,cluster_rows = F,
#         column_names_side = "top",row_names_side = "right"
#         #col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))
#         ) 

?Heatmap

# inter -------------------------------------------------------------------
rrinter_ori <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) 
str(rrinter_ori)
rrinter <- rrinter_ori %>% mutate(coverage= (`coverage.Yeast-polyA-8-cycle-repl-2`+`coverage.Yeast-polyA-8-cycle-repl-1`)/2) %>%
  select(rna1,group1,rna2,group2,coverage)
str(rrinter)

str(mRNA)
str(rrinter)
f1_rrinter <- merge(rrinter,mRNA,by.x = "rna1",by.y = "gene_id",all.x = T) %>% mutate(inter_val_1 = coverage/mRNA_level) %>% 
  rename(mRNA_level_1 = mRNA_level)
f2_rrinter <- merge(f1_rrinter,mRNA,by.x = "rna2",by.y = "gene_id",all.x = T) %>% mutate(inter_val_2 = coverage/mRNA_level) %>%
  rename(mRNA_level_2 = mRNA_level)
# 因为 inter 的 rna1 和 rna2 一样，不一样的是 位置信息，即 id，所以
str(f2_rrinter)
"e" -> f2_rrinter[is.na(f2_rrinter$group1),3]
"e" -> f2_rrinter[is.na(f2_rrinter$group2),4]
na.omit(f2_rrinter) -> f2_rrinter
final_inter <- select(f2_rrinter,rna1,group1,rna2,group2,inter_val_1,inter_val_2) %>% arrange(group1,group2,rna1) %>%
  mutate(id = seq(1,nrow(f2_rrinter))) 
str(final_inter)
head(final_inter)
dat1 <- select(final_inter,inter_val_1,inter_val_2)
rownames(dat1) <- final_inter$id
head(dat1)
summary(log2(dat1))
summary(dat1)

anno_row1 <- data.frame(group1 = final_inter$group1)
rownames(anno_row1) <- rownames(dat1)
pheatmap(log2(dat1), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row1,cellheight = 1,
         color = colorRampPalette(colors = c("blue","yellow","red"))(1000),
         main = "RNA-RNA inter"
)# 587 个基因
?pheatmap

# 按照 group2 排序
# 因为前面对 f2_rrinter 进行了 控制，所以这里就不用了
final_inter2 <- select(f2_rrinter,rna1,group1,rna2,group2,inter_val_1,inter_val_2) %>% arrange(group2,group1,rna2) %>%
  mutate(id = seq(1,nrow(f2_rrinter))) 
str(final_inter2)
head(final_inter2)
dat2 <- select(final_inter2,inter_val_1,inter_val_2)
rownames(dat2) <- final_inter2$id
head(dat2)
summary(log2(dat2))
summary(dat2)

anno_row2 <- data.frame(group2 = final_inter2$group2)
rownames(anno_row2) <- rownames(dat2)
pheatmap(log2(dat2), cluster_row = F, cluster_cols = F,
         display_numbers = F,show_rownames = F
         ,cellwidth = 40,annotation_row = anno_row2,cellheight = 1,
         color = colorRampPalette(colors = c("blue","yellow","red"))(1000),
         main = "RNA-RNA inter"
)# 587 个基因


# Heatmap(log2(dat1),show_row_names = T,cluster_rows = F,
#         column_names_side = "top",row_names_side = "right"
#         #col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))
#         ) 

?Heatmap



# cicro score -------------------------------------------------------------
mRNA <- read_csv("~/rstudio/mrna_looped/Factors/ho/mRNA_level.csv",col_names = T) %>%  
  mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3) %>%
  filter(mRNA_level != 0 , mRNA_level != Inf) %>% select(gene_id,mRNA_level)
str(mRNA)
TE <- read_csv("~/rstudio/mrna_looped/Factors/OE/TE.csv",col_names = T) %>%  
  filter(TE != 0 , TE != Inf) %>% select(gene_id,TE)
str(TE)

rrintra_ori <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIintra.xlsx",colNames = T) 
str(rrintra_ori)
rrintra <- rrintra_ori %>% mutate(cover= (`coverage.Yeast-polyA-8-cycle-repl-2`+`coverage.Yeast-polyA-8-cycle-repl-1`)/2,
                                  span = start2-start1) %>% mutate(multip = cover*span) %>%
  select(rna1,cover,span,multip,len1)
str(rrintra)
f1_rrintra <- group_by(rrintra,rna1,len1) %>% summarise(add_mul = sum(multip),add_cover = sum(cover))

write_csv(f1_rrintra,"~/rstudio/mrna_looped/check1.csv")
write_csv(rrintra,"~/rstudio/mrna_looped/check2.csv")

str(f1_rrintra)
str(mRNA)
############ 
f2_rrintra <- merge(TE,f1_rrintra,by.x="gene_id",by.y = "rna1") %>% mutate(circul = add_mul/(add_cover*len1))
str(f2_rrintra)
f3_rrintra <- merge(mRNA,f2_rrintra,by="gene_id") %>% mutate(circul2 = circul/mRNA_level)
str(f3_rrintra)
cor1 <- cor.test(f2_rrintra$TE,f2_rrintra$circul,method = "s",exact = F)
cor1 <- cor.test(f3_rrintra$TE,f3_rrintra$circul,method = "s",exact = F)

sig <- paste(signif(cor1$estimate,3),signif(cor1$p.value,4),sep = " ")
sig
# mean_val_33 <- summarise(group_by(d_both,id),mean_te=mean(log2(TE)))

(p100a <- pp%+%f2_rrintra + aes(x = log2(TE), y=circul)
  + geom_point()
  + labs(x ="log2(TE)",y="circular_score")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
             label = sig, size = 10))
