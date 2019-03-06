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

# #<2> 加载 exon length -----------------------------------------------------
!!!!!这里错了。不过自己也没有改了
dexonlen <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length")
head(dexonlen);nrow(dexonlen) # 6584

# #<1> 加载 data ------------------------------------------------------------
dmrna1 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_1_mrna_final.csv", col_names = T )
dribo1 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_1_ribo_final.csv", col_names = T )

dmrna2 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_2_mrna_final.csv", col_names = T )
dribo2 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_2_ribo_final.csv", col_names = T )

dmrna3 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_3_mrna_final.csv", col_names = T )
dribo3 <- read_csv("~/rstudio/mrna_looped/Hh/last_Hh_3_ribo_final.csv", col_names = T )

# 单个样本
dmrna1 <- read_csv("~/rstudio/mrna_looped/F1/last_F1_mrna_final.csv", col_names = T )
dribo1 <- read_csv("~/rstudio/mrna_looped/F1/last_F1_ribo_final.csv", col_names = T )

# #<3> 计算 rpkm ------------------------------------------------------------
getRPKM <- function(data,dexonlen){
  d2 <- group_by(data,gene_id) %>% summarise(counts = n()) %>% merge(dexonlen, by = "gene_id")
  allreads <- sum(d2$counts)/1000000
  drpkm <- mutate(d2, mRNA_level = counts/(allreads*exon_length)) %>% filter(grepl("^Y",gene_id))
  return(drpkm)
}

drpkm1 <- getRPKM(dmrna1,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_1 = mRNA_level)
drpkm2 <- getRPKM(dmrna2,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_2 = mRNA_level)
drpkm3 <- getRPKM(dmrna3,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_3 = mRNA_level)
drpkm123 <- merge(drpkm1,drpkm2, by = "gene_id", all = T) %>% merge(drpkm3,by = "gene_id", all = T); drpkm123[is.na(drpkm123)] <- 0

write_csv(drpkm123,"~/rstudio/mrna_looped/gal1/mRNA_level.csv",col_names = T)
write_csv(drpkm1,"~/rstudio/mrna_looped/F1/mRNA_level.csv",col_names = T)

# # TRU seq 建库 根据 htseq count 结果求fpkm ------------------------------------------------------------
?read_delim
dmrna1 <- read_tsv("~/mrna_looped/mho/mho_1_mrna/report-count-mho_1_mrna.txt", col_names = F ) %>%
  select(-X2) %>% rename(gene_id = X1, counts = X3)
dmrna2 <- read_tsv("~/mrna_looped/mho/mho_2_mrna/report-count-mho_2_mrna.txt", col_names = F ) %>%
  select(-X2) %>% rename(gene_id = X1, counts = X3)
dmrna3 <- read_tsv("~/mrna_looped/mho/mho_3_mrna/report-count-mho_3_mrna.txt", col_names = F ) %>%
  select(-X2) %>% rename(gene_id = X1, counts = X3)

head(dexonlen) # 6584
getRPKM <- function(data,dexonlen){
  d2 <- merge(data,dexonlen, by = "gene_id")
  allreads <- sum(d2$counts)/1000000
  drpkm <- mutate(d2, mRNA_level = counts/(allreads*exon_length)) %>% filter(grepl("^Y",gene_id))
  return(drpkm)
}

drpkm1 <- getRPKM(dmrna1,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_1 = mRNA_level)
drpkm2 <- getRPKM(dmrna2,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_2 = mRNA_level)
drpkm3 <- getRPKM(dmrna3,dexonlen) %>% select(c("gene_id","mRNA_level")) %>% rename(mRNA_level_3 = mRNA_level)
drpkm123 <- merge(drpkm1,drpkm2, by = "gene_id", all = T) %>% merge(drpkm3,by = "gene_id", all = T)
drpkm123[is.na(drpkm123)] <- 0

write_csv(drpkm123,"~/rstudio/mrna_looped/Factors/mho/mRNA_level.csv",col_names = T)


# # 简单看一下各种分布 -------------------------------------------------------------
# s_drpkm123 <- gather(drpkm123,exper,mRNA_level,-gene_id)
# (p00 <- p%+%s_drpkm123 + aes(x= factor(exper), y= log(mRNA_level)) + geom_boxplot(notch = T))
# (p00 <- p%+%s_drpkm123 + aes(color = factor(exper), log(mRNA_level)) + geom_line(stat = "density"))
# (plot(density(log(drpkm1$mRNA_level_1))))
# (lines(density(log(drpkm2$mRNA_level_2))))
# (lines(density(log(drpkm3$mRNA_level_3))))

# # 重复性比较 -----------------------------------------------------------------
pd <- merge(drpkm1,drpkm2, by = "gene_id")
(p11 <- p%+%pd + aes(log(mRNA_level_1),log(mRNA_level_2)) + geom_point() + 
    geom_smooth(method = "lm", se=T, color="red", formula = y ~ x) )
cor.test(pd$mRNA_level_1,pd$mRNA_level_2,method = "s",exact = F)

pd <- merge(drpkm1,drpkm3, by = "gene_id")
(p11 <- p%+%pd + aes(log(mRNA_level_1),log(mRNA_level_3)) + geom_point() + 
    geom_smooth(method = "lm", se=T, color="red", formula = y ~ x) )
cor.test(pd$mRNA_level_1,pd$mRNA_level_3,method = "s",exact = F)

pd <- merge(drpkm2,drpkm3, by = "gene_id")
(p11 <- p%+%pd + aes(log(mRNA_level_2),log(mRNA_level_3)) + geom_point() + 
    geom_smooth(method = "lm", se=T, color="red", formula = y ~ x) )
cor.test(pd$mRNA_level_2,pd$mRNA_level_3,method = "s",exact = F)


# # <4> 计算 rpf ------------------------------------------------------------
getRPF <- function(data,dexonlen){
  d2 <- group_by(data,gene_id) %>% summarise(counts = n()) %>% merge(dexonlen, by = "gene_id")
  allreads <- sum(d2$counts)/1000000
  drpf <- mutate(d2, RPF = counts/(allreads*exon_length)) %>% filter(grepl("^Y",gene_id))
  return(drpf)
}
drpf1 <- getRPF(dribo1,dexonlen) %>% select(c("gene_id","RPF")) %>% rename(RPF_1 = RPF)
drpf2 <- getRPF(dribo2,dexonlen) %>% select(c("gene_id","RPF")) %>% rename(RPF_2 = RPF)
drpf3 <- getRPF(dribo3,dexonlen) %>% select(c("gene_id","RPF")) %>% rename(RPF_3 = RPF)
drpf123 <- merge(drpf1,drpf2, by = "gene_id", all = T) %>% merge(drpf3,by = "gene_id", all = T); drpf123[is.na(drpf123)] <- 0

write_csv(drpf123,"~/rstudio/mrna_looped/Hh/RPF.csv",col_names = T)
write_csv(drpf1,"~/rstudio/mrna_looped/F1/RPF.csv",col_names = T)

# # 简单看一下各种分布 -------------------------------------------------------------
# s_drpf123 <- gather(drpf123,exper,RPF,-gene_id)
# (p00 <- p%+%s_drpf123 + aes(x= factor(exper), y= log(RPF)) + geom_boxplot(notch = T,varwidth = T))
# (p00 <- p%+%s_drpf123 + aes(color = factor(exper), log(RPF)) + geom_line(stat = "density"))
# (plot(density(log(drpf1$RPF_1))))
# (lines(density(log(drpf2$RPF_2))))
# (lines(density(log(drpf3$RPF_3))))

# # 重复性比较 -----------------------------------------------------------------
pd2 <- merge(drpf1,drpf2, by = "gene_id")
(p11 <- p%+%pd2 + aes(log(RPF_1),log(RPF_2)) + geom_point() + 
    geom_smooth(method = "lm", se=T, color="red", formula = y ~ x) )
cor.test(pd2$RPF_1,pd2$RPF_2,method = "s",exact = F)

pd2 <- merge(drpf1,drpf3, by = "gene_id")
(p11 <- p%+%pd2 + aes(log(RPF_1),log(RPF_3)) + geom_point() + 
    geom_smooth(method = "lm", se=T, color="red", formula = y ~ x) )
cor.test(pd2$RPF_1,pd2$RPF_3,method = "s",exact = F)

pd2 <- merge(drpf2,drpf3, by = "gene_id")
(p11 <- p%+%pd2 + aes(log(RPF_2),log(RPF_3)) + geom_point() + 
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x) )
cor.test(pd2$RPF_2,pd2$RPF_3,method = "s",exact = F)


# #<5> 计算 TE --------------------------------------------------------------
mRNA_level <- read_csv("~/rstudio/mrna_looped/Factors/Hh/mRNA_level.csv",col_names = T) %>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)
# mRNA_level <- read_csv("~/rstudio/mrna_looped/Hh/mRNA_level.csv",col_names = T) %>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_3)/2)
mRNA_level <- read_csv("~/rstudio/mrna_looped/Factors/F1/mRNA_level.csv",col_names = T) %>% mutate(mRNA_level = mRNA_level_1)

RPF <- read_csv("~/rstudio/mrna_looped/Factors/Hh/RPF.csv",col_names = T) %>% mutate(RPF = (RPF_1+RPF_2+RPF_3)/3)
# RPF <- read_csv("~/rstudio/mrna_looped/Hh/RPF.csv",col_names = T) %>% mutate(RPF = (RPF_1+RPF_2)/2)
RPF <- read_csv("~/rstudio/mrna_looped/Factors/F1/RPF.csv",col_names = T) %>% mutate(RPF = RPF_1)

TE <- merge(select(mRNA_level,c("gene_id","mRNA_level")), select(RPF,c("gene_id","RPF")), by = "gene_id") %>% 
  mutate(TE = RPF/mRNA_level) %>% filter(TE > 0)
write_csv(TE,"~/rstudio/mrna_looped/F1/TE.csv")
#check 一下
#check 一下有没有小于0 的
plot(density(log(TE$TE)))


# #<6> 从MS 数据计算 TE ratio -----------------------------------------------------
# install.packages("openxlsx")
library("openxlsx")
d4g1 <- read.xlsx("~/rstudio/mrna_looped/protein groups2.xlsx",sheet = 1) %>% as.tibble() %>% select(1:7)
dmfs <- read.xlsx("~/rstudio/mrna_looped/protein groups2.xlsx",sheet = 2) %>% as.tibble() %>% select(1:7)

# 4g1 ---------------------------------------------------------------------
# 获得蛋白水平 ratio mean_eIF4g1/mean_HO ## 需要取倒数
str(d4g1)
# 绝对水平的
# dat_4g1 <- rename(d4g1,protein_id = Protein.IDs,eIF4g1_1=`Ratio.H/L.eIF4G1_1`,
#                   eIF4g1_2=`Ratio.H/L.eIF4G1_2`,eIF4g1_3=`Ratio.H/L.eIF4G1_3`,
#                   HO_1= `Ratio.H/L.HO_1`,HO_2= `Ratio.H/L.HO_2`,HO_3= `Ratio.H/L.HO_3`) %>%
#   filter(grepl("^Y",protein_id))
# normalized 相对水平的
dat_4g1 <- rename(d4g1,protein_id = Protein.IDs,eIF4g1_1=`Ratio.H/L.normalized.eIF4G1_1`,
                  eIF4g1_2=`Ratio.H/L.normalized.eIF4G1_2`,eIF4g1_3=`Ratio.H/L.normalized.eIF4G1_3`,
                  HO_1= `Ratio.H/L.normalized.HO_1`,HO_2= `Ratio.H/L.normalized.HO_2`,HO_3= `Ratio.H/L.normalized.HO_3`) %>%
  filter(grepl("^Y",protein_id))
str(dat_4g1)
dat_4g1$eIF4g1_1 <- as.numeric(dat_4g1$eIF4g1_1);dat_4g1$eIF4g1_2 <- as.numeric(dat_4g1$eIF4g1_2);
dat_4g1$eIF4g1_3 <- as.numeric(dat_4g1$eIF4g1_3);dat_4g1$HO_1 <- as.numeric(dat_4g1$HO_1);
dat_4g1$HO_2 <- as.numeric(dat_4g1$HO_2);dat_4g1$HO_3 <- as.numeric(dat_4g1$HO_3);

str(dat_4g1)
dat_4g1[is.na(dat_4g1)] <- 0
#先保守操作，只选择 3次重复都能测到的 gene 
dat1 <- filter(dat_4g1,(eIF4g1_1 & eIF4g1_2 & eIF4g1_3 & HO_1 &HO_2 & HO_3) != 0) %>%
  mutate(mean_eIF4g1 = (eIF4g1_1+eIF4g1_2+eIF4g1_3)/3, mean_HO= (HO_1+HO_2+HO_3)/3) %>%
  mutate(pro_ratio = 1/(mean_eIF4g1/mean_HO))

# 保守操作效果不理想，采用尽量算全部的操作
str(dat_4g1)
dat001 <- dat_4g1 %>% filter(eIF4g1_1>0 &eIF4g1_2>0 &eIF4g1_3>0) %>% mutate(mean_eIF4g1 = (eIF4g1_1+eIF4g1_2+eIF4g1_3)/3)
dat002 <- dat_4g1 %>% filter((eIF4g1_1>0 & (eIF4g1_2>0) & eIF4g1_3==0) | (eIF4g1_1>0 & (eIF4g1_2==0) & eIF4g1_3>0) | 
                               (eIF4g1_1==0 & (eIF4g1_2>0) & eIF4g1_3>0)) %>% mutate(mean_eIF4g1 = (eIF4g1_1+eIF4g1_2+eIF4g1_3)/2)
dat003 <- dat_4g1 %>% filter((eIF4g1_1>0 & (eIF4g1_2==0) & eIF4g1_3==0) | (eIF4g1_1==0 & (eIF4g1_2==0) & eIF4g1_3>0) | 
                               (eIF4g1_1==0 & (eIF4g1_2>0) & eIF4g1_3==0)) %>% mutate(mean_eIF4g1 = (eIF4g1_1+eIF4g1_2+eIF4g1_3)/1)
dat_4g1_1 <- rbind(dat001,dat002,dat003)

dat001 <- dat_4g1_1 %>% filter(HO_1>0 &HO_2>0 &HO_3>0) %>% mutate(mean_HO = (HO_1+HO_2+HO_3)/3)
dat002 <- dat_4g1_1 %>% filter((HO_1>0 & (HO_2>0) & HO_3==0) | (HO_1>0 & (HO_2==0) & HO_3>0) | 
                              (HO_1==0 & (HO_2>0) & HO_3>0)) %>% mutate(mean_HO = (HO_1+HO_2+HO_3)/2)
dat003 <- dat_4g1_1 %>% filter((HO_1>0 & (HO_2==0) & HO_3==0) | (HO_1==0 & (HO_2==0) & HO_3>0) | 
                              (HO_1==0 & (HO_2>0) & HO_3==0)) %>% mutate(mean_HO = (HO_1+HO_2+HO_3)/1)
dat_4g1_2 <- rbind(dat001,dat002,dat003)

dat1 <- dat_4g1_2 %>% mutate(pro_ratio = 1/(mean_eIF4g1/mean_HO))

write_csv(dat1,"~/rstudio/mrna_looped/MS_m4g1_mho_all_normalized.csv")

# 获得 mRNA level
mRNA_mut <- read_csv("~/rstudio/mrna_looped/Factors/m4g1/mRNA_level.csv",col_names = T) %>% mutate(mRNA_mut = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)
mRNA_wt <- read_csv("~/rstudio/mrna_looped/Factors/mho/mRNA_level.csv",col_names = T) %>% mutate(mRNA_wt = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)
mRNA_ratio <- merge(select(mRNA_mut,c("gene_id","mRNA_mut")), select(mRNA_wt,c("gene_id","mRNA_wt")), by = "gene_id") %>%
  mutate(mRNA_ratio = mRNA_mut/mRNA_wt) %>% filter(mRNA_ratio > 0,mRNA_ratio != Inf )
str(mRNA_ratio)

#计算 TEfold change
dat1 <- read_csv("~/rstudio/mrna_looped/MS_m4g1_mho_all_normalized.csv")
FC <- merge(select(mRNA_ratio,c("gene_id","mRNA_ratio")), select(dat1,c("protein_id","pro_ratio")), by.x = "gene_id",by.y = "protein_id") %>% 
  mutate(FC = pro_ratio/mRNA_ratio) %>% filter(FC > 0, FC != Inf )

# write_csv(FC,"~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho.csv")
write_csv(FC,"~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized.csv")
# 去掉那一个
FC2 <- filter(FC,gene_id != "YGR162W")
write_csv(FC2,"~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized2.csv")

# mfs ---------------------------------------------------------------------
# 获得蛋白水平 ratio mean_mfs/mean_gal1  ## 需要取倒数
str(dmfs)
# 绝对水平的
# dat_mfs <- rename(dmfs,protein_id = Protein.IDs,mfs_1=`Ratio.H/L.Fusion_1`,
#                   mfs_2=`Ratio.H/L.Fusion_2`,mfs_3=`Ratio.H/L.Fusion_3`,
#                   gal1_1= `Ratio.H/L.Gal1_1`,gal1_2= `Ratio.H/L.Gal1_2`,gal1_3= `Ratio.H/L.Gal1_3`) %>%
#   filter(grepl("^Y",protein_id))
# normalized 相对水平的
dat_mfs <- rename(dmfs,protein_id = Protein.IDs,mfs_1=`Ratio.H/L.normalized.Fusion_1`,
                  mfs_2=`Ratio.H/L.normalized.Fusion_2`,mfs_3=`Ratio.H/L.normalized.Fusion_3`,
                  gal1_1= `Ratio.H/L.normalized.Gal1_1`,gal1_2= `Ratio.H/L.normalized.Gal1_2`,gal1_3= `Ratio.H/L.normalized.Gal1_3`) %>%
  filter(grepl("^Y",protein_id))
str(dat_mfs)
dat_mfs$mfs_1 <- as.numeric(dat_mfs$mfs_1);dat_mfs$mfs_2 <- as.numeric(dat_mfs$mfs_2);
dat_mfs$mfs_3 <- as.numeric(dat_mfs$mfs_3);dat_mfs$gal1_1 <- as.numeric(dat_mfs$gal1_1);
dat_mfs$gal1_2 <- as.numeric(dat_mfs$gal1_2);dat_mfs$gal1_3 <- as.numeric(dat_mfs$gal1_3);
str(dat_mfs)
dat_mfs[is.na(dat_mfs)] <- 0

#先保守操作，只选择 3次重复都能测到的 gene 
head(dat_mfs)
dat2 <- filter(dat_mfs,(mfs_1 & mfs_2 &mfs_3 &gal1_1& gal1_2 &gal1_3) !=0) %>% 
  mutate(mean_mfs = (mfs_1+mfs_2+mfs_3)/3, mean_gal1= (gal1_1+gal1_2+gal1_3)/3) %>%
  mutate(pro_ratio = 1/(mean_mfs/mean_gal1))

# 保守操作效果不理想，采用尽量算全部的操作
str(dat_mfs)
dat001 <- dat_mfs %>% filter(mfs_1>0 &mfs_2>0 &mfs_3>0) %>% mutate(mean_mfs = (mfs_1+mfs_2+mfs_3)/3)
dat002 <- dat_mfs %>% filter((mfs_1>0 & (mfs_2>0) & mfs_3==0) | (mfs_1>0 & (mfs_2==0) & mfs_3>0) | 
                               (mfs_1==0 & (mfs_2>0) & mfs_3>0)) %>% mutate(mean_mfs = (mfs_1+mfs_2+mfs_3)/2)
dat003 <- dat_mfs %>% filter((mfs_1>0 & (mfs_2==0) & mfs_3==0) | (mfs_1==0 & (mfs_2==0) & mfs_3>0) | 
                               (mfs_1==0 & (mfs_2>0) & mfs_3==0)) %>% mutate(mean_mfs = (mfs_1+mfs_2+mfs_3)/1)
dat_mfs_1 <- rbind(dat001,dat002,dat003)

dat001 <- dat_mfs_1 %>% filter(gal1_1>0 &gal1_2>0 &gal1_3>0) %>% mutate(mean_gal1 = (gal1_1+gal1_2+gal1_3)/3)
dat002 <- dat_mfs_1 %>% filter((gal1_1>0 & (gal1_2>0) & gal1_3==0) | (gal1_1>0 & (gal1_2==0) & gal1_3>0) | 
                                 (gal1_1==0 & (gal1_2>0) & gal1_3>0)) %>% mutate(mean_gal1 = (gal1_1+gal1_2+gal1_3)/2)
dat003 <- dat_mfs_1 %>% filter((gal1_1>0 & (gal1_2==0) & gal1_3==0) | (gal1_1==0 & (gal1_2==0) & gal1_3>0) | 
                                 (gal1_1==0 & (gal1_2>0) & gal1_3==0)) %>% mutate(mean_gal1 = (gal1_1+gal1_2+gal1_3)/1)
dat_mfs_2 <- rbind(dat001,dat002,dat003)

dat2 <- dat_mfs_2 %>% mutate(pro_ratio = 1/(mean_mfs/mean_gal1))

write_csv(dat2,"~/rstudio/mrna_looped/MS_mfs_gal1_all_normalized.csv")

# 获得 mRNA level
mRNA_mut <- read_csv("~/rstudio/mrna_looped/Factors/mfs/mRNA_level.csv",col_names = T) %>% mutate(mRNA_mut = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)
mRNA_wt <- read_csv("~/rstudio/mrna_looped/Factors/gal1/mRNA_level.csv",col_names = T) %>% mutate(mRNA_wt = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)
mRNA_ratio <- merge(select(mRNA_mut,c("gene_id","mRNA_mut")), select(mRNA_wt,c("gene_id","mRNA_wt")), by = "gene_id") %>%
  mutate(mRNA_ratio = mRNA_mut/mRNA_wt) %>% filter(mRNA_ratio > 0,mRNA_ratio != Inf )
str(mRNA_ratio)
#计算 TEfold change
FC <- merge(select(mRNA_ratio,c("gene_id","mRNA_ratio")), select(dat2,c("protein_id","pro_ratio")), by.x = "gene_id", by.y = "protein_id") %>% 
  mutate(FC = pro_ratio/mRNA_ratio) %>% filter(FC > 0, FC != Inf )

# write_csv(FC,"~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1.csv")
# write_csv(FC,"~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all.csv")
write_csv(FC,"~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all_normalized.csv")




# get reads (试一下跟report -count 那里对比一下，结果是一样的，但是只看了一个样本）---------------------------------------------------------------
# check 一下单个样本的 -结果 某一个样本的 htseq count 的结果 和我手提的结果是一样的-----------------------------------------------------------------
# FS CAF20d
dmrna1 <- read_csv("~/rstudio/mrna_looped/Factors/CAF20d/last_CAF20d_1_mrna_final.csv", col_names = T )
dribo1 <- read_csv("~/rstudio/mrna_looped/Factors/CAF20d/last_CAF20d_1_ribo_final.csv", col_names = T )
check1 <- read_delim("~/mrna_looped/CAF20d/CAF20d_1_mrna/report-count-CAF20d_1_mrna.txt",delim = "\t",col_names = F) %>%
  rename(gene_id=X1,reads=X3) %>% select(gene_id,reads) %>% filter(reads>0, grepl("^Y",gene_id))
check2 <- read_delim("~/mrna_looped/CAF20d/CAF20d_1_ribo/report-count-CAF20d_1_ribo.txt",delim = "\t",col_names = F) %>%
  rename(gene_id=X1,reads=X3) %>% select(gene_id,reads) %>% filter(reads>0, grepl("^Y",gene_id))

getreads <- function(data,dexonlen){
  d2 <- group_by(data,gene_id) %>% summarise(counts = n()) %>% merge(dexonlen, by = "gene_id")
  return(d2)
}
dreads1 <- getreads(dmrna1,dexonlen) %>% select(c("gene_id","counts")) %>% rename(counts_1 = counts)
dreads2 <- getreads(dribo1,dexonlen) %>% select(c("gene_id","counts")) %>% rename(counts_2 = counts)
same <- merge(dreads1,check1,by="gene_id",all=T)
same2 <- merge(dreads2,check2,by="gene_id",all=T)
sum(dreads1$counts_1);sum(check1$reads)
sum(dreads2$counts_2);sum(check2$reads)



