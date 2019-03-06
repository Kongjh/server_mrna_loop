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

##<1> 看一下  riboseq 和 mrna 数据 ############
#<1.1> 加载 big data （只选取了 “Y。。” 基因）#####
dmrna1 <- read_csv("~/rstudio/mrna_looped/Factors/mfs/last_mfs_1_mrna_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))
dribo1 <- read_csv("~/rstudio/mrna_looped/mfs/last_mfs_1_ribo_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))

dmrna2 <- read_csv("~/rstudio/mrna_looped/Factors/mfs/last_mfs_2_mrna_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))
dribo2 <- read_csv("~/rstudio/mrna_looped/mfs/last_mfs_2_ribo_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))

dmrna3 <- read_csv("~/rstudio/mrna_looped/Factors/mfs/last_mfs_3_mrna_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))
dribo3 <- read_csv("~/rstudio/mrna_looped/mfs/last_mfs_3_ribo_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))

#<1.2> 分析 ribo mrna 周期性 长度 p5位置 #####
# n1 <- "ribo1";d1 <- dribo1
# n1 <- "mrna1";d1 <- dmrna1
# n1 <- "ribo2";d1 <- dribo2
# n1 <- "mrna2";d1 <- dmrna2
# n1 <- "ribo3";d1 <- dribo3
# n1 <- "mrna3";d1 <- dmrna3

nlist <- c("ribo1","mrna1","ribo2","mrna2","ribo3","mrna3")
dlist <- c("dribo1","dmrna1","dribo2","dmrna2","dribo3","dmrna3")

nlist <- c("mrna1","mrna2","mrna3")
dlist <- c("dmrna1","dmrna2","dmrna3")

dmrna1 <- read_csv("~/rstudio/mrna_looped/F1/last_F1_mrna_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))
dribo1 <- read_csv("~/rstudio/mrna_looped/F1/last_F1_ribo_final.csv", col_names = T ) %>% filter(grepl("^Y", gene_id))
nlist <- c("ribo1","mrna1")
dlist <- c("dribo1","dmrna1")

# for (i in seq(1,6)){
for (i in seq(1,3)){
  n1 <- nlist[i]
  d1 <- get(dlist[i])
  p_title1 <- paste(n1,"length distribution"); p_title2 <- paste(n1,"frame distribution");
  p_title3 <- paste(n1,"5' end"); p_title4 <-paste(n1,"(28nt) 5' end");p_title5 <-paste(n1," 5' end")
  ### 这里使用 facet 
  # tidyverse::tidyverse_conflicts
  print((p_rs <- p%+%d1 + 
      aes(length) + geom_bar() + 
      labs(title = p_title1)))
  print((p_rs <- p%+%d1 + 
      aes(length, fill = factor(frame)) + geom_bar(position = "dodge2") + 
      labs(title = p_title2)))
  print((p_p5_rs <- p%+%(filter(d1, length == 28)) + aes(p5) + geom_bar()+ xlim(c(-30,300))
    + labs(title = p_title4, x = "5' end position")))
  print((p_p5_rs <- p%+%d1 + aes(p5) + geom_bar()+ xlim(c(-30,300))
    + labs(title = p_title5, x = "5' end position")))
}
a

# # 以下是原版 没有循环 --------------------------------------------------------
p_title1 <- paste(n1,"length distribution"); p_title2 <- paste(n1,"frame distribution");
p_title3 <- paste(n1,"5' end"); p_title4 <-paste(n1,"(28nt) 5' end");p_title5 <-paste(n1," 5' end")
### 这里使用 facet 
# tidyverse::tidyverse_conflicts
(p_rs <- p%+%d1 + 
    aes(length) + geom_bar() + 
    labs(title = p_title1))
(p_rs <- p%+%d1 + 
    aes(length, fill = factor(frame)) + geom_bar(position = "dodge2") + 
    labs(title = p_title2))

# (p_p5_rs <- p%+%(filter(d1, id == name)) + aes(p5) + geom_bar()
#   + labs(title = p_title3, x = "5' end position"))
# (p_p5_rs <- p%+%(filter(d1, id == name)) + aes(p5) + geom_bar()+ xlim(c(-1,500))
#   + labs(title = p_title3, x = "5' end position"))
# (p_p5_rs <- p%+%(filter(d1, id == name)) + aes(p5) + geom_bar()+ xlim(c(-1,2000))
#   + labs(title = p_title3, x = "5' end position"))

# (p_p5_rs <- p%+%(filter(d1, id == name, length == 28)) + aes(p5) + geom_bar()
#   + labs(title = p_title4, x = "5' end position"))
(p_p5_rs <- p%+%(filter(d1, length == 28)) + aes(p5) + geom_bar()+ xlim(c(-30,300))
  + labs(title = p_title4, x = "5' end position"))
(p_p5_rs <- p%+%d1 + aes(p5) + geom_bar()+ xlim(c(-30,300))
  + labs(title = p_title5, x = "5' end position"))
# (p_p5_rs <- p%+%(filter(d1, id == n1)) + aes(p5) + geom_bar()+ xlim(c(-30,30))
  # + labs(title = p_title4, x = "5' end position"))
# (p_p5_rs <- p%+%(filter(d1, id == name, length == 20)) + aes(p5) + geom_bar()
#   + labs(title = "ribo-seq 20nt 5' end", x = "5' end position"))
# (p_p5_rs <- p%+%(filter(d1, id == name, length == 20)) + aes(p5) + geom_bar() + xlim(c(-1,400))
#   + labs(title = "ribo-seq 20nt 5' end", x = "5' end position"))




# 之后用get rpkm 的代码获得 rpkm --------------------------------------------------


# <2> 计算 fold change -----------------------------------------------------
# <2.1>加载 looped gene， 及得到 unlooped gene 用的是 3A+3B ----------------------------------------------------------
# 加载gene ------------------------------------------------------------------
# 这里同时获得 那篇 paper 里的 3000 来个基因
# gene <- read_csv("~/rstudio/mrna_looped/paper_gene.csv", col_names = T) %>% select("gene_id") # 2767
gene <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
g_3A <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3A.txt",delim = " ", col_names = T ) # 246
g_3B <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3B.txt",delim = " ", col_names = T ) # 149
looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395 395
unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) %>%select(-"exon_length")
looped <- mutate(looped_gene,id = "looped")
unlooped <- mutate(unlooped_gene,id = "unlooped")
loop_unloop <- rbind(looped,unlooped)
# <2.2>calculate   TE FC 并 save 运行一次就够 ------------------------------------------------------------
#get looped FC
wt_list <- c("~/rstudio/mrna_looped/Factors/KLRK/TE.csv","~/rstudio/mrna_looped/Factors/ho/TE.csv","~/rstudio/mrna_looped/Factors/Hd2/TE.csv",
             "~/rstudio/mrna_looped/Factors/ho/TE.csv","~/rstudio/mrna_looped/Factors/Hd2/TE.csv","~/rstudio/mrna_looped/Factors/ho/TE.csv",
             "~/rstudio/mrna_looped/Factors/HD/TE.csv","~/rstudio/mrna_looped/Factors/Hh/TE.csv","~/rstudio/mrna_looped/Factors/ho/TE.csv",
             "~/rstudio/mrna_looped/Factors/H1s/TE.csv","~/rstudio/mrna_looped/Factors/H1/TE.csv","~/rstudio/mrna_looped/Factors/WT/TE.csv",
             "~/rstudio/mrna_looped/Factors/WT/TE.csv","~/rstudio/mrna_looped/Factors/OE/TE.csv","~/rstudio/mrna_looped/Factors/AD/TE.csv",
             "~/rstudio/mrna_looped/Factors/ho/TE.csv")

mut_list <- c("~/rstudio/mrna_looped/Factors/4A/TE.csv","~/rstudio/mrna_looped/Factors/CAF20d/TE.csv","~/rstudio/mrna_looped/Factors/5Bd/TE.csv",
             "~/rstudio/mrna_looped/Factors/4G1d/TE.csv","~/rstudio/mrna_looped/Factors/4G1dchx/TE.csv","~/rstudio/mrna_looped/Factors/PAB1D/TE.csv",
             "~/rstudio/mrna_looped/Factors/PAB1Dchx/TE.csv","~/rstudio/mrna_looped/Factors/PAB1d2/TE.csv","~/rstudio/mrna_looped/Factors/4ED/TE.csv",
             "~/rstudio/mrna_looped/Factors/F1s/TE.csv","~/rstudio/mrna_looped/Factors/F1/TE.csv","~/rstudio/mrna_looped/Factors/FS/TE.csv",
             "~/rstudio/mrna_looped/Factors/OE/TE.csv","~/rstudio/mrna_looped/Factors/FS/TE.csv","~/rstudio/mrna_looped/Factors/HD/TE.csv",
             "~/rstudio/mrna_looped/Factors/EAP1/TE.csv")

save_list <- c("~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC.csv",
               "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC.csv",
               "~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC.csv",
               "~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC.csv")

save_list2 <- c("~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC2.csv",
                "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC2.csv",
                "~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC2.csv",
                "~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC2.csv","~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC2.csv")

# 那篇 paper 的
save_list <- c("~/rstudio/mrna_looped/paper_Fold_Change/4A_KLRK_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/CAF20d_ho_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/5Bd_Hd_final_FC.csv",
               "~/rstudio/mrna_looped/paper_Fold_Change/4G1d_ho_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/4G1dchx_Hd_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1D_ho_final_FC.csv",
               "~/rstudio/mrna_looped/paper_Fold_Change/PAB1Dchx_HD_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1d_Hh_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/4ED_ho_final_FC.csv",
               "~/rstudio/mrna_looped/paper_Fold_Change/F1s_H1s_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/F1_H1_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_WT_final_FC.csv",
               "~/rstudio/mrna_looped/paper_Fold_Change/OE_WT_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_OE_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/AD_HD_final_FC.csv",
               "~/rstudio/mrna_looped/paper_Fold_Change/EAP1_ho_final_FC.csv")

save_list2 <- c("~/rstudio/mrna_looped/paper_Fold_Change/4A_KLRK_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/CAF20d_ho_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/5Bd_Hd_final_FC2.csv",
                "~/rstudio/mrna_looped/paper_Fold_Change/4G1d_ho_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/4G1dchx_Hd_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1D_ho_final_FC2.csv",
                "~/rstudio/mrna_looped/paper_Fold_Change/PAB1Dchx_HD_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1d_Hh_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/4ED_ho_final_FC2.csv",
                "~/rstudio/mrna_looped/paper_Fold_Change/F1s_H1s_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/F1_H1_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_WT_final_FC2.csv",
                "~/rstudio/mrna_looped/paper_Fold_Change/OE_WT_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_OE_final_FC2.csv","~/rstudio/mrna_looped/paper_Fold_Change/AD_HD_final_FC2.csv",
                "~/rstudio/mrna_looped/paper_Fold_Change/EAP1_ho_final_FC2.csv")
num2 <- seq(1,11)
num2 <- seq(1,16)
num2
for (i in num2){
  wt_file <- wt_list[i]
  mut_file <- mut_list[i]
  save_FCfile <- save_list[i]
  save_FCfile2 <- save_list2[i]
  # wt_file <- "~/rstudio/mrna_looped/KLRK/TE.csv"
  # mut_file <- "~/rstudio/mrna_looped/4A/TE.csv"
  # save_FCfile <- "~/rstudio/mrna_looped/4A_KLRK_final_FC.csv"
  # save_FCfile2 <- "~/rstudio/mrna_looped/4A_KLRK_final_FC2.csv"
  #这是 3A+3B
  looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395 395
  unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) #%>%select(-"exon_length") # 6189  2372
  #length(unique(looped_gene$gene_id));length(unique(unlooped_gene$gene_id))
  
  looped_wt_TE <- read_csv(wt_file) %>% merge(looped_gene, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(wt_TE = TE)# a+b 393  a 246
  looped_mut_TE <- read_csv(mut_file) %>% merge(looped_gene, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(mut_TE = TE)# 394  246
  
  looped_FC <- merge(looped_wt_TE,looped_mut_TE,by = "gene_id") %>%
    mutate(FC = mut_TE/wt_TE) # 393 246
  looped_FC <- cbind(looped_FC,id = rep("looped",nrow(looped_FC)))
  # head(looped_FC)
  
  #get unlooped FC
  unlooped_wt_TE <- read_csv(wt_file) %>% merge(unlooped_gene, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(wt_TE = TE)# 5101
  unlooped_mut_TE <- read_csv(mut_file) %>% merge(unlooped_gene, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(mut_TE = TE)# 5135
  
  unlooped_FC <- merge(unlooped_wt_TE,unlooped_mut_TE,by = "gene_id") %>%
    mutate(FC = mut_TE/wt_TE) # 4912
  unlooped_FC <- cbind(unlooped_FC,id = rep("unlooped",nrow(unlooped_FC)))
  
  FC <- rbind(looped_FC,unlooped_FC)
  # head(FC)
  
  write_csv(FC,save_FCfile)
  
  
  # 下面是只用 3A的基因，代码与上面一样，只是 所有命名后面都有一个 2  --------------------------------------------------------------
  # <2.1>加载 looped gene， 及得到 unlooped gene 用的是只有 3A ----------------------------------------------------------
  # gene <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
  # g_3A <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3A.txt",delim = " ", col_names = T ) # 246
  # g_3B <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3B.txt",delim = " ", col_names = T ) # 149
  looped_gene2 <- g_3A$group_3A %>% as.tibble() %>% rename(gene_id = value) # 246
  unlooped_gene2 <- filter(gene, !(gene_id %in% looped_gene2$gene_id)) #%>%select(-"exon_length") # 6189 
  #length(unique(looped_gene$gene_id));length(unique(unlooped_gene$gene_id))
  # <2.2>calculate   TE FC 并 save------------------------------------------------------------
  #get looped FC
  looped_wt_TE2 <- read_csv(wt_file) %>% merge(looped_gene2, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(wt_TE2 = TE)# a+b 393  a 246
  looped_mut_TE2 <- read_csv(mut_file) %>% merge(looped_gene2, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(mut_TE2 = TE)# 394  246
  
  looped_FC2 <- merge(looped_wt_TE2,looped_mut_TE2,by = "gene_id") %>%
    mutate(FC2 = mut_TE2/wt_TE2) # 393 246
  looped_FC2 <- cbind(looped_FC2,id = rep("looped",nrow(looped_FC2)))
  # head(looped_FC2)
  
  #get unlooped FC
  unlooped_wt_TE2 <- read_csv(wt_file) %>% merge(unlooped_gene2, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(wt_TE2 = TE)# 5101
  unlooped_mut_TE2 <- read_csv(mut_file) %>% merge(unlooped_gene2, by = "gene_id") %>%select(c("gene_id","TE")) %>%
    rename(mut_TE2 = TE)# 5135
  
  unlooped_FC2 <- merge(unlooped_wt_TE2,unlooped_mut_TE2,by = "gene_id") %>%
    mutate(FC2 = mut_TE2/wt_TE2) # 4912
  unlooped_FC2 <- cbind(unlooped_FC2,id = rep("unlooped",nrow(unlooped_FC2)))
  
  FC2 <- rbind(looped_FC2,unlooped_FC2)
  # head(FC2)
  
  write_csv(FC2,save_FCfile2)
  print(i)
}


# 特例 计算MS 数据得来的 TEfc  ------------------------------------------------------
# # 在 get_rpkm_rpf.r 那个脚本那里 ;但是在这里加上id looped 和 unlooped  运行一次就够 ----------
# 获得 looped 或者 unlooped gene ----------------------------------------------
looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395 395
unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) %>%select(-"exon_length")
looped <- mutate(looped_gene,id = "looped")
unlooped <- mutate(unlooped_gene,id = "unlooped")
loop_unloop <- rbind(looped,unlooped)

# 有两次 分别是保守的 以及 _all 的 以及 all_normalized 的
ms_4g1 <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized2.csv") %>% merge(loop_unloop,by = "gene_id")
nrow(ms_4g1[ms_4g1$id == "looped",]) # 170 203
nrow(ms_4g1[ms_4g1$id == "unlooped",]) # 778 1123
write_csv(ms_4g1,"~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized2.csv")

ms_mfs <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all_normalized.csv") %>% merge(loop_unloop,by = "gene_id")
nrow(ms_mfs[ms_mfs$id == "looped",]) # 134 175
nrow(ms_mfs[ms_mfs$id == "unlooped",]) # 538 801
write_csv(ms_mfs,"~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all_normalized.csv")

# # 分析 FC 3A 以及 3S+3B -----------------------------------------------------------------
# save_list <- c("~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC.csv",
#                "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC.csv",
#                "~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC.csv",
#                "~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC.csv")


# "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC.csv"
# "~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC.csv"

# fc <- "~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC.csv"
# fc2 <- "~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC2.csv"
# 
# fc <- "~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC.csv"
# fc2 <- "~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC2.csv"
# 
# fc <- "~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC.csv"
# fc2 <- "~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC2.csv"
# 
# fc <- "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv"
# fc2 <- "~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC2.csv"

fc <- c("~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/AD_HD_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/OE_WT_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/FS_WT_final_FC.csv",
        "~/rstudio/mrna_looped/Fold_Change/EAP1_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/FS_OE_final_FC.csv")
# paper 的 gene
fc <- c("~/rstudio/mrna_looped/paper_Fold_Change/4A_KLRK_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/CAF20d_ho_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/5Bd_Hd_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/4G1d_ho_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/4G1dchx_Hd_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1D_ho_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/PAB1Dchx_HD_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/PAB1d_Hh_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/4ED_ho_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/F1s_H1s_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/F1_H1_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/AD_HD_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/OE_WT_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_WT_final_FC.csv",
        "~/rstudio/mrna_looped/paper_Fold_Change/EAP1_ho_final_FC.csv","~/rstudio/mrna_looped/paper_Fold_Change/FS_OE_final_FC.csv")

x_title <- c("log2(FC=4A_KLRK)","log2(FC=CAF20d_ho)","log2(FC=5Bd_Hd)","log2(FC=4G1d_ho)",
             "log2(FC=4G1dchx_Hd)","log2(FC=PAB1D_ho)","log2(FC=PAB1Dchx_HD)","log2(FC=PAB1d_Hh)",
             "log2(FC=4ED_ho)","log2(FC=F1s_H1s)","log2(FC=F1_H1)","log2(FC=AD_HD)",
             "log2(FC=OE_WT)","log2(FC=FS_WT)","log2(FC=EAP1_ho)","log2(FC=FS_OE)")

p_num <- c("p001","p002","p003","p004","p005","p006","p007","p008",
           "p009","p010","p011","p012","p013","p014","p015","p016")

num3 <- seq(1,16) 

### 针对一次的
fc <- c("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all.csv") 
x_title <- c("log2(FC=MS_m4g1_mho_all)")
p_num <- c("p001")
num3 <- seq(1)
FC <- read_csv(fc,col_names = T) %>% filter(FC > 0, FC != Inf)
str(FC)
mean_val <- summarise(group_by(FC,id),mean_fc=mean(log2(FC)))
### 针对2次的
fc <- c("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all.csv",
        "~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all.csv") 
x_title <- c("log2(FC=MS_m4g1_mho_all)","log2(FC=MS_mfs_gal1_all)")
p_num <- c("p001","p002")
num3 <- seq(2)
### 针对3次的
fc <- c("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized.csv",
        "~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all_normalized.csv",
        "~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized2.csv") 
x_title <- c("log2(FC=MS_m4g1_mho_all_normalized)","log2(FC=MS_mfs_gal1_all_normalized)","log2(FC=MS_m4g1_mho_all_normalized2")
p_num <- c("p001","p002","p003")
num3 <- seq(3)
num3
for (i in num3){
  FC <- read_csv(fc[i],col_names = T) %>% filter(FC > 0, FC != Inf)
  # head(FC)
  pvalue <- wilcox.test(FC[FC$id == "looped",]$FC,FC[FC$id == "unlooped",]$FC) #p-value = 0.15
  mean_val <- summarise(group_by(FC,id),mean_fc=mean(log2(FC)))
  # pvalue$p.value
  (p00 <- p%+%FC + aes(x = log2(FC), fill = factor(id),color= factor(id)) 
    + geom_density(alpha = 0.5) 
    + labs(x = x_title[i])
    + geom_vline(aes(xintercept=mean_fc,color = factor(id)),data=mean_val,size=1,lty="dashed")
    + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
               label = signif(pvalue$p.value,4) , size = 10))
  # + scale_fill_manual(values=c("#F8766D", "#00BFC4"), 
  #                     breaks=c("unlooped", "looped"),
  #                     labels=c("unlooped", "looped")) )
  assign(p_num[i],p00)
  print(i)
  # # 2表示 只有3A
  # FC2 <- read_csv(fc2,col_names = T)
  # head(FC2)
  # wilcox.test(FC2[FC2$id == "looped",]$FC2,FC2[FC2$id == "unlooped",]$FC2) #p-value = 0.15
  # (p01 <- p%+%FC2 + aes(x = log2(FC2), fill = factor(id),color= factor(id))
  #   + geom_density(alpha = 0.5)
  #   + labs(x = "log2(FC=FS/OE) 3A")
  #   + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
  #              label = "p-value =" , size = 6))
  # # + scale_fill_manual(values=c("#F8766D", "#00BFC4"),
  # #                     breaks=c("unlooped", "looped"),
  # #                     labels=c("unlooped", "looped")) )
}
print(p002)
# # 保存 --------------------------------------------------------------------
#计算了 3A+ 3B
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(16,16))) ####将页面分成2*2矩阵
# p_num <- c("p001","p002","p003","p004","p005","p006","p007","p008",
#            "p009","p010","p011","p012","p013","p014","p015","p016")

print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
print(p003, vp = vplayout(2:4,10:12)) 
print(p004, vp = vplayout(2:4,14:16)) 

print(p005, vp = vplayout(6:8,2:4))   ###将(2,1)的位置画图b
print(p006, vp = vplayout(6:8,6:8)) 
print(p007, vp = vplayout(6:8,10:12)) 
print(p008, vp = vplayout(6:8,14:16))

print(p009, vp = vplayout(10:12,2:4))   ###将(2,1)的位置画图b
print(p010, vp = vplayout(10:12,6:8)) 
print(p011, vp = vplayout(10:12,10:12)) 
print(p012, vp = vplayout(10:12,14:16))

print(p013, vp = vplayout(14:16,2:4))   ###将(2,1)的位置画图b
print(p014, vp = vplayout(14:16,6:8)) 
print(p015, vp = vplayout(14:16,10:12)) 
print(p016, vp = vplayout(14:16,14:16))

dev.off()
# 22 9 
# 64 46

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
dev.off()
# 45 20inch

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,13))) ####将页面分成2*2矩阵
print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
print(p003, vp = vplayout(2:4,10:12)) 
dev.off()
# 60 20inch


# 用MS的基因算其他factors的FC 因为MS 的结果不好 看看是不是MS 的实验做得不好 --------------------------
file <- c("~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv","~/rstudio/mrna_looped/Fold_Change/FS_WT_final_FC.csv")
x_text <- c("log2(FC=4G1d_ho_msgene)","log2(FC=FS_WT_msgene)")
p_n <- c("p100","p101")
MS_flie <- c("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all.csv","~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all.csv")
# str(MS_genelist) # 1326
for (i in seq(1,2)){
  MS_genelist <- read_csv(MS_flie[i],col_names = T) %>% select(gene_id)
  new_fc <- read_csv(file[i],col_names = T) %>% select(gene_id,FC,id) %>%
    merge(MS_genelist,by="gene_id")
  # str(new_fc)
  mean_un <- summarise(group_by(new_fc,id),mean_fc=mean(log2(FC)))
  pvalue <- wilcox.test(new_fc[new_fc$id == "looped",]$FC,new_fc[new_fc$id == "unlooped",]$FC) #p-value = 0.15
  (p00 <- p%+%new_fc + aes(x = log2(FC), fill = factor(id),color= factor(id)) 
    + geom_density(alpha = 0.5,kernel="gaussian") 
    + geom_vline(aes(xintercept=mean_fc,color = factor(id)),data=mean_un,size=1,lty="dashed")
    + labs(x = x_text[i])
    + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
               label = signif(pvalue$p.value,4) , size = 10))
  # + scale_fill_manual(values=c("#F8766D", "#00BFC4"), 
  #                     breaks=c("unlooped", "looped"),
  #                     labels=c("unlooped", "looped")) )
  assign(p_n[i],p00)
}
print(p100)
# check
check1 <- p101$data
length(check1[check1$id == "looped",]$gene_id); 
length(check1[check1$id == "unlooped",]$gene_id);

# # 保存 --------------------------------------------------------------------
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
print(p100, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p101, vp = vplayout(2:4,6:8)) 
dev.off()
# 45 20inch









# 聚类 ----------------------------------------------------------------------
# 获得 big data 运行并save 一次就行-------------------------------------------------------------
o4A <- read_csv("~/rstudio/mrna_looped/Fold_Change/4A_KLRK_final_FC.csv") %>% select(gene_id,FC,id) %>% rename(o4A = FC) %>%
  filter(o4A != Inf)
CAF20d <- read_csv("~/rstudio/mrna_looped/Fold_Change/CAF20d_ho_final_FC.csv")%>% select(gene_id,FC,id) %>% rename( CAF20d= FC)
o5Bd <- read_csv("~/rstudio/mrna_looped/Fold_Change/5Bd_Hd_final_FC.csv")%>% select(gene_id,FC,id) %>% rename( o5Bd= FC)
o4G1d <- read_csv("~/rstudio/mrna_looped/Fold_Change/4G1d_ho_final_FC.csv")%>% select(gene_id,FC,id) %>% rename( o4G1d= FC)
o4G1dchx <- read_csv("~/rstudio/mrna_looped/Fold_Change/4G1dchx_Hd_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(o4G1dchx = FC)
PAB1D <- read_csv("~/rstudio/mrna_looped/Fold_Change/PAB1D_ho_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(PAB1D = FC)
PAB1Dchx <- read_csv("~/rstudio/mrna_looped/Fold_Change/PAB1Dchx_HD_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(PAB1Dchx = FC)
PAB1d <- read_csv("~/rstudio/mrna_looped/Fold_Change/PAB1d_Hh_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(PAB1d = FC)
PAB1d$PAB1d <- as.numeric(PAB1d$PAB1d)
PAB1d <- filter(PAB1d, PAB1d > 0, PAB1d != Inf )
o4ED <- read_csv("~/rstudio/mrna_looped/Fold_Change/4ED_ho_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(o4ED = FC)
F1s <- read_csv("~/rstudio/mrna_looped/Fold_Change/F1s_H1s_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(F1s = FC)
F1 <- read_csv("~/rstudio/mrna_looped/Fold_Change/F1_H1_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(F1 = FC)
EAP1 <- read_csv("~/rstudio/mrna_looped/Fold_Change/EAP1_ho_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(EAP1 = FC)
OE <- read_csv("~/rstudio/mrna_looped/Fold_Change/OE_WT_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(OE = FC)
FS <- read_csv("~/rstudio/mrna_looped/Fold_Change/FS_WT_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(FS = FC)
AD <- read_csv("~/rstudio/mrna_looped/Fold_Change/AD_HD_final_FC.csv")%>% select(gene_id,FC,id) %>% rename(AD = FC)
# <- read_csv("")
# <- read_csv("")
head(CAF20d)
dfactors <- merge(o4A,CAF20d,by=c("gene_id","id")) %>% merge(o5Bd,by=c("gene_id","id")) %>% merge(o4G1d,by=c("gene_id","id")) %>%
  merge(o4G1dchx,by=c("gene_id","id"))%>% merge(PAB1D,by=c("gene_id","id"))%>% merge(PAB1Dchx,by=c("gene_id","id"))%>% 
  merge(PAB1d,by=c("gene_id","id"))%>% merge(o4ED,by=c("gene_id","id"))%>% merge(F1s,by=c("gene_id","id"))%>% 
  merge(F1,by=c("gene_id","id"))%>% merge(EAP1,by=c("gene_id","id"))%>% merge(OE,by=c("gene_id","id"))%>% 
  merge(FS,by=c("gene_id","id")) %>% merge(AD,by=c("gene_id","id")) 
head(dfactors)
write_csv(dfactors,"~/rstudio/mrna_looped/Fold_Change/dfactors_all_FC.csv")

# length(factor(dfactors[dfactors$id=="looped",]$id)) # 389

# 开始分析 --------------------------------------------------------------------
dfactors <- read_csv("~/rstudio/mrna_looped/Fold_Change/dfactors_all_FC.csv")

# 不加药 ---------------------------------------------------------------------
head(dfactors);colnames(dfactors)
tmp1 <- dfactors %>% select("gene_id","id","CAF20d","PAB1D","PAB1d","o4ED","o4G1d","EAP1","OE","FS")
data1 <- tmp1[,-c(1,2)]
rownames(data1) <- tmp1$gene_id
head(data1)
length(tmp1[tmp1$id == "unlooped",]$id)
length(tmp1[tmp1$id == "looped",]$id)

head(data1)
# cor2 <- rcorr(as.matrix(data1),type = "s")
# cor2
cor1 <- corr.test(data1,method = "s")
corrplot(cor1$r, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
corrplot(cor1$p, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
chart.Correlation(log2(data1),histogram = TRUE,pch=19)

pheatmap(cor1$r, cluster_row = T, fontsize_number = 16,fontsize = 14,
         display_numbers = TRUE,clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",cutree_rows = 3)
?pheatmap



# 相关性聚类 所有 factor 尝试分成两组 与成环相关和无关 -----------------------------------------------
str(dfactors)
colnames(dfactors)
# "gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","OE","FS","AD"
tmp1 <- dfactors %>% select("gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","AD")
data1 <- tmp1[,-c(1,2)]
rownames(data1) <- tmp1$gene_id
str(data1)
length(tmp1[tmp1$id == "unlooped",]$id) # 4596
length(tmp1[tmp1$id == "looped",]$id) # 389

# cor2 <- rcorr(as.matrix(data1),type = "s")
# cor2
cor1 <- corr.test(data1,method = "s")
c1 <- corrplot(cor1$r, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
c2 <- corrplot(cor1$p, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
# chart.Correlation(log2(data1),histogram = TRUE,pch=19)
p1 <- pheatmap(cor1$r, cluster_row = T, fontsize_number = 16,fontsize = 14,
         display_numbers = TRUE,clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",cutree_rows = 3)

Heatmap(cor1$r,show_row_names = F,clustering_distance_rows = "euclidean",clustering_distance_columns ="euclidean",
        row_dend_side = "right",column_dend_side = "top",column_names_side = "top",
        split = 5,show_row_dend = T,column_dend_reorder = T,
        column_dend_height = unit(1.8,"cm"),
        col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))) 

summary(log2(data1))
range(log2(data1)) #-4.6903  4.7768
bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
bk
length(bk)
colors()
# color = colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)
str(tmp1)
box1 <- tmp1 %>% select(-"id") %>% gather(key = "dat_factors",value = "fc",-"gene_id")
str(box1)
table(box1$dat_factors)
num <- seq(1,13)
(pp00 <- pp%+%box1 + aes(x= factor(dat_factors) ,y=log2(fc),fill=factor(dat_factors)) + geom_boxplot()
  + theme(axis.text.x  = element_text(size = 12,face = "bold"))
  + labs(x=NULL)
  #+ facet_wrap(~dat_factors,scales = "free", nrow = 4) 
  ) 
p2 <- pheatmap(log2(data1), cluster_row = T, fontsize_number = 16,fontsize = 14,
               display_numbers = F,clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",cutree_rows = 7,
               show_rownames = F,scale = "none",
               color = colorRampPalette(colors = c("purple","blue","white","red","yellow"))(1200))
?pheatmap
summary(data1)
summary(log2(data1))
summary(log10(data1))
summary(scale(data1))
summary(scale(log2(data1)))
###
str(data1)
cc <- scale(data1)
str(cc)
mean(data1$CAF20d) # 0.99083 1.0444
sd(data1$CAF20d) # 0.27215 0.4879
(0.83063-0.99083)/0.27215
(1.16561-1.0444)/0.4879

?scale
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")
library(ComplexHeatmap)
require(circlize)
str(data1)
?Heatmap
summary(log2(data1))
# col = colorRamp2(c(-5, 0, 5), c("green", "white", "red")))
Heatmap(log2(data1),show_row_names = F,clustering_distance_rows = "euclidean",clustering_distance_columns ="euclidean",
        row_dend_side = "right",column_dend_side = "top",column_names_side = "top",
        split = 5,show_row_dend = T,column_dend_reorder = T,
        column_dend_height = unit(1.8,"cm"),col = colorRamp2(breaks = c(-4,-2,0,2,4),colors = c("purple","blue","white","red","yellow"))) 
breaks = c(-2,-1,0,1,2),
?Heatmap
colors()
str(a)
# save
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
pushViewport(viewport(layout = grid.layout(5,5))) ####将页面分成2*2矩阵
print(p1, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p101, vp = vplayout(2:4,6:8)) 
dev.off()
# 45 20inch



# 用 TE foldchange 存在显著性的来做 聚类 分出基因 ---------------------------------------------
str(dfactors)
colnames(dfactors)
# "gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","OE","FS" 
tmp1 <- dfactors %>% select("gene_id","id","EAP1","OE","FS","CAF20d","PAB1D","o4ED","o4G1d")
data1 <- tmp1[,-c(1,2)]
rownames(data1) <- tmp1$gene_id
str(data1)
length(tmp1[tmp1$id == "unlooped",]$id) # 4596
length(tmp1[tmp1$id == "looped",]$id) # 389

# cor2 <- rcorr(as.matrix(data1),type = "s")
# cor2
cor1 <- corr.test(data1,method = "s")
corrplot(cor1$r, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
corrplot(cor1$p, type = "lower", method = "number", number.cex = 1.4,bg = "black",tl.cex = 1.4,number.digits =3)
# chart.Correlation(log2(data1),histogram = TRUE,pch=19)

p1 <- pheatmap(cor1$r, cluster_row = T, fontsize_number = 16,fontsize = 14,
               display_numbers = TRUE,clustering_distance_cols = "correlation",
               clustering_distance_rows = "correlation",cutree_rows = 3)
Heatmap(log2(data1),show_row_names = F,clustering_distance_rows = "euclidean",clustering_distance_columns ="euclidean",
        row_dend_side = "right",column_dend_side = "top",column_names_side = "top",
        split = 5,show_row_dend = T,column_dend_reorder = T,
        column_dend_height = unit(1.8,"cm"),col = colorRamp2(breaks = c(-2,-1,0,1,2),
                                                             colors = c("purple","blue","white","red","yellow"))) 

?colorRampPalette
range(log2(data1)) #-4.6903  4.7768
summary(log2(data1))
?pheatmap




# k-means 聚类 --------------------------------------------------------------
str(dfactors)
colnames(dfactors)
# "gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","OE","FS" 
tmp1 <- dfactors %>% select("gene_id","id","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1d","F1s","EAP1","OE","FS")
data1 <- tmp1[,-c(1,2)]
rownames(data1) <- tmp1$gene_id
str(data1)
length(tmp1[tmp1$id == "unlooped",]$id) # 4596
length(tmp1[tmp1$id == "looped",]$id) # 389
aa <- kmeans(log2(data1),4)
aa
plot(log2(data1$o5Bd),log2(data1$OE), col=aa$cluster,pch="*")
?kmeans

iris2<-iris[,1:4]
str(iris2)
iris.kmeans<-kmeans(iris2,3)
iris.kmeans
table(iris$Species,iris.kmeans$cluster)
plot(iris2$Sepal.Length,iris2$Sepal.Width,col=iris.kmeans$cluster,pch="*")
points(iris.kmeans$centers,pch="X",cex=1.5,col=4)


# 从我们数据获得基因 "CAF20d","PAB1D","PAB1d","o4ED","o4G1d","EAP1","OE","FS"--试试sample 1000个基因------------------------------------------------------------
sam_ple <- sample(dfactors$gene_id,1000)
head(dfactors);colnames(dfactors)
tmp2 <- dfactors %>% select("gene_id","id","CAF20d","PAB1D","PAB1d","o4ED","o4G1d","EAP1","OE","FS") %>%
  filter(gene_id %in% sam_ple)

data2 <- tmp2[,-c(1,2)]
rownames(data2) <- tmp2$gene_id
head(data2)
length(tmp2[tmp2$id == "unlooped",]$id)
length(tmp2[tmp2$id == "looped",]$id) 

c <- pheatmap(log2(data2), cluster_row = T, fontsize_number = 16,fontsize = 10,
         display_numbers = F,clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",cutree_rows = 8,show_rownames = F)
cc <- c$tree_row
cc$labels
cc$order

row_cluster <- cutree(c$tree_row,k=8)
row_cluster2 <- as.data.frame(row_cluster)
num <- group_by(row_cluster2,row_cluster) %>% summarise(count= n())

?corrplot


# overlap gene MS 和 paper 的gene -------------------------------------------
paper_gene <- read_csv("~/rstudio/mrna_looped/paper_gene.csv") %>% select(gene_id) # 2767
MS_flie1 <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all.csv") %>% select(gene_id,id) # 1326
MS_flie2 <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all.csv") %>% select(gene_id,id) # 976
s1 <- merge(paper_gene,MS_flie1,by="gene_id") # 972
s2 <- merge(paper_gene,MS_flie2,by="gene_id") # 734

length(s1[s1$id=="looped",]$id);length(s1[s1$id=="unlooped",]$id) # 203 769



# 插曲 看看有没异常值被忽略了 mRNA level -----------------------------------------------
# 加载 data 运行一次 save就行 -----------------------------------------------------
# 注意 这里 几个 factor 没用3个重复的 mrna 比如 Hh 也有一些只有一个样本 比如 F1这类
o4A <- read_csv("~/rstudio/mrna_looped/Factors/4A/mRNA_level.csv") %>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4A = mRNA_level) 
CAF20d <- read_csv("~/rstudio/mrna_looped/Factors/CAF20d/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(CAF20d = mRNA_level) 
o5Bd <- read_csv("~/rstudio/mrna_looped/Factors/5Bd/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o5Bd = mRNA_level) 
o4G1d <- read_csv("~/rstudio/mrna_looped/Factors/4G1d/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4G1d = mRNA_level) 
o4G1dchx <- read_csv("~/rstudio/mrna_looped/Factors/4G1dchx/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4G1dchx = mRNA_level) 
PAB1D <- read_csv("~/rstudio/mrna_looped/Factors/PAB1D/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1D = mRNA_level) 
PAB1Dchx <- read_csv("~/rstudio/mrna_looped/Factors/PAB1Dchx/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1Dchx = mRNA_level) 
PAB1d <- read_csv("~/rstudio/mrna_looped/Factors/PAB1d2/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1d = mRNA_level) 
# PAB1d$PAB1d <- as.numeric(PAB1d$PAB1d)
# PAB1d <- filter(PAB1d, PAB1d > 0 )
o4ED <- read_csv("~/rstudio/mrna_looped/Factors/4ED/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4ED = mRNA_level) 
AD <- read_csv("~/rstudio/mrna_looped/Factors/AD/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(AD = mRNA_level) 
HD <- read_csv("~/rstudio/mrna_looped/Factors/HD/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(HD = mRNA_level) 
Hd <- read_csv("~/rstudio/mrna_looped/Factors/Hd2/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(Hd = mRNA_level) 
ho <- read_csv("~/rstudio/mrna_looped/Factors/ho/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(ho = mRNA_level) 
wt <- read_csv("~/rstudio/mrna_looped/Factors/WT/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(wt = mRNA_level) 
Hh <- read_csv("~/rstudio/mrna_looped/Factors/Hh/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_3)/2)  %>% 
  select(gene_id,mRNA_level) %>% rename(Hh = mRNA_level) 
H1 <- read_csv("~/rstudio/mrna_looped/Factors/H1/mRNA_level.csv") %>% select(gene_id,mRNA_level_1) %>% rename(H1 = mRNA_level_1) 
H1s <- read_csv("~/rstudio/mrna_looped/Factors/H1s/mRNA_level.csv")%>% select(gene_id,mRNA_level_1) %>% rename(H1s = mRNA_level_1) 
F1 <- read_csv("~/rstudio/mrna_looped/Factors/F1/mRNA_level.csv") %>% select(gene_id,mRNA_level_1) %>% rename(F1 = mRNA_level_1) 
F1s <- read_csv("~/rstudio/mrna_looped/Factors/F1s/mRNA_level.csv")%>% select(gene_id,mRNA_level_1) %>% rename(F1s = mRNA_level_1) 
mfs <- read_csv("~/rstudio/mrna_looped/Factors/mfs/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(mfs = mRNA_level) 
gal1 <- read_csv("~/rstudio/mrna_looped/Factors/gal1/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(gal1 = mRNA_level) 
m4g1 <- read_csv("~/rstudio/mrna_looped/Factors/m4g1/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(m4g1 = mRNA_level) 
mho <- read_csv("~/rstudio/mrna_looped/Factors/mho/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(mho = mRNA_level) 
KLRK <- read_csv("~/rstudio/mrna_looped/Factors/KLRK/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(KLRK = mRNA_level) 
OE <- read_csv("~/rstudio/mrna_looped/Factors/OE/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(OE = mRNA_level) 
FS <- read_csv("~/rstudio/mrna_looped/Factors/FS/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(FS = mRNA_level) 

head(CAF20d)
dfac_mrna_t <- merge(o4A,CAF20d,by=c("gene_id"),all = T) %>% merge(o5Bd,by=c("gene_id"),all = T) %>% merge(o4G1d,by=c("gene_id"),all = T) %>%
  merge(o4G1dchx,by=c("gene_id"))%>% merge(PAB1D,by=c("gene_id"))%>% merge(PAB1Dchx,by=c("gene_id"),all = T)%>% 
  merge(PAB1d,by=c("gene_id"))%>% merge(o4ED,by=c("gene_id"))%>% merge(AD,by=c("gene_id"),all = T)%>% 
  merge(HD,by=c("gene_id"),all = T)%>% merge(Hd,by=c("gene_id"),all = T)%>% merge(ho,by=c("gene_id"),all = T)%>% 
  merge(wt,by=c("gene_id"),all = T)%>% merge(Hh,by=c("gene_id"),all = T)%>% merge(H1,by=c("gene_id"),all = T)%>% 
  merge(H1s,by=c("gene_id"),all = T)%>% merge(F1,by=c("gene_id"),all = T)%>% merge(F1s,by=c("gene_id"),all = T)%>% 
  merge(mfs,by=c("gene_id"),all = T)%>% merge(gal1,by=c("gene_id"),all = T)%>% merge(m4g1,by=c("gene_id"),all = T)%>% 
  merge(mho,by=c("gene_id"),all = T)%>% merge(KLRK,by=c("gene_id"),all = T)%>% merge(OE,by=c("gene_id"),all = T)%>% 
  merge(FS,by=c("gene_id"),all = T)
head(dfac_mrna_t)
write_csv(dfac_mrna_t,"~/rstudio/mrna_looped/Factors/dfac_mrna_allT.csv")


# 开始分析 --------------------------------------------------------------------
dfac_mrna_t <- read_csv("~/rstudio/mrna_looped/Factors/dfac_mrna_allT.csv",col_names = T)
na_dfac <- dfac_mrna_t[which(rowSums(is.na(dfac_mrna_t)) > 0 ),]

a <- select(na_dfac,gene_id,FS,wt)
a <- select(na_dfac,gene_id,OE,wt)
a <- select(na_dfac,gene_id,PAB1d,Hh) # 有几个 100 忽略了
a <- select(na_dfac,gene_id,CAF20d,ho) #有个 200 几的 
a <- select(na_dfac,gene_id,o4A,KLRK) #有个 1400 几的 

# 有空



# 插曲 引入 score值 ------------------------------------------------------------



# ORF transcript looped gene short ORF length ----------------------------------------------------------
# 加载 looped gene， 及得到 unlooped gene 用的是 3A+3B ----------------------------------------------------------
gene <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
g_3A <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3A.txt",delim = " ", col_names = T ) # 246
g_3B <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3B.txt",delim = " ", col_names = T ) # 149
#这是 3A+3B
looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395
unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) %>%select(-"exon_length") # 6189 
#length(unique(looped_gene$gene_id));length(unique(unlooped_gene$gene_id))

looped_len <- read_csv("~/rstudio/gene_exonlength.csv") %>% merge(looped_gene, by = "gene_id") %>%select(c("gene_id","exon_length")) %>%
  rename(long = exon_length) %>% mutate(id = "looped")

unlooped_len <- read_csv("~/rstudio/gene_exonlength.csv") %>% merge(unlooped_gene, by = "gene_id") %>%select(c("gene_id","exon_length")) %>%
  rename(long = exon_length) %>% mutate(id = "unlooped")

orf_FC <- rbind(looped_len,unlooped_len)

wilcox.test(log2(looped_len$long),log2(unlooped_len$long))
(pp05 <- p%+%orf_FC + aes(x = log2(long), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4,bw = 0.4) + labs(x = "log2(ORF length)")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
             label = "p-value <2e-16" , size = 6)
)



# <3> ORF length 与 TEfc 的相关性 --------------------------------------------------
# # <3.1>  in OE/WT FC ----------------------------------------------------
OE_WT <- read_csv("~/rstudio/mrna_looped/Fold_Change/OE_WT_final_FC.csv",col_names = T)
str(OE_WT)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
OE_WT <- merge(OE_WT, gene_len,by = "gene_id", all.x = T) # 5314

# 不分开，all gene 
head(OE_WT)
# 分开 这是 looped
OE_WT_loop <- OE_WT %>% filter(id == "looped")
# 分开 这是 unlooped
OE_WT_unloop <- OE_WT %>% filter(id == "unlooped")

cor.test(OE_WT$exon_length,OE_WT$FC, method = "s",exact = F)
(p006a <- p%+%OE_WT + 
   aes(x = log2(exon_length), y = log2(FC)) +
   labs(x="log2(ORF length)",y="log2(TE fold change in OE)") +
   geom_point(alpha = .4)+
   geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.45287\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

#尝试使用cut bin
?cut2
as.numeric(cut2(log2(OE_WT$exon_length),g=10)) -> OE_WT$bin_len
str(OE_WT)
table(OE_WT$bin_len)
tapply(log2(OE_WT$FC), INDEX=OE_WT$bin_len, FUN = mean) -> mean_fc_log2
bar1 <- tapply(log2(OE_WT$FC),INDEX=OE_WT$bin_len, FUN = sd)/sqrt(
  tapply(log2(OE_WT$FC),INDEX=OE_WT$bin_len, FUN = length)
)

(p006a1 <- pp + 
    aes(x = seq(1,10), y = mean_fc_log2) +
    labs(x="(ORF length)bins",y="log2(TE fold change in OE)") +
    geom_point(size=2.5)
  + geom_errorbar(aes(ymax=mean_fc_log2+bar1,ymin=mean_fc_log2-bar1),width=.1,size=.5)
  # + geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.45287\n P<2e-16" , size = 12)
  #+ scale_x_continuous(labels = seq(1,10))
)

cor.test(OE_WT_loop$exon_length,OE_WT_loop$FC, method = "s",exact = F)
(p006b <- p%+%OE_WT_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in OE)") +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.53544\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor.test(OE_WT_unloop$exon_length,OE_WT_unloop$FC, method = "s",exact = F)
(p006c <- p%+%OE_WT_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in OE)") +
    geom_point(alpha=.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.41758\n P<2e-16" , size = 12)
  
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'



# # <3.1>  in FS/WT FC ----------------------------------------------------
FS_WT <- read_csv("~/rstudio/mrna_looped/FS_WT_final_FC.csv",col_names = T)
head(FS_WT)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
FS_WT <- merge(FS_WT, gene_len,by = "gene_id", all.x = T) # 5314

# 不分开，all gene 
head(FS_WT)
# 分开 这是 looped
FS_WT_loop <- FS_WT %>% filter(id == "looped")
# 分开 这是 unlooped
FS_WT_unloop <- FS_WT %>% filter(id == "unlooped")

cor.test(FS_WT$exon_length,FS_WT$FC, method = "s",exact = F)
(p007a <- p%+%FS_WT + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS)") +
    geom_point(alpha = .4)+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.47932\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor.test(FS_WT_loop$exon_length,FS_WT_loop$FC, method = "s",exact = F)
(p007b <- p%+%FS_WT_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS)") +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.2967\n P= 2e-09" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor.test(FS_WT_unloop$exon_length,FS_WT_unloop$FC, method = "s",exact = F)
(p007c <- p%+%FS_WT_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS)") +
    geom_point(alpha=.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.45525\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# # <3.1>  in FS/OE FC ----------------------------------------------------
FS_OE <- read_csv("~/rstudio/mrna_looped/FS_OE_final_FC.csv",col_names = T)
head(FS_OE)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
FS_OE <- merge(FS_OE, gene_len,by = "gene_id", all.x = T) # 5314

# 不分开，all gene 
head(FS_OE)
# 分开 这是 looped
FS_OE_loop <- FS_OE %>% filter(id == "looped")
# 分开 这是 unlooped
FS_OE_unloop <- FS_OE %>% filter(id == "unlooped")

cor.test(FS_OE$exon_length,FS_OE$FC, method = "s",exact = F)
(p008a <- p%+%FS_OE + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS/OE)") +
    geom_point(alpha = .4)+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.20272\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor.test(FS_OE_loop$exon_length,FS_OE_loop$FC, method = "s",exact = F)
(p008b <- p%+%FS_OE_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS/OE)") +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=-0.070479\n P=0.16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor.test(FS_OE_unloop$exon_length,FS_OE_unloop$FC, method = "s",exact = F)
(p008c <- p%+%FS_OE_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in FS/OE)") +
    geom_point(alpha=0.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  # + geom_fill_
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = "rho=0.20268\n P<2e-16" , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# save --------------------------------------------------------------------
###新建画图页面
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(18,14))) ####将页面分成2*2矩阵
print(p006a, vp = vplayout(2:4,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p006b, vp = vplayout(2:4,6:8))   ###将(2,1)的位置画图b
print(p006c, vp = vplayout(2:4,10:12))  ###将（2,2)的位置画图a

print(p007a, vp = vplayout(6:8,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p007b, vp = vplayout(6:8,6:8))   ###将(2,1)的位置画图b
print(p007c, vp = vplayout(6:8,10:12))  ###将（2,2)的位置画图a

print(p008a, vp = vplayout(10:12,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p008b, vp = vplayout(10:12,6:8))   ###将(2,1)的位置画图b
print(p008c, vp = vplayout(10:12,10:12))  ###将（2,2)的位置画图a

print(pp05, vp = vplayout(14:16,2:4))   ###将（1,1)和(1,2)的位置画图c

dev.off() ##画下一幅图，记得关闭窗口
# 50 57

# color = "#00BFC4"





# #  in MS数据的 m4g1/mho FC ----------------------------------------------------
m4g1_mho <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_m4g1_mho_all_normalized.csv",col_names = T)
str(m4g1_mho)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
m4g1_mho <- merge(m4g1_mho, gene_len,by = "gene_id", all.x = T) # 5314

# 不分开，all gene 
head(m4g1_mho)
# 分开 这是 looped
m4g1_mho_loop <- m4g1_mho %>% filter(id == "looped")
# 分开 这是 unlooped
m4g1_mho_unloop <- m4g1_mho %>% filter(id == "unlooped")

m4g1_a <- cor.test(m4g1_mho$exon_length,m4g1_mho$FC, method = "s",exact = F)
m4g1_a1 <-paste(round(m4g1_a$estimate,2),signif(m4g1_a$p.value,5),sep = "  ") 
m4g1_a1

(p009a <- p%+%m4g1_mho + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in m4g1_mho)") +
    geom_point(alpha = .4)+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = m4g1_a1 , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

m4g1_b <- cor.test(m4g1_mho_loop$exon_length,m4g1_mho_loop$FC, method = "s",exact = F)
m4g1_b1 <-paste(round(m4g1_b$estimate,2),signif(m4g1_b$p.value,5),sep = "  ") 
m4g1_b1
(p009b <- p%+%m4g1_mho_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in m4g1_mho)") +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label =m4g1_b1  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

m4g1_c <- cor.test(m4g1_mho_unloop$exon_length,m4g1_mho_unloop$FC, method = "s",exact = F)
m4g1_c1 <-paste(round(m4g1_c$estimate,2),signif(m4g1_c$p.value,5),sep = "  ") 
m4g1_c1

(p009c <- p%+%m4g1_mho_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in m4g1_mho)") +
    geom_point(alpha=0.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  # + geom_fill_
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = m4g1_c1 , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# #  in MS数据的 mfs/gal1 FC ----------------------------------------------------
mfs_gal1 <- read_csv("~/rstudio/mrna_looped/Fold_Change/MS_mfs_gal1_all_normalized.csv",col_names = T)
str(mfs_gal1)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
mfs_gal1 <- merge(mfs_gal1, gene_len,by = "gene_id", all.x = T) # 5314

# 不分开，all gene 
head(mfs_gal1)
# 分开 这是 looped
mfs_gal1_loop <- mfs_gal1 %>% filter(id == "looped")
# 分开 这是 unlooped
mfs_gal1_unloop <- mfs_gal1 %>% filter(id == "unlooped")

mfs_gal1_a <- cor.test(mfs_gal1$exon_length,mfs_gal1$FC, method = "s",exact = F)
mfs_gal1_a1 <-paste(round(mfs_gal1_a$estimate,2),signif(mfs_gal1_a$p.value,5),sep = "  ") 
mfs_gal1_a1

(p010a <- p%+%mfs_gal1 + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in mfs_gal1)") +
    geom_point(alpha = .4)+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = mfs_gal1_a1 , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

mfs_gal1_b <- cor.test(mfs_gal1_loop$exon_length,mfs_gal1_loop$FC, method = "s",exact = F)
mfs_gal1_b1 <-paste(round(mfs_gal1_b$estimate,2),signif(mfs_gal1_b$p.value,5),sep = "  ") 
mfs_gal1_b1
(p010b <- p%+%mfs_gal1_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in mfs_gal1)") +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label =mfs_gal1_b1  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

mfs_gal1_c <- cor.test(mfs_gal1_unloop$exon_length,mfs_gal1_unloop$FC, method = "s",exact = F)
mfs_gal1_c1 <-paste(round(mfs_gal1_c$estimate,2),signif(mfs_gal1_c$p.value,5),sep = "  ") 
mfs_gal1_c1

(p010c <- p%+%mfs_gal1_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y="log2(TE fold change in mfs_gal1)") +
    geom_point(alpha=0.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  # + geom_fill_
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = mfs_gal1_c1 , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# MS 数据的 save --------------------------------------------------------------------
###新建画图页面
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(9,13))) ####将页面分成2*2矩阵
print(p009a, vp = vplayout(2:4,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p009b, vp = vplayout(2:4,6:8))   ###将(2,1)的位置画图b
print(p009c, vp = vplayout(2:4,10:12))  ###将（2,2)的位置画图a

print(p010a, vp = vplayout(6:8,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p010b, vp = vplayout(6:8,6:8))   ###将(2,1)的位置画图b
print(p010c, vp = vplayout(6:8,10:12))  ###将（2,2)的位置画图a

dev.off() ##画下一幅图，记得关闭窗口
# 30 18

# color = "#00BFC4"


# <4> TE in WT OE FS 只看 3A+3B 了----------------------------------------------------------
# <4.1>加载 looped gene， 及得到 unlooped gene 用的是 3A+3B ----------------------------------------------------------
gene <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
g_3A <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3A.txt",delim = " ", col_names = T ) # 246
g_3B <- read_delim("~/rstudio/mrna_looped/looped_gene/Group3B.txt",delim = " ", col_names = T ) # 149
#这是 3A+3B
looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395
unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) %>%select(-"exon_length") # 6189 
#length(unique(looped_gene$gene_id));length(unique(unlooped_gene$gene_id))

# <4.1> TE in WT
looped_WT_TE <- read_csv("~/rstudio/mrna_looped/WT/TE.csv") %>% merge(looped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "looped")
unlooped_WT_TE <- read_csv("~/rstudio/mrna_looped/WT/TE.csv") %>% merge(unlooped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "unlooped")

WT_TE <- rbind(looped_WT_TE,unlooped_WT_TE)

wilcox.test(log2(looped_WT_TE$TE),log2(unlooped_WT_TE$TE))
(pp00 <- p%+%WT_TE + aes(x = log2(TE), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4) + labs(x = "log2(TE) in WT")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
             label = "p-value <2e-16" , size = 6))

(plot(density(log2(looped_WT_TE$WT_TE)),col =2))
(lines(density(log2(unlooped_WT_TE$WT_TE)), col = 3))
legend("topleft",c("looped WT TE","unlooped WT TE"), col = 2:3, lty = 1,lwd = 3,bty = "n",cex=1.7)
# ?geom_density

# <4.2>TE in OE
looped_OE_TE <- read_csv("~/rstudio/mrna_looped/OE/TE.csv") %>% merge(looped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "looped")
unlooped_OE_TE <- read_csv("~/rstudio/mrna_looped/OE/TE.csv") %>% merge(unlooped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "unlooped")

OE_TE <- rbind(looped_OE_TE,unlooped_OE_TE)

wilcox.test(log2(looped_OE_TE$TE),log2(unlooped_OE_TE$TE))
(pp01 <- p%+%OE_TE + aes(x = log2(TE), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4) + labs(x = "log2(TE) in OE")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
             label = "p-value = 4.2e-12" , size = 6))

(plot(density(log2(looped_OE_TE$OE_TE)),col =2))
(lines(density(log2(unlooped_OE_TE$OE_TE)), col = 3))
legend("topleft",c("looped OE TE","unlooped OE TE"), col = 2:3, lty = 1,lwd = 3,bty = "n",cex=1.7)

# <4.3>TE in FS
looped_FS_TE <- read_csv("~/rstudio/mrna_looped/FS/TE.csv") %>% merge(looped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "looped")
unlooped_FS_TE <- read_csv("~/rstudio/mrna_looped/FS/TE.csv") %>% merge(unlooped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "unlooped")

FS_TE <- rbind(looped_FS_TE,unlooped_FS_TE)

wilcox.test(log2(looped_FS_TE$TE),log2(unlooped_FS_TE$TE))
(pp02 <- p%+%FS_TE + aes(x = log2(TE), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4) + labs(x = "log2(TE) in FS")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
             label = "p-value = 0.11" , size = 6))

(plot(density(log2(looped_FS_TE$FS_TE)),col =2))
(lines(density(log2(unlooped_FS_TE$FS_TE)), col = 3))
legend("topleft",c("looped FS TE","unlooped FS TE"), col = 2:3, lty = 1,lwd = 3,bty = "n",cex=1.7)

# <4.3>TE in AD
looped_AD_TE <- read_csv("~/rstudio/mrna_looped/AD/TE.csv") %>% merge(looped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "looped")
unlooped_AD_TE <- read_csv("~/rstudio/mrna_looped/AD/TE.csv") %>% merge(unlooped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "unlooped")

AD_TE <- rbind(looped_AD_TE,unlooped_AD_TE)

wilcox.test(log2(looped_AD_TE$TE),log2(unlooped_AD_TE$TE))
(pp03 <- p%+%AD_TE + aes(x = log2(TE), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4) + labs(x = "log2(TE) in AD")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
             label = "p-value <2e-16" , size = 6))

(plot(density(log2(looped_FS_TE$FS_TE)),col =2))
(lines(density(log2(unlooped_FS_TE$FS_TE)), col = 3))
legend("topleft",c("looped FS TE","unlooped FS TE"), col = 2:3, lty = 1,lwd = 3,bty = "n",cex=1.7)

# <4.3>TE in EAP1
looped_EAP1_TE <- read_csv("~/rstudio/mrna_looped/EAP1/TE.csv") %>% merge(looped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "looped")
unlooped_EAP1_TE <- read_csv("~/rstudio/mrna_looped/EAP1/TE.csv") %>% merge(unlooped_gene, by = "gene_id") %>% 
  select(c("gene_id","TE")) %>% mutate(id = "unlooped")

EAP1_TE <- rbind(looped_EAP1_TE,unlooped_EAP1_TE)

wilcox.test(log2(looped_EAP1_TE$TE),log2(unlooped_EAP1_TE$TE))
(pp04 <- p%+%EAP1_TE + aes(x = log2(TE), fill = factor(id),color= factor(id)) 
  + geom_density(alpha = 0.4) + labs(x = "log2(TE) in EAP1")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0, 
             label = "p-value =3.7e-09" , size = 6))

(plot(density(log2(looped_FS_TE$FS_TE)),col =2))
(lines(density(log2(unlooped_FS_TE$FS_TE)), col = 3))
legend("topleft",c("looped FS TE","unlooped FS TE"), col = 2:3, lty = 1,lwd = 3,bty = "n",cex=1.7)


# 保存 ----------------------------------------------------------------------
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(10,14))) ####将页面分成2*2矩阵
print(pp00, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(pp01, vp = vplayout(2:4,6:8))
print(pp02, vp = vplayout(2:4,10:12))
print(pp03, vp = vplayout(6:8,2:4))
print(pp04, vp = vplayout(6:8,6:8))   

dev.off()



## 对 orf 进行cutoff 看是否长度才是##############






# 另一种展示方式  boxplot --------------------------------------------------------
dfactors <- read_csv("~/rstudio/mrna_looped/Fold_Change/dfactors_all_FC.csv")
str(dfactors)
colnames(dfactors)
# "gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","OE","FS","AD"
tmp1 <- dfactors %>% select("gene_id","id","o4A","CAF20d","o5Bd","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4ED","F1s","F1","EAP1","OE","FS","AD")
# data1 <- tmp1[,-c(1,2)]
# rownames(data1) <- tmp1$gene_id
# str(data1)
length(tmp1[tmp1$id == "unlooped",]$id) # 4596
length(tmp1[tmp1$id == "looped",]$id) # 389
str(tmp1)
dat1 <- tmp1 %>% gather(key = "factors",value = "TE_fc",-c("gene_id","id"))
str(dat1)
order <- c("FS","OE","o4ED","o4G1d","o4G1dchx","PAB1D","PAB1Dchx","PAB1d","o4A","AD","CAF20d","EAP1","o5Bd","F1s","F1")
(p100a <- pp%+%dat1 + aes(x= factor(factors),y=log2(TE_fc),fill=factor(id))
  + geom_boxplot(notch = T)
  #+ geom_jitter(alpha = 0.3,fill=factor(id))
  + labs(x= NULL)
  + scale_x_discrete(limits=order)
  + theme(axis.text.x = element_text(size = 21,face = "bold"))
  )


# DEseq2 --------------------------------------------------------------------
library(DESeq2)
??DEseq2








# 师姐 ----------------------------------------------------------------------
先做了 FS，再换成 OE 直接 FS OE 替换，或者 fs oe
之后再把基因 换成 inter 的 重新跑一遍 命名我没改了，直接 ai 编辑里面改了（此时 还是 OE 的）
再把 oe 换成 fs
# rri new_OE TEfc ---------------------------------------------------------
fs_tefc <- read_csv("~/rstudio/mrna_looped/Fold_Change/FS_WT_final_FC.csv",col_names = T) # 5336 5405
length(fs_tefc[fs_tefc$id == "looped",]$id) # 393 393
length(fs_tefc[fs_tefc$id == "unlooped",]$id) # 4943 5012
str(fs_tefc)
library(openxlsx)
# 下面是 intra 的
rri <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIintra.xlsx",colNames = T) %>% select(rna1) %>%
  rename(gene_id = rna1)
# 下面是 inter 的
rri_1 <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) %>% select(rna1) %>%
  rename(gene_id = rna1) # 603
rri_2 <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) %>% select(rna2) %>%
  rename(gene_id = rna2) # 603
rri <- rbind(rri_1,rri_2)

# 下面开始是共同的
rri <- rri[!duplicated(rri$gene_id),]
rri <- data.frame(gene_id = rri) # 264

str(rri)

rri_fc <- merge(rri,fs_tefc,by="gene_id") # 
new_fs <- filter(fs_tefc,!(gene_id %in% rri_fc$gene_id)) # 5282 5351

write_csv(rri_fc,"~/rstudio/mrna_looped/Fold_Change/rri_fc_from_fs_inter.csv")
write_csv(new_fs,"~/rstudio/mrna_looped/Fold_Change/new_fs_inter.csv")

### 针对2次的
fc <- c("~/rstudio/mrna_looped/Fold_Change/rri_fc_from_fs_inter.csv",
        "~/rstudio/mrna_looped/Fold_Change/new_fs_inter.csv") 
x_title <- c("log2(FC=rri_fs_WT)","log2(FC=newfs_WT)")
p_num <- c("p001","p002")
num3 <- seq(2)
num3
for (i in num3){
  FC <- read_csv(fc[i],col_names = T) %>% filter(FC > 0, FC != Inf)
  # head(FC)
  pvalue <- wilcox.test(FC[FC$id == "looped",]$FC,FC[FC$id == "unlooped",]$FC) #p-value = 0.15
  mean_val <- summarise(group_by(FC,id),mean_fc=mean(log2(FC)))
  # pvalue$p.value
  (p00 <- p%+%FC + aes(x = log2(FC), fill = factor(id),color= factor(id))
    + geom_density(alpha = 0.5)
    + labs(x = x_title[i])
    + geom_vline(aes(xintercept=mean_fc,color = factor(id)),data=mean_val,size=1,lty="dashed")
    + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
               label = signif(pvalue$p.value,4) , size = 10))
  # (p00 <- pp%+%FC + aes(x= factor(id) ,y=log2(FC),fill=factor(id)) + geom_boxplot(notch = T)
  #   + labs(x=NULL)
    #+ facet_wrap(~dat_factors,scales = "free", nrow = 4) 
  # ) 
  # + scale_fill_manual(values=c("#F8766D", "#00BFC4"), 
  #                     breaks=c("unlooped", "looped"),
  #                     labels=c("unlooped", "looped")) )
  assign(p_num[i],p00)
  print(i)
  # # 2表示 只有3A
  # FC2 <- read_csv(fc2,col_names = T)
  # head(FC2)
  # wilcox.test(FC2[FC2$id == "looped",]$FC2,FC2[FC2$id == "unlooped",]$FC2) #p-value = 0.15
  # (p01 <- p%+%FC2 + aes(x = log2(FC2), fill = factor(id),color= factor(id))
  #   + geom_density(alpha = 0.5)
  #   + labs(x = "log2(FC=FS/OE) 3A")
  #   + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
  #              label = "p-value =" , size = 6))
  # # + scale_fill_manual(values=c("#F8766D", "#00BFC4"),
  # #                     breaks=c("unlooped", "looped"),
  # #                     labels=c("unlooped", "looped")) )
}
# # 保存 --------------------------------------------------------------------
#计算了 3A+ 3B
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(16,16))) ####将页面分成2*2矩阵
# p_num <- c("p001","p002","p003","p004","p005","p006","p007","p008",
#            "p009","p010","p011","p012","p013","p014","p015","p016")

print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
print(p003, vp = vplayout(2:4,10:12)) 
print(p004, vp = vplayout(2:4,14:16)) 

print(p005, vp = vplayout(6:8,2:4))   ###将(2,1)的位置画图b
print(p006, vp = vplayout(6:8,6:8)) 
print(p007, vp = vplayout(6:8,10:12)) 
print(p008, vp = vplayout(6:8,14:16))

print(p009, vp = vplayout(10:12,2:4))   ###将(2,1)的位置画图b
print(p010, vp = vplayout(10:12,6:8)) 
print(p011, vp = vplayout(10:12,10:12)) 
print(p012, vp = vplayout(10:12,14:16))

print(p013, vp = vplayout(14:16,2:4))   ###将(2,1)的位置画图b
print(p014, vp = vplayout(14:16,6:8)) 
print(p015, vp = vplayout(14:16,10:12)) 
print(p016, vp = vplayout(14:16,14:16))

dev.off()
# 22 9 
# 64 46

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
dev.off()
# 45 20inch

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,13))) ####将页面分成2*2矩阵
print(p001, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p002, vp = vplayout(2:4,6:8)) 
print(p003, vp = vplayout(2:4,10:12)) 
dev.off()
# 60 20inch


# 画 3者的 TE fc -------------------------------------------------------------
# rri_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc.csv") %>% select(-id) %>% mutate(id = "RRi")
rri_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc_from_fs_inter.csv") %>% select(-id) %>% mutate(id = "RRi")
new_fs <- read_csv("~/rstudio/mrna_looped/Fold_Change/new_fs_inter.csv")
dboth <- rbind(new_fs,rri_fc)
str(dboth)

pvalue_3 <- wilcox.test(dboth[dboth$id == "looped",]$FC,dboth[dboth$id == "unlooped",]$FC) #p-value = 0.15
pvalue_4 <- wilcox.test(dboth[dboth$id == "looped",]$FC,dboth[dboth$id == "RRi",]$FC) #p-value = 0.15
pvalue_5 <- wilcox.test(dboth[dboth$id == "unlooped",]$FC,dboth[dboth$id == "RRi",]$FC) #p-value = 0.15

sig <- paste(signif(pvalue_3$p.value,4),signif(pvalue_4$p.value,4),signif(pvalue_5$p.value,4),sep = " ")
sig
mean_val_3 <- summarise(group_by(dboth,id),mean_fc=mean(log2(FC)))
(p100a <- p%+%dboth + aes(x = log2(FC), fill = factor(id),color= factor(id))
  + geom_density(alpha = 0.4)
  + labs(x ="log2(FC=rri_newfs)" )
  + geom_vline(aes(xintercept=mean_fc,color = factor(id)),data=mean_val_3,size=1,lty="dashed")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
             label = sig, size = 10))
# 10 8 
# # ORF length 与 TEfc 的相关性 ------------------------------------------------
# 前面已经去过重了，所以这里没有去重
# in 上面的 rri --------------------------------------------------------------
# rri_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc.csv")
rri_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc_from_fs_inter.csv")
rri_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc.csv")
str(rri_fc)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
rri_fc <- merge(rri_fc, gene_len,by = "gene_id", all.x = T) # 54

# 不分开，all gene 
head(rri_fc)
# 分开 这是 looped
rri_fc_loop <- rri_fc %>% filter(id == "looped") # 30  30
# 分开 这是 unlooped
rri_fc_unloop <- rri_fc %>% filter(id == "unlooped") # 24 24
ytitle_e <- "log2(FC=rri_fs_WT)"

cor1e <- cor.test(rri_fc$exon_length,rri_fc$FC, method = "s",exact = F)
(p01e <- p%+%rri_fc + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y=ytitle_e ) +
    geom_point()+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = paste(round(cor1e$estimate,2),signif(cor1e$p.value,5),sep = "  ") , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor2e <- cor.test(rri_fc_loop$exon_length,rri_fc_loop$FC, method = "s",exact = F)
(p02e <- p%+%rri_fc_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y=ytitle_e) +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = paste(round(cor2e$estimate,2),signif(cor2e$p.value,5),sep = "  ")  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor3e <- cor.test(rri_fc_unloop$exon_length,rri_fc_unloop$FC, method = "s",exact = F)
(p03e <- p%+%rri_fc_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y=ytitle_e) +
    geom_point(alpha=0.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  # + geom_fill_
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label =paste(round(cor3e$estimate,2),signif(cor3e$p.value,5),sep = "  ")  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# in 上面的 new_fc -----------------------------------------------------------
new_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/new_fs.csv")
str(new_fc)
gene_len <- read_csv("~/rstudio/gene_exonlength.csv", col_names = T) %>% select("gene_id","exon_length") # 6584
new_fc <- merge(new_fc, gene_len,by = "gene_id", all.x = T) # 5282

# 不分开，all gene 
head(new_fc)
# 分开 这是 looped
new_fc_loop <- new_fc %>% filter(id == "looped") # 363 
# 分开 这是 unlooped
new_fc_unloop <- new_fc %>% filter(id == "unlooped") # 4919 4988
ytitle_f <- "log2(FC=newfs_WT)"

cor1f <- cor.test(new_fc$exon_length,new_fc$FC, method = "s",exact = F)
(p01f <- p%+%new_fc + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y=ytitle_f ) +
    geom_point()+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = paste(round(cor1f$estimate,2),signif(cor1f$p.value,5),sep = "  ") , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

#试一下画那个密度图
(p01f <- p%+%new_fc + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y=ytitle_f ) +
    geom_point()+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = paste(round(cor1f$estimate,2),signif(cor1f$p.value,5),sep = "  ") , size = 12)
)

cor2f <- cor.test(new_fc_loop$exon_length,new_fc_loop$FC, method = "s",exact = F)
(p02f <- p%+%new_fc_loop + 
    aes(x = log2(exon_length), y = log2(FC)) +
    labs(x="log2(ORF length)",y=ytitle_f) +
    geom_point(color="#F8766D")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label = paste(round(cor2f$estimate,2),signif(cor2f$p.value,5),sep = "  ")  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

cor3f <- cor.test(new_fc_unloop$exon_length,new_fc_unloop$FC, method = "s",exact = F)
(p03f <- p%+%new_fc_unloop + 
    aes(x = log2(exon_length), y = log2(FC),color=factor(id)) +
    labs(x="log2(ORF length)",y=ytitle_f) +
    geom_point(alpha=0.4,color="#00BFC4")+
    geom_smooth(method = "auto", se=T, color="red", formula = y ~ x)
  # + geom_fill_
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.4, vjust =1.3, 
             label =paste(round(cor3f$estimate,2),signif(cor3f$p.value,5),sep = "  ")  , size = 12)
)
#`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'


# save --------------------------------------------------------------------
###新建画图页面
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(9,13))) ####将页面分成2*2矩阵
print(p01e, vp = vplayout(2:4,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p02e, vp = vplayout(2:4,6:8))   ###将(2,1)的位置画图b
print(p03e, vp = vplayout(2:4,10:12))  ###将（2,2)的位置画图a

print(p01f, vp = vplayout(6:8,2:4))   ###将（1,1)和(1,2)的位置画图c
print(p02f, vp = vplayout(6:8,6:8))   ###将(2,1)的位置画图b
print(p03f, vp = vplayout(6:8,10:12))  ###将（2,2)的位置画图a

dev.off() ##画下一幅图，记得关闭窗口
# 50 30
# color = "#00BFC4"

# TE distribution in ho---------------------------------------------------------
looped_gene <- c(g_3A$group_3A,g_3B$group_3B) %>% as.tibble() %>% rename(gene_id = value)# 395 395
unlooped_gene <- filter(gene, !(gene_id %in% looped_gene$gene_id)) %>%select(-"exon_length")
looped <- mutate(looped_gene,id = "looped")
unlooped <- mutate(unlooped_gene,id = "unlooped")
loop_unloop <- rbind(looped,unlooped)
# rri new_ho TEfc ---------------------------------------------------------
ho_te <- read_csv("~/rstudio/mrna_looped/Factors/ho/TE.csv",col_names = T) %>%select(-c(mRNA_level,RPF))# 5442
str(ho_te)
final_ho <- merge(ho_te,loop_unloop,by="gene_id") # 5442
str(final_ho)

length(final_ho[final_ho$id == "looped",]$id) # 394
length(final_ho[final_ho$id == "unlooped",]$id) # 5048
str(final_ho)

library(openxlsx)
# 下面是 intra 的
rri <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIintra.xlsx",colNames = T) %>% select(rna1) %>%
  rename(gene_id = rna1)
# 下面是 inter 的
rri_1 <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) %>% select(rna1) %>%
  rename(gene_id = rna1) # 603
rri_2 <- read.xlsx("~/rstudio/mrna_looped/rna_rna_interation/RRIinter.xlsx",colNames = T) %>% select(rna2) %>%
  rename(gene_id = rna2) # 603
rri <- rbind(rri_1,rri_2)

# 下面开始是通用的
rri <- rri[!duplicated(rri$gene_id),]
rri <- data.frame(gene_id = rri)
str(rri)

rri_te <- merge(rri,final_ho,by="gene_id") # 54
new_ho <- filter(final_ho,!(gene_id %in% rri_te$gene_id)) # 5388

write_csv(rri_te,"~/rstudio/mrna_looped/Factors/rri_te_from_ho_inter.csv")
write_csv(new_ho,"~/rstudio/mrna_looped/Factors/new_ho_te_inter.csv")

# 画 3者的 distribution ----------------------------------------------------------
rri_te <- read_csv("~/rstudio/mrna_looped/Factors/rri_te_from_ho_inter.csv") %>% select(-id) %>% mutate(id = "RRi")# 54
new_ho <- read_csv("~/rstudio/mrna_looped/Factors/new_ho_te_inter.csv") # 5388
length(new_ho[new_ho$id == "looped",]$id) # 364
length(new_ho[new_ho$id == "unlooped",]$id) # 5024
d_both <- rbind(rri_te,new_ho) # 5405
str(d_both)

pvalue_a <- wilcox.test(d_both[d_both$id == "looped",]$TE,d_both[d_both$id == "unlooped",]$TE) #p-value = 0.15
pvalue_b <- wilcox.test(d_both[d_both$id == "looped",]$TE,d_both[d_both$id == "RRi",]$TE) #p-value = 0.15
pvalue_c <- wilcox.test(d_both[d_both$id == "unlooped",]$TE,d_both[d_both$id == "RRi",]$TE) #p-value = 0.15

sig <- paste(signif(pvalue_a$p.value,4),signif(pvalue_b$p.value,4),signif(pvalue_c$p.value,4),sep = " ")
sig
mean_val_33 <- summarise(group_by(d_both,id),mean_te=mean(log2(TE)))

(p100a <- p%+%d_both + aes(x = log2(TE), fill = factor(id),color= factor(id))
  + geom_density(alpha = 0.4)
  + labs(x ="log2(TE)" )
  + geom_vline(aes(xintercept=mean_te,color = factor(id)),data=mean_val_33,size=1,lty="dashed")
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
             label = sig, size = 10))
# 10 8







# 4 组 ---------------------------------------------------------------------

inter_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc_from_fs_inter.csv",col_names = T) %>% 
  filter(FC > 0, FC != Inf) %>% mutate(kind="inter")
intra_fc <- read_csv("~/rstudio/mrna_looped/Fold_Change/rri_fc.csv",col_names = T) %>% 
  filter(FC > 0, FC != Inf) %>% mutate(kind="intra")

da_both <- rbind(inter_fc,intra_fc)

str(da_both)
# head(FC)
pvalue1 <- wilcox.test(da_both[(da_both$id == "looped" & da_both$kind == "inter"),]$FC,
                               da_both[(da_both$id == "looped" & da_both$kind == "intra"),]$FC)  #p-value = 0.15
pvalue1

pvalue2 <- wilcox.test(da_both[(da_both$id == "unlooped" & da_both$kind == "inter"),]$FC,
                      da_both[(da_both$id == "unlooped" & da_both$kind == "intra"),]$FC)  #p-value = 0.15
pvalue2
?wilcox.test
a <- da_both[(da_both$id == "unlooped" & da_both$kind == "inter"),]

pvalue <- wilcox.test(FC[(FC$id == "looped"),]$FC,
                      FC[(FC$id == "unlooped"),]$FC)  #p-value = 0.15
pvalue
mean_val <- summarise(group_by(FC,id),mean_fc=mean(log2(FC)))

# pvalue$p.value
(p300 <- pp%+%da_both + aes(x = factor(id), y=log2(FC), fill = factor(kind),color= factor(kind))
  + geom_boxplot(notch = T,varwidth = T)
  + labs(x = NULL)
  # + geom_vline(aes(xintercept=mean_fc,color = factor(id)),data=mean_val,size=1,lty="dashed")
  # + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust =1.0,
             # label = signif(pvalue$p.value,4) , size = 10)
  )
