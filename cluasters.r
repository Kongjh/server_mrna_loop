# enviroment --------------------------------------------------------------
setwd("~/rstudio/mrna_looped/")
getwd()
library(psych)
library(corrplot)#载入两个包
library("grid")
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
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.20,"cm"),
        axis.line = element_line(size = 1.5),
        axis.text = element_text(size = 30,face = "bold"),  
        axis.title = element_text(size = 33, face = "bold"))
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
a


# Hcluster ----------------------------------------------------------------
iris <- iris
head(iris)
dim(iris)#返回行列数 150 5
idx<-sample(1:150,40)
iris3<-iris[idx,-5]
iris3
hc<-hclust(dist(iris3),method = "ave")  #注意hcluster里边传入的是dist返回值对象
plot(hc,hang=-1) 
plot(hc,hang=-1,labels=iris$Species[idx])  #这里的hang=-1使得树的节点在下方对齐
#将树分为3块
rect.hclust(hc,k=3)  
groups<-cutree(hc,k=3)


# pheatmap -----------------------------------------------------------------
install.packages("pheatmap")
library(pheatmap)

