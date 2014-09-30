setwd("genetic_corr/")###一定要修改这个目录，推荐使用Rstudio来运行
library(gtools)
source("gcorr_function.R")
read.table("TN_N34.txt",header=TRUE,na.string=".")->TN
#gcor为计算遗传相关的函数
gcor(TN,c("Block","Genotype","Environment"),c("SP","SN","PW","PN"),"Block + Genotype * Environment","Genotype")
#gcor的参数:
######## 数据集（TN） #########
#TN: 读入的数据集合，格式见附件的TN_N32.txt

######## 因变量(c("Block","Genotype","Environment")) #########
#c("Block","Genotype","Environment"):处理，即方差分析的因变量，包括各水平，其中：
#Block：为重复/区组
#Genotype：株系/品系/品种等
#Environment：为环境

######## 计算变量(c(c("SP","SN","PW","PN"))) ######
#SP,SN,PW,PN等:为计算的各个性状值trait

######## 统计模型("Block + Genotype * Environment") ######
#这里给出方差分析的模型

######## 统计模型("Genotype") ######
#Genotype:要返回的因变量的遗传相关系数