# litterature
library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)
library(psych)
attach(mtcars)

#################################################################
# Set directory/Import files 
#################################################################

setwd(paste0("/Users/guillaume/Desktop/File_Script/Fig_1_2/"))

############
# Figure 1 #
############

Gene_DELDUP_all_score = fread(paste0("Sum_IID_CNV_Gene_TYPE_LOEUF_score.tsv"), header = TRUE, sep = "\t", na.strings = c("NA", "#NULL!", ".", ""), dec = ".")
LOEUF_scores_autosomal_genes = read.csv("LOEUF_scores_autosomal_genes.tsv", header = TRUE, sep = "\t", na.strings = c("NA", "#NULL!", ".", ""), dec = ".", fileEncoding = "UTF-8-BOM")

Gene_DELDUP_all_score$Score_oe_lof_upper = ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.35,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.5,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.65,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.8,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.95,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.1,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.25,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.4,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.55,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.7,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.85,1,0)


Gene_DELDUP_all_score$Score_oe_lof_upper_plus = 
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.05,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.1,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.15,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.2,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.25,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.3,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.35,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.4,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.45,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.5,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.55,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.6,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.65,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.7,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.75,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.8,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.85,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.9,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 0.95,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.0,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.05,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.1,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.15,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.2,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.25,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.3,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.35,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.4,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.45,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.5,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.55,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.6,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.65,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.7,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.75,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.8,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.85,1,0) +
  
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.9,1,0) +
  ifelse(Gene_DELDUP_all_score$oe_lof_upper >= 1.95,1,0)


Gene_DELDUP_all_score_1 = subset(Gene_DELDUP_all_score, All_persons == 1)
Gene_DELDUP_all_score_10 = subset(Gene_DELDUP_all_score, All_persons >=30)

table_Gene_DELDUP_all_score = data.frame(cbind(table(Gene_DELDUP_all_score$Score_oe_lof_upper,Gene_DELDUP_all_score$TYPE),
                                               table(Gene_DELDUP_all_score_1$Score_oe_lof_upper,Gene_DELDUP_all_score_1$TYPE),
                                               table(Gene_DELDUP_all_score_10$Score_oe_lof_upper,Gene_DELDUP_all_score_10$TYPE),
                                               c("]0,0.35]","]0.35,0.5]","]0.5,0.65]","]0.65,0.8]","]0.8,0.95]","]0.95,1.1]","]1.1,1.25]","]1.25,1.4]","]1.4,1.55]","]1.55,1.7]","]1.7,1.85]","]1.85,2]"),
                                               c(0.175, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35, 1.5, 1.65, 1.8, 1.95)))

colnames(table_Gene_DELDUP_all_score) = c("DEL","DUP","DEL_1","DUP_1","DEL_10","DUP_10","Interval","Size")

table_Gene_DELDUP_all_score$ALL_gene_LOEUF = 0
table_Gene_DELDUP_all_score$DEL = as.numeric(table_Gene_DELDUP_all_score$DEL)
table_Gene_DELDUP_all_score$DUP =as.numeric(table_Gene_DELDUP_all_score$DUP)
table_Gene_DELDUP_all_score$DEL_1 = as.numeric(table_Gene_DELDUP_all_score$DEL_1)
table_Gene_DELDUP_all_score$DUP_1 =as.numeric(table_Gene_DELDUP_all_score$DUP_1)
table_Gene_DELDUP_all_score$DEL_10 = as.numeric(table_Gene_DELDUP_all_score$DEL_10)
table_Gene_DELDUP_all_score$DUP_10 =as.numeric(table_Gene_DELDUP_all_score$DUP_10)
table_Gene_DELDUP_all_score$Size =as.numeric(table_Gene_DELDUP_all_score$Size)

wind_c = c(0,0.35, 0.5, 0.65, 0.8, 0.95, 1.1, 1.25, 1.4, 1.55, 1.7, 1.85,2)

col_size = c(0.175, 0.450, 0.600, 0.750, 0.900, 1.050, 1.200, 1.350, 1.500, 1.650, 1.800, 1.950)

for (size in 1:length(col_size)) {
  table_Gene_DELDUP_all_score$ALL_gene_LOEUF[size] = dim(subset(LOEUF_scores_autosomal_genes, oe_lof_upper>wind_c[size] & oe_lof_upper<= wind_c[size+1]))[1]
}

table_Gene_DELDUP_all_score$DEL_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DEL / table_Gene_DELDUP_all_score$ALL_gene_LOEUF
table_Gene_DELDUP_all_score$DUP_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DUP / table_Gene_DELDUP_all_score$ALL_gene_LOEUF
table_Gene_DELDUP_all_score$DEL_1_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DEL_1 / table_Gene_DELDUP_all_score$ALL_gene_LOEUF
table_Gene_DELDUP_all_score$DUP_1_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DUP_1 / table_Gene_DELDUP_all_score$ALL_gene_LOEUF
table_Gene_DELDUP_all_score$DEL_10_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DEL_10 / table_Gene_DELDUP_all_score$ALL_gene_LOEUF
table_Gene_DELDUP_all_score$DUP_10_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score$DUP_10 / table_Gene_DELDUP_all_score$ALL_gene_LOEUF

###########################
#### Del proportion genes - step 1
table_Gene_DELDUP_all_score$DEL_on_ALL_gene_LOEUF_select = table_Gene_DELDUP_all_score$DEL_on_ALL_gene_LOEUF - table_Gene_DELDUP_all_score$DEL_1_on_ALL_gene_LOEUF - table_Gene_DELDUP_all_score$DEL_10_on_ALL_gene_LOEUF

DEL_on_ALL_gene_LOEUF = cbind(subset(table_Gene_DELDUP_all_score, select = c("Size","DEL_on_ALL_gene_LOEUF_select")), "1_DEL_on_ALL_gene_LOEUF")
DEL_1_on_ALL_gene_LOEUF_select = cbind(subset(table_Gene_DELDUP_all_score, select = c("Size","DEL_1_on_ALL_gene_LOEUF")), "0_DEL_1_on_ALL_gene_LOEUF_select")
DEL_10_on_ALL_gene_LOEUF = cbind(subset(table_Gene_DELDUP_all_score, select = c("Size","DEL_10_on_ALL_gene_LOEUF")), "2_DEL_10_on_ALL_gene_LOEUF")

names(DEL_on_ALL_gene_LOEUF) = c("Size","Freq","Freq_model")
names(DEL_1_on_ALL_gene_LOEUF_select) = c("Size","Freq","Freq_model")
names(DEL_10_on_ALL_gene_LOEUF) = c("Size","Freq","Freq_model")

DEL_on_ALL_gene_LOEUF_all = rbind(DEL_on_ALL_gene_LOEUF,DEL_1_on_ALL_gene_LOEUF_select,DEL_10_on_ALL_gene_LOEUF)

#####################

table_Gene_DELDUP_all_score$ND_DEL_gene_LOEUF = table_Gene_DELDUP_all_score$ALL_gene_LOEUF - table_Gene_DELDUP_all_score$DEL
table_Gene_DELDUP_all_score$DEL_2_9_gene_LOEUF = table_Gene_DELDUP_all_score$DEL - table_Gene_DELDUP_all_score$DEL_1 - table_Gene_DELDUP_all_score$DEL_10


ND_DEL_gene_LOEUF = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","ND_DEL_gene_LOEUF")), "DEL_0")
DEL_1 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DEL_1")), "DEL_01")
DEL_2_9 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DEL_2_9_gene_LOEUF")), "DEL_02_09")
DEL_10 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DEL_10")), "DEL_10")


names(ND_DEL_gene_LOEUF) = c("Interval","NB","Freq_model")
names(DEL_1) = c("Interval","NB","Freq_model")
names(DEL_2_9) = c("Interval","NB","Freq_model")
names(DEL_10) = c("Interval","NB","Freq_model")

DEL_on_ALL_gene_LOEUF_catg = rbind(ND_DEL_gene_LOEUF, DEL_1, DEL_2_9, DEL_10)

##################################

table_Gene_DELDUP_all_score$ND_DUP_gene_LOEUF = table_Gene_DELDUP_all_score$ALL_gene_LOEUF - table_Gene_DELDUP_all_score$DUP
table_Gene_DELDUP_all_score$DUP_2_9_gene_LOEUF = table_Gene_DELDUP_all_score$DUP - table_Gene_DELDUP_all_score$DUP_1 - table_Gene_DELDUP_all_score$DUP_10

ND_DUP_gene_LOEUF = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","ND_DUP_gene_LOEUF")), "DUP_0")
DUP_1 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DUP_1")), "DUP_01")
DUP_2_9 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DUP_2_9_gene_LOEUF")), "DUP_02_09")
DUP_10 = cbind(subset(table_Gene_DELDUP_all_score, select = c("Interval","DUP_10")), "DUP_10")


names(ND_DUP_gene_LOEUF) = c("Interval","NB","Freq_model")
names(DUP_1) = c("Interval","NB","Freq_model")
names(DUP_2_9) = c("Interval","NB","Freq_model")
names(DUP_10) = c("Interval","NB","Freq_model")

DUP_on_ALL_gene_LOEUF_catg = rbind(ND_DUP_gene_LOEUF, DUP_1, DUP_2_9, DUP_10)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

table_Gene_DELDUP_all_score_plus = data.frame(cbind(table(Gene_DELDUP_all_score$Score_oe_lof_upper_plus,Gene_DELDUP_all_score$TYPE),
                                                    table(Gene_DELDUP_all_score_1$Score_oe_lof_upper_plus,Gene_DELDUP_all_score_1$TYPE),
                                                    rbind(c(0,0),table(Gene_DELDUP_all_score_10$Score_oe_lof_upper_plus,Gene_DELDUP_all_score_10$TYPE)),
                                                    c("]0,0.05]", "]0.05,0.1]", "]0.1,0.15]", "]0.15,0.2]", "]0.2,0.25]", "]0.25,0.3]", "]0.3,0.35]", "]0.35,0.4]", "]0.4,0.45]", "]0.45,0.5]", "]0.5,0.55]", "]0.55,0.6]", "]0.6,0.65]", "]0.65,0.7]", "]0.7,0.75]", "]0.75,0.8]", "]0.8,0.85]", "]0.85,0.9]", "]0.9,0.95]", "]0.95,1]", "]1,1.05]", "]1.05,1.1]", "]1.1,1.15]", "]1.15,1.2]", "]1.2,1.25]", "]1.25,1.3]", "]1.3,1.35]", "]1.35,1.4]", "]1.4,1.45]", "]1.45,1.5]", "]1.5,1.55]", "]1.55,1.6]", "]1.6,1.65]", "]1.65,1.7]", "]1.7,1.75]", "]1.75,1.8]", "]1.8,1.85]", "]1.85,1.9]", "]1.9,1.95]", "]1.95,2]"),
                                                    c(0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975)))

Gene_DEL_all_score_only = subset(Gene_DELDUP_all_score, TYPE == "DEL")

score_40_nb = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39)
score_40_details = c("]0,0.05]", "]0.05,0.1]", "]0.1,0.15]", "]0.15,0.2]", "]0.2,0.25]", "]0.25,0.3]", "]0.3,0.35]", "]0.35,0.4]", "]0.4,0.45]", "]0.45,0.5]", "]0.5,0.55]", "]0.55,0.6]", "]0.6,0.65]", "]0.65,0.7]", "]0.7,0.75]", "]0.75,0.8]", "]0.8,0.85]", "]0.85,0.9]", "]0.9,0.95]", "]0.95,1]", "]1,1.05]", "]1.05,1.1]", "]1.1,1.15]", "]1.15,1.2]", "]1.2,1.25]", "]1.25,1.3]", "]1.3,1.35]", "]1.35,1.4]", "]1.4,1.45]", "]1.45,1.5]", "]1.5,1.55]", "]1.55,1.6]", "]1.6,1.65]", "]1.65,1.7]", "]1.7,1.75]", "]1.75,1.8]", "]1.8,1.85]", "]1.85,1.9]", "]1.9,1.95]", "]1.95,2]")

Carrier_1_DEL = c()
Carrier_2_9_DEL = c()
Carrier_10_DEL = c()

for (Nb_score_40 in 1:length(score_40_nb)) {
  Carrier_1_DEL = c(Carrier_1_DEL, sum((subset(Gene_DEL_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons == 1))$All_persons))
  Carrier_2_9_DEL = c(Carrier_2_9_DEL, sum((subset(Gene_DEL_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons > 1 & All_persons < 30))$All_persons))
  Carrier_10_DEL = c(Carrier_10_DEL, sum((subset(Gene_DEL_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons >= 30))$All_persons))
}

Carrier_1_DEL
Carrier_2_9_DEL
Carrier_10_DEL

Carrier_DEL = data.frame(cbind(c(Carrier_1_DEL,Carrier_2_9_DEL,Carrier_10_DEL),
                               c(score_40_details,score_40_details,score_40_details),
                               c("01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL", "01_Carrier_DEL",
                                 "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL", "02_09_Carrier_DEL",
                                 "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL", "10_Carrier_DEL")))

names(Carrier_DEL) = c("Carrier","Interval","Freq_model")

Carrier_DEL$Carrier = as.numeric(Carrier_DEL$Carrier)

##########################################################################################################################################################################

Gene_DUP_all_score_only = subset(Gene_DELDUP_all_score, TYPE == "DUP")

Carrier_1_DUP = c()
Carrier_2_9_DUP = c()
Carrier_10_DUP = c()

for (Nb_score_40 in 1:length(score_40_nb)) {
  Carrier_1_DUP = c(Carrier_1_DUP, sum((subset(Gene_DUP_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons == 1))$All_persons))
  Carrier_2_9_DUP = c(Carrier_2_9_DUP, sum((subset(Gene_DUP_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons > 1 & All_persons < 30))$All_persons))
  Carrier_10_DUP = c(Carrier_10_DUP, sum((subset(Gene_DUP_all_score_only, Score_oe_lof_upper_plus == score_40_nb[Nb_score_40] & All_persons >= 30))$All_persons))
}

Carrier_1_DUP
Carrier_2_9_DUP
Carrier_10_DUP

Carrier_DUP = data.frame(cbind(c(Carrier_1_DUP,Carrier_2_9_DUP,Carrier_10_DUP),
                               c(score_40_details,score_40_details,score_40_details),
                               c("01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP", "01_Carrier_DUP",
                                 "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP", "02_09_Carrier_DUP",
                                 "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP", "10_Carrier_DUP")))

names(Carrier_DUP) = c("Carrier","Interval","Freq_model")

Carrier_DUP$Carrier = as.numeric(Carrier_DUP$Carrier)

##########################################################################################################################################################################

table_Gene_DELDUP_all_score_plus_carriers = data.frame(cbind(table(Gene_DEL_all_score_only$Score_oe_lof_upper_plus, Gene_DEL_all_score_only$TYPE),
                                                             table(Gene_DEL_all_score_only$Score_oe_lof_upper_plus, Gene_DEL_all_score_only$TYPE),
                                                             table(Gene_DEL_all_score_only$Score_oe_lof_upper_plus, Gene_DEL_all_score_only$TYPE),
                                                             c("]0,0.05]", "]0.05,0.1]", "]0.1,0.15]", "]0.15,0.2]", "]0.2,0.25]", "]0.25,0.3]", "]0.3,0.35]", "]0.35,0.4]", "]0.4,0.45]", "]0.45,0.5]", "]0.5,0.55]", "]0.55,0.6]", "]0.6,0.65]", "]0.65,0.7]", "]0.7,0.75]", "]0.75,0.8]", "]0.8,0.85]", "]0.85,0.9]", "]0.9,0.95]", "]0.95,1]", "]1,1.05]", "]1.05,1.1]", "]1.1,1.15]", "]1.15,1.2]", "]1.2,1.25]", "]1.25,1.3]", "]1.3,1.35]", "]1.35,1.4]", "]1.4,1.45]", "]1.45,1.5]", "]1.5,1.55]", "]1.55,1.6]", "]1.6,1.65]", "]1.65,1.7]", "]1.7,1.75]", "]1.75,1.8]", "]1.8,1.85]", "]1.85,1.9]", "]1.9,1.95]", "]1.95,2]"),
                                                             c(0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975)))


colnames(table_Gene_DELDUP_all_score_plus) = c("DEL","DUP","DEL_1","DUP_1","DEL_10","DUP_10","Interval","Size")

table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF = 0
table_Gene_DELDUP_all_score_plus$DEL = as.numeric(table_Gene_DELDUP_all_score_plus$DEL)
table_Gene_DELDUP_all_score_plus$DUP =as.numeric(table_Gene_DELDUP_all_score_plus$DUP)
table_Gene_DELDUP_all_score_plus$DEL_1 = as.numeric(table_Gene_DELDUP_all_score_plus$DEL_1)
table_Gene_DELDUP_all_score_plus$DUP_1 =as.numeric(table_Gene_DELDUP_all_score_plus$DUP_1)
table_Gene_DELDUP_all_score_plus$DEL_10 = as.numeric(table_Gene_DELDUP_all_score_plus$DEL_10)
table_Gene_DELDUP_all_score_plus$DUP_10 =as.numeric(table_Gene_DELDUP_all_score_plus$DUP_10)
table_Gene_DELDUP_all_score_plus$Size =as.numeric(table_Gene_DELDUP_all_score_plus$Size)

wind_c = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2)

for (size in 1:40) {
  table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF[size] = dim(subset(LOEUF_scores_autosomal_genes, oe_lof_upper>wind_c[size] & oe_lof_upper<= wind_c[size+1]))[1]
}

table_Gene_DELDUP_all_score_plus$DEL_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DEL / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF
table_Gene_DELDUP_all_score_plus$DUP_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DUP / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF
table_Gene_DELDUP_all_score_plus$DEL_1_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DEL_1 / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF
table_Gene_DELDUP_all_score_plus$DUP_1_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DUP_1 / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF
table_Gene_DELDUP_all_score_plus$DEL_10_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DEL_10 / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF
table_Gene_DELDUP_all_score_plus$DUP_10_on_ALL_gene_LOEUF =  table_Gene_DELDUP_all_score_plus$DUP_10 / table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF

###########################
#### Del proportion genes - step 1
table_Gene_DELDUP_all_score_plus$DEL_on_ALL_gene_LOEUF_select = table_Gene_DELDUP_all_score_plus$DEL_on_ALL_gene_LOEUF - table_Gene_DELDUP_all_score_plus$DEL_1_on_ALL_gene_LOEUF - table_Gene_DELDUP_all_score_plus$DEL_10_on_ALL_gene_LOEUF

DEL_on_ALL_gene_LOEUF_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Size","DEL_on_ALL_gene_LOEUF_select")), "1_DEL_on_ALL_gene_LOEUF")
DEL_1_on_ALL_gene_LOEUF_select_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Size","DEL_1_on_ALL_gene_LOEUF")), "0_DEL_1_on_ALL_gene_LOEUF_select")
DEL_10_on_ALL_gene_LOEUF_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Size","DEL_10_on_ALL_gene_LOEUF")), "2_DEL_10_on_ALL_gene_LOEUF")

names(DEL_on_ALL_gene_LOEUF_plus) = c("Size","Freq","Freq_model")
names(DEL_1_on_ALL_gene_LOEUF_select_plus) = c("Size","Freq","Freq_model")
names(DEL_10_on_ALL_gene_LOEUF_plus) = c("Size","Freq","Freq_model")

DEL_on_ALL_gene_LOEUF_all_plus = rbind(DEL_on_ALL_gene_LOEUF_plus, DEL_1_on_ALL_gene_LOEUF_select_plus, DEL_10_on_ALL_gene_LOEUF_plus)

#####################

table_Gene_DELDUP_all_score_plus$ND_DEL_gene_LOEUF = table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF - table_Gene_DELDUP_all_score_plus$DEL
table_Gene_DELDUP_all_score_plus$DEL_2_9_gene_LOEUF = table_Gene_DELDUP_all_score_plus$DEL - table_Gene_DELDUP_all_score_plus$DEL_1 - table_Gene_DELDUP_all_score_plus$DEL_10
table_Gene_DELDUP_all_score_plus$DEL_1_9_gene_LOEUF = table_Gene_DELDUP_all_score_plus$DEL_1 + table_Gene_DELDUP_all_score_plus$DEL_2_9_gene_LOEUF


ND_DEL_gene_LOEUF_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","ND_DEL_gene_LOEUF")), "DEL_0")
DEL_1_9_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","DEL_1_9_gene_LOEUF")), "DEL_01")
DEL_10_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","DEL_10")), "DEL_10")


names(ND_DEL_gene_LOEUF_plus) = c("Interval","NB","Freq_model")
names(DEL_1_9_plus) = c("Interval","NB","Freq_model")
names(DEL_10_plus) = c("Interval","NB","Freq_model")

DEL_on_ALL_gene_LOEUF_catg_plus = rbind(ND_DEL_gene_LOEUF_plus, DEL_1_9_plus, DEL_10_plus)


Sum_DEL_0_ND_gene_LOEUF_plus = sum(table_Gene_DELDUP_all_score_plus$ND_DEL_gene_LOEUF)
Sum_DEL_01_09_plus = sum(table_Gene_DELDUP_all_score_plus$DEL_1_9_gene_LOEUF)
Sum_DEL_10_plus = sum(table_Gene_DELDUP_all_score_plus$DEL_10)
Sum_DEL_all = sum(table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF)

Sum_DEL_all_catg = rbind(Sum_DEL_0_ND_gene_LOEUF_plus,
                         Sum_DEL_01_09_plus,
                         Sum_DEL_10_plus)

names <- rownames(Sum_DEL_all_catg)
rownames(Sum_DEL_all_catg) <- NULL
Sum_DEL_all_catg <- data.frame(cbind(names,Sum_DEL_all_catg))
Sum_DEL_all_catg$V2 = as.numeric(Sum_DEL_all_catg$V2)
Sum_DEL_all_catg$Freq_model = 100 * (Sum_DEL_all_catg$V2 / Sum_DEL_all)
names(Sum_DEL_all_catg) = c("Name","NB","Freq_model")
Sum_DEL_all_catg$Interval = "All genes with LOEUF"

#####################

table_Gene_DELDUP_all_score_plus$ND_DUP_gene_LOEUF = table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF - table_Gene_DELDUP_all_score_plus$DUP
table_Gene_DELDUP_all_score_plus$DUP_2_9_gene_LOEUF = table_Gene_DELDUP_all_score_plus$DUP - table_Gene_DELDUP_all_score_plus$DUP_1 - table_Gene_DELDUP_all_score_plus$DUP_10
table_Gene_DELDUP_all_score_plus$DUP_1_9_gene_LOEUF = table_Gene_DELDUP_all_score_plus$DUP_1 + table_Gene_DELDUP_all_score_plus$DUP_2_9_gene_LOEUF


ND_DUP_gene_LOEUF_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","ND_DUP_gene_LOEUF")), "DUP_0")
DUP_1_9_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","DUP_1_9_gene_LOEUF")), "DUP_01")
DUP_10_plus = cbind(subset(table_Gene_DELDUP_all_score_plus, select = c("Interval","DUP_10")), "DUP_10")


names(ND_DUP_gene_LOEUF_plus) = c("Interval","NB","Freq_model")
names(DUP_1_9_plus) = c("Interval","NB","Freq_model")
names(DUP_10_plus) = c("Interval","NB","Freq_model")

DUP_on_ALL_gene_LOEUF_catg_plus = rbind(ND_DUP_gene_LOEUF_plus, DUP_1_9_plus, DUP_10_plus)


Sum_DUP_0_ND_gene_LOEUF_plus = sum(table_Gene_DELDUP_all_score_plus$ND_DUP_gene_LOEUF)
Sum_DUP_01_09_plus = sum(table_Gene_DELDUP_all_score_plus$DUP_1_9_gene_LOEUF)
Sum_DUP_10_plus = sum(table_Gene_DELDUP_all_score_plus$DUP_10)
Sum_DUP_all = sum(table_Gene_DELDUP_all_score_plus$ALL_gene_LOEUF)

Sum_DUP_all_catg = rbind(Sum_DUP_0_ND_gene_LOEUF_plus,
                         Sum_DUP_01_09_plus,
                         Sum_DUP_10_plus)


names <- rownames(Sum_DUP_all_catg)
rownames(Sum_DUP_all_catg) <- NULL
Sum_DUP_all_catg <- data.frame(cbind(names,Sum_DUP_all_catg))
Sum_DUP_all_catg$V2 = as.numeric(Sum_DUP_all_catg$V2)
Sum_DUP_all_catg$Freq_model = 100 * (Sum_DUP_all_catg$V2 / Sum_DUP_all)
names(Sum_DUP_all_catg) = c("Name","NB","Freq_model")
Sum_DUP_all_catg$Interval = "All genes with LOEUF"

##################################

#### Figure 1 del et dup

LOEUF_DEL_1_detail = ggplot(data=DEL_on_ALL_gene_LOEUF_catg_plus, aes(x=Interval, y=NB, fill=Freq_model)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("grey","red", "darkred")) + theme(legend.position = c(0.8, 0.8)) +# theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

LOEUF_DEL_1_carrier = ggplot(data=Carrier_DEL, aes(x=Interval, y=Carrier, fill=Freq_model)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("red", "red", "darkred")) + theme(legend.position = c(0.8, 0.8)) +# theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))


LOEUF_DUP_1_detail = ggplot(data=DUP_on_ALL_gene_LOEUF_catg_plus, aes(x=Interval, y=NB, fill=Freq_model)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("grey","dodgerblue", "darkblue")) + theme(legend.position = c(0.8, 0.8)) +# theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

LOEUF_DUP_1_carrier = ggplot(data=Carrier_DUP, aes(x=Interval, y=Carrier, fill=Freq_model)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("dodgerblue","dodgerblue", "darkblue")) + theme(legend.position = c(0.8, 0.8)) +# theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))




LOEUF_DEL_1_detail_Ok = ggplot(data=DEL_on_ALL_gene_LOEUF_catg, aes(x=Interval, y=NB, fill=Freq_model)) +
  ylab("Numbre of unique gene") + xlab("Interval of LOEUF") +
  geom_bar(stat="identity") + scale_fill_manual(values=c("grey","red", "red", "darkred")) + theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

LOEUF_DEL_1_carrier_Ok = ggplot(data=Carrier_DEL, aes(x=Interval, y=Carrier, fill=Freq_model)) + ylim(0,7000) +
  ylab("Numbre of gene deleted") + xlab("Interval of LOEUF") +
  geom_bar(stat="identity") + scale_fill_manual(values=c("red", "red", "darkred")) + theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

LOEUF_DUP_1_carrier_Ok = ggplot(data=Carrier_DUP, aes(x=Interval, y=Carrier, fill=Freq_model)) + ylim(0,7000) +
  ylab("Numbre of gene deleted") + xlab("Interval of LOEUF") +
  geom_bar(stat="identity") + scale_fill_manual(values=c("dodgerblue", "dodgerblue", "darkblue")) + theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

############################################################

DEL_on_ALL_gene_LOEUF_catg$TYPE = "Deleted gene"
DUP_on_ALL_gene_LOEUF_catg$TYPE = "Duplicated gene"

DELDUP_on_ALL_gene_LOEUF_catg = rbind(DEL_on_ALL_gene_LOEUF_catg,DUP_on_ALL_gene_LOEUF_catg)

DEL_on_ALL_gene_LOEUF_catg_plus$TYPE = "Deleted gene"
DUP_on_ALL_gene_LOEUF_catg_plus$TYPE = "Duplicated gene"

DELDUP_on_ALL_gene_LOEUF_catg = rbind(DEL_on_ALL_gene_LOEUF_catg_plus,DUP_on_ALL_gene_LOEUF_catg_plus)

DELDUP_on_ALL_gene_LOEUF_catg$Categories = ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DEL_0" , "No carrier",
                                                  ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DEL_01" , "Number of deletion < 10",
                                                         ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DEL_02_09" , "Number of deletion < 10",
                                                                ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DEL_10" , "Number of deletion >= 10",
                                                                       ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DUP_0", "No carrier",
                                                                              ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DUP_01", "Number of duplication < 10",
                                                                                     ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DUP_02_09", "Number of duplication < 10",
                                                                                            ifelse(DELDUP_on_ALL_gene_LOEUF_catg$Freq_model == "DUP_10", "Number of duplication >= 10","ND"))))))))


Fig_DELDUP_on_ALL_gene_LOEUF_catg = ggplot(data=DELDUP_on_ALL_gene_LOEUF_catg, aes(x=TYPE, y=NB, fill=Categories)) + geom_col() +
  facet_grid(~Interval) + ylab("Numbre of unique gene") + xlab("Interval of LOEUF") + 
  theme(legend.position='none', axis.text.x=element_blank() , axis.ticks.x=element_blank()) + 
  scale_fill_manual(values=c("grey","red", "darkred", "dodgerblue", "darkblue"))

############################################################

Carrier_DEL$TYPE = "Deleted gene"
Carrier_DUP$TYPE = "Duplicated gene"

Carrier_DELDUP = rbind(Carrier_DEL,Carrier_DUP)

Carrier_DELDUP$Categories = ifelse(Carrier_DELDUP$Freq_model == "01_Carrier_DEL" , "Carriers of deletion < 10",
                                   ifelse(Carrier_DELDUP$Freq_model == "02_09_Carrier_DEL" , "Carriers of deletion < 10",
                                          ifelse(Carrier_DELDUP$Freq_model == "10_Carrier_DEL" , "Carriers of deletion >= 10",
                                                 ifelse(Carrier_DELDUP$Freq_model == "01_Carrier_DUP", "Carriers of duplication < 10",
                                                        ifelse(Carrier_DELDUP$Freq_model == "02_09_Carrier_DUP", "Carriers of duplication < 10",
                                                               ifelse(Carrier_DELDUP$Freq_model == "10_Carrier_DUP", "Carriers of duplication >= 10","ND"))))))




LOEUF_DELDUP_1_carrier_Ok = ggplot(data=Carrier_DELDUP, aes(x=TYPE, y=Carrier, fill=Categories)) + geom_col() +
  facet_grid(~Interval) + ylab("Numbre of carriers") + xlab("Interval of LOEUF") + 
  theme(legend.position='none', axis.text.x=element_blank() , axis.ticks.x=element_blank()) + 
  scale_fill_manual(values=c("red", "darkred", "dodgerblue", "darkblue"))

##############################

Sum_DEL_all_catg$TYPE = "Del."
Sum_DUP_all_catg$TYPE = "Dup."

Sum_DELDUP_all_catg = rbind(Sum_DEL_all_catg,Sum_DUP_all_catg)

Sum_DELDUP_all_catg$Categories = c("No carrier","Number of deletion < 10","Number of deletion >= 10",
                                   "No carrier","Number of duplication < 10","Number of duplication >= 10")

Sum_DEL_all_catg_size <- ddply(subset(Sum_DELDUP_all_catg, TYPE == "Del."), .(Interval), transform, pos = 100 - (cumsum(Freq_model) - (0.5 * Freq_model)))
Sum_DUP_all_catg_size <- ddply(subset(Sum_DELDUP_all_catg, TYPE == "Dup."), .(Interval), transform, pos = 100 - (cumsum(Freq_model) - (0.5 * Freq_model)))

Sum_DELDUP_all_catg = rbind(Sum_DEL_all_catg_size,Sum_DUP_all_catg_size)


LOEUF_DELDUP_1 = ggplot() + 
  geom_bar(aes(x=TYPE, y=Freq_model, fill=Categories) , data = Sum_DELDUP_all_catg , stat="identity") + 
  geom_label(data=Sum_DELDUP_all_catg, aes(x=TYPE, y = pos, label = paste0(round(Freq_model, digits=1),"%"))) +
  scale_fill_manual(values=c("grey","coral1", "darkred", "cornflowerblue", "darkblue")) + ylab("Percentage (%)")  + theme_minimal() + theme(legend.position='none', axis.title.x=element_blank())

##############################

Sum_carriers = c(sum(subset(Carrier_DELDUP,Carrier_DELDUP$Categories == "Carriers of deletion < 10")$Carrier),
                 sum(subset(Carrier_DELDUP,Carrier_DELDUP$Categories == "Carriers of deletion >= 10")$Carrier),
                 sum(subset(Carrier_DELDUP,Carrier_DELDUP$Categories == "Carriers of duplication < 10")$Carrier),
                 sum(subset(Carrier_DELDUP,Carrier_DELDUP$Categories == "Carriers of duplication >= 10")$Carrier))

Sum_carriers_TYPE = c("Del.","Del.","Dup.","Dup.")

Sum_carriers_all_type = c(sum(subset(Carrier_DELDUP,TYPE == "Deleted gene")$Carrier),sum(subset(Carrier_DELDUP,TYPE == "Deleted gene")$Carrier),
                          sum(subset(Carrier_DELDUP,TYPE == "Duplicated gene")$Carrier),sum(subset(Carrier_DELDUP,TYPE == "Duplicated gene")$Carrier))

Carrier_DELDUP_all_catg_perentage = data.frame(cbind(Sum_carriers,Sum_carriers_TYPE,Sum_carriers_all_type))

names(Carrier_DELDUP_all_catg_perentage) = c("NB","TYPE","total")

Carrier_DELDUP_all_catg_perentage$NB = as.numeric(Carrier_DELDUP_all_catg_perentage$NB)
Carrier_DELDUP_all_catg_perentage$total = as.numeric(Carrier_DELDUP_all_catg_perentage$total)

Carrier_DELDUP_all_catg_perentage$Freq_model = 100 * (Carrier_DELDUP_all_catg_perentage$NB / Carrier_DELDUP_all_catg_perentage$total)

Carrier_DELDUP_all_catg_perentage <- ddply(Carrier_DELDUP_all_catg_perentage, .(TYPE), transform, pos = 100 - (cumsum(Freq_model) - (0.5 * Freq_model)))

Carrier_DELDUP_all_catg_perentage$Categories = c("Number of deletion < 10","Number of deletion >= 10",
                                                 "Number of duplication < 10","Number of duplication >= 10")

LOEUF_DELDUP_1_carrier = ggplot() + 
  geom_bar(aes(x=TYPE, y=Freq_model, fill=Categories) , data = Carrier_DELDUP_all_catg_perentage , stat="identity") +
  geom_label(data=Carrier_DELDUP_all_catg_perentage, aes(x=TYPE, y = pos, label = paste0(round(Freq_model, digits=1),"%"))) +
  scale_fill_manual(values=c("coral1", "darkred", "cornflowerblue", "darkblue")) + ylab("Percentage (%)")  + theme_minimal() + theme(legend.position='none', axis.title.x=element_blank()) 

LOEUF_DUP_1_detail_Ok = ggplot(data=DUP_on_ALL_gene_LOEUF_catg, aes(x=Interval, y=NB, fill=Freq_model)) +
  ylab("Numbre of unique gene") + xlab("Interval of LOEUF") +
  geom_bar(stat="identity") + scale_fill_manual(values=c("grey","blue", "blue", "darkblue")) + theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

LOEUF_DUP_1_carrier_Ok = ggplot(data=Carrier_DUP, aes(x=Interval, y=Carrier, fill=Freq_model)) + ylim(0,16000) +
  ylab("Numbre of gene duplicated") + xlab("Interval of LOEUF") +
  geom_bar(stat="identity") + scale_fill_manual(values=c("blue", "blue", "darkblue")) + theme_classic() +
  theme(legend.position='none',axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))


Sum_DUP_all_catg <- ddply(Sum_DUP_all_catg, .(Interval), transform, pos = 100 -(cumsum(Freq_model) - (0.5 * Freq_model)))


fill <- c("grey","blue", "dodgerblue", "darkblue")

LOEUF_DUP_1 = ggplot() + 
  geom_bar(aes(x=Interval, y=Freq_model, fill=Name) , data = Sum_DUP_all_catg , stat="identity") + #scale_fill_manual(values=c("grey","red", "coral1", "darkred")) +
  geom_label(data=Sum_DUP_all_catg, aes(x=Interval, y = pos, label = paste0(round(Freq_model, digits=1),"%"))) +
  scale_fill_manual(values=fill) + ylab("Percentage (%)")  + theme_minimal() + theme(legend.position='none') #+ theme_minimal(base_size = 24)

##########################################################################################################################################################################
##########################################################################################################################################################################
#library("ggVennDiagram")
library(ggvenn)
library("eulerr")

DEL_1_9 = as.character((subset(Gene_DELDUP_all_score, TYPE == "DEL" & All_persons <30))$gene_id)
DUP_1_9 = as.character((subset(Gene_DELDUP_all_score, TYPE == "DUP" & All_persons <30))$gene_id)

DEL_10 = as.character((subset(Gene_DELDUP_all_score, TYPE == "DEL" & All_persons >=30))$gene_id)
DUP_10 = as.character((subset(Gene_DELDUP_all_score, TYPE == "DUP" & All_persons >=30))$gene_id)

DEL_DUP_venn = list(
  A = DEL_1_9, 
  B = DUP_1_9, 
  C = DEL_10,
  D = DUP_10
)


# write.csv(DEL_1_9,"/home/guillaf/projects/rrg-jacquese/guillaf/DEL_DUP_Intell_article_all_Cohorts/Figure_Table_Article_and_Supp/File_del_less_30.txt", row.names = F)
# write.csv(DUP_1_9,"/home/guillaf/projects/rrg-jacquese/guillaf/DEL_DUP_Intell_article_all_Cohorts/Figure_Table_Article_and_Supp/File_dup_less_30.txt" , row.names = F)
# write.csv(DEL_10,"/home/guillaf/projects/rrg-jacquese/guillaf/DEL_DUP_Intell_article_all_Cohorts/Figure_Table_Article_and_Supp/File_del_more_30.txt" , row.names = F)
# write.csv(DUP_10,"/home/guillaf/projects/rrg-jacquese/guillaf/DEL_DUP_Intell_article_all_Cohorts/Figure_Table_Article_and_Supp/File_dup_more_30.txt" , row.names = F)

ggvenn(
  DEL_DUP_venn, 
  fill_color = c("coral1", "darkred", "cornflowerblue", "darkblue"),
  stroke_size = 0.5, set_name_size = 4
)


col <- c("coral1", "cornflowerblue", "darkred","darkblue")

# 30 carriers
Venn_del_dup = plot(euler(c("A" = 1294, "B" = 6843, "C" = 18, "D" = 108,
                            "A&B" = 4852,"B&C" = 90, "C&D" = 170,
                            "A&B&C" = 0, "B&C&D" = 0,
                            "A&B&C&D" = 0,
                            "A&C" = 0,"A&C&D" = 0,"A&B&D" = 0,"B&D" = 0,
                            "A&D" = 372), shape = "ellipse"),             
                    key = TRUE, counts = TRUE,
                    quantities = list(type = c("counts", "percent"), font=1, round=2, col = "white"),
                    fills = list(fill = col, alpha = 0.9),
                    edges=list(lty = 0),
                    factor_names = T)

#################################################################
# Meta analysis on cognitive function
#################################################################

model_IQ_DelDup = fread(
  "META_Model_with_DDD_ClinGen_DelDup_2024.tsv",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".")

library(grid)
library(metafor)
library(meta)
library(stringr)
library(tidyr)
library("cowplot")

names_meta = c("Model_2_DDD")

names_list = c("DEL_0","DEL_1", "DEL_2", "DEL_3", "DEL_4", "DEL_5", "DEL_6", "DEL_7", "DEL_8", "DEL_9", "DEL_10", "DEL_11", "DEL_12", "DEL_13", "DEL_14", "DEL_15", "DEL_16", "DEL_17", "DEL_18", "DEL_19", "DEL_20", "DEL_21", "DEL_22", "DEL_23", "DEL_24", "DEL_25", "DEL_26", "DEL_27", "DEL_28", "DEL_29", "DEL_30", "DEL_31", "DEL_32", "DEL_33", "DEL_34", "DEL_35", "DEL_36", "DEL_37", "DEL_38",
               "DUP_0","DUP_1", "DUP_2", "DUP_3", "DUP_4", "DUP_5", "DUP_6", "DUP_7", "DUP_8", "DUP_9", "DUP_10", "DUP_11", "DUP_12", "DUP_13", "DUP_14", "DUP_15", "DUP_16", "DUP_17", "DUP_18", "DUP_19", "DUP_20", "DUP_21", "DUP_22", "DUP_23", "DUP_24", "DUP_25", "DUP_26", "DUP_27", "DUP_28", "DUP_29", "DUP_30", "DUP_31", "DUP_32", "DUP_33", "DUP_34", "DUP_35", "DUP_36", "DUP_37", "DUP_38")

all_info_meta_DelDup_wind_modif_final = c()

model_IQ_DelDup$Estimate = model_IQ_DelDup$Value

model_IQ_DelDup = subset(model_IQ_DelDup, model_IQ_DelDup$Pheno == "ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final" & Cohorts != "ALL" & !(Cohorts == "MSSNG" | Cohorts == "SSC" | Cohorts == "SPARK"))

model_IQ_DelDup$Estimate = model_IQ_DelDup$Value

### with ASD

for (j in 1:length(names_meta)) {
  #j=1
  list_info_Hetero_est_fix_DelDup = c()
  list_info_Hetero_SD_fix_DelDup = c()
  list_info_Hetero_pvalue_fix_DelDup = c()
  
  list_info_Hetero_est_random_DelDup = c()
  list_info_Hetero_SD_random_DelDup = c()
  list_info_Hetero_pvalue_random_DelDup = c()
  
  list_info_Hetero_DelDup = c()
  list_info_Hetero_I_DelDup = c()
  list_info_Hetero_t2_DelDup = c()
  
  for (i in 1:length(names_list)) {
    #i=1
    File_Analyses = subset(model_IQ_DelDup, names == names_list[i] & Model == names_meta[j] & (catg == 13 | catg == 14) & !(is.na(Estimate)))# & !(Cohorts == "SSC" | Cohorts == "SPARK"))
    
    studies = c(File_Analyses$Cohorts)
    beta = File_Analyses$Estimate
    beta.se = File_Analyses$Std.Error
    
    # model with meta packgae fixed and random effect
    me = metagen(beta, beta.se, studlab = studies)
    
    list_info_Hetero_est_fix_DelDup = c(list_info_Hetero_est_fix_DelDup, me$TE.fixed)
    list_info_Hetero_SD_fix_DelDup = c(list_info_Hetero_SD_fix_DelDup, me$seTE.fixed)
    list_info_Hetero_pvalue_fix_DelDup = c(list_info_Hetero_pvalue_fix_DelDup, me$pval.fixed)
    
    list_info_Hetero_est_random_DelDup = c(list_info_Hetero_est_random_DelDup, me$TE.random)
    list_info_Hetero_SD_random_DelDup = c(list_info_Hetero_SD_random_DelDup, me$seTE.random)
    list_info_Hetero_pvalue_random_DelDup = c(list_info_Hetero_pvalue_random_DelDup, me$pval.random)
    
    list_info_Hetero_DelDup = c(list_info_Hetero_DelDup, me$pval.Q)
    list_info_Hetero_I_DelDup = c(list_info_Hetero_I_DelDup, me$I2)
    list_info_Hetero_t2_DelDup = c(list_info_Hetero_t2_DelDup, me$tau2)
    #forest(me, xlim = c(-0.3,0.1))
  }
  
  names_list_size = c(0, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975,
                      1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925,
                      0, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975,
                      1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925)
  
  names_list_TYPE = c("DEL","DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DUP","DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", 
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP",
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP",
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP")
  
  all_info_meta_DelDup_wind_modif = cbind(names_list, names_list_TYPE, names_list_size,
                                          list_info_Hetero_DelDup, list_info_Hetero_I_DelDup, list_info_Hetero_t2_DelDup,
                                          list_info_Hetero_est_fix_DelDup,list_info_Hetero_SD_fix_DelDup,list_info_Hetero_pvalue_fix_DelDup,
                                          list_info_Hetero_est_random_DelDup,list_info_Hetero_SD_random_DelDup,list_info_Hetero_pvalue_random_DelDup,names_meta[j],"with_ASD")
  
  colnames(all_info_meta_DelDup_wind_modif) = c("Names", "TYPE", "Size",
                                                "Heterog_pvalue","I", "T2",
                                                "Estimate_fix","Std.Error_fix","p.value_fix",
                                                "Estimate_random","Std.Error_random","p.value_random","Model","Cohorts")
  
  all_info_meta_DelDup_wind_modif = data.frame(all_info_meta_DelDup_wind_modif)
  
  all_info_meta_DelDup_wind_modif$Heterog_pvalue_adj = p.adjust(all_info_meta_DelDup_wind_modif$Heterog_pvalue, method ="fdr")
  
  all_info_meta_DelDup_wind_modif_final = rbind(all_info_meta_DelDup_wind_modif_final,-999, all_info_meta_DelDup_wind_modif)
  
}

for (j in 1:length(names_meta)) {
  
  list_info_Hetero_est_fix_DelDup = c()
  list_info_Hetero_SD_fix_DelDup = c()
  list_info_Hetero_pvalue_fix_DelDup = c()
  
  list_info_Hetero_est_random_DelDup = c()
  list_info_Hetero_SD_random_DelDup = c()
  list_info_Hetero_pvalue_random_DelDup = c()
  
  list_info_Hetero_DelDup = c()
  list_info_Hetero_I_DelDup = c()
  list_info_Hetero_t2_DelDup = c()
  
  for (i in 1:length(names_list)) {
    
    File_Analyses = subset(model_IQ_DelDup, names == names_list[i] & Model == names_meta[j] & (catg == 13 | catg == 14) & !(Cohorts == "SSC" | Cohorts == "SPARK"| Cohorts == "MSSNG"))
    
    studies = c(File_Analyses$Cohort_final)
    #studies = c("Imagen","SYS-child","SYS-parent","LBC","CaG-GSA","CaG-Omni2.5","G-Scot","SSC-1MV1","SSC-1MV3","SSC-Omni2.5","MSSNG")
    beta = File_Analyses$Estimate
    beta.se = File_Analyses$Std.Error
    
    # model with meta packgae fixed and random effect
    me = metagen(beta, beta.se, studlab = studies)
    
    list_info_Hetero_est_fix_DelDup = c(list_info_Hetero_est_fix_DelDup, me$TE.fixed)
    list_info_Hetero_SD_fix_DelDup = c(list_info_Hetero_SD_fix_DelDup, me$seTE.fixed)
    list_info_Hetero_pvalue_fix_DelDup = c(list_info_Hetero_pvalue_fix_DelDup, me$pval.fixed)
    
    list_info_Hetero_est_random_DelDup = c(list_info_Hetero_est_random_DelDup, me$TE.random)
    list_info_Hetero_SD_random_DelDup = c(list_info_Hetero_SD_random_DelDup, me$seTE.random)
    list_info_Hetero_pvalue_random_DelDup = c(list_info_Hetero_pvalue_random_DelDup, me$pval.random)
    
    list_info_Hetero_DelDup = c(list_info_Hetero_DelDup, me$pval.Q)
    list_info_Hetero_I_DelDup = c(list_info_Hetero_I_DelDup, me$I2)
    list_info_Hetero_t2_DelDup = c(list_info_Hetero_t2_DelDup, me$tau2)
    #forest(me, xlim = c(-0.3,0.1))
  }
  
  names_list_size = c(0, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975,
                      1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925,
                      0, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975,
                      1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925)
  
  names_list_TYPE = c("DEL","DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
                      "DUP","DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", 
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP",
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP",
                      "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP", "DUP")
  
  all_info_meta_DelDup_wind_modif = cbind(names_list, names_list_TYPE, names_list_size,
                                          list_info_Hetero_DelDup, list_info_Hetero_I_DelDup, list_info_Hetero_t2_DelDup,
                                          list_info_Hetero_est_fix_DelDup,list_info_Hetero_SD_fix_DelDup,list_info_Hetero_pvalue_fix_DelDup,
                                          list_info_Hetero_est_random_DelDup,list_info_Hetero_SD_random_DelDup,list_info_Hetero_pvalue_random_DelDup,names_meta[j],"without_ASD")
  
  colnames(all_info_meta_DelDup_wind_modif) = c("Names", "TYPE", "Size",
                                                "Heterog_pvalue","I", "T2",
                                                "Estimate_fix","Std.Error_fix","p.value_fix",
                                                "Estimate_random","Std.Error_random","p.value_random","Model","Cohorts")
  
  all_info_meta_DelDup_wind_modif = data.frame(all_info_meta_DelDup_wind_modif)
  
  all_info_meta_DelDup_wind_modif$Heterog_pvalue_adj = p.adjust(all_info_meta_DelDup_wind_modif$Heterog_pvalue, method ="fdr")
  
  all_info_meta_DelDup_wind_modif_final = rbind(all_info_meta_DelDup_wind_modif_final,-999, all_info_meta_DelDup_wind_modif)
  
}

Meta_all_info_meta_DelDup_wind_modif_final = all_info_meta_DelDup_wind_modif_final

test_1 = subset(all_info_meta_DelDup_wind_modif_final,
                all_info_meta_DelDup_wind_modif_final$Estimate_fix != -999 & 
                  all_info_meta_DelDup_wind_modif_final$Cohorts == "with_ASD") 


#################################################################
# Pool analysis on cognitive function
#################################################################

Model_C0_DelDup_all_final = read.csv(
  paste0("POOL_Model_with_DDD_ClinGen_DelDup_2024.tsv"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)

#################################################################

score_cutoff_LOEUF = c(10000,60,40,20,10)
age_value = c(100*120,75*12,70*12,65*12,60*12,55*12)
ancestry = c("ALL","EUR")
Diagnostic_asd = c("with_ASD","without_ASD")
Model_stat = c("Model_02")

#################################################################

Model_C0_DelDup_select_all_final = subset(Model_C0_DelDup_all_final, (catg == 13 | catg == 14) & PC == "with_PC")

Model_C0_DelDup_select_all_final_pFDR = c()

for (anc_code in 1:length(ancestry)) {
  print(ancestry[anc_code])
  
  for (i in 1:length(Diagnostic_asd)) {
    print(Diagnostic_asd[i])
    
    for (l in 1:length(age_value)) {
      print(age_value[l])
      
      for (k in 1:length(score_cutoff_LOEUF)) {
        print(score_cutoff_LOEUF[k])
        
        for (j in 1:length(Model_stat)) {
          
          ### Deletion
          Fichier_model_DEL = subset(Model_C0_DelDup_select_all_final, TYPE_CNV == "DEL" & scores_catg == score_cutoff_LOEUF[k] & Model == Model_stat[j] & Age == age_value[l] & ASD == Diagnostic_asd[i] & Anc == ancestry[anc_code])
          print(paste0("DEL_Model_", dim(Fichier_model_DEL)[1]))
          
          Fichier_model_DEL$p.value_FDR = p.adjust(Fichier_model_DEL$p.value, method ="fdr")
          
          ### Duplication
          Fichier_model_DUP = subset(Model_C0_DelDup_select_all_final, TYPE_CNV == "DUP" & scores_catg == score_cutoff_LOEUF[k] & Model == Model_stat[j] & Age == age_value[l] & ASD == Diagnostic_asd[i] & Anc == ancestry[anc_code])
          print(paste0("DUP_Model_", dim(Fichier_model_DUP)[1]))
          
          Fichier_model_DUP$p.value_FDR = p.adjust(Fichier_model_DUP$p.value, method ="fdr")
          ### Merge information
          
          Model_C0_DelDup_select_all_final_pFDR = rbind(Model_C0_DelDup_select_all_final_pFDR, Fichier_model_DEL, Fichier_model_DUP)
          
        }}}}}

#################################################################

Model_C0_DelDup_2 = subset(Model_C0_DelDup_select_all_final_pFDR, catg == 13 | catg == 14)

Model_C0_DelDup_2$Size = as.numeric(Model_C0_DelDup_2$Size)
Model_C0_DelDup_2$Estimate = as.numeric(Model_C0_DelDup_2$Estimate)
Model_C0_DelDup_2$Std.Error = as.numeric(Model_C0_DelDup_2$Std.Error)
Model_C0_DelDup_2$Age = as.numeric(Model_C0_DelDup_2$Age)
Model_C0_DelDup_2$scores_catg = as.numeric(Model_C0_DelDup_2$scores_catg)
Model_C0_DelDup_2$p.value = as.numeric(Model_C0_DelDup_2$p.value)
Model_C0_DelDup_2$p.value_FDR = as.numeric(Model_C0_DelDup_2$p.value_FDR)

##################################################################################################################################

Model_C0_Del_Model_All = subset(Model_C0_DelDup_2, scores_catg >= 10000 & Age >= 100*120 & Anc == "ALL" & ASD == "without_ASD" & Model == "Model_02" )

Model_C0_Del_Model_All$TYPE = Model_C0_Del_Model_All$TYPE_CNV

Model_C0_Del_Model_All$Size = ifelse(Model_C0_Del_Model_All$Size == 0.025, 0 , Model_C0_Del_Model_All$Size)

Info_Meta_Mega_Merge = merge(test_1, Model_C0_Del_Model_All, by = c("TYPE", "Size"))

Info_Meta_Mega_Merge$Estimate_Het = ifelse(Info_Meta_Mega_Merge$Heterog_pvalue >= 0.1, Info_Meta_Mega_Merge$Estimate_fix, Info_Meta_Mega_Merge$Estimate_random)

##################################################################################################################################

all_info_meta_DelDup_wind_modif_final$Estimate = ifelse(all_info_meta_DelDup_wind_modif_final$Heterog_pvalue >= 0.1,
                                                        all_info_meta_DelDup_wind_modif_final$Estimate_fix, all_info_meta_DelDup_wind_modif_final$Estimate_random)

all_info_meta_DelDup_wind_modif_final$p.value = ifelse(all_info_meta_DelDup_wind_modif_final$Heterog_pvalue >= 0.1,
                                                       all_info_meta_DelDup_wind_modif_final$p.value_fix, all_info_meta_DelDup_wind_modif_final$p.value_random)

all_info_meta_DelDup_wind_modif_final$Std.Error = ifelse(all_info_meta_DelDup_wind_modif_final$Heterog_pvalue >= 0.1,
                                                         all_info_meta_DelDup_wind_modif_final$Std.Error_fix, all_info_meta_DelDup_wind_modif_final$Std.Error_random)

############################################

test_1_del = subset(subset(all_info_meta_DelDup_wind_modif_final, all_info_meta_DelDup_wind_modif_final$TYPE == "DEL" &
                             all_info_meta_DelDup_wind_modif_final$Estimate_fix != -999 & 
                             all_info_meta_DelDup_wind_modif_final$Cohorts == "without_ASD"), select = c(Estimate, Size, p.value, Std.Error))
test_1_del$Model = "Meta"

test_1_del$p.value_FDR = p.adjust(test_1_del$p.value, method ="fdr")

Model_C0_Del_Model_All_selection = subset(subset(Model_C0_Del_Model_All, TYPE == "DEL"), select = c(Estimate, Size, p.value_FDR, Std.Error))
Model_C0_Del_Model_All_selection$Model = "Pool"

Meta_Pool_DEL = rbind(subset(test_1_del, select = c(Estimate, Size, p.value_FDR, Std.Error, Model)), Model_C0_Del_Model_All_selection)

Meta_Pool_DEL$Estimate = as.numeric(Meta_Pool_DEL$Estimate)
Meta_Pool_DEL$Size = as.numeric(Meta_Pool_DEL$Size)
Meta_Pool_DEL$p.value_FDR = as.numeric(Meta_Pool_DEL$p.value_FDR)
Meta_Pool_DEL$Std.Error = as.numeric(Meta_Pool_DEL$Std.Error)

############################################

test_1_dup = subset(subset(all_info_meta_DelDup_wind_modif_final, all_info_meta_DelDup_wind_modif_final$TYPE == "DUP" &
                             all_info_meta_DelDup_wind_modif_final$Estimate_fix != -999 & 
                             all_info_meta_DelDup_wind_modif_final$Cohorts == "without_ASD"), select = c(Estimate, Size, p.value, Std.Error))
test_1_dup$Model = "Meta"
test_1_dup$p.value_FDR = p.adjust(test_1_dup$p.value, method ="fdr")

Model_C0_Dup_Model_All_selection = subset(subset(Model_C0_Del_Model_All, TYPE == "DUP"), select = c(Estimate, Size, p.value_FDR, Std.Error))
Model_C0_Dup_Model_All_selection$Model = "Pool"

Meta_Pool_DUP = rbind(subset(test_1_dup, select = c(Estimate, Size, p.value_FDR, Std.Error, Model)) , Model_C0_Dup_Model_All_selection)

Meta_Pool_DUP$Estimate = as.numeric(Meta_Pool_DUP$Estimate)
Meta_Pool_DUP$Size = as.numeric(Meta_Pool_DUP$Size)
Meta_Pool_DUP$p.value_FDR = as.numeric(Meta_Pool_DUP$p.value_FDR)
Meta_Pool_DUP$Std.Error = as.numeric(Meta_Pool_DUP$Std.Error)

#################################################################################################################################################
# For Deletion

pd <- position_dodge(0.02)
colors_2_fig = c("red", "red")

DEL_Model_sensitivity_p2 = 
  ggplot(Meta_Pool_DEL, aes(x=Size, y=Estimate, colour=Model)) + ylim(-0.25,0.1) + xlim(0.01,2) + labs(x ="LOEUF") + 
  geom_errorbar(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=Model),position=pd) +
  scale_color_manual(values=colors_2_fig) +
  scale_fill_manual(values=colors_2_fig) +
  geom_line(aes(linetype=Model), position=pd, size=1) + scale_linetype_manual(values=c("dashed", "solid")) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(Meta_Pool_DEL$p.value_FDR < 0.05, 5, 2),
             shape=ifelse(Meta_Pool_DEL$Model == "Meta" & Meta_Pool_DEL$p.value_FDR < 0.05 , 21,
                          ifelse(Meta_Pool_DEL$Model == "Meta" & Meta_Pool_DEL$p.value_FDR >= 0.05 , 1,
                                 ifelse(Meta_Pool_DEL$Model != "Meta" & Meta_Pool_DEL$p.value_FDR < 0.05 , 22, 0))),
             fill="white") + theme(legend.position = c(0.8, 0.2), axis.title.y = element_blank())


####

pd <- position_dodge(0.005)

DEL_Model_sensitivity_p1 = 
  ggplot(Meta_Pool_DEL, aes(x=Size, y=Estimate, colour=Model)) + ylim(-0.8, 0.15) + xlim(-0.01,0.01) + labs(x ="LOEUF") + 
  geom_errorbar(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=Model),position=pd) +
  scale_color_manual(values=colors_2_fig) +
  scale_fill_manual(values=colors_2_fig) +
  geom_line(position=pd, size=2) +
  geom_hline(yintercept=0) +
  geom_point(position=pd, size=ifelse(Meta_Pool_DEL$p.value_FDR < 0.05, 5, 2),
             shape=ifelse(Meta_Pool_DEL$Model == "Meta" & Meta_Pool_DEL$p.value_FDR < 0.05 , 21,
                          ifelse(Meta_Pool_DEL$Model == "Meta" & Meta_Pool_DEL$p.value_FDR >= 0.05 , 1,
                                 ifelse(Meta_Pool_DEL$Model != "Meta" & Meta_Pool_DEL$p.value_FDR < 0.05 , 22, 0))),
             fill="white") + theme(legend.position='none')

#################################################################################################################################################
# For Duplication

pd <- position_dodge(0.02)
colors_2_fig = c("blue", "blue")

DUP_Model_sensitivity_p2 = 
  ggplot(Meta_Pool_DUP, aes(x=Size, y=Estimate, colour=Model)) + ylim(-0.15,0.06) + xlim(0.01,2) + labs(x ="LOEUF") + 
  geom_errorbar(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=Model),position=pd) +
  scale_color_manual(values=colors_2_fig) +
  scale_fill_manual(values=colors_2_fig) +
  geom_line(aes(linetype=Model), position=pd, size=1) + scale_linetype_manual(values=c("dashed", "solid")) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(Meta_Pool_DUP$p.value_FDR < 0.05, 5, 2),
             shape=ifelse(Meta_Pool_DUP$Model == "Meta" & Meta_Pool_DUP$p.value_FDR < 0.05 , 21,
                          ifelse(Meta_Pool_DUP$Model == "Meta" & Meta_Pool_DUP$p.value_FDR >= 0.05 , 1,
                                 ifelse(Meta_Pool_DUP$Model != "Meta" & Meta_Pool_DUP$p.value_FDR < 0.05 , 22, 0))),
             fill="white") + theme(legend.position = c(0.8, 0.2), axis.title.y = element_blank())

####

pd <- position_dodge(0.005)

DUP_Model_sensitivity_p1 = 
  ggplot(Meta_Pool_DUP, aes(x=Size, y=Estimate, colour=Model)) + ylim(-0.6, 0) + xlim(-0.01,0.01) + labs(x ="LOEUF") + 
  geom_errorbar(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=Model),position=pd) +
  scale_color_manual(values=colors_2_fig) +
  scale_fill_manual(values=colors_2_fig) +
  geom_line(position=pd, size=2) +
  geom_hline(yintercept=0) +
  geom_point(position=pd, size=ifelse(Meta_Pool_DUP$p.value_FDR < 0.05, 5, 2),
             shape=ifelse(Meta_Pool_DUP$Model == "Meta" & Meta_Pool_DUP$p.value_FDR < 0.05 , 21,
                          ifelse(Meta_Pool_DUP$Model == "Meta" & Meta_Pool_DUP$p.value_FDR >= 0.05 , 1,
                                 ifelse(Meta_Pool_DUP$Model != "Meta" & Meta_Pool_DUP$p.value_FDR < 0.05 , 22, 0))),
             fill="white") + theme(legend.position='none')


########################################################################
# Figure 1 A,B,C
ggdraw() +
  
  draw_plot(LOEUF_DELDUP_1, 0.02, .5, 0.25, .5) +
  draw_plot(LOEUF_DELDUP_1_carrier, 0.02, 0, 0.25, .5) +
  draw_plot(Venn_del_dup, 0.35, 0, 0.55, 1) +
  
  draw_plot_label(c("A", "B", "C"),
                  c(0, 0, 0.3),
                  c(1, 0.5, 1), size = 15)

# Figure 2B
ggdraw() +

  draw_plot(DEL_Model_sensitivity_p1, 0, .50, 0.15, .50) +
  draw_plot(DEL_Model_sensitivity_p2, 0.20, .50, 0.70, .50) +
  draw_plot(DUP_Model_sensitivity_p1, 0, 0, 0.15, .50) +
  draw_plot(DUP_Model_sensitivity_p2, 0.20, 0, 0.70, .50) +
  
  draw_plot_label(c("A", "B"),
                  c(0, 0),
                  c(1, 0.5), size = 15)
