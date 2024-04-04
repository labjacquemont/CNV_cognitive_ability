library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)
'%!in%' <- Negate('%in%')

###################################################################################################
setwd(paste0("/Users/guillaume/Desktop/File_Script/Fig_1_2"))
###################################################################################################

List_pheno_IID_CNV_gene = read.csv("List_pheno_IID_CNV_gene_ENSG00000198142.tsv", header = TRUE, sep = "\t", na.strings = c("NA", "#NULL!", ".", ""), dec = ".", fileEncoding = "UTF-8-BOM")

###################################################################################################

List_pheno_IID_CNV_gene_0 = subset(List_pheno_IID_CNV_gene, Cohort == "UKBB" & DUP.loeuf_inv == 0 & DEL.loeuf_inv == 0 & Merge_final_ancestry == "EUR" & !(is.na(z_score_FI_all_adj_age_sex_PC1_10)))
List_pheno_IID_CNV_gene_carriers = subset(List_pheno_IID_CNV_gene, Cohort == "UKBB" & (DUP.loeuf_inv > 0 | DEL.loeuf_inv > 0) & gene_id_carriers == 1 & Merge_final_ancestry == "EUR" & !(is.na(z_score_FI_all_adj_age_sex_PC1_10)))
List_pheno_IID_CNV_gene_no_carriers = subset(List_pheno_IID_CNV_gene, Cohort == "UKBB" & (DUP.loeuf_inv > 0 | DEL.loeuf_inv > 0) & gene_id_carriers == 0 & Merge_final_ancestry == "EUR" & !(is.na(z_score_FI_all_adj_age_sex_PC1_10)))

######
List_pheno_IID_CNV_gene_0$Grp = "Control"
List_pheno_IID_CNV_gene_carriers$Grp = "CNV"
List_pheno_IID_CNV_gene_no_carriers$Grp = "CNV_other"

##################

List_pheno_IID_CNV_gene_all = rbind(List_pheno_IID_CNV_gene_0, List_pheno_IID_CNV_gene_carriers, List_pheno_IID_CNV_gene_no_carriers)

##################

test_1 = t.test( List_pheno_IID_CNV_gene_0$z_score_FI_all_adj_age_sex_PC1_10,
                 List_pheno_IID_CNV_gene_carriers$z_score_FI_all_adj_age_sex_PC1_10,
                 alternative = c("two.sided"),
                 conf.level = 0.95)

test_2 = t.test(List_pheno_IID_CNV_gene_0$z_score_FI_all_adj_age_sex_PC1_10,
                List_pheno_IID_CNV_gene_no_carriers$z_score_FI_all_adj_age_sex_PC1_10,
                alternative = c("two.sided"),
                conf.level = 0.95)

test_3 = t.test(List_pheno_IID_CNV_gene_carriers$z_score_FI_all_adj_age_sex_PC1_10,
                List_pheno_IID_CNV_gene_no_carriers$z_score_FI_all_adj_age_sex_PC1_10,
                alternative = c("two.sided"),
                conf.level = 0.95)

##################

F1_chr2 = 
  ggplot(List_pheno_IID_CNV_gene_all, aes(x=Grp, y=ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final, fill=Grp)) + 
  #geom_violin() +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(width=0.5,outlier.shape = NA) + 
  scale_fill_manual(values=c("#4DAC26", "#F4A582", "#56B4E9")) +
  ylim(-3, 4) + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), legend.position = "none", text = element_text(size = 16)) +
  labs( y = "Z-score cog. ability adj.")

##################

List_pheno_IID_CNV_gene_UKBB = subset(List_pheno_IID_CNV_gene, Merge_final_ancestry == "EUR" & !(is.na(z_score_FI_all_adj_age_sex_PC1_10)))

m2 <- lm(z_score_FI_all_adj_age_sex_PC1_10 ~ DEL.loeuf_inv + DUP.loeuf_inv, List_pheno_IID_CNV_gene_UKBB , na.action = na.omit)
pred_mean <- predict(m2, newdata = data.frame(DEL.loeuf_inv = mean(List_pheno_IID_CNV_gene_UKBB$DEL.loeuf_inv), DUP.loeuf_inv = mean(List_pheno_IID_CNV_gene_UKBB$DUP.loeuf_inv)))
List_pheno_IID_CNV_gene_UKBB$z_score_FI_all_adj_age_sex_PC1_10_CNV <- pred_mean + m2$residuals

t_test = t.test(subset(List_pheno_IID_CNV_gene_UKBB, gene_id_carriers == 0 )$z_score_FI_all_adj_age_sex_PC1_10_CNV,
                subset(List_pheno_IID_CNV_gene_UKBB, gene_id_carriers == 1 )$z_score_FI_all_adj_age_sex_PC1_10_CNV,
                alternative = c("two.sided"),
                conf.level = 0.95)

F2_chr2 = 
  ggplot(List_pheno_IID_CNV_gene_UKBB, aes(x=as.factor(gene_id_carriers), y=z_score_FI_all_adj_age_sex_PC1_10_CNV, fill=as.factor(gene_id_carriers))) + 
  #geom_violin() +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(width=0.5,outlier.shape = NA) + 
  scale_fill_manual(values=c("#80CDC1", "#D01C8B")) +
  ylim(-3, 4) + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none", text = element_text(size = 16)) #+

###################################################################################################
###################################################################################################

model_stat = c("1_0")
Model_Gene_DelDup_select_all_final_pFDR_score_Ok = read.csv("Model_Gene_DelDup_select_all_final_pFDR_score_Ok_filtre_2024.tsv", header = TRUE, sep = "\t", na.strings = c("NA", "#NULL!", ".", ""), dec = ".", fileEncoding = "UTF-8-BOM")

###################################################################################################

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL = subset(Model_Gene_DelDup_select_all_final_pFDR_score_Ok,# total_carriers >= 10 &
                                                                   TYPE == "DEL" & model_type == model_stat & !(is.na(chromosome))) #& p.value_FDR < 0.05


Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$CHR = as.numeric(Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$chromosome)

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$BP = ((Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$end_position - Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$start_position)/2) + Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$start_position

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$P = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$pvalue

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$P_log10 = -log10(Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$pvalue)

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$SNP = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL$gene_list_name


#########################

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP = subset(Model_Gene_DelDup_select_all_final_pFDR_score_Ok, #total_carriers >= 10 &
                                                                   TYPE == "DUP" & model_type == model_stat & !(is.na(chromosome))) #& p.value_FDR < 0.05

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$CHR = as.numeric(Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$chromosome)

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$BP = ((Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$end_position - Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$start_position)/2) + Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$start_position

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$P = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$pvalue

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$P_log10 = -log10(Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$pvalue)

Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$SNP = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP$gene_list_name

#########################

CNV_del = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DEL
CNV_dup = Model_Gene_DelDup_select_all_final_pFDR_score_Ok_test_DUP

#########################

# ajouter la taille max des chromosomes et calculer les position sur l'axe x du miami plot
chromosomesize = read.delim2(paste0("chromosomesize.txt"))

chromosomesize_bonus_start = chromosomesize
chromosomesize_bonus_stop = chromosomesize

chromosomesize_bonus_stop$P = 1
chromosomesize_bonus_stop$p.value_FDR = 1
chromosomesize_bonus_stop$SNP = "Position_Stop"
chromosomesize_bonus_stop$total_carriers = 0
chromosomesize_bonus_stop$Estimate = 0
names(chromosomesize_bonus_stop) = c("CHR","BP","P","p.value_FDR","SNP","total_carriers","Estimate")

chromosomesize_bonus_start$Basepairs = 1
chromosomesize_bonus_start$P = 1
chromosomesize_bonus_start$p.value_FDR = 1
chromosomesize_bonus_start$SNP = "Position_Start"
chromosomesize_bonus_start$total_carriers = 0
chromosomesize_bonus_start$Estimate = 0
names(chromosomesize_bonus_start) = c("CHR","BP","P","p.value_FDR","SNP","total_carriers","Estimate")

CNV_del_ok = rbind(subset(CNV_del, select = c("CHR","BP","P","p.value_FDR","SNP","total_carriers","Estimate")), chromosomesize_bonus_start, chromosomesize_bonus_stop)
CNV_dup_ok = rbind(subset(CNV_dup, select = c("CHR","BP","P","p.value_FDR","SNP","total_carriers","Estimate")), chromosomesize_bonus_start, chromosomesize_bonus_stop)

CNV_del_ok$TYPE = "DEL"
CNV_dup_ok$TYPE = "DUP"

CNV_deldup_ok = rbind(CNV_del_ok,CNV_dup_ok)

# ajouter la taille max des chromosomes et calculer les position sur l'axe x du miami plot
chromosomesize = read.delim2(paste0("chromosomesize.txt"))

Gene.info_deldup = merge(CNV_deldup_ok, chromosomesize, by ="CHR", all.x=T)

# Prepare the dataset
CNV_deldup_ok_1 <- Gene.info_deldup %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(Basepairs)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(Gene.info_deldup, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

CNV_del_ok_1 = subset(CNV_deldup_ok_1, TYPE == "DEL")
CNV_dup_ok_1 = subset(CNV_deldup_ok_1, TYPE == "DUP")

# Prepare X axis
all = rbind(CNV_del_ok_1[,c("CHR","BPcum")],CNV_dup_ok_1[,c("CHR","BPcum")])

# Prepare X axis
axisdf = all %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#Max and min positions of genes
min_pos=min(CNV_del_ok_1$BPcum, CNV_dup_ok_1$BPcum)
max_pos=max(CNV_del_ok_1$BPcum, CNV_dup_ok_1$BPcum)




GWAS_CNV_DEL = ggplot() +
  
  #All CNV
  # Show all points
  
  geom_point(data=CNV_del_ok_1, aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha=0.8, size=4,
             #shape= ifelse(CNV_del_ok_1$p.value_FDR < 0.05 & CNV_del_ok_1$Estimate < 0 & CNV_del_ok_1$total_carriers >= 30, 25, 16), fill="white") +
             shape= ifelse(CNV_del_ok_1$p.value_FDR < 0.05 & CNV_del_ok_1$Estimate < 0 & CNV_del_ok_1$total_carriers >= 30, 25,
                           ifelse(CNV_del_ok_1$p.value_FDR < 0.05 & CNV_del_ok_1$Estimate > 0 & CNV_del_ok_1$total_carriers >= 30, 24,16)), fill="white") +
  scale_color_manual(values = rep(c("red", "darkred"), 21 )) +
  
  ## All CNV
  #geom_point(data=subset(CNV_del, p.value_FDR < 0.05 & Estimate > 0 & total_carriers >= 30), aes(x=BPcum, y=-log10(pvalue), color="Significant DEL"), size=2.5, shape= 24) +
  #geom_point(data=subset(CNV_del, p.value_FDR < 0.05 & Estimate < 0 & total_carriers >= 30), aes(x=BPcum, y=-log10(pvalue), color="Significant DEL"), size=2.5, shape= 25, fill="white") +
  
  geom_segment(data=CNV_del_ok_1,aes(x = min_pos, y = -log10(0.05), xend = max_pos, yend = -log10(0.05)), linetype = "dotted", color = "black") +
  
  # custom X axis:
  geom_label(data=axisdf, aes( y=-1, x= center, label=CHR), fill="white", fontface = "bold", alpha = 0.5, hjust=0.5) +
  
  scale_y_continuous(breaks=seq(0, 20, 10)) + #trans="sqrt",
  
  # Custom the theme:
  theme_bw() +
  theme(
    plot.title = element_text(hjust=0.5, size = 15),
    axis.title.x = element_blank(),
    axis.title = element_text(size=18),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 13),
    #legend.key.size = unit(2, 'cm'),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
  )

GWAS_CNV_DUP = ggplot() +
  
  #All CNV
  # Show all points
  
  geom_point(data=CNV_dup_ok_1, aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha=0.8, size=4,
             shape= ifelse(CNV_dup_ok_1$p.value_FDR < 0.05 & CNV_dup_ok_1$Estimate < 0 & CNV_dup_ok_1$total_carriers >= 30, 25,
                           ifelse(CNV_dup_ok_1$p.value_FDR < 0.05 & CNV_dup_ok_1$Estimate > 0 & CNV_dup_ok_1$total_carriers >= 30, 24,16)), fill="white") + 
  scale_color_manual(values = rep(c("darkblue", "cyan3"), 21 )) +
  
  ## All CNV
  #geom_point(data=subset(CNV_del, p.value_FDR < 0.05 & Estimate > 0 & total_carriers >= 30), aes(x=BPcum, y=-log10(pvalue), color="Significant DEL"), size=2.5, shape= 24) +
  #geom_point(data=subset(CNV_del, p.value_FDR < 0.05 & Estimate < 0 & total_carriers >= 30), aes(x=BPcum, y=-log10(pvalue), color="Significant DEL"), size=2.5, shape= 25, fill="white") +
  
  geom_segment(data=CNV_dup_ok_1,aes(x = min_pos, y = -log10(0.05), xend = max_pos, yend = -log10(0.05)), linetype = "dotted", color = "black") +
  
  # custom X axis:
  #geom_label(data=axisdf, aes( y=-0.5, x= center, label=CHR), fill="white", fontface = "bold", alpha = 0.5, hjust=0.5) +
  
  scale_y_continuous(breaks=seq(0, 20, 10)) + scale_y_reverse() + #trans="sqrt",
  
  # Custom the theme:
  theme_bw() +
  theme(
    plot.title = element_text(hjust=0.5, size = 15),
    axis.title.x = element_blank(),
    axis.title = element_text(size=18),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 13),
    #legend.key.size = unit(2, 'cm'),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
  )

########################################################################


library("cowplot")

########

ggdraw() +
  
  draw_plot(GWAS_CNV_DEL , 0, 0.5, 0.6, .40) +
  draw_plot(GWAS_CNV_DUP , 0, 0.1, 0.6, .40) +
  
  draw_plot(F1_chr2 , 0.62, 0, 0.20, 1) +
  draw_plot(F2_chr2 , 0.82, 0, 0.13, 1) +
  
  draw_plot_label(c("A", "B"),
                  c(0, 0.6),
                  c(1, 1), size = 15)
########
