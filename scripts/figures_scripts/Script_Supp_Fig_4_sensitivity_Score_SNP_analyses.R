library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)
'%!in%' <- Negate('%in%')

#################################################################

list_file = c("Score_30_SNP_10", "Score_30_SNP_15", "Score_30_SNP_30", "Score_30_SNP_50", "Score_30_SNP_100",
              "Score_40_SNP_10", "Score_50_SNP_10", "Score_75_SNP_10", "Score_100_SNP_10")

setwd(paste0("/Users/guillaume/Desktop/File_Script/Supp_Fig_4/"))

Model_C0_DelDup_all_final = read.csv(
  paste0("Model_with_DDD_ClinGen_DelDup_sensitivity_SCORE_SNP_2024.tsv"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)


#################################################################

score_cutoff_LOEUF = c(10000)
age_value = c(100*120)
ancestry = c("ALL","EUR")
Diagnostic_asd = c("without_ASD")
Model_stat = c("Model_02")
Files_pheno_score = list_file

summary(as.factor(Model_C0_DelDup_all_final$file_name_filtre))

#################################################################

Model_C0_DelDup_select_all_final = subset(Model_C0_DelDup_all_final, catg == 13 | catg == 14)

Model_C0_DelDup_select_all_final_pFDR = c()

for (Code in 1:length(Files_pheno_score)) {
  print(Files_pheno_score[Code])
  
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
            Fichier_model_DEL = subset(Model_C0_DelDup_select_all_final, TYPE_CNV == "DEL" & file_name_filtre == Files_pheno_score[Code] &
                                         scores_catg == score_cutoff_LOEUF[k] & Age == age_value[l] & ASD == Diagnostic_asd[i] & Anc == ancestry[anc_code] & Model == Model_stat[j])
            #print(paste0("DEL_Model_", dim(Fichier_model_DEL)[1]))
            
            Fichier_model_DEL$p.value_FDR = p.adjust(Fichier_model_DEL$p.value, method ="fdr")
            
            ### Duplication
            Fichier_model_DUP = subset(Model_C0_DelDup_select_all_final, TYPE_CNV == "DUP" & file_name_filtre == Files_pheno_score[Code] & scores_catg == score_cutoff_LOEUF[k] & Age == age_value[l] & ASD == Diagnostic_asd[i] & Anc == ancestry[anc_code] & Model == Model_stat[j])
            #print(paste0("DUP_Model_", dim(Fichier_model_DUP)[1]))
            
            Fichier_model_DUP$p.value_FDR = p.adjust(Fichier_model_DUP$p.value, method ="fdr")
            ### Merge information
            
            Model_C0_DelDup_select_all_final_pFDR = rbind(Model_C0_DelDup_select_all_final_pFDR, Fichier_model_DEL, Fichier_model_DUP)
            
          }}}}}}

#################################################################

Model_C0_DelDup_2 = subset(Model_C0_DelDup_select_all_final_pFDR, catg == 13 | catg == 14)

Model_C0_DelDup_2$Size = as.numeric(Model_C0_DelDup_2$Size)
Model_C0_DelDup_2$Estimate = as.numeric(Model_C0_DelDup_2$Estimate)
Model_C0_DelDup_2$Std.Error = as.numeric(Model_C0_DelDup_2$Std.Error)
Model_C0_DelDup_2$Age = as.numeric(Model_C0_DelDup_2$Age)
Model_C0_DelDup_2$scores_catg = as.numeric(Model_C0_DelDup_2$scores_catg)
Model_C0_DelDup_2$p.value = as.numeric(Model_C0_DelDup_2$p.value)
Model_C0_DelDup_2$p.value_FDR = as.numeric(Model_C0_DelDup_2$p.value_FDR)

#################################################################

summary(as.factor(Model_C0_DelDup_2$file_name_filtre))

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_30_SNP_10", "Score_030_SNP_010",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_30_SNP_15", "Score_030_SNP_015",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_30_SNP_30", "Score_030_SNP_030",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_30_SNP_50", "Score_030_SNP_050",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_30_SNP_100", "Score_030_SNP_100",
                                            Model_C0_DelDup_2$file_name_filtre)

###

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_40_SNP_10", "Score_040_SNP_010",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_50_SNP_10", "Score_050_SNP_010",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_75_SNP_10", "Score_075_SNP_010",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$file_name_filtre = ifelse(Model_C0_DelDup_2$file_name_filtre == "Score_100_SNP_10", "Score_100_SNP_010",
                                            Model_C0_DelDup_2$file_name_filtre)

Model_C0_DelDup_2$Size = ifelse(Model_C0_DelDup_2$catg == "14", 0,
                                            Model_C0_DelDup_2$Size)

###

Test = subset(Model_C0_DelDup_2, TYPE_CNV == "DEL" & Anc == "ALL" & 
                (file_name_filtre == "Score_030_SNP_010"
                 | file_name_filtre == "Score_030_SNP_100"
                 | file_name_filtre == "Score_030_SNP_015"
                 | file_name_filtre == "Score_030_SNP_030"
                 | file_name_filtre == "Score_030_SNP_050"))

pd <- position_dodge(0) # move them .05 to the left and right

DEL_All_SNP_2 = ggplot(subset(Test, catg != 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg != 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position = c(0.75, 0.3), axis.title.y = element_blank())

pd <- position_dodge(0.01)

DEL_All_SNP_1 =   ggplot(subset(Test, catg == 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #ylim(-0.5, 0) +
  xlim(-0.01,0.01) + 
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg == 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position='none')
  
#################################################################

Test = subset(Model_C0_DelDup_2, TYPE_CNV == "DUP" & Anc == "ALL" & 
                (file_name_filtre == "Score_030_SNP_010"
                 | file_name_filtre == "Score_030_SNP_100"
                 | file_name_filtre == "Score_030_SNP_015"
                 | file_name_filtre == "Score_030_SNP_030"
                 | file_name_filtre == "Score_030_SNP_050"))

pd <- position_dodge(0) # move them .05 to the left and right

DUP_All_SNP_2 = ggplot(subset(Test, catg != 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg != 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position = c(0.75, 0.3), axis.title.y = element_blank())

pd <- position_dodge(0.01)

DUP_All_SNP_1 =   ggplot(subset(Test, catg == 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #ylim(-0.5, 0) +
  xlim(-0.01,0.01) + 
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg == 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position='none')

###

Test = subset(Model_C0_DelDup_2, TYPE_CNV == "DEL" & Anc == "ALL" & 
                (file_name_filtre == "Score_030_SNP_010"
                 | file_name_filtre == "Score_040_SNP_010"
                 | file_name_filtre == "Score_050_SNP_010"
                 | file_name_filtre == "Score_075_SNP_010"
                 | file_name_filtre == "Score_100_SNP_010"))

pd <- position_dodge(0) # move them .05 to the left and right

DEL_All_Score_2 = ggplot(subset(Test, catg != 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg != 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position = c(0.75, 0.3), axis.title.y = element_blank())

pd <- position_dodge(0.01)

DEL_All_Score_1 =   ggplot(subset(Test, catg == 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #ylim(-0.5, 0) +
  xlim(-0.01,0.01) + 
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg == 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position='none')

#################################################################

Test = subset(Model_C0_DelDup_2, TYPE_CNV == "DUP" & Anc == "ALL" & 
                (file_name_filtre == "Score_030_SNP_010"
                 | file_name_filtre == "Score_040_SNP_010"
                 | file_name_filtre == "Score_050_SNP_010"
                 | file_name_filtre == "Score_075_SNP_010"
                 | file_name_filtre == "Score_100_SNP_010"))

pd <- position_dodge(0) # move them .05 to the left and right

DUP_All_Score_2 = ggplot(subset(Test, catg != 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
  #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
  #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
  scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg != 14)$file_name_filtre))), palette = "Temps")) +
  geom_line(position=pd, size=3) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(subset(Test, catg != 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position = c(0.75, 0.3), axis.title.y = element_blank())

pd <- position_dodge(0.01)

DUP_All_Score_1 =   ggplot(subset(Test, catg == 14), aes(x=Size, y=Estimate, colour=as.factor(file_name_filtre))) +
    #ylim(-0.5, 0) +
    xlim(-0.01,0.01) + 
    #geom_errorbar(aes(ymin=Est-(1.96*SE), ymax=Est+(1.96*SE)), position=pd, size=1) +
    #geom_ribbon(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error),fill=as.factor(Age)),alpha=0.3) +
    scale_color_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
    scale_fill_manual(values=hcl.colors(length(summary(as.factor(subset(Test, catg == 14)$file_name_filtre))), palette = "Temps")) +
    geom_line(position=pd, size=3) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
    geom_point(position=pd, size=ifelse(subset(Test, catg == 14)$p.value_FDR < 0.05, 5, 1) , shape=21, fill="white") + 
  theme(legend.position='none')


######

library("cowplot")

ggdraw() +
  draw_plot(DEL_All_Score_1, 0, 0.5, 0.08, .5) +  #A
  draw_plot(DEL_All_Score_2, 0.08, 0.5, 0.42, .5) +  #A
  
  draw_plot(DUP_All_Score_1, 0.5, 0.5, 0.08, .5) +  #B
  draw_plot(DUP_All_Score_2, 0.58, 0.5, 0.42, .5) +  #B
  
  draw_plot(DEL_All_SNP_1, 0, 0, 0.08, .5) +  #C
  draw_plot(DEL_All_SNP_2, 0.08, 0, 0.42, .5) +  #C
  
  draw_plot(DUP_All_SNP_1, 0.5, 0, 0.08, .5) +  #D
  draw_plot(DUP_All_SNP_2, 0.58, 0, 0.42, .5) +  #D
  
  draw_plot_label(c("A", "B", "C", "D"),
                  c(0, 0.5, 0, 0.5),
                  c(1, 1, 0.5, 0.5), size = 15)
