library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)

'%!in%' <- Negate('%in%')

#################################################################

all_data_info_pheno_scores_SNP_Artefact = fread( "/Users/guillaume/Desktop/File_Script/all_data_info_pheno_scores_SNP_Artefact_select_cohorts_20230505.tsv",
  #"all_data_info_pheno_scores_SNP_Artefact_select_cohorts_20220909.tsv",
  header = TRUE,sep = "\t", na.strings = c("NA", "#NULL!", ".", ""), dec = "." )

#########################################################################################################################
#########################################################################################################################

C0_Del_2 = c()
ancestry = c("ALL","EUR")
type_list = c("DEL","DUP")
list_DDD = "DDD_ClinGen"
score_tested = "LOEUF"
DDD_selection = list_DDD
DDD_model = paste0("without",list_DDD)
cohorts = c("ALL","CaG","G-Scot","Imagen","LBC","SPARK","SSC","SYS","Z_Gfact_5_Donald_et_al_2016","Z_Gfact_5_online","Z_score_FI","Z_score_FI_online","MSSNG")
model_IQ_DelDup = c()
score_cutoff_LOEUF = 10000


pheno_info = c("ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final","ZScore_IQ_adj_test_age_sex_PC1_10_with_online_UKBB_Final","ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB_Final","ZScore_IQ_adj_age_sex_with_online_UKBB_Final")

# With PC in formula
for (i_pheno in 1:length(pheno_info)) {

  all_data_info_pheno_scores_SNP_Artefact$pheno_test = all_data_info_pheno_scores_SNP_Artefact[,pheno_info[i_pheno]]

for (i in 1:length(cohorts)) {

  if (cohorts[i] != "ALL") {Fichier = subset(all_data_info_pheno_scores_SNP_Artefact, Cohort_final == cohorts[i] & !(is.na(pheno_test)))}
  if (cohorts[i] == "ALL") {Fichier = all_data_info_pheno_scores_SNP_Artefact}
  
  for (type in type_list) {
    
    scores_catg = paste0("catg_", DDD_model, "_", score_tested, "_", 1:40)
    scores_catg_bis = paste0("catg_", score_tested, "_", 1:40)
    
    scores_catg_DDD = DDD_selection
    
    LOEUF_Sum_1 = paste0("catg_", score_tested, "_", 1:20 )
    LOEUF_Sum_2 = paste0("catg_", score_tested, "_", 21:40)
    
    LOEUF_Sum_1_DDD = paste0("catg_", DDD_model, "_", score_tested, "_", 1:20 )
    LOEUF_Sum_2_DDD = paste0("catg_", DDD_model, "_", score_tested, "_", 21:40)
    
    List_scores_catg = c(paste0(score_tested, "_Sum_", 1:2))
    List_scores_catg_DDD = c(paste0(score_tested, "_Sum_", 1:2, "_DDD"))
    
    size = c(0, 13)
    
    ######################################
    ####   Boucle pour les analysis   ####
    ######################################
    
    
    for (ll in 1:(length(scores_catg_DDD))) { #- 3)
      
      for (catg_win in List_scores_catg_DDD) {
        catg_win_without_win_13 = paste0(type, ".", get(catg_win))
        Fichier_w = subset(Fichier, select = catg_win_without_win_13)
        Fichier[paste0(type, ".", catg_win)] = rowSums(Fichier_w)
      }
      
      # On fait tournée le model
      list_outside_wind = paste0(type, ".", List_scores_catg_DDD)
      formula_13_win = as.formula(paste(
        "pheno_test ~",
        paste0(type, ".", scores_catg_DDD),
        "+",paste(list_outside_wind, collapse = " + ") ,"+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
      
      model_IQ_Del = lm(formula_13_win,
                        Fichier,
                        na.action = na.omit)
      
      A1 = summary(model_IQ_Del)$coefficients[1:4,1:4]
      
      
      Sum_13 = sum(Fichier[paste0(type, ".", scores_catg_DDD)])
      Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1_DDD")])
      Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2_DDD")])
      
      Sum_ABC = data.frame(
        rbind(
          0,
          Sum_13,
          Sum_1,
          Sum_2
        )
      )
      win_ABC = data.frame(rbind(0, 14, 1, 2))
      
      if (Sum_13 > 0) { C1 = cbind(A1, Sum_ABC, 0, win_ABC,
                                   type,scores_catg_DDD,cohorts[i],paste0(type,"_0"),"Model_2_DDD",pheno_info[i_pheno],"Yes")
      
      }
      
      if (Sum_13 == 0) { C1 = cbind(rbind(A1[1,1:4],c(NA,NA,NA,NA),A1[2:3,1:4]), Sum_ABC, 0, win_ABC,
                                    type,scores_catg_DDD,cohorts[i],paste0(type,"_0"),"Model_2_DDD",pheno_info[i_pheno],"Yes")
      
      }

      colnames(C1) = c("Value",
                       "Std.Error",
                       "t.value",
                       "p.value",
                       "rbind.Sum.",
                       "ll",
                       "catg",
                       "TYPE_CNV","Model_DDD","Cohorts","names","Model","Pheno","PC")
      
      C0_Del_2 = rbind(C0_Del_2, C1)
      
    }
    
    ###################################################################################################################        
    
    for (ll in 1:(length(scores_catg) - 2)) {

      # On créat la win 13 et sa somme
      Fichier[paste0(type, ".", score_tested, "_Sum_13")] = Fichier[, paste0(type, ".", scores_catg[ll])] + Fichier[, paste0(type, ".", scores_catg[ll + 1])] + Fichier[, paste0(type, ".", scores_catg[ll + 2])]
      sum(Fichier[paste0(type, ".", score_tested, "_Sum_13")])
      scores_catg_13 = c(
        paste0(type, ".", scores_catg[ll]),
        paste0(type, ".", scores_catg[ll + 1]),
        paste0(type, ".", scores_catg[ll + 2]))
      
      scores_catg_13_bis = c(
        paste0(type, ".", scores_catg_bis[ll]),
        paste0(type, ".", scores_catg_bis[ll + 1]),
        paste0(type, ".", scores_catg_bis[ll + 2]))
      
      
      # On calcule les autre win sans les valeur utilisé dans la win 13
      
      for (catg_win in List_scores_catg) {
        
        catg_win_without_win_13 = paste0(type, ".", get(catg_win)[paste0(type, ".", get(catg_win)) %!in% scores_catg_13_bis])
        
        if (catg_win_without_win_13[1] != paste0(type, ".")) {
          Fichier_w = subset(Fichier, select = catg_win_without_win_13)
          Fichier[paste0(type, ".", catg_win)] = rowSums(Fichier_w)
        }
      }
      
      # On fait tournée le model
      list_outside_wind = paste0(type, ".", List_scores_catg)
      formula_13_win = as.formula(paste(
        "pheno_test ~",
        paste0(type, ".", score_tested, "_Sum_13"),"+", paste(list_outside_wind, collapse = " + "),"+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
      
      model_IQ_Del = lm(formula_13_win,
                        Fichier,
                        na.action = na.omit)
      
      A1 = summary(model_IQ_Del)$coefficients[1:4,1:4]
      
      
      Sum_13 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_13")])
      Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1")])
      Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2")])
      
      Sum_ABC = data.frame(
        rbind(
          0,
          Sum_13,
          Sum_1,
          Sum_2
        )
      )
      win_ABC = data.frame(rbind(0, 13, 1, 2))
      C1 = cbind(A1, Sum_ABC, ll, win_ABC,
                 type,scores_catg_DDD,cohorts[i],paste0(type,"_",ll),"Model_2_DDD",pheno_info[i_pheno],"Yes")
      
      colnames(C1) = c("Value",
                       "Std.Error",
                       "t.value",
                       "p.value",
                       "rbind.Sum.",
                       "ll",
                       "catg",
                       "TYPE_CNV","Model_DDD","Cohorts","names","Model","Pheno","PC")
      
      C0_Del_2 = rbind(C0_Del_2, C1)
      
    }
    
    Model_C0_Del_2 = data.frame(C0_Del_2)
    Model_C0_Del_2$Names = 0
    Model_C0_Del_2$Names = rownames(Model_C0_Del_2)
    
  }}}

# Without PC in formula
for (i_pheno in 1:length(pheno_info)) {
  
  all_data_info_pheno_scores_SNP_Artefact$pheno_test = all_data_info_pheno_scores_SNP_Artefact[,pheno_info[i_pheno]]
  
  for (i in 1:length(cohorts)) {
    if (cohorts[i] != "ALL") {Fichier = subset(all_data_info_pheno_scores_SNP_Artefact, Cohort_final == cohorts[i] & !(is.na(pheno_test)))}
    if (cohorts[i] == "ALL") {Fichier = all_data_info_pheno_scores_SNP_Artefact}
    
    for (type in type_list) {
      
      scores_catg = paste0("catg_", DDD_model, "_", score_tested, "_", 1:40)
      scores_catg_bis = paste0("catg_", score_tested, "_", 1:40)
      
      scores_catg_DDD = DDD_selection
      
      LOEUF_Sum_1 = paste0("catg_", score_tested, "_", 1:20 )
      LOEUF_Sum_2 = paste0("catg_", score_tested, "_", 21:40)
      
      LOEUF_Sum_1_DDD = paste0("catg_", DDD_model, "_", score_tested, "_", 1:20 )
      LOEUF_Sum_2_DDD = paste0("catg_", DDD_model, "_", score_tested, "_", 21:40)
      
      List_scores_catg = c(paste0(score_tested, "_Sum_", 1:2))
      List_scores_catg_DDD = c(paste0(score_tested, "_Sum_", 1:2, "_DDD"))
      
      size = c(0, 13)
      
      ######################################
      ####   Boucle pour les analysis   ####
      ######################################
      
      
      for (ll in 1:(length(scores_catg_DDD))) {
        
        for (catg_win in List_scores_catg_DDD) {
          catg_win_without_win_13 = paste0(type, ".", get(catg_win))
          Fichier_w = subset(Fichier, select = catg_win_without_win_13)
          Fichier[paste0(type, ".", catg_win)] = rowSums(Fichier_w)
        }
        
        # On fait tournée le model
        list_outside_wind = paste0(type, ".", List_scores_catg_DDD)
        formula_13_win = as.formula(paste(
          "pheno_test ~",
          paste0(type, ".", scores_catg_DDD),
          "+",paste(list_outside_wind, collapse = " + ")))
        
        model_IQ_Del = lm(formula_13_win,
                          Fichier,
                          na.action = na.omit)
        
        A1 = summary(model_IQ_Del)$coefficients
        
        Sum_13 = sum(Fichier[paste0(type, ".", scores_catg_DDD)])
        Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1_DDD")])
        Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2_DDD")])
        
        Sum_ABC = data.frame(
          rbind(
            0,
            Sum_13,
            Sum_1,
            Sum_2
          )
        )
        win_ABC = data.frame(rbind(0, 14, 1, 2))
        
        if (Sum_13 > 0) { C1 = cbind(A1, Sum_ABC, 0, win_ABC,
                                     type,scores_catg_DDD,cohorts[i],paste0(type,"_0"),"Model_2_DDD",pheno_info[i_pheno],"No")
        
        }
        
        if (Sum_13 == 0) { C1 = cbind(rbind(A1[1,1:4],c(NA,NA,NA,NA),A1[2:3,1:4]), Sum_ABC, 0, win_ABC,
                                      type,scores_catg_DDD,cohorts[i],paste0(type,"_0"),"Model_2_DDD",pheno_info[i_pheno],"No")
        
        }
        colnames(C1) = c("Value",
                         "Std.Error",
                         "t.value",
                         "p.value",
                         "rbind.Sum.",
                         "ll",
                         "catg",
                         "TYPE_CNV","Model_DDD","Cohorts","names","Model","Pheno","PC")
        
        C0_Del_2 = rbind(C0_Del_2, C1)
        
      }
      
      ###################################################################################################################        
      ###################################################################################################################
      ###################################################################################################################   
      
      for (ll in 1:(length(scores_catg) - 2)) { #- 3))

        # On créat la win 13 et sa somme
        Fichier[paste0(type, ".", score_tested, "_Sum_13")] = Fichier[, paste0(type, ".", scores_catg[ll])] + Fichier[, paste0(type, ".", scores_catg[ll + 1])] + Fichier[, paste0(type, ".", scores_catg[ll + 2])]
        sum(Fichier[paste0(type, ".", score_tested, "_Sum_13")])
        scores_catg_13 = c(
          paste0(type, ".", scores_catg[ll]),
          paste0(type, ".", scores_catg[ll + 1]),
          paste0(type, ".", scores_catg[ll + 2]))
        
        scores_catg_13_bis = c(
          paste0(type, ".", scores_catg_bis[ll]),
          paste0(type, ".", scores_catg_bis[ll + 1]),
          paste0(type, ".", scores_catg_bis[ll + 2]))
        
        
        # On calcule les autre win sans les valeur utilisé dans la win 13
        
        for (catg_win in List_scores_catg) {
          
          catg_win_without_win_13 = paste0(type, ".", get(catg_win)[paste0(type, ".", get(catg_win)) %!in% scores_catg_13_bis])
          
          if (catg_win_without_win_13[1] != paste0(type, ".")) {
            Fichier_w = subset(Fichier, select = catg_win_without_win_13)
            Fichier[paste0(type, ".", catg_win)] = rowSums(Fichier_w)
          }
        }
        
        # On fait tournée le model
        list_outside_wind = paste0(type, ".", List_scores_catg)
        formula_13_win = as.formula(paste(
          "pheno_test ~",
          paste0(type, ".", score_tested, "_Sum_13"),
          "+",paste(list_outside_wind, collapse = " + ")))
        
        model_IQ_Del = lm(formula_13_win,
                          Fichier,
                          na.action = na.omit)
        
        A1 = summary(model_IQ_Del)$coefficients
        
        Sum_13 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_13")])
        Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1")])
        Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2")])
        
        Sum_ABC = data.frame(
          rbind(
            0,
            Sum_13,

            Sum_1,
            Sum_2
          )
        )
        
        win_ABC = data.frame(rbind(0, 13, 1, 2))
        C1 = cbind(A1, Sum_ABC, ll, win_ABC,
                   type,scores_catg_DDD,cohorts[i],paste0(type,"_",ll),"Model_2_DDD",pheno_info[i_pheno],"No")
        
        colnames(C1) = c("Value",
                         "Std.Error",
                         "t.value",
                         "p.value",
                         "rbind.Sum.",
                         "ll",
                         "catg",
                         "TYPE_CNV","Model_DDD","Cohorts","names","Model","Pheno","PC")
        
        C0_Del_2 = rbind(C0_Del_2, C1)
        
      }
      
      Model_C0_Del_2 = data.frame(C0_Del_2)
      Model_C0_Del_2$Names = 0
      Model_C0_Del_2$Names = rownames(Model_C0_Del_2)
      
    }}}

###########################################

model_IQ_DelDup = Model_C0_Del_2

model_IQ_DelDup_ALL = subset(model_IQ_DelDup, Cohorts == "ALL" & (catg == "13" | catg == "14") & TYPE_CNV == "DEL")#& Anc == "ALL")

model_IQ_DelDup_ALL$Pheno_TYPE_CNV = paste(model_IQ_DelDup_ALL$Pheno, model_IQ_DelDup_ALL$TYPE_CNV,"PC", model_IQ_DelDup_ALL$PC, sep = "_")

pd = position_dodge(0)

ggplot(model_IQ_DelDup_ALL, aes(x=ll, y=as.numeric(Value), colour=Pheno_TYPE_CNV)) +labs(x ="LOEUF") + xlim(1,40) + ylim(-0.3,0.07) +

  geom_line(position=pd, size=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0.35, linetype="dashed", color = "red") +
  geom_point(position=pd, size=ifelse(model_IQ_DelDup_ALL$p.value < 0.05, 5, 3),
             shape=ifelse(model_IQ_DelDup_ALL$p.value < 0.05 , 23,
                          ifelse(model_IQ_DelDup_ALL$p.value >= 0.05 , 9, 7)), fill="white") +
  theme(legend.position = c(0.8, 0.2))

write.table(model_IQ_DelDup, "META_Model_with_DDD_ClinGen_DelDup_2024.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

#################################################################