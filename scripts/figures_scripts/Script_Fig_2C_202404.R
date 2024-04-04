####################################################################
####################################################################
# Fig 2C
####################################################################
####################################################################

library(ggplot2)
library(ggnewscale)

#################################################################
# Set directory/Import files 
#################################################################

setwd(paste0("/Users/guillaume/Desktop/File_Script/Fig_1_2/"))

#################################################################

Model_C0_DelDup_all_final = read.csv(
  paste0("Model_with_DDD_ClinGen_DelDup_sensitivity_catg_2024.tsv"),
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

Model_C0_DelDup_select_all_final = subset(Model_C0_DelDup_all_final, (catg == 13 | catg == 14))

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

#################################################################

Model_C0_DelDup_2$catg_info = paste(Model_C0_DelDup_2$ASD, Model_C0_DelDup_2$Anc, Model_C0_DelDup_2$Age, Model_C0_DelDup_2$scores_catg, sep = "_")

####################################

Model_C0_DelDup_2_test = data.frame(subset(Model_C0_DelDup_2, (catg_info == "with_ASD_ALL_12000_10000"
                                                                           | catg_info == "with_ASD_EUR_12000_10000"
                                                                           | catg_info == "with_ASD_ALL_12000_60"
                                                                           | catg_info == "with_ASD_ALL_12000_40"
                                                                           | catg_info == "with_ASD_ALL_12000_20"
                                                                           | catg_info == "with_ASD_ALL_900_10000"
                                                                           | catg_info == "with_ASD_ALL_840_10000"
                                                                           | catg_info == "with_ASD_ALL_780_10000"
                                                                           | catg_info == "with_ASD_ALL_720_10000"
                                                                           | catg_info == "with_ASD_ALL_660_10000"
                                                                           
                                                                           | catg_info == "without_ASD_ALL_12000_10000"
                                                                           | catg_info == "without_ASD_EUR_12000_10000"
                                                                           | catg_info == "without_ASD_ALL_12000_60"
                                                                           | catg_info == "without_ASD_ALL_12000_40"
                                                                           | catg_info == "without_ASD_ALL_12000_20"
                                                                           | catg_info == "without_ASD_ALL_900_10000"
                                                                           | catg_info == "without_ASD_ALL_840_10000"
                                                                           | catg_info == "without_ASD_ALL_780_10000"
                                                                           | catg_info == "without_ASD_ALL_720_10000"
                                                                           | catg_info == "without_ASD_ALL_660_10000")))



Model_C0_DelDup_2_test$catg_info_plus = ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_12000_10000", "29_with_ASD_ALL_12000_10000" ,
                                          ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_EUR_12000_10000", "28_with_ASD_EUR_12000_10000" ,
                                                 
                                                 ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_12000_60", "22_with_ASD_ALL_12000_60" ,
                                                        ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_12000_40", "21_with_ASD_ALL_12000_40" ,
                                                               ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_12000_20", "20_with_ASD_ALL_12000_20" ,
                                                                      
                                                                      ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_900_10000", "27_with_ASD_ALL_900_10000" ,
                                                                             ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_840_10000", "26_with_ASD_ALL_840_10000" ,
                                                                                    ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_780_10000", "25_with_ASD_ALL_780_10000" ,
                                                                                           ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_720_10000", "24_with_ASD_ALL_720_10000" ,
                                                                                                  ifelse(Model_C0_DelDup_2_test$catg_info == "with_ASD_ALL_660_10000", "23_with_ASD_ALL_660_10000" ,
                                            
                                                                                                         ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_12000_10000", "19_without_ASD_ALL_12000_10000" ,
                                                                                                                ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_EUR_12000_10000", "18_without_ASD_EUR_12000_10000" ,
                                                                                                                       
                                                                                                                       ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_12000_60", "12_without_ASD_ALL_12000_60" ,
                                                                                                                              ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_12000_40", "11_without_ASD_ALL_12000_40" ,
                                                                                                                                     ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_12000_20", "10_without_ASD_ALL_12000_20" ,
                                                                                                                                            
                                                                                                                                            ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_900_10000", "17_without_ASD_ALL_900_10000" ,
                                                                                                                                                   ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_840_10000", "16_without_ASD_ALL_840_10000" ,
                                                                                                                                                          ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_780_10000", "15_without_ASD_ALL_780_10000" ,
                                                                                                                                                                 ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_720_10000", "14_without_ASD_ALL_720_10000" ,
                                                                                                                                                                        ifelse(Model_C0_DelDup_2_test$catg_info == "without_ASD_ALL_660_10000", "13_without_ASD_ALL_660_10000", "0"))))))))))))))))))))


###############################################

ggplot(mapping = aes(x = ll, y = catg_info_plus, size = log10(carriers.Sum.))) + facet_wrap(~ TYPE_CNV) + theme(panel.background = element_rect(fill="transparent")) +
  # This bit is for making scales
  scale_radius(range = c(2, 9)) +
  geom_point(data=subset(Model_C0_DelDup_2_test, Estimate <= 0), aes(fill=Estimate), 
             shape = 21,
             col = paste(ifelse((subset(Model_C0_DelDup_2_test, Estimate <= 0)$p.value_FDR) < 0.05, "black", "grey")),
             stroke = as.numeric(paste(ifelse((subset(Model_C0_DelDup_2_test, Estimate <= 0)$p.value_FDR) < 0.05, 1, 0.2)))) +
  #scale_fill_gradient2(low="darkblue", mid="skyblue3", high="white", midpoint = -1) +
  scale_fill_gradient2(low="purple4", mid="purple", high="white", midpoint = -1) +
  
  new_scale_fill() +
  geom_point(data=subset(Model_C0_DelDup_2_test, Estimate > 0) , aes(fill = Estimate),
             shape = 21,
             col = paste(ifelse((subset(Model_C0_DelDup_2_test, Estimate > 0)$p.value_FDR) < 0.05, "black", "grey")),
             stroke = as.numeric(paste(ifelse((subset(Model_C0_DelDup_2_test, Estimate > 0)$p.value_FDR) < 0.05 , 1, 0.2)))) +
  scale_fill_gradient2(low = "white", mid="darkorange", high = "red", midpoint = 0.05) +
  new_scale_fill() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1 , size =12)) +
  geom_hline(yintercept = seq(0.5, 20, 1), size = .2) +
  
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

###############################################