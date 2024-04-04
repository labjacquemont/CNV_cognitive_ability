library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)

#################################################################

'%!in%' <- Negate('%in%')
user = "guillaf"

#################################################################
# Set directory/Import files
#################################################################

setwd("/project/rrg-jacquese/All_user_common_folder/Cohorts_sum_of_scores_CNV/Cohorts_sum_of_scores_CNV/DEL_DUP_Intell_article_all_Cohorts/All_cohorts_filtre_artefact_SNP/sum_score_gene_overlap")

List_info_FID_IID_largeCNV_0 = read.csv(
  "Summary_Mega_MSSNG_SPARK_UKBB_Pheno_final_with_largeCNV_T21_ChrX_20210619.dat",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)

Ancestry_info_0 = read.csv(
  "Mega_MSSNG_SPARK_UKBB_MSSNG_Pheno_final_20220606_Geno_ancestry_catg_details.txt",
  #"Mega_MSSNG_SPARK_UKBB_Pheno_final_20210924_Geno_ancestry_catg_details.txt",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)

#summary(as.factor(Ancestry_info_0$Ancestry))

Ancestry_info_0$Merge_final_ancestry = Ancestry_info_0$Ancestry
Ancestry_info_0$Catg_ancestry = "KING"

List_info_FID_IID_largeCNV = merge(
  List_info_FID_IID_largeCNV_0,
  subset(
    Ancestry_info_0,
    select = c(
      "individual",
      "PC1",
      "PC2",
      "PC3",
      "PC4",
      "PC5",
      "PC6",
      "PC7",
      "PC8",
      "PC9",
      "PC10",
      "Merge_final_ancestry",
      "Catg_ancestry"
    )
  ),
  by = c("individual")
)

###########################################################
# new pheno ###############################################

New_pheno_FID_IID = read.csv(
  "/project/rrg-jacquese/All_user_common_folder/Cohorts_sum_of_scores_CNV/Cohorts_sum_of_scores_CNV/Mega_MSSNG_SPARK_UKBB_Pheno_final_PC_20221226.tsv",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)

New_pheno_FID_IID = subset(New_pheno_FID_IID, select =c("individual", "ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final"))

###########################################################

List_info_FID_IID_largeCNV = merge(List_info_FID_IID_largeCNV, New_pheno_FID_IID, by = c("individual"))

#################################################################

Model_C0_DelDup_all_final = data.frame()
Model_C0_DelDup_all = data.frame()

list_file = c("Score_30_SNP_10", "Score_30_SNP_15", "Score_30_SNP_30", "Score_30_SNP_50", "Score_30_SNP_100",
              "Score_40_SNP_10", "Score_50_SNP_10", "Score_75_SNP_10", "Score_100_SNP_10")


list_DDD_complete = c("DDD_ClinGen")


#for (list_DDD in list_DDD_complete) {
list_DDD = list_DDD_complete[1]
print(list_DDD)

for (file_name in list_file) {

    print(file_name)

all_data_info_pheno_scores_SNP_Artefact = fread(
  paste0("/home/",user,"/projects/rrg-jacquese/guillaf/0_scratch_2022/UKBB-MEGA-SPARK-MSSNG_Sum_scores_overlap_norecip_50_1freq_no_pool_CNV_",file_name,"_20230419_catg_LOEUF.tsv"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = "."
)

type_list = c("DEL", "DUP")
# Create list of score columns
score_tested = "LOEUF"
#gene_overlap_cat = "full"
list_score = c(
  paste0("DUP.catg_without",list_DDD,"_",score_tested,"_", 1:40),
  paste0("DEL.catg_without",list_DDD,"_",score_tested,"_", 1:40),
  paste0("DUP.catg_",score_tested,"_", 1:40),
  paste0("DEL.catg_",score_tested,"_", 1:40)
)
type_list = c("DEL", "DUP")

window_list = c(0.175,
                0.45,
                0.6,
                0.75,
                0.9,
                1.05,
                1.2,
                1.35,
                1.5,
                1.65,
                1.8,
                1.95)

all_data_info_pheno_scores_SNP_Artefact = merge(
  List_info_FID_IID_largeCNV,
  subset(
    all_data_info_pheno_scores_SNP_Artefact,
    select = c(
      "individual",
      "ZScore_IQ_adj_test_age_sex_PC1_10_with_online_UKBB_Final",
      "Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB",
      "ZScore_IQ_adj_age",
      "ZScore_IQ_adj_age_sex",
      "ZScore_IQ_adj_age_PC",
      "ZScore_IQ_adj_age_sex_PC",
      "z_score_FI_all",
      "Age_of_FI_all",
      "z_score_FI_online_all",
      "Age_of_FI_all_online",
      "Z_Gfact_5_Donald_et_al_2016",
      "Age_of_Gfac5_all_Donald_et_al_2016",
      "Z_Gfact_4_Cox_et_al_2019",
      "f.53.Gfactor_2and3.Age_of_all_Cox_et_al_2019",
      "Z_Gfact_5_online_clinic",
      "Age_of_Gfac5_2_all_online_clinic",
      "Z_Gfact_5_online",
      "z_score_FI_all_adj_age_sex",
      "z_score_FI_all_adj_age",
      "z_score_FI_all_adj_age_sex_PC1_10",
      "z_score_FI_all_adj_age_PC1_10",
      "z_score_FI_online_all_adj_age_sex",
      "z_score_FI_online_all_adj_age",
      "z_score_FI_online_all_adj_age_sex_PC1_10",
      "z_score_FI_online_all_adj_age_PC1_10",
      "Z_Gfact_5_Donald_et_al_2016_adj_age_sex",
      "Z_Gfact_5_Donald_et_al_2016_adj_age",
      "Z_Gfact_5_Donald_et_al_2016_adj_age_sex_PC1_10",
      "Z_Gfact_5_Donald_et_al_2016_adj_age_PC1_10",
      "Z_Gfact_4_Cox_et_al_2019_adj_age_sex",
      "Z_Gfact_4_Cox_et_al_2019_adj_age",
      "Z_Gfact_4_Cox_et_al_2019_adj_age_sex_PC1_10",
      "Z_Gfact_4_Cox_et_al_2019_adj_age_PC1_10",
      "Z_Gfact_5_online_clinic_adj_age_sex",
      "Z_Gfact_5_online_clinic_adj_age",
      "Z_Gfact_5_online_clinic_adj_age_sex_PC1_10",
      "Z_Gfact_5_online_clinic_adj_age_PC1_10",
      "Z_Gfact_5_online_adj_age_sex",
      "Z_Gfact_5_online_adj_age",
      "Z_Gfact_5_online_adj_age_sex_PC1_10",
      "Z_Gfact_5_online_adj_age_PC1_10",
      "ZScore_IQ_adj_age_with_online_UKBB",
      "Type_ZScore_IQ_adj_age_with_online_UKBB",
      "ZScore_IQ_adj_age_sex_with_online_UKBB",
      "Type_ZScore_IQ_adj_age_sex_with_online_UKBB",
      "ZScore_IQ_adj_age_PC1_10_with_online_UKBB",
      "Type_ZScore_IQ_adj_age_PC1_10_with_online_UKBB",
      "ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB",
      "Type_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB",
      "ZScore_IQ_adj_age_sex_with_online_UKBB_Final",
      "Type_ZScore_IQ_adj_age_sex_with_online_UKBB_Final",
      "ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB_Final",
      "Type_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB_Final",
      #"Gene_freq",
      "DUP.loeuf_inv",
      "DEL.loeuf_inv",
      list_score, paste0("DEL.",list_DDD), paste0("DUP.",list_DDD)
      #"PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Merge_final_ancestry", "Catg_ancestry"
    )
  )
  ,
  by = "individual"
)


all_data_info_pheno_scores_SNP_Artefact = subset(
  all_data_info_pheno_scores_SNP_Artefact,
  DUP_10Mb == 0 &
    DEL_10Mb == 0 &
    Dup_chrX == 0 &
    Large_dup_X == 0 & large_Dup_chr21 == 0 & T21 == 0
)

#########################################################################################################################
#########################################################################################################################
removed_individuals = c(
  "SP0107612",
  "SP0080521",
  "2284951_Batch_b058",
  "3952562_UKBiLEVEAX_b3",
  "5377011_Batch_b057",
  "11100837",
  "11101241",
  "11101860",
  "11102913",
  "11106217",
  "11106812",
  "11108960",
  "97178300763",
  "97178301491",
  "97178301690",
  "97178302178",
  "2904393-610Kq",
  "4584699124_R01C01",
  "5370019131_R01C01",
  "6042279056_R02C01",
  "6057809079_R02C01",
  "6312937069_R02C01",
  "7340542076_R02C01",
  "7956321114_R04C01"
)

removed_sub_cohorts = c(
  "MSSNG-father",
  "MSSNG-mother",
  "MSSNG-sibling",
  "SPARK-father",
  "SPARK-mother",
  "SPARK-sibling",
  "SSC-father",
  "SSC-mother",
  "SSC-sibling"
)

all_data_info_pheno_scores_SNP_Artefact = subset(
  all_data_info_pheno_scores_SNP_Artefact,!(
    is.na(ZScore_IQ_adj_test_age_sex_PC1_10_with_online_UKBB_Final)
  ) &
    individual %!in% removed_individuals &
    Sub_cohort %!in% removed_sub_cohorts
)

##################################################################################################################################

list_individual_duplicat = subset(
  fread(
    paste0("/home/",user,"/projects/rrg-jacquese/All_user_common_folder/Cohorts_sum_of_scores_CNV/Update_file_pheno_CNV_duplicate_20220609/Duplication_Total_all_cohorts_final_2022_IQ.txt"),
    header = TRUE,
    sep = "\t",
    na.strings = c("NA", "#NULL!", ".", ""),
    dec = "."
  ) ,
  select = c(individual_remove)
)

list_individual_duplicat$remove = "remove"

all_data_info_pheno_scores_SNP_Artefact = subset(
  merge(
    all_data_info_pheno_scores_SNP_Artefact,
    list_individual_duplicat,
    by.x = "individual",
    by.y = "individual_remove",
    all.x = T
  ),
  is.na(remove)
)

###################################################################################################
###################################################################################################
###################################################################################################

all_data_info_pheno_scores_SNP_Artefact$Cohort_final = ifelse(
  all_data_info_pheno_scores_SNP_Artefact$Cohort == "UKBB",
  all_data_info_pheno_scores_SNP_Artefact$Type_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB_Final,
  all_data_info_pheno_scores_SNP_Artefact$Cohort
)

all_data_info_pheno_scores_SNP_Artefact$Cohort_final_bis = ifelse(
  all_data_info_pheno_scores_SNP_Artefact$Cohort_final == "Z_score_FI" |
    all_data_info_pheno_scores_SNP_Artefact$Cohort_final == "Z_score_FI_online" ,
  "Z_score_FI_all",
  all_data_info_pheno_scores_SNP_Artefact$Cohort_final
)

summary(as.factor(all_data_info_pheno_scores_SNP_Artefact$Cohort_final))
summary(as.factor(all_data_info_pheno_scores_SNP_Artefact$Cohort_final_bis))

all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB = ifelse( is.na(all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB) , -9 ,
                                                                                                    all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB)


Fichier_0 = all_data_info_pheno_scores_SNP_Artefact
Fichier_0$pheno_test = Fichier_0$ZScore_IQ_adj_test_age_sex_PC1_10_with_online_UKBB_Final
summary(as.factor(Fichier_0$Cohort))

Fichier_0$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB = ifelse( is.na(Fichier_0$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB) , -9 , Fichier_0$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB)

#fwrite(all_data_info_pheno_scores_SNP_Artefact,"all_data_info_pheno_scores_SNP_Artefact_select_cohorts_20220406.tsv", quote=F, sep="\t", row.names=T, col.names=T)

# POOL wind and Cohorts

#Fichier_test = subset(Fichier, !(is.na(pheno_test)))
#windows_data_frame = data.frame()

#################################################################

score_cutoff_LOEUF = c(10000)#,60,40,20,10)#,4,2.85)

age_value = c(100*120)#,75*12,70*12,65*12,60*12,55*12)

Files_pheno_score = list(all_data_info_pheno_scores_SNP_Artefact)

Name_files_pheno_score = c("all_data_info_pheno_scores_SNP_Artefact")

Pheno_clinic_online = c("ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final") #,"ZScore_IQ_adj_test_age_sex_PC1_10_with_online_UKBB_Final")

ancestry = c("ALL", "EUR")#,"No_EUR")

Diagnostic_asd = c("without_ASD")#,"with_ASD")#,"only_ASD")

#Fichier_0 = Fichier
#Fichier_0 = subset(Fichier, Cohort == "UKBB")

Fichier_0$DEL.loeuf_inv_full = Fichier_0$DEL.loeuf_inv
Fichier_0$DUP.loeuf_inv_full = Fichier_0$DUP.loeuf_inv

# # Create list of score columns
# score_tested = "LOEUF"
#gene_overlap_cat = "full"
# list_score_ok = c(
#   paste0("DUP.catg_",score_tested,"_", 1:40,"_",gene_overlap_cat),
#   paste0("DEL.catg_",score_tested,"_", 1:40,"_",gene_overlap_cat)
# )
#
#
#
# for (i in 1:length(list_score)) {
#   #i=1
#   Fichier_0[list_score_ok[i]] = Fichier_0[list_score[i]]
# }

DDD_model = paste0("without",list_DDD)
DDD_selection = list_DDD
C0_Del_2 = c()
Model_C0_Del_2 = c()
Model_C0_DelDup_all_final = c()
#Model_C0_DelDup_all = c()

#summary(Fichier_0$PC2)

#Model_C0_DelDup_2
#################################################################

for (anc_code in 1:length(ancestry)) {
  #anc_code = 1
  print(ancestry[anc_code])

  if (ancestry[anc_code] == "ALL") { Fichier_1 = Fichier_0}

  if (ancestry[anc_code] == "EUR") { Fichier_1 = subset(Fichier_0, Merge_final_ancestry == ancestry[anc_code]) }

  if (ancestry[anc_code] == "No_EUR") { Fichier_1 = subset(Fichier_0, Merge_final_ancestry != "EUR") }

  for (i in 1:length(Diagnostic_asd)) {
    print(Diagnostic_asd[i])
    #i=1
    if (Diagnostic_asd[i] == "without_ASD") { Fichier_2 =  subset(Fichier_1, !(Cohort_final == "MSSNG" | Cohort_final == "SSC" | Cohort_final == "SPARK"))}

    if (Diagnostic_asd[i] == "only_ASD") { Fichier_2 =  subset(Fichier_1, (Cohort_final == "MSSNG" | Cohort_final == "SSC" | Cohort_final == "SPARK"))}

    if (Diagnostic_asd[i] == "with_ASD") { Fichier_2 =  Fichier_1}

    for (l in 1:length(age_value)) {
      print(age_value[l])
      #l=1
      Fichier_3 = subset(Fichier_2, Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB <= age_value[l])

      for (type in type_list) {
        #type = c("DEL")
        Fichier_3 = subset(Fichier_3,!(is.na(pheno_test)))

        for (j in 1:length(score_cutoff_LOEUF)) {
          #j=1
          if (type == "DEL") {
            Fichier = subset(Fichier_3, DEL.loeuf_inv_full < score_cutoff_LOEUF[j])
          }

          if (type == "DUP") {
            Fichier = subset(Fichier_3, DUP.loeuf_inv_full < score_cutoff_LOEUF[j])
          }

          print(dim(Fichier)[1])

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
              "+",# paste0(type, ".DDD")  ,"+",
              paste(list_outside_wind, collapse = " + ") ,
              "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
            ))

            model_IQ_Del = lm(formula_13_win,
                              Fichier,
                              na.action = na.omit)

            A1 = summary(model_IQ_Del)$coefficients[1:4,]


            Sum_13 = sum(Fichier[paste0(type, ".", scores_catg_DDD)])
            Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1_DDD")])
            Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2_DDD")])

            Sum_13_IID = dim( subset(Fichier, Fichier[paste0(type, ".", scores_catg_DDD)] >=1))[1]
            Sum_1_IID = dim( subset(Fichier, Fichier[paste0(type, ".", score_tested, "_Sum_1_DDD")] >=1))[1]
            Sum_2_IID = dim( subset(Fichier, Fichier[paste0(type, ".", score_tested, "_Sum_2_DDD")] >=1))[1]

            Sum_ABC = data.frame(
              rbind(
                0,
                Sum_13,
                #Sum_DDD,
                Sum_1,
                Sum_2#,
                #Sum_3#,
                #Sum_4,
                #Sum_5,
                #Sum_6#,
                # Sum_7,
                # Sum_8,
                # Sum_9,
                # Sum_10,
                # Sum_11,
                # Sum_12
              )
            )

            Sum_ABC_IID = data.frame(
              rbind(
                0,
                Sum_13_IID,
                Sum_1_IID,
                Sum_2_IID
              )
            )

            win_ABC = data.frame(rbind(0, 14, 1, 2))# ,14, 3))#, 4, 5, 6))#, 7, 8, 9, 10, 11, 12))

            if (Sum_13 > 0) { C1 = cbind(A1, Sum_ABC, Sum_ABC_IID, 0, win_ABC,
                                         score_cutoff_LOEUF[j],age_value[l],ancestry[anc_code],Diagnostic_asd[i],dim(Fichier)[1],type,scores_catg_DDD)

            }

            if (Sum_13 == 0) { C1 = cbind(rbind(A1[1,1:4],c(NA,NA,NA,NA),A1[2:3,1:4]), Sum_ABC, Sum_ABC_IID, 0, win_ABC,
                                          score_cutoff_LOEUF[j],age_value[l],ancestry[anc_code],Diagnostic_asd[i],dim(Fichier)[1],type,scores_catg_DDD)

            }
            #row.names(C1) = c(paste0("DEL.",scores_catg[ll+1],"_",ll+2,"_",ll+3),paste0("DEL._40_without_",ll+1,"_",ll+2,"_",ll+3))
            colnames(C1) = c("Value",
                             "Std.Error",
                             "t.value",
                             "p.value",
                             "rbind.Sum.",
                             "carriers.Sum.",
                             "ll",
                             "catg",
                             "scores_catg","Age","Anc","ASD","N","TYPE_CNV","Model_DDD")

            C0_Del_2 = rbind(C0_Del_2, C1)

          }

          ###################################################################################################################
          ###################################################################################################################
          ###################################################################################################################

          for (ll in 1:(length(scores_catg) - 2)) { #- 3))
            #ll=1
            #print(scores_catg[ll])
            #type = "DEL"
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

              # if (catg_win_without_win_13 ==paste0(type,".")) {
              #   Fichier[paste0("DEL.", catg_win)] = 0
              # }

              #catg_win_without_win_13 = paste0("DEL.",scores_catg_1[paste0("DEL.",scores_catg_1) %!in% scores_catg_13])
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
              "+",# paste0(type, ".DDD")  ,"+",
              paste(list_outside_wind, collapse = " + "),
              "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
            ))

            model_IQ_Del = lm(formula_13_win,
                              Fichier,
                              na.action = na.omit)

            A1 = summary(model_IQ_Del)$coefficients[1:4,]


            Sum_13 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_13")])
            #Sum_DDD = sum(Fichier[paste0(type, ".DDD")])
            Sum_1 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_1")])
            Sum_2 = sum(Fichier[paste0(type, ".", score_tested, "_Sum_2")])

            Sum_13_IID = dim( subset(Fichier, Fichier[paste0(type, ".", score_tested, "_Sum_13")] >=1))[1]
            Sum_1_IID = dim( subset(Fichier, Fichier[paste0(type, ".", score_tested, "_Sum_1_DDD")] >=1))[1]
            Sum_2_IID = dim( subset(Fichier, Fichier[paste0(type, ".", score_tested, "_Sum_2_DDD")] >=1))[1]

            Sum_ABC = data.frame(
              rbind(
                0,
                Sum_13,
                #Sum_DDD,
                Sum_1,
                Sum_2#,
                #Sum_3#,
                #Sum_4,
                #Sum_5,
                #Sum_6#,
                # Sum_7,
                # Sum_8,
                # Sum_9,
                # Sum_10,
                # Sum_11,
                # Sum_12
              )
            )

            Sum_ABC_IID = data.frame(
              rbind(
                0,
                Sum_13_IID,
                Sum_1_IID,
                Sum_2_IID
              )
            )

            win_ABC = data.frame(rbind(0, 13, 1, 2))# ,14, 3))#, 4, 5, 6))#, 7, 8, 9, 10, 11, 12))

            C1 = cbind(A1, Sum_ABC,Sum_ABC_IID, ll, win_ABC,
                       score_cutoff_LOEUF[j],age_value[l],ancestry[anc_code],Diagnostic_asd[i],dim(Fichier)[1],type,scores_catg_DDD)

            #row.names(C1) = c(paste0("DEL.",scores_catg[ll+1],"_",ll+2,"_",ll+3),paste0("DEL._40_without_",ll+1,"_",ll+2,"_",ll+3))
            colnames(C1) = c("Value",
                             "Std.Error",
                             "t.value",
                             "p.value",
                             "rbind.Sum.",
                             "carriers.Sum.",
                             "ll",
                             "catg",
                             "scores_catg","Age","Anc","ASD","N","TYPE_CNV","Model_DDD")

            C0_Del_2 = rbind(C0_Del_2, C1)

          }

          Model_C0_Del_2 = data.frame(C0_Del_2)
          Model_C0_Del_2$Names = 0
          Model_C0_Del_2$Names = rownames(Model_C0_Del_2)

          # pdf(paste0(type, "13_win.pdf"))
          # plot((subset(Model_C0_Del_2, catg == 13))$ll, (subset(Model_C0_Del_2, catg == 13))$Value, ylim = c(-0.37, 0.15))
          # dev.off()
          # fwrite(
          #   Model_C0_Del_2,
          #   "C0_Del_Dup_step1.tsv",
          #   quote = F,
          #   sep = "\t",
          #   row.names = T,
          #   col.names = T
          # )

        }}
    }}}


Model_C0_DelDup_2 = Model_C0_Del_2

#################################################################

Model_C0_DelDup_2$Model =  "Model_02"
Model_C0_DelDup_2$file_name_filtre = file_name

Model_C0_DelDup_all = rbind(Model_C0_DelDup_all, Model_C0_DelDup_2)

}

Model_C0_DelDup_all$Estimate = as.numeric(Model_C0_DelDup_all$Value)
Model_C0_DelDup_all$Std.Error = as.numeric(Model_C0_DelDup_all$Std.Error)
Model_C0_DelDup_all$Age = as.numeric(Model_C0_DelDup_all$Age)
Model_C0_DelDup_all$scores_catg = as.numeric(Model_C0_DelDup_all$scores_catg)
Model_C0_DelDup_all$p.value = as.numeric(Model_C0_DelDup_all$p.value)
Model_C0_DelDup_all$Size = Model_C0_DelDup_all$ll*0.05 + 0.025


Model_C0_DelDup_all_final = rbind(Model_C0_DelDup_all_final, Model_C0_DelDup_all)

fwrite(Model_C0_DelDup_all_final,paste0("Model_with_DDD_ClinGen_DelDup_sensitivity_SCORE_SNP_2024.tsv"), quote=F, sep="\t", row.names=T, col.names=T)