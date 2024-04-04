library(stringr)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(nlme)
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

list_file = c(#"CNVfiltered_overlapartifact_Freq1_noreciprocitie_50_UKBBMEGASPARKMSSNG_Score_30_SNP_10_20220527_catg_LOEUF_DDD"
  # "CNVfiltered_overlapartifact_Freq1_noreciprocitie_50_UKBBMEGASPARKMSSNG_Score_30_SNP_10_20220629_catg_LOEUF_DDD"
  # "CNVfiltered_overlapartifact_Freq1_noreciprocitie_50_UKBBMEGASPARKMSSNG_Score_30_SNP_10_20220707_catg_LOEUF_DDD_LOEUF_INV"
  "UKBB-MEGA-SPARK-MSSNG_Sum_scores_overlap_norecip_50_1freq_no_pool_CNV_Score_30_SNP_10_20230419_catg_LOEUF"
)

# list_DDD_complete = c("DDD_old_35", "DDD_old_ID_Satt_35", "DDD_old_ASD_ID_Satt_35", "DDD_2022_Jade_35", "DDD_old_ClinGen_35", "DDD_old_ID_Satt_ClinGen_35", "DDD_old_ASD_ID_Satt_ClinGen_35", "DDD_ClinGen_35",
#                       "DDD_PanelApp_35", "DDD_Fulgent_35",
#                       "DDD_old", "DDD_old_ID_Satt", "DDD_old_ASD_ID_Satt", "DDD_2022_Jade", "DDD_old_ClinGen", "DDD_old_ID_Satt_ClinGen", "DDD_old_ASD_ID_Satt_ClinGen", "DDD_ClinGen", "DDD_PanelApp", "DDD_Fulgent")


list_DDD_complete = c("DDD_ClinGen"#"DDD_old_35", "DDD_old_ID_Satt_35", "DDD_old_ASD_ID_Satt_35", "DDD_2022_Jade_35", "DDD_old_ClinGen_35", "DDD_old_ID_Satt_ClinGen_35", "DDD_old_ASD_ID_Satt_ClinGen_35", "DDD_ClinGen_35", "DDD_PanelApp_35", "DDD_Fulgent_35",
                      #"DDD_old", "DDD_old_ID_Satt", "DDD_old_ASD_ID_Satt", "DDD_2022_Jade", "DDD_old_ClinGen", "DDD_old_ID_Satt_ClinGen", "DDD_old_ASD_ID_Satt_ClinGen", "DDD_ClinGen", "DDD_PanelApp", "DDD_Fulgent",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_35", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_35",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_35_Falling", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_35_Rising", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_35_Non_transitional",
                      #"DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_35_Falling", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_35_Rising", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_35_Non_transitional",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_Falling", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_Rising", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_Non_transitional",
                      #"DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_Falling", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_Rising", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_Non_transitional",
                      #"Gene_DDD_2022_Jade_35_Falling", "Gene_DDD_2022_Jade_35_Rising", "Gene_DDD_2022_Jade_35_Non_transitional",
                      #"Gene_DDD_2022_Jade_Falling", "Gene_DDD_2022_Jade_Rising", "Gene_DDD_2022_Jade_Non_transitional",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_35_NA", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_35_NA", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Fulgent_NA", "DDD_old_ID_Satt_ClinGen_PanelApp_Fulgent_NA",
                      #"Gene_DDD_2022_Jade_NA", "Gene_DDD_2022_Jade_35_NA",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp", "DDD_old_ID_Satt_ClinGen_PanelApp",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Falling", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Rising", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_Non_transitional", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_NA",
                      #"DDD_old_ID_Satt_ClinGen_PanelApp_Falling", "DDD_old_ID_Satt_ClinGen_PanelApp_Rising", "DDD_old_ID_Satt_ClinGen_PanelApp_Non_transitional", "DDD_old_ID_Satt_ClinGen_PanelApp_NA",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_35", "DDD_old_ID_Satt_ClinGen_PanelApp_35",
                      #"DDD_old_ASD_ID_Satt_ClinGen_PanelApp_35_Falling", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_35_Rising", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_35_Non_transitional", "DDD_old_ASD_ID_Satt_ClinGen_PanelApp_35_NA",
                      #"DDD_old_ID_Satt_ClinGen_PanelApp_35_Falling", "DDD_old_ID_Satt_ClinGen_PanelApp_35_Rising", "DDD_old_ID_Satt_ClinGen_PanelApp_35_Non_transitional", "DDD_old_ID_Satt_ClinGen_PanelApp_35_NA",
                      #"DDD_ID_Satt","DDD_ID_Satt_35","DDD_ASD_ID_Satt","DDD_ASD_ID_Satt_35","DDD_ClinGen_ID_Satt", "DDD_ClinGen_ID_Satt_35","DDD_ClinGen_ASD_ID_Satt","DDD_ClinGen_ASD_ID_Satt_35"
)

#for (list_DDD in list_DDD_complete) {
list_DDD = list_DDD_complete[1]
print(list_DDD)

#  for (file in 1:length(list_file)) {
file = 1
file_name = list_file[file]

all_data_info_pheno_scores_SNP_Artefact = fread(
  paste0("/home/guillaf/projects/rrg-jacquese/guillaf/0_scratch_2022/",file_name,".tsv"),
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


summary(as.factor(all_data_info_pheno_scores_SNP_Artefact$Cohort))

all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB = ifelse( is.na(all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB) , -9 , all_data_info_pheno_scores_SNP_Artefact$Age_ZScore_IQ_adj_age_sex_PC1_10_with_online_UKBB)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

dim(all_data_info_pheno_scores_SNP_Artefact)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

gene.score.info = fread(
  paste0("/home/",user,"/projects/rrg-jacquese/guillaf/0_scratch_2022/Test_SNP_cover_all_Cohorts/geneuniquefiltered_overlapartifactFreq_UKBBMEGASPARKMSSNG_DEL_0_DUP_0_20230419.tsv"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = "."
)

gene.score.info = subset(subset(gene.score.info, proportion_gene_overlap == 1), select = c("CHR", "START", "STOP", "TYPE", "gene_id", "gene", "oe_lof_upper", "proportion_gene_overlap") )

ALL_IID_CNV_filtre = fread(
  paste0("/home/",user,"/projects/rrg-jacquese/guillaf/0_scratch_2022/Test_SNP_cover_all_Cohorts/CNVfiltered_overlapartifact_Freq1_noreciprocitie_50_UKBBMEGASPARKMSSNG_Score_30_SNP_10_20230419.tsv"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = "."
)

###################################################################################################
###################################################################################################

# List_IID_CNV_gene_chr2_0 = merge(ALL_IID_CNV_filtre, subset(gene.score.info, CHR == "chr2" & TYPE == "DUP" & START <= 118732006 & STOP >= 104142057), by = c("CHR", "START", "STOP", "TYPE"))
# 
# List_IID_gpop = subset(subset(List_pheno_IID_CNV_gene, Cohort != "SSC" & Cohort != "SPARK" & Cohort != "MSSNG"), select = c("individual"))
#   
# List_IID_CNV_gene_chr2 = merge(List_IID_CNV_gene_chr2_0, List_IID_gpop, by = c("individual"))
# 
# List_IID_CNV_gene_chr2_gene_interect = unique(subset(subset(List_IID_CNV_gene_chr2, gene_id == "ENSG00000135960" | gene_id == "ENSG00000172985" | gene_id == "ENSG00000186522" | gene_id == "ENSG00000198142"),
#                                               select = c("CHR", "START", "STOP", "individual", "TYPE")))
# 
# List_IID_CNV_gene_chr2_gene_interect$Important = "YES"
# 
# List_IID_CNV_gene_chr2_gene_all = unique(subset(List_IID_CNV_gene_chr2, select = c("CHR", "START", "STOP", "individual", "TYPE")))
# 
# List_IID_CNV_gene_chr2_gene_all_OK = merge(List_IID_CNV_gene_chr2_gene_all, List_IID_CNV_gene_chr2_gene_interect, by = c("CHR", "START", "STOP", "individual", "TYPE"), all = T)
# 
# summary(as.factor(List_IID_CNV_gene_chr2_gene_all_OK$Important))
# 
# List_IID_CNV_gene_chr2_gene_nointerect = unique(subset(subset(List_IID_CNV_gene_chr2_gene_all_OK, is.na(Important)), select = c("CHR", "START", "STOP", "individual", "TYPE")))
# 
# List_IID_CNV_gene_chr2_gene_interect = unique(subset(subset(List_IID_CNV_gene_chr2, gene_id == "ENSG00000135960" | gene_id == "ENSG00000172985" | gene_id == "ENSG00000186522" | gene_id == "ENSG00000198142"),
#                                                      select = c("CHR", "START", "STOP", "individual", "TYPE")))
# 
# fwrite(List_IID_CNV_gene_chr2_gene_interect,paste0("/home/guillaf/projects/rrg-jacquese/guillaf/List_IID_CNV_gene_chr2_gene_interect.tsv"), quote=F, sep="\t", row.names=F, col.names=F)
# 
# fwrite(List_IID_CNV_gene_chr2_gene_nointerect, paste0("/home/guillaf/projects/rrg-jacquese/guillaf/List_IID_CNV_gene_chr2_gene_nointerect.tsv"), quote=F, sep="\t", row.names=F, col.names=F)
# 
# length(unique(List_IID_CNV_gene_chr2$gene_id))
# 
# List_gene_observation_chr2 = unique(List_IID_CNV_gene_chr2$gene_id)
# 
# ###################################################################################################
# ###################################################################################################
# 
# List_IID_CNV_gene = merge(ALL_IID_CNV_filtre, subset(gene.score.info, gene_id == "ENSG00000170412"), by = c("CHR", "START", "STOP", "TYPE")) # chr17
# List_IID_CNV_gene = merge(ALL_IID_CNV_filtre, subset(gene.score.info, gene_id == "ENSG00000198142"), by = c("CHR", "START", "STOP", "TYPE")) # chr2
# #List_IID_CNV_gene = merge(ALL_IID_CNV_filtre, subset(gene.score.info, gene_id == "ENSG00000174943"), by = c("CHR", "START", "STOP", "TYPE")) # chr16
# List_IID_CNV_gene = merge(ALL_IID_CNV_filtre, subset(gene.score.info, gene_id == "ENSG00000170113"), by = c("CHR", "START", "STOP", "TYPE")) # chr15
# 
# ###################################################################################################

Gene_CNV_GWAS = fread(
  paste0("/home/",user,"/projects/rrg-jacquese/guillaf/DEL_DUP_Intell_article_all_Cohorts/Gene_CNV_GWAS.txt"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = "."
)

###################################################################################################
###################################################################################################

Infor_CNV_GWAS_Gene = c()
List_Cohort = c("All_wo_ASD", "CaG", "G-Scot", "Imagen", "LBC", "MSSNG", "SPARK", "SSC", "SYS", "Z_Gfact_5_Donald_et_al_2016", "Z_Gfact_5_online", "Z_score_FI", "Z_score_FI_online")

for (Gene_line in 1:dim(Gene_CNV_GWAS)[1]) {
  
  #Gene_line = 1
  
  TYPE_gene = (Gene_CNV_GWAS$TYPE)[Gene_line]
  Gene_info = (Gene_CNV_GWAS$unique_genes)[Gene_line]
  
  print(paste0(Gene_line/dim(Gene_CNV_GWAS)[1]))
  
  List_IID_CNV_gene = merge(ALL_IID_CNV_filtre, subset(gene.score.info, gene_id == Gene_info & TYPE == TYPE_gene), by = c("CHR", "START", "STOP", "TYPE")) # chr2
  
  List_pheno_IID_CNV_gene = merge(all_data_info_pheno_scores_SNP_Artefact, List_IID_CNV_gene, by = c("individual"), all.x = T)
  
  List_pheno_IID_CNV_gene$gene_id_carriers = ifelse( !(is.na(List_pheno_IID_CNV_gene$gene_id)), 1, 0)
  
  ###################################################################################################
  
  for (Cohort_name in List_Cohort) {
    
    C1 = c()
    
    if (Cohort_name == "All_wo_ASD") {
      
      model_IQ_DelDup = lm(ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final ~ as.factor(gene_id_carriers) +
                             PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 ,
                           subset(List_pheno_IID_CNV_gene, Cohort_final != "MSSNG" & Cohort_final !="SPARK" & Cohort_final !="SSC"), na.action = na.omit)
      
      #print(summary(List_pheno_IID_CNV_gene_test))
      
      A1 = summary(model_IQ_DelDup)$coefficients[1:2,1:4]
      
      Sum_1 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final != "MSSNG" & Cohort_final !="SPARK" & Cohort_final !="SSC" & gene_id_carriers > 0))[1]
      Sum_2 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final != "MSSNG" & Cohort_final !="SPARK" & Cohort_final !="SSC" & gene_id_carriers == 0))[1]
      Sum_3 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final != "MSSNG" & Cohort_final !="SPARK" & Cohort_final !="SSC"))[1]
      
      Sum_ABC1 = data.frame(
        rbind(
          0,
          Sum_1
        )
      )
      
      Sum_ABC2 = data.frame(
        rbind(
          0,
          Sum_2
        )
      )
      
      Sum_ABC3 = data.frame(
        rbind(
          0,
          Sum_3
        )
      )
      
      C1 = cbind(A1, Sum_ABC1, Sum_ABC2, Sum_ABC3, TYPE_gene,Gene_info,Cohort_name)
      
      colnames(C1) = c("Value",
                       "Std.Error",
                       "t.value",
                       "p.value",
                       "rbind.Sum.carrier",
                       "rbind.Sum.nocarrier",
                       "All_ind",
                       "TYPE_CNV","Gene_name","Cohorts")
      
      Infor_CNV_GWAS_Gene = rbind(Infor_CNV_GWAS_Gene, C1)
      
    }
    
    
    if (dim(subset(List_pheno_IID_CNV_gene, Cohort_final == Cohort_name & gene_id_carriers == 1))[1] > 0) {
      
      model_IQ_DelDup = lm(ZScore_IQ_adj_test2_age_sex_with_online_UKBB_Final ~ as.factor(gene_id_carriers) +
                             PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 ,
                           subset(List_pheno_IID_CNV_gene, Cohort_final == Cohort_name), na.action = na.omit)
      
      #print(summary(List_pheno_IID_CNV_gene_test))
      
      A1 = summary(model_IQ_DelDup)$coefficients[1:2,1:4]
      
      Sum_1 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final == Cohort_name & gene_id_carriers > 0))[1]
      Sum_2 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final == Cohort_name & gene_id_carriers == 0))[1]
      Sum_3 = dim(subset(List_pheno_IID_CNV_gene, Cohort_final == Cohort_name))[1]
      
      Sum_ABC1 = data.frame(
        rbind(
          0,
          Sum_1
        )
      )
      
      Sum_ABC2 = data.frame(
        rbind(
          0,
          Sum_2
        )
      )
      
      Sum_ABC3 = data.frame(
        rbind(
          0,
          Sum_3
        )
      )
      
      C1 = cbind(A1, Sum_ABC1, Sum_ABC2, Sum_ABC3, TYPE_gene,Gene_info,Cohort_name)
      
      colnames(C1) = c("Value",
                       "Std.Error",
                       "t.value",
                       "p.value",
                       "rbind.Sum.carrier",
                       "rbind.Sum.nocarrier",
                       "All_ind",
                       "TYPE_CNV","Gene_name","Cohorts")
      
      Infor_CNV_GWAS_Gene = rbind(Infor_CNV_GWAS_Gene, C1)
      
    }
  }
}

Infor_CNV_GWAS_Gene_Ok = subset(Infor_CNV_GWAS_Gene,All_ind !=0)

fwrite(Infor_CNV_GWAS_Gene_Ok,paste0("Meta_CNV_GWAS_Gene.tsv"), quote=F, sep="\t", row.names=F, col.names=T)
