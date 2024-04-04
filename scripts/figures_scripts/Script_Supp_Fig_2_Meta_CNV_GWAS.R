#################################################################

library(grid)
library(metafor)
library(meta)
library(stringr)
library(tidyr)

###################################################################################################

setwd(paste0("/Users/guillaume/Desktop/File_Script/Supp_Fig_2/"))

Infor_CNV_GWAS_Gene_Ok_del_dup = read.csv(
  "Meta_CNV_GWAS_Gene.tsv",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = ".",
  fileEncoding = "UTF-8-BOM"
)

Gene_CNV_GWAS = fread(
  paste0("Gene_CNV_GWAS.txt"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "#NULL!", ".", ""),
  dec = "."
)

Gene_CNV_GWAS$TYPE_Gene_name = paste(Gene_CNV_GWAS$TYPE, Gene_CNV_GWAS$unique_genes, sep = "_")

###################################################################################################

Infor_CNV_GWAS_Gene_Ok_del_dup$TYPE_Gene_name = paste(Infor_CNV_GWAS_Gene_Ok_del_dup$TYPE_CNV, Infor_CNV_GWAS_Gene_Ok_del_dup$Gene_name, sep = "_")

Infor_CNV_GWAS_Gene_Ok_del_dup_final = merge(Infor_CNV_GWAS_Gene_Ok_del_dup, Gene_CNV_GWAS, by.x = c("Gene_name","TYPE_CNV","TYPE_Gene_name"), by.y = c("unique_genes","TYPE","TYPE_Gene_name"))

all_info_meta_DelDup_wind_modif_final = c()

Infor_CNV_GWAS_Gene_Ok_del_dup_final$Estimate = Infor_CNV_GWAS_Gene_Ok_del_dup_final$Value

model_IQ_DelDup = subset(Infor_CNV_GWAS_Gene_Ok_del_dup_final, Cohorts != "All_wo_ASD" & !(Cohorts == "MSSNG" | Cohorts == "SSC" | Cohorts == "SPARK"))

model_IQ_DelDup$Estimate = model_IQ_DelDup$Value

names_list = c("DEL_ENSG00000013364", "DEL_ENSG00000058600", "DEL_ENSG00000072864", "DEL_ENSG00000079616", "DEL_ENSG00000085721", "DEL_ENSG00000090238", "DEL_ENSG00000102879", "DEL_ENSG00000102882", "DEL_ENSG00000102886", "DEL_ENSG00000103222", "DEL_ENSG00000103319", "DEL_ENSG00000103485", "DEL_ENSG00000103495", "DEL_ENSG00000103502", "DEL_ENSG00000116128", "DEL_ENSG00000121634", "DEL_ENSG00000131778", "DEL_ENSG00000131781", "DEL_ENSG00000131791", "DEL_ENSG00000133392", "DEL_ENSG00000133393", "DEL_ENSG00000140157", "DEL_ENSG00000140740", "DEL_ENSG00000140743", "DEL_ENSG00000140749", "DEL_ENSG00000144152", "DEL_ENSG00000144161", "DEL_ENSG00000149922", "DEL_ENSG00000149923", "DEL_ENSG00000149925", "DEL_ENSG00000149926", "DEL_ENSG00000149927", "DEL_ENSG00000149929", "DEL_ENSG00000149930", "DEL_ENSG00000149932", "DEL_ENSG00000153093", "DEL_ENSG00000153094", "DEL_ENSG00000153107", "DEL_ENSG00000153208", "DEL_ENSG00000153214", "DEL_ENSG00000153666", "DEL_ENSG00000155714", "DEL_ENSG00000156968", "DEL_ENSG00000157045", "DEL_ENSG00000162836", "DEL_ENSG00000166780", "DEL_ENSG00000166783", "DEL_ENSG00000167194", "DEL_ENSG00000167371", "DEL_ENSG00000169592", "DEL_ENSG00000170113", "DEL_ENSG00000174938", "DEL_ENSG00000174939", "DEL_ENSG00000174943", "DEL_ENSG00000174992", "DEL_ENSG00000175267", "DEL_ENSG00000182974", "DEL_ENSG00000183706", "DEL_ENSG00000183793", "DEL_ENSG00000183921", "DEL_ENSG00000185716", "DEL_ENSG00000185905", "DEL_ENSG00000188092", "DEL_ENSG00000197006", "DEL_ENSG00000197414", "DEL_ENSG00000197471", "DEL_ENSG00000203836", "DEL_ENSG00000265107", "DEL_ENSG00000273749",
               "DUP_ENSG00000013364", "DUP_ENSG00000040608", "DUP_ENSG00000063515", "DUP_ENSG00000070010", "DUP_ENSG00000070371", "DUP_ENSG00000070413", "DUP_ENSG00000072864", "DUP_ENSG00000079616", "DUP_ENSG00000085721", "DUP_ENSG00000090238", "DUP_ENSG00000091262", "DUP_ENSG00000093009", "DUP_ENSG00000093010", "DUP_ENSG00000099889", "DUP_ENSG00000099899", "DUP_ENSG00000099901", "DUP_ENSG00000099904", "DUP_ENSG00000099910", "DUP_ENSG00000099917", "DUP_ENSG00000099937", "DUP_ENSG00000099940", "DUP_ENSG00000099942", "DUP_ENSG00000099949", "DUP_ENSG00000099957", "DUP_ENSG00000099960", "DUP_ENSG00000100033", "DUP_ENSG00000100056", "DUP_ENSG00000100075", "DUP_ENSG00000100084", "DUP_ENSG00000102879", "DUP_ENSG00000102882", "DUP_ENSG00000102886", "DUP_ENSG00000103222", "DUP_ENSG00000103485", "DUP_ENSG00000103495", "DUP_ENSG00000103502", "DUP_ENSG00000117262", "DUP_ENSG00000117281", "DUP_ENSG00000121851", "DUP_ENSG00000128185", "DUP_ENSG00000128191", "DUP_ENSG00000131779", "DUP_ENSG00000131788", "DUP_ENSG00000133392", "DUP_ENSG00000133393", "DUP_ENSG00000134160", "DUP_ENSG00000135960", "DUP_ENSG00000143127", "DUP_ENSG00000149922", "DUP_ENSG00000149923", "DUP_ENSG00000149925", "DUP_ENSG00000149926", "DUP_ENSG00000149927", "DUP_ENSG00000149929", "DUP_ENSG00000149930", "DUP_ENSG00000149932", "DUP_ENSG00000156968", "DUP_ENSG00000157045", "DUP_ENSG00000166780", "DUP_ENSG00000166783", "DUP_ENSG00000166912", "DUP_ENSG00000167194", "DUP_ENSG00000167371", "DUP_ENSG00000168509", "DUP_ENSG00000169592", "DUP_ENSG00000169918", "DUP_ENSG00000169926", "DUP_ENSG00000172985", "DUP_ENSG00000174827",
               "DUP_ENSG00000174938", "DUP_ENSG00000174939", "DUP_ENSG00000174943", "DUP_ENSG00000174992", "DUP_ENSG00000179889", "DUP_ENSG00000183426", "DUP_ENSG00000183597", "DUP_ENSG00000183628", "DUP_ENSG00000183773", "DUP_ENSG00000183793", "DUP_ENSG00000184058", "DUP_ENSG00000184113", "DUP_ENSG00000184436", "DUP_ENSG00000184470", "DUP_ENSG00000184702", "DUP_ENSG00000185252", "DUP_ENSG00000185608", "DUP_ENSG00000185838", "DUP_ENSG00000185905", "DUP_ENSG00000186141", "DUP_ENSG00000186364", "DUP_ENSG00000186522", "DUP_ENSG00000187905", "DUP_ENSG00000188280", "DUP_ENSG00000197471", "DUP_ENSG00000198142", "DUP_ENSG00000198483", "DUP_ENSG00000198690", "DUP_ENSG00000203618", "DUP_ENSG00000206203", "DUP_ENSG00000234409", "DUP_ENSG00000241973", "DUP_ENSG00000242259", "DUP_ENSG00000244486", "DUP_ENSG00000265491", "DUP_ENSG00000265972", "DUP_ENSG00000271601", "DUP_ENSG00000272031", "DUP_ENSG00000273611", "DUP_ENSG00000273706", "DUP_ENSG00000275066", "DUP_ENSG00000275410", "DUP_ENSG00000275700", "DUP_ENSG00000275793", "DUP_ENSG00000276023", "DUP_ENSG00000276234", "DUP_ENSG00000277161", "DUP_ENSG00000278053", "DUP_ENSG00000278259", "DUP_ENSG00000278311", "DUP_ENSG00000278505", "DUP_ENSG00000278535", "DUP_ENSG00000278540", "DUP_ENSG00000278619")

### with ASD

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

  File_Analyses = subset(model_IQ_DelDup, TYPE_Gene_name == names_list[i] & !(is.na(Estimate)))
  
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

}

all_info_meta_DelDup_wind_modif = cbind(names_list,
                                        list_info_Hetero_DelDup, list_info_Hetero_I_DelDup, list_info_Hetero_t2_DelDup,
                                        list_info_Hetero_est_fix_DelDup,list_info_Hetero_SD_fix_DelDup,list_info_Hetero_pvalue_fix_DelDup,
                                        list_info_Hetero_est_random_DelDup,list_info_Hetero_SD_random_DelDup,list_info_Hetero_pvalue_random_DelDup,"without_ASD")

colnames(all_info_meta_DelDup_wind_modif) = c("Names", 
                                              "Heterog_pvalue","I", "T2",
                                              "Estimate_fix","Std.Error_fix","p.value_fix",
                                              "Estimate_random","Std.Error_random","p.value_random","Cohorts")

all_info_meta_DelDup_wind_modif = data.frame(all_info_meta_DelDup_wind_modif)

all_info_meta_DelDup_wind_modif$Heterog_pvalue_adj = p.adjust(all_info_meta_DelDup_wind_modif$Heterog_pvalue, method ="fdr")

all_info_meta_DelDup_wind_modif_final = rbind(all_info_meta_DelDup_wind_modif_final, all_info_meta_DelDup_wind_modif)

all_info_meta_DelDup_wind_modif_final$TYPE = substr(all_info_meta_DelDup_wind_modif_final$Names, 1, 3)

all_info_meta_DelDup_wind_modif_final_Ok = merge(Gene_CNV_GWAS, all_info_meta_DelDup_wind_modif_final, by.x = c("TYPE_Gene_name"), by.y = c("Names"))

table(all_info_meta_DelDup_wind_modif_final_Ok$TYPE.x, all_info_meta_DelDup_wind_modif_final_Ok$TYPE.y)

all_info_meta_DelDup_wind_modif_final_Ok$Heterog_pvalue = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$Heterog_pvalue)


all_info_meta_DelDup_wind_modif_final_Ok$Estimate_fix = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$Estimate_fix)
all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_fix = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_fix)
all_info_meta_DelDup_wind_modif_final_Ok$p.value_fix = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$p.value_fix)

all_info_meta_DelDup_wind_modif_final_Ok$Estimate_random = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$Estimate_random)
all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_random = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_random)
all_info_meta_DelDup_wind_modif_final_Ok$p.value_random = as.numeric(all_info_meta_DelDup_wind_modif_final_Ok$p.value_random)

all_info_meta_DelDup_wind_modif_final_Ok$Estimate_final = ifelse(all_info_meta_DelDup_wind_modif_final_Ok$Heterog_pvalue > 0.1, all_info_meta_DelDup_wind_modif_final_Ok$Estimate_fix, all_info_meta_DelDup_wind_modif_final_Ok$Estimate_random)
all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_final = ifelse(all_info_meta_DelDup_wind_modif_final_Ok$Heterog_pvalue > 0.1, all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_fix, all_info_meta_DelDup_wind_modif_final_Ok$Std.Error_random)
all_info_meta_DelDup_wind_modif_final_Ok$p.value_final = ifelse(all_info_meta_DelDup_wind_modif_final_Ok$Heterog_pvalue > 0.1, all_info_meta_DelDup_wind_modif_final_Ok$p.value_fix, all_info_meta_DelDup_wind_modif_final_Ok$p.value_random)


###

Info_CNV_meta = subset(all_info_meta_DelDup_wind_modif_final_Ok, select = c("CNV", "TYPE.x", "TYPE_Gene_name", "Estimate_final","Std.Error_final", "p.value_final"))
colnames(Info_CNV_meta) = c("CNV", "TYPE_CNV", "TYPE_Gene_name", "Estimate","Std.Error", "p.value")
Info_CNV_meta$Cohorts = c("Meta_All_wo_ASD")

Info_CNV_meta_pool = rbind(Info_CNV_meta, subset(subset(Infor_CNV_GWAS_Gene_Ok_del_dup_final,Cohorts == "All_wo_ASD") , select = c("CNV", "TYPE_CNV", "TYPE_Gene_name", "Estimate","Std.Error", "p.value", "Cohorts")))

####

List_max_carrier_CNV = c()

for (i in c("DEL-1_1q21.1_distal", "DEL-2_2q13", "DEL-3_15q11.2", "DEL-4_16p13.11", "DEL-5_16p12.1", "DEL-6_16p11.2_proximal",
            "DUP-1_1q21.1_TAR", "DUP-2", "DUP-3_15q13.3_BP4-BP5", "DUP-4_16p13.11", "DUP-5_16p11.2_proximal", "DUP-6_17q12", "DUP-7_22q11.2_proximal")){
  
  List_i = c()
  name_i = c()
  nb_name_i = c()
  
  List_i = max(subset(Infor_CNV_GWAS_Gene_Ok_del_dup_final, CNV == i & Cohorts == "All_wo_ASD")$rbind.Sum.carrier)
  name_i = c(unique(subset(Infor_CNV_GWAS_Gene_Ok_del_dup_final, CNV == i & rbind.Sum.carrier == List_i & Cohorts == "All_wo_ASD")$TYPE_Gene_name))
  nb_name_i = length(subset(Infor_CNV_GWAS_Gene_Ok_del_dup_final, CNV == i & rbind.Sum.carrier == List_i & Cohorts == "All_wo_ASD")$TYPE_Gene_name)
  List_max_carrier_CNV = rbind(List_max_carrier_CNV, data.frame(i,List_i,name_i,nb_name_i ))
  
}

data.frame(List_max_carrier_CNV)
unique(subset(List_max_carrier_CNV, select = c("i","List_i")))

filtered_data <- List_max_carrier_CNV %>%
  group_by(i) %>%
  slice(1)

data.frame(filtered_data)

Info_CNV_meta_pool_filtre = merge(data.frame(filtered_data), Info_CNV_meta_pool, by.x = c("name_i"), by.y = c("TYPE_Gene_name"), all.x = T)

pd <- position_dodge(0.9)

library(forcats)

# Convertir 'category' en facteur et inverser l'ordre
Info_CNV_meta_pool_filtre$CNV <- factor(Info_CNV_meta_pool_filtre$CNV)
Info_CNV_meta_pool_filtre$CNV <- fct_rev(Info_CNV_meta_pool_filtre$CNV)

test = ggplot(Info_CNV_meta_pool_filtre, aes(y=Estimate, x=CNV, colour=Cohorts, grp = CNV)) + 
  
  geom_errorbar(aes(ymin=Estimate-(1.96*Std.Error), ymax=Estimate+(1.96*Std.Error)), position=pd) +
  geom_hline(yintercept=0) +
  #scale_color_manual(values=c("purple", "orange")) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(position=pd, size=ifelse(Info_CNV_meta_pool_filtre$p.value < 0.05, 5, 3),
             shape=23, fill="white") +
  theme(legend.position='bottom')

test + coord_flip()
