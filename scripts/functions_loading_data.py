import pandas as pd
import numpy as np
def prepare_Megadataset_dataframes(score_names):
	'''
	Prepare dataframes for the Megadataset analysis.

	Args:
		score_names (list): A list of score names.

	Returns:
		clean_individuals (DataFrame): Dataframe containing clean individual information.
		cnv_indd (DataFrame): Dataframe containing CNV individual information.
		gene_score_info_overlap_split (DataFrame): Dataframe containing gene score information with overlapping genes.

	'''
	# Different files used for analysis
	clean_data = "<Path to the clean list of individuals passing QC filters and having information of their ancestry>"
	ind_info_file = "<Path to the complete list of individuals with various phenotypic information (age, sex, conditions, cognitive ability etc...)>"
	cnv_ind_file = "<Path to the list of CNVs passing CNV QC filters and larger than 50kb>"
	gene_info_file = "<Path to the list of all genes overlapped by a CNV with its scores and the protpotion of overlap with the CNV>"

	# Open files
	cnv_indd = pd.read_csv(cnv_ind_file, sep='\t')
	individual_info = pd.read_csv(ind_info_file, sep='\t')     
	gene_score_info = pd.read_csv(gene_info_file, sep='\t', usecols=['CHR', 'START', 'STOP','TYPE', 'proportion_gene_overlap', 'pLI_2019', 'NB_Genes' , 'loeuf_inv', 'gene', 'gene_id'] + score_names)
	clean_individual_list = pd.read_csv(clean_data, sep="\t")
	gene_score_info_overlap_split = gene_score_info[gene_score_info["proportion_gene_overlap"] >= 1].copy()
	gene_score_info_overlap_split = gene_score_info_overlap_split.reset_index()
	clean_individuals = pd.merge(clean_individual_list.loc[:,['individual', 'Merge_final_ancestry']], individual_info, on="individual", how='inner')

	return clean_individuals, cnv_indd, gene_score_info_overlap_split

def get_DDD_genes():
	# List of DDD genes
	DDD_genes = pd.read_csv("<Path to the list of clingen DDD gees>", sep='\t')
	DDD_genes.rename(columns={"gene_gnomad":"gene"}, inplace=True)
	return DDD_genes