import numpy as np
import pandas as pd
import time
import statsmodels.formula.api as sm


def compute_loeuf_cat_specific_genes(
        data: pd.DataFrame,
        gene_column_name: str,
        gene_list: list,
        ddd_genes: list) -> pd.DataFrame:
    """
    Create LOEUF categories for specific gene-set in a given dataset.

    Args:
        data (pandas.DataFrame): The input dataset.
        gene_column_name (str): The name of the column containing gene names.
        gene_list (list): A list of genes to consider.
        ddd_genes (list): A list of DDD (Deciphering Developmental Disorders) genes.

    Returns:
        pandas.DataFrame: The input dataset with additional columns for LOEUF categories and inverse of oe_lof_upper.

    The function assigns LOEUF categories to genes based on specific conditions and creates additional columns in the dataset
    to indicate the category for each gene. It also calculates the inverse of oe_lof_upper and assigns the values to a new column.

    The LOEUF categories are defined as follows:
    - '1_gene_list': Genes in gene_list with oe_lof_upper < 0.35
    - '2_gene_list': Genes in gene_list with 0.35 <= oe_lof_upper < 1
    - '3_gene_list': Genes in gene_list with oe_lof_upper >= 1
    - '1_outside': Genes outside gene_list with oe_lof_upper < 1 and not in ddd_genes
    - '2_outside': Genes outside gene_list with oe_lof_upper >= 1 and not in ddd_genes
    - 'ddd': Genes in ddd_genes but not in gene_list

    Note: The function modifies the input dataset in-place by adding new columns and dropping a temporary column.
    """
    # List of gene categories used as covariates in regression models
    cat_list = [
        '1_gene_list',
        '2_gene_list',
        '3_gene_list',
        '1_outside',
        '2_outside',
        'ddd']
    # Assigning genes to the given categories based on following conditions
    data["tmp_cat"] = np.select([(data[gene_column_name].isin(gene_list)) & (data.oe_lof_upper < 0.35),
                                 (data[gene_column_name].isin(gene_list)) & (data.oe_lof_upper < 1) & (
        data.oe_lof_upper >= 0.35),
        (data[gene_column_name].isin(gene_list)) & (data.oe_lof_upper >= 1),
        (~data[gene_column_name].isin(gene_list)) & (data.oe_lof_upper < 1) & (
                                     ~data[gene_column_name].isin(ddd_genes)),
        (~data[gene_column_name].isin(gene_list)) & (data.oe_lof_upper >= 1) & (
                                     ~data[gene_column_name].isin(ddd_genes)),
        (data[gene_column_name].isin(ddd_genes)) & (~data[gene_column_name].isin(gene_list))
    ], cat_list)
    # For each cateory will create the associated column with binary values
    # indicating if gene within category
    for cat in cat_list:
        # Create column fill with 0
        data.loc[:, "LOEUF_cat_{}".format(cat)] = 0
        # Get indexes of genes matching the category
        indexes = data[(data.tmp_cat == cat)].index
        # Replace 0 by 1 for selected indexes
        data.loc[indexes, "LOEUF_cat_{}".format(cat)] = 1
    # Create 1 over LOEUF score
    data["loeuf_inv_full"] = 1 / data.oe_lof_upper
    data.drop(columns={'tmp_cat'})
    return data


def aggregate_ind_cnv_genes(cnv_indd, gene_score_info, individual_info):
    SUM_IND_type_SCORE = {}
    for cnv_type in ["DEL", 'DUP']:
        # Seperate data based on CNV type to compute DEL and DUP separetely
        gene_score_info_type = gene_score_info[gene_score_info["TYPE"] == cnv_type]
        cnv_indd_type = cnv_indd[cnv_indd["TYPE"] == cnv_type]
        # Create gene file with gene coordinates as index using chromosome start
        # and stop info
        df2 = pd.DataFrame(gene_score_info_type.loc[:, ["CHR", "START", "STOP"]])
        df2["position"] = df2.apply(
            lambda row: "-".join([str((val)) for val in row]), axis=1)
        print(df2["position"].nunique())
        gene_score_info_type.index = df2["position"]
    
        # Create CNV file with gene coordinates as index using chromosome start
        # and stop info
        df3 = pd.DataFrame(cnv_indd_type.loc[:, ["CHR", "START", "STOP"]])
        df3["position"] = df3.apply(
            lambda row: "-".join([str((val)) for val in row]), axis=1)
        cnv_indd_type.index = cnv_indd_type["individual"]
    
        
        cnv_indd_type_INDVID = df3.loc[:, ["position"]]
        
        cnv_indd_type_INDVID0 = np.column_stack(
            (cnv_indd_type_INDVID, cnv_indd_type["individual"])
        )
        
        cnv_indd_type_IND_POS = pd.DataFrame(
            cnv_indd_type_INDVID0, columns=["position", "individual"]
        )
    
        # Merge two datasets
        POS_IND_score_type = pd.merge(
            cnv_indd_type_IND_POS, gene_score_info_type, on="position", how="outer"
        )
        print(POS_IND_score_type.shape)
        # sum over individuals
        SUM_IND_type_SCORE[cnv_type] = POS_IND_score_type 
        
    
    # FINAL OUTPUT-Merge DUP and DEL

    POS_IND_score_merge = pd.concat([
            SUM_IND_type_SCORE['DUP'],
            SUM_IND_type_SCORE['DEL']], axis = 0)
    # add
    
    return POS_IND_score_merge

def computeUniqueGenes(
        cnv_indd: pd.DataFrame,
        gene_score_info_overlap_split: pd.DataFrame,
        cnv_type: str,
        cat: str) -> pd.DataFrame:
    """
    Computes unique genes lists by individual based on the given copy number variant (CNV) information and gene score information.

    Args:
        cnv_indd (pandas.DataFrame): DataFrame containing CNV information by individual.
        gene_score_info_overlap_split (pandas.DataFrame): DataFrame containing gene score information.
        cnv_type (str): Type of CNV to consider.
        cat (str): Category of gene score information to consider.

    Returns:
        pandas.DataFrame: DataFrame with individual and corresponding unique genes.
    """
    start_time = time.time()
    # Filter cnv_indd DataFrame based on cnv_type
    cnv_indd_tmp = cnv_indd[cnv_indd.TYPE == cnv_type].copy()
    # Create CNV file with gene coordinates as index using chromosome start
    # and stop info
    cnv_indd_tmp["POS"] = cnv_indd_tmp[['CHR', 'START', 'STOP']].apply(
        lambda row: "--".join([str((val)) for val in row]), axis=1)
    # Extract only genes matching the LOEUF category given by the cat
    # argument. If argument is all then whole dataset is copied
    if cat != 'all':
        gene_score_info_overlap_split_cat = gene_score_info_overlap_split[gene_score_info_overlap_split.tmp_cat == cat].copy()
    else:
        gene_score_info_overlap_split_cat = gene_score_info_overlap_split.copy()
    # Create gene file with gene coordinates as index using chromosome start
    # and stop info
    gene_score_info_overlap_split_cat["POS"] = gene_score_info_overlap_split_cat[[
        'CHR', 'START', 'STOP']].apply(lambda row: "--".join([str((val)) for val in row]), axis=1)
    # Filter genes_clean DataFrame by removing rows where oe_lof_upper is NaN
    genes_clean = gene_score_info_overlap_split_cat[~gene_score_info_overlap_split_cat.oe_lof_upper.isna()]
    # Merge cnv_indd_tmp and genes_clean DataFrames based on genes defined by
    # their POS value
    tmp2 = pd.merge(cnv_indd_tmp, genes_clean, on="POS")
    # Create a dataframe containing for each individual the list of genes
    # carried joined by a comma
    ind_genes = tmp2.groupby(['individual'],
                             as_index=False).agg({'gene_id': ','.join})
    print("-- %s seconds --" % (time.time() - start_time), flush=True)

    return ind_genes


def compute_unique_gene_list(list_uniques, gene_list, carriers, cnv_type):
    list_uniques_type = {}
    for (key, value) in list_uniques.items():
        if key.startswith(cnv_type):
            list_uniques_type[key] = value

    uniques = {}
    if len(carriers) == 2:
        loeuf_list = ['all']
        uniques['{}_all'.format(
            cnv_type)] = list_uniques_type['{}_all'.format(cnv_type)]

    else:
        loeuf_list = ['intol', 'inter', 'tol']
        for (key, value) in list_uniques_type.items():
            if not key.endswith('all'):
                uniques[key] = value

    unique_genes = [[]]
    for i, c in enumerate(carriers[1:]):
        tmp_genes = uniques['{}_{}'.format(cnv_type, loeuf_list[i])][uniques['{}_{}'.format(
            cnv_type, loeuf_list[i])].individual.isin(c)].copy()
        all_genes = [genes.split(',') for genes in tmp_genes.gene_id]
        all_g = [g for sublist in all_genes for g in sublist]
        unique_genes.append([gl_g for gl_g in list(
            set(all_g)) if gl_g in gene_list.gene_id.tolist()])
    return unique_genes


def runModel(
        data: pd.DataFrame,
        model_name: str,
        cnv_type: str,
        pheno_score: str,
        list_scores: list,
        gene_list: list,
        list_uniques: list,
        other_covariates: str) -> pd.DataFrame:
    # Get number of covariates added to the model.
    len_covariates = len(other_covariates.split(sep='+')) -1
    # Documentation for each model:
    #   - tmp_scores will be the list of scores (sum of genes by categories) used in the formula as covariate.
    #   - formula will aggregate the phenotype score name, the score covariates defined and the other optional covariates
    #   - cat will contains the category description, None for intercept values, g for score covariable of the gene_list and o for score covariate outside the list of interest
    #   - win will describe the position of each window, useful to distinguish covariates within a stritified model.
    #   - col_sum = sum of genes counted in each score covariate
    #   - carriers = number of individuals carrying at least one gene of the category.
    #   - unique_genes will be the list of unique genes identified in our dataset for each score covariate.

    if model_name == "3_3":
        tmp_scores = ['{}_{}'.format(cnv_type, val) for val in list_scores]
        formula = "{} ~ {}".format(pheno_score,
                                   ' + '.join(tmp_scores)) + other_covariates
        cat = [None, 'g', 'g', 'g', 'o', 'o', 'o']
        win = [0, 1, 2, 3, 1, 2, 0]
        col_sum = [0] + [data[val].sum() for val in tmp_scores]
        carriers = [[]] + [data.loc[data[val] > 0].individual.tolist()
                           for val in tmp_scores]
        unique_genes = compute_unique_gene_list(
            list_uniques, gene_list, carriers[0:4], cnv_type) + [[], [], []]

    elif model_name == "1_3":
        tmp_scores = ['{}_{}'.format(cnv_type, val) for val in list_scores if not (val.endswith("gene_list"))]
        data.eval('{0}_LOEUF_gene_list = {0}_LOEUF_cat_1_gene_list + {0}_LOEUF_cat_2_gene_list + {0}_LOEUF_cat_3_gene_list'.format(cnv_type),
            inplace=True)
        formula = "{} ~ {}_LOEUF_gene_list + {}".format(pheno_score, cnv_type, ' + '.join(tmp_scores)) + other_covariates
        cat = [None, 'g', 'o', 'o', 'o']
        win = [0, 1, 0, 1, 2]
        col_sum = [0, data['{0}_LOEUF_gene_list'.format(cnv_type)].sum()] + [data[val].sum() for val in tmp_scores]
        carriers = [[]] + [data[data['{0}_LOEUF_gene_list'.format(cnv_type)] > 0].individual.unique().tolist()] + [data[data[val] > 0].individual.unique().tolist() for val in tmp_scores]
        unique_genes = compute_unique_gene_list(
            list_uniques, gene_list, carriers[0:2], cnv_type) + [[], [], []]

    elif model_name == "1_1":
        data.eval(
            '{0}_LOEUF_gene_list = {0}_LOEUF_cat_1_gene_list + {0}_LOEUF_cat_2_gene_list + {0}_LOEUF_cat_3_gene_list'.format(cnv_type),
            inplace=True)
        data.eval(
            '{0}_LOEUF_outside = {0}_LOEUF_cat_1_outside + {0}_LOEUF_cat_2_outside + {0}_LOEUF_cat_ddd'.format(cnv_type),
            inplace=True)
        formula = "{0} ~ {1}_LOEUF_gene_list +  {1}_LOEUF_outside ".format(
            pheno_score, cnv_type) + other_covariates
        cat = [None, 'g', 'o']
        win = [0, 1, 1]
        col_sum = [0, data['{0}_LOEUF_gene_list'.format(cnv_type)].sum(
        ), data['{0}_LOEUF_outside'.format(cnv_type)].sum()]
        carriers = [[]] + [data.loc[data['{0}_LOEUF_gene_list'.format(cnv_type)] > 0].individual.unique(
        ).tolist(), data.loc[data['{0}_LOEUF_outside'.format(cnv_type)] > 0].individual.unique().tolist()]
        unique_genes = compute_unique_gene_list(
            list_uniques, gene_list, carriers[0:2], cnv_type) + [[]]

    elif model_name == '1_0':
        data.eval(
            '{0}_LOEUF_gene_list = {0}_LOEUF_cat_1_gene_list + {0}_LOEUF_cat_2_gene_list + {0}_LOEUF_cat_3_gene_list'.format(cnv_type),
            inplace=True)
        formula = "{0} ~ {1}_LOEUF_gene_list".format(
            pheno_score, cnv_type) + other_covariates
        cat = [None, 'g']
        win = [0, 1]
        col_sum = [0, data['{0}_LOEUF_gene_list'.format(cnv_type)].sum()]
        carriers = [
            []] + [data.loc[data['{0}_LOEUF_gene_list'.format(cnv_type)] > 0].individual.unique().tolist()]
        unique_genes = compute_unique_gene_list(
            list_uniques, gene_list, carriers[0:2], cnv_type)

    elif model_name == "3_0":
        tmp_scores = [
            '{}_{}'.format(
                cnv_type,
                val) for val in list_scores if (
                val.endswith("gene_list"))]
        formula = "{} ~ + {}".format(pheno_score,
                                     ' + '.join(tmp_scores)) + other_covariates
        cat = [None, 'g', 'g', 'g']
        win = [0, 1, 2, 3]
        col_sum = [0] + [data[val].sum() for val in tmp_scores]
        carriers = [[]] + [data[data[val] > 0].individual.unique().tolist()
                           for val in tmp_scores]
        unique_genes = compute_unique_gene_list(
            list_uniques, gene_list, carriers[0:4], cnv_type)

    else:
        print("Unknwon model")
        exit()

    # Run the linear regression following formula and based on data. Capture error and print exception
    try:
        reg = sm.gls(formula, data=data).fit()
    except Exception as e:
        print("Error: ", e)
        return pd.DataFrame()

    if len_covariates == 0:
        return pd.DataFrame({'Estimate': reg.params,
							 'SE': reg.bse,
							 'tvalues': reg.tvalues,
							 'pvalue': reg.pvalues,
							 'cat': cat,
							 'win': win,
							 'aic': reg.aic,
							 'bic': reg.bic,
							 'model_type': model_name,
							 'TYPE': cnv_type,
							 'N_total_genes': col_sum,
							 'total_carriers': [len(c) for c in carriers],
							#   'carriers': [','.join(c) for c in carriers],
							 'unique_genes': [','.join(ug) for ug in unique_genes],
							 'n_unique_genes': [len(ug) for ug in unique_genes]
							 }
							)
    else:
        return pd.DataFrame({'Estimate': reg.params[:-len_covariates],
							 'SE': reg.bse[:-len_covariates],
							 'tvalues': reg.tvalues[:-len_covariates],
							 'pvalue': reg.pvalues[:-len_covariates],
							 'cat': cat,
							 'win': win,
							 'aic': reg.aic,
							 'bic': reg.bic,
							 'model_type': model_name,
							 'TYPE': cnv_type,
							 'N_total_genes': col_sum,
							 'total_carriers': [len(c) for c in carriers],
							#   'carriers': [','.join(c) for c in carriers],
							 'unique_genes': [','.join(ug) for ug in unique_genes],
							 'n_unique_genes': [len(ug) for ug in unique_genes]
							 }
							)

def create_dict_uniques_genes_by_loeuf_cat(gene_score_info_overlap_split, cnv_indd):
    loeuf_list = ['intol', 'inter', 'tol']

    gene_score_info_overlap_split["tmp_cat"] = np.select([gene_score_info_overlap_split.oe_lof_upper<0.35,
              (gene_score_info_overlap_split.oe_lof_upper<1) & (gene_score_info_overlap_split.oe_lof_upper>=0.35),
               (gene_score_info_overlap_split.oe_lof_upper>=1)],
              loeuf_list)

    list_uniques = {}
    for cnv_type in ['DEL', 'DUP']:
        for cat in loeuf_list + ['all']:
            list_uniques['{}_{}'.format(cnv_type, cat)] = computeUniqueGenes(cnv_indd, gene_score_info_overlap_split, cnv_type, cat)
    return list_uniques


def compute_linear_models(data, cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates, dict_filters):
	tmp_summary = pd.DataFrame()
	model_summary_3_3 = runModel(data, "3_3", cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates)
	# model_summary_3_0 = runModel(data, "3_0", cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates)

	model_summary_1_3 = runModel(data, "1_3", cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates)
	model_summary_1_1 = runModel(data, "1_1", cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates)
	# model_summary_1_0 = runModel(data, "1_0", cnv_type, pheno_score, list_scores, gene_list, list_uniques, other_covariates)


	tmp_summary = pd.concat([model_summary_3_3, model_summary_1_3, model_summary_1_1])

	tmp_summary["gene_list_n"] = gene_list.shape[0]

	tmp_summary["gene_list_name"] = dict_filters['list_name']
	tmp_summary["ancestry"] = dict_filters['ancestry']
	tmp_summary["diagnosis"] = dict_filters['diag']
	tmp_summary["max_age"] = dict_filters['max_age']
	tmp_summary["max_score"] = dict_filters['max_s']
	tmp_summary.insert(0, "variable", tmp_summary.index)
	return tmp_summary