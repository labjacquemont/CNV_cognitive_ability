# Effects of gene dosage on cognitive ability: A function-based association study across brain and non-brain processes

### Authors:

Guillaume Huguet**; Thomas Renne**; Cécile Poulain; Alma Dubuc; Kuldeep Kumar; Sayeh Kazem; Worrawat Engchuan; Omar Shanta;  Elise Douard; Catherine Proulx; Martineau Jean-Louis; Zohra Saci; Josephine Mollon; Laura M Schultz; Emma E M Knowles; Simon R. Cox; David Porteous; Gail Davies; Paul Redmond; Sarah E. Harris; Gunter Schumann; Guillaume Dumas; Aurélie Labbe;  Zdenka Pausova; Tomas Paus; Stephen W Scherer; Jonathan Sebat; Laura Almasy; David C Glahn; Sébastien Jacquemont

** shared first authorship

### Corresponding authors:

**Guillaume Huguet**\
Sainte Justine University Hospital\
3175 chemin de la Côte-Sainte-Catherine,\
Montréal, QC H3T 1C5\
guillaumeaf.huguet@gmail.com
 
**Sébastien Jacquemont**\
Sainte Justine University Hospital\
3175 chemin de la Côte-Sainte-Catherine,\
Montréal, QC H3T 1C5\
sebastien.jacquemont@umontreal.ca

### Data availability
All general population data are available to other investigators online; Imagen (https://www.cataloguementalhealth.ac.uk), LBC (https://lothian-birth-cohorts.ed.ac.uk/), SYS (contact: tomas.paus@umontreal.ca), CaG (https://portal.canpath.ca/), Generation Scotland (https://www.ed.ac.uk/generation-scotland), and UK Biobank (ukbiobank.ac.uk). All ASD population data are available to other investigators online; SSC (https://www.sfari.org/), SPARK (https://www.sfari.org/), and MSSNG (https://research.mss.ng/). All derived measures used in this study are available upon request (S.J., sebastien.jacquemont@umontreal.ca). The rest of the CNV carriers’ data cannot be shared as participants did not provide consent. Source data are provided with this paper. Summary statistics and gene-sets studied are available through figshare: XXX

## Abstract
Genomic Copy Number Variants (CNVs) that increase risk for neurodevelopmental disorders are also associated with lower cognitive ability in general population cohorts. Studies have focussed on a small set of recurrent CNVs, but burden analyses suggested that the vast majority of CNVs affecting cognitive ability are too rare to reach variant-level association. As a result, the full range of gene-dosage-sensitive biological processes linked to cognitive ability remains unknown.

To investigate this issue, we identified all CNVs >50 kilobases in 258k individuals from 6 general population cohorts with assessments of general cognitive abilities. We performed a CNV-GWAS and functional burden analyses, which tested 6502 gene-sets defined by tissue and cell-type transcriptomics as well as gene ontology disrupted by all rare coding CNVs.

CNV-GWAS identified a novel duplication at 2q12.3 associated with higher performance in cognitive ability. Among the 864 gene-sets associated with cognitive ability, only 11% showed significant effects for both deletions and duplication. Accordingly, we systematically observed negative correlations between deletion and duplication effect sizes across all levels of biological observations. We quantified the preferential effects of deletions versus duplication using tagDS, a new normalized metric. Cognitive ability was preferentially affected by cortical, presynaptic, and negative-regulation gene-sets when duplicated. In contrast, preferential effects of deletions were observed for subcortical, post-synaptic, and positive-regulation gene-sets. A large proportion of gene-sets assigned to non-brain organs were associated with cognitive ability due to low tissue specificity genes, which were associated with higher sensitive to haploinsufficiency. Overall, most biological functions associated with cognitive ability are divided into those sensitive to either deletion or duplications.

## Repository architecture
<pre>
└── scripts
    ├── Compute_effect_sizes_for_genesets.ipynb
    ├── Compute_effect_sizes_permutations.ipynb
    ├── Compute_effect_sizes_random_gene_sets.ipynb
    ├── figures_scripts
        ├── Compare_tagDS_Collins.ipynb
        ├── Correlation_tagDS.ipynb
        ├── Correlation_tagDS_LOEUF_strat.ipynb
        ├── create_HPA_matrix.ipynb
        ├── create_matrix_celltype.ipynb
        ├── Matrix_GTEx.ipynb
        └── represent_significant_GO_terms.ipynb
    ├── functions_compute_effect_sizes.py
    ├── functions_loading_data.py
    └── Random_tagDS_distribution.ipynbpute_effect_sizes_for_genesets.ipynb

</pre>
