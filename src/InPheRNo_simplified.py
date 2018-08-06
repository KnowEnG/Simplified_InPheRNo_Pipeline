# -*- coding: utf-8 -*-
"""
@author: emad2

The script is an implementation of InPheRNo-simplified. This is a method for 
obtaining phenotype-relevant co-expression networks. This method works by 
combining the p-value of association of (transcription factors, genes) and the
p-value of association between (gene-phenotype) using Fisher's method.


This script uses the normalized transcriptomic data to generate P-value of 
gene-TF correlation (Pearson's correlation).

As input, this script takes in 3 files: 1) A list of all transcription factors 
(TFs), 2) a file containing p-values of gene-phenotype associations for genes 
of interest, and 3) the expression profiles (properly normalized) of genes and 
TFs. If there are shared gene name sbetween list of TFs and gene-phenotype pvalue
file, the script drops those. First, the p-value of Pearson's correlation between
expresison profiles of genes and TFs are calculated; then these p-values are
combined using Fisher's method to obtain a final p-value. 

The script generates two outputs: a (gene x TF) file containing the combined p-values
and a (gene x TF) file containign the statistic corresponding to those p-values. 
Any of these files can be thresholded in order to obtain a final phenotype-relevant
co-expression network. 


"""

import numpy as np
import pandas as pd
import os
import scipy.stats as ss
import argparse


###############################################################################
# Parse command line options/arguments
parser = argparse.ArgumentParser()

parser.add_argument('-id', '--input_directory', default = './Data', help = 'Address of directory containing input files')
parser.add_argument('-od', '--output_directory', default = './Results', help = 'output directory adddress')
parser.add_argument('-it', '--input_tf', default = 'TF_Ensemble.csv', help = 'Name of the file containing list of TFs in a csv file. The file should not have a header.')
parser.add_argument('-ie', '--input_expression', default = 'expr_sample.csv', help = 'A file containing gene and TF expression data (gene x samples). The file has a header (sample names).')
parser.add_argument('-igp', '--input_gene_phenotype_interest', default = 'Pvalue_gene_phenotype_interest.csv', help = 'A file (gene x pvalue) containing p-values of gene-phenotype only for genes of interest (and not all genes), sorted in an ascending order based on the p-value (smallest p-values appear first). Only include genes of interest to reduce computation time. The file has a header.')
parser.add_argument('-onp', '--output_network_pvalue', default = 'Network_pvalue.csv', help = 'A file (gene x TF) representing the network in which the values correspond to the combined p-values of (gene-phenotype) and (TF-gene) assoiations.')
parser.add_argument('-ons', '--output_network_statistic', default = 'Network_statistic.csv', help = 'A file (gene x TF) representing the network in which the values correspond to the statistic of the combined p-values of (gene-phenotype) and (TF-gene) assoiations.')

args = parser.parse_args()

###############################################################################
delim_tl = ','
delim_ex = ','
delim_gp = ','
if args.input_tf[-3:] in ['tsv', 'txt']:
    delim_tl = '\t'
if args.input_expression[-3:] in ['tsv', 'txt']:
    delim_ex = '\t'
if args.input_gene_phenotype_interest[-3:] in ['tsv', 'txt']:
    delim_gp = '\t'


address_TF = os.path.join(args.input_directory, args.input_tf)
address_out_dir = args.output_directory
if not os.path.exists(address_out_dir):
    os.makedirs(address_out_dir)

address_in_expr = os.path.join(args.input_directory, args.input_expression) 
address_in_gene_pheno = os.path.join(args.input_directory, args.input_gene_phenotype_interest) 

address_out_pvalue = os.path.join(args.output_directory, args.output_network_pvalue)
address_out_statistic = os.path.join(args.output_directory, args.output_network_statistic)




##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
TF_list = list(pd.read_csv(address_TF, sep=delim_tl, header=None, index_col=0).index.values)
expr_all = pd.read_csv(address_in_expr, sep=delim_ex, index_col=0)    
gene_pheno_pval = pd.read_csv(address_in_gene_pheno, sep=delim_gp, index_col=0)
gene_pheno_pval = gene_pheno_pval.dropna()
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########

eps = 3e-308

gene_pheno_list = list(gene_pheno_pval.index.values)
gene_TF_list = list(expr_all.index.values)


only_TF_list = list(set(TF_list).intersection(set(gene_TF_list)))
only_TF_list.sort()
only_gene_list = list(set(gene_TF_list) - set(TF_list))
only_gene_list = list(set(only_gene_list).intersection(gene_pheno_list))
only_gene_list.sort()


#find intersetion of gene_pheno_pval file and genes with expression values
gene_pheno_pval = gene_pheno_pval.loc[only_gene_list]
pvalue_gp_df = gene_pheno_pval.sort_values(by=gene_pheno_pval.columns[0])    #sort
only_gene_list = list(pvalue_gp_df.index.values)

expr_gene = expr_all.loc[only_gene_list]    #this contains expression of all samples
expr_TF = expr_all.loc[only_TF_list]


###############################################################################
####TF-gene correlation
tf_gene_cor_pval = np.zeros((len(only_gene_list),len(only_TF_list))) 
for i_g in range(len(only_gene_list)):
    print('calculating gene-tf correlation for gene number', i_g)
    gene = only_gene_list[i_g]
    for i_t in range(len(only_TF_list)):
        tf = only_TF_list[i_t]
        mask = ~np.isnan(expr_gene.loc[gene]) & ~np.isnan(expr_TF.loc[tf])
        expr_gene_tmp = expr_gene.loc[gene][mask]
        expr_TF_tmp = expr_TF.loc[tf][mask]
        _, pval_tmp = ss.pearsonr(expr_gene_tmp, expr_TF_tmp)
        if pval_tmp < eps:
            pval_tmp = eps
        tf_gene_cor_pval[i_g, i_t] = pval_tmp

tf_gene_cor_pval_DF = pd.DataFrame(tf_gene_cor_pval, index=only_gene_list, columns=only_TF_list)

###############################################################################
####Combinging p-values
statistic = np.zeros((len(only_gene_list), len(only_TF_list)))
pval_comb = np.zeros((len(only_gene_list), len(only_TF_list)))
i_g = 0
for gene in only_gene_list:
    print('Combining p-values for gene number', i_g)
    i_t = 0
    for tf in only_TF_list:
        statistic[i_g][i_t], pval_comb[i_g][i_t] = ss.combine_pvalues([tf_gene_cor_pval_DF.loc[gene][tf], pvalue_gp_df.loc[gene]['PValue']], method='fisher')
        i_t += 1
    i_g += 1
statistic_DF = pd.DataFrame(statistic, index=only_gene_list, columns=only_TF_list)
pval_comb_DF = pd.DataFrame(pval_comb, index=only_gene_list, columns=only_TF_list)


##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
statistic_DF.to_csv(address_out_statistic) #sorted pvalues of gene-pheno
pval_comb_DF.to_csv(address_out_pvalue) #sorted pvalues of gene-TF (sorted based on gene-pheno)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########




