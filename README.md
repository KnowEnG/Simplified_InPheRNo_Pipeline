# Simplified_InPheRNo_Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Simplified InPheRNo Pipeline.
**** 

# Input Files
| **Example Data File** | **Requirements** |
| --------------------------------------- | ---------------------------------------- |
| /TF_Ensemble.csv | csv/tsv, no-header - names of regulators (TFs) |
| /Pvalue_gene_phenotype_interest.csv | csv/tsv, genes x 1-p-value (with header) |
| /expr_sample.csv | csv/tsv, gene/TF x samples |

# Output Files
| **Example Data File** | **Format** |
| --------------------------------------- | ---------------------------------------- |
| /Network_statistic.csv | csv, genes x genes |
| /Network_pvalue.csv | csv genes x genes |

****
# How to install and run this pipeline with the example data.
1) Clone this repository to your directory with all KnowEnG python3 libraries installed.

```git clone https://github.com/dlanier/Simplified_InPheRNo_Pipeline.git```

2) Change to the Simplified_InPheRNo_Pipeline/test directory.

```cd Simplified_InPheRNo_Pipeline/test```

3) Set up the environment.

```make env_setup```

4) Run the default data to test the installation.

```make run_InPheRNo_simplified```

****
# How to run with your data:
1) Use steps 1 - 3 above to setup the environment and place the template yaml file in the run_directory.

2) Edit the yaml file to reflect your input and output directory and input file names.

3) Run step 4 above.
****

# Run with docker

[Simplified InPheRNo Pipeline on Docker Hub](https://hub.docker.com/r/knowengdev/simplified_inpherno_pipeline/)
