# Simplified_InPheRNo_Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Simplified InPheRNo Pipeline.
Simplified InPheRNo (Inference of Phenotype-relevant Regulatory Networks) is a computational tool to reconstruct phenotype-relevant transcriptional regulatory networks (TRNs) using transcriptomic data.
**** 

| **Example Data File*** | **Requirements** |
| --------------------------------------- | ---------------------------------------- |
| /TF_Ensemble.csv | csv/tsv, no-header - names of regulators (TFs) |
| /Pvalue_gene_phenotype_interest.csv | csv/tsv, genes x 1-p-value (with header) |
| /expr_sample.csv | csv/tsv, gene/TF x samples |

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

# Run with docker

[Simplified InPheRNo Pipeline on Docker Hub](https://hub.docker.com/r/knowengdev/simplified_inpherno_pipeline/)
