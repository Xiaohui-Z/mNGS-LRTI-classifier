# mNGS-LRTI-classifier
Host response from the clinical metagenomic sequencing data enable accurate diagnosis of lower respiratory infections


This document contains all the scripts required to performed DEG analysis and wgcna analysis, and python script to establish a classifier to distinguish LRTI from non-LRTIs


Hardware requirements

This code requires only a standard computer with enough RAM to support the in-memory operations.




OS Requirements

This package is supported for macOS and Linux. The package has been tested on the following systems:

macOS: Big Sur (11.7.10)

Python Dependencies

numpy (1.19.2)

scikit-learn (1.2.0)

pandas (1.1.3)




The python scripts were written in Python Python 3.8.9.

R scripts were written using Rstudio


Step 1: Utilize edgeR for transcriptomic analysis to identify differentially expressed genes between LRTI and non-LRTI patients.
Step 2: Perform KEGG and GO enrichment analyses on the differentially expressed genes.
Step 3: Conduct WGCNA to explore the association between macrotranscriptomic genes and clinical features of patients.
Step 4: Carry out sensitivity analysis after excluding non-bacterial LRTI cases.
Step 5: Perform feature selection to identify the optimal combination of differential genes.
Step 6: Construct a machine learning model using the selected gene features.
