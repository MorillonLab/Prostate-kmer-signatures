

Scripts for reproducing results of "Reference-free transcriptome exploration reveals novel RNAs for prostate cancer diagnosis"


 1. **ROC-AUC_for_contigs.R**: Infers a signature from a set of RNA contigs obtained in the Discovery Set (PAIR cohort) and tests this signature in the Validation Set (TCGA-PRAD cohort). Results are shown as AUC.

 2. **ROC-AUC_for_genes.R**: Infers a signature using conventional gene expression counting. The script starts from a gene expression table in the Discovery Set, then performs selection using DE-seq and stability selection using penalized logistic regression. The signature is then used to build a predictor in the Validation Set.

 3. **ROC-AUC_for_random_contigs.R**: This script does the same as ROC-AUC_for_contigs.R but each RNA contig is quantified in the Validation Set using a set of randomly sampled kmers.

