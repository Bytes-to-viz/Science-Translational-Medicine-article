# Science-Translational-Medicine-article
Code used in some of the analyses from the "Science Translational Medicine" article titled "Loss of TGFβ signaling increases alternative end-joining DNA repair that sensitizes to genotoxic therapies across cancer types" (link: https://stm.sciencemag.org/content/13/580/eabc4465). 

The article is about the developement of a score -- termed βalt-- to predict cancer patient's prognosis based on their expression phenotype of the genes from two signatures (the TGFβ pathway and a DNA damage repair pathway called alternative end-joining). 

Each folder (TCGA-OV, TCGA-LUSC, TCGA_GBM) contains the code that was used to analyze each cancer dataset from The Cancer Genome Atlas (TCGA), consisting mainly of the following analytic processes:
1. Single sample gene set enrichment analysis (ssGSEA) scores calculation. 
2. βalt score calculation. 
3. Analysis of the correlation between the TGFβ and the altEJ gene expression signatures. 
4. Analysis of the association between tumors' genomic alterations and their βalt score.
5. Analysis of the association of the βalt score with patients' overall survival and progression-free survival.



