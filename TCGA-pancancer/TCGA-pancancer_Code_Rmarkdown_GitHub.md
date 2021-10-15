---
title: "TCGA-pancancer"
Author: "Bytes-to-viz"
Date: "August 2020"
output:
  pdf_document: 
    latex_engine: xelatex
  word_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Code from the Science Translational Medicine article about predicting patients' prognosis based on their transcriptomic phenotype

This R Markdown document is part of a series containing the code that was used to perform some of the analyses from the article "Loss of TGFB signaling increases alternative end-joining DNA repair that sensitizes to genotoxic therapies across cancer types", published on Science Translational Medicine on 2021 (link: https://stm.sciencemag.org/content/13/580/eabc4465).

Specifically, this is the code that was used to analyze the pancancer dataset from The Cancer Genome Atlas (TCGA). Many of the results from this analysis can be seen in Figure 6 of the article.



## PART 1: PREPARE THE ENVIRONMENT

#Step 1: Load (+/- install) the necessary packages.
```{r message = FALSE, warning = FALSE}
library("GSVA")
library("GSEABase")
library("methods")
library("edgeR")
library("geneplotter")
library("genefilter")
library("BiocGenerics")
library("Biobase")
library("graph")
library("XML")
library("lattice")
library("limma")
library("shinythemes")
library("shiny")
library("RColorBrewer")
library("parallel")
library("cluster")
library("Matrix")
library("locfit")
library("snow")
library("dplyr")
library("ggplot2")
library ("remotes")
library("OIsurv") 
library("survival")
library("KMsurv")
library("splines")
library("survminer")
library("readxl")
library("reshape2")
library("data.table")
```



## PART 2: IMPORT THE FILES FOR THE ANALYSIS

#Step 1: Import the files with the ssGSEA scores.
```{r message = FALSE, warning = FALSE}
sBLCA = fread("Input/ssgsea scores/BLCA.ssgseas.tsv", data.table=FALSE)
sBRCA = fread("Input/ssgsea scores/BRCA.ssgseas.tsv", data.table=FALSE)
sCOAD = fread("Input/ssgsea scores/COAD.ssgseas.tsv", data.table=FALSE)
sESCA = fread("Input/ssgsea scores/ESCA.ssgseas.tsv", data.table=FALSE)
sGBM = fread("Input/ssgsea scores/GBM.ssgseas.tsv", data.table=FALSE)
sHNSC = fread("Input/ssgsea scores/HNSC.ssgseas.tsv", data.table=FALSE)
sKIRC = fread("Input/ssgsea scores/KIRC.ssgseas.tsv", data.table=FALSE)
sLIHC = fread("Input/ssgsea scores/LIHC.ssgseas.tsv", data.table=FALSE)
sLUAD = fread("Input/ssgsea scores/LUAD.ssgseas.tsv", data.table=FALSE)
sLUSC = fread("Input/ssgsea scores/LUSC.ssgseas.tsv", data.table=FALSE)
sOV = fread("Input/ssgsea scores/OV.ssgseas.tsv", data.table=FALSE)
sPAAD = fread("Input/ssgsea scores/PAAD.ssgseas.tsv", data.table=FALSE)
sPRAD = fread("Input/ssgsea scores/PRAD.ssgseas.tsv", data.table=FALSE)
sSKCM = fread("Input/ssgsea scores/SKCM.ssgseas.tsv", data.table=FALSE)
sTGCT = fread("Input/ssgsea scores/TGCT.ssgseas.tsv", data.table=FALSE)
sTHCA = fread("Input/ssgsea scores/THCA.ssgseas.tsv", data.table=FALSE)
sUCEC = fread("Input/ssgsea scores/UCEC.ssgseas.tsv", data.table=FALSE)
##These files were sent by collaborators (Miquel Angel, Roderic and Luis) on March 2020.  
##Contain the ssSGSEA scores from TCGA-pancancer patients calculated by them. 
```

#Step 2: Merge all the ssGSEA scores in one file. 
```{r message = FALSE, warning = FALSE}
ssgsea <- rbind(sBLCA,sBRCA,sCOAD,sESCA,sGBM,sHNSC,sKIRC,sLIHC,sLUAD,sLUSC,sOV,sPAAD,sPRAD,sSKCM,sTGCT,sTHCA,sUCEC)
ssgsea$sampleID <- ssgsea$V1
ssgsea$sampleID <- chartr(".", "-", ssgsea$sampleID) #turn "."s into "-"s
ssgsea$sampleID <- substring(ssgsea$sampleID,1,15) #keep only characters 1-15
##Dimensions: 7115 samples x 5 variables.
```

#Step 3: Import the file with clinical information. 
```{r message = FALSE, warning = FALSE}
PanCancer_ClinicalInformation <- read_excel("Input/PanCancer_ClinicalInformation.xlsx")
##This file contatins clinical information from TCGA-pancancer patients with primary solid tumors. 
##It was sent by a collaborator (Mao) on March 2020, who generated it from the from the integrated TCGA pan-cancer clinical dataset from the Genomic Data Commons (GDC).
##Dimensions: 10967 samples x 13 variables.
```



## PART 3: PUT ALL THE INFORMATION THAT WE WILL NEED FOR FURTHER ANALYSES IN ONE SINGLE DATAFRAME 

#Step 1: Merge the "ssgsea" and the "PanCancer_ClinicalInformation" dataframes.
```{r message = FALSE, warning = FALSE}
ssgsea$ID2 <- substr(ssgsea$sampleID, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(ssgsea, PanCancer_ClinicalInformation, 
             by.x = "ID2", by.y = "Sample ID",
             all.x=FALSE, all.y=FALSE) #7115 patients with ssGSEA scores, 6949 of them with clinical data.
dim(ALL) #[1] 6949  18
```

#Step 2: Prepare the survival variables so that they are in the appropriate format.
```{r message = FALSE, warning = FALSE, results='hide'}
table(ALL$`Overall Survival Status`) 
ALL$OSstatus <- ifelse(ALL$`Overall Survival Status` == "DECEASED", 1, 
                       ifelse(ALL$`Overall Survival Status` == "LIVING", 0, NA)) 
class(ALL$OSstatus) #[1] "numeric"
class(ALL$`Overall Survival (Months)`) #[1] "character"
ALL$`Overall Survival (Months)` <- as.numeric(ALL$`Overall Survival (Months)`)
```



## PART 4: ELIMINATE RECURRENT AND NORMAL TISSUE SAMPLES 

#Step 1: Create a variable that shows the sample type. 
```{r message = FALSE, warning = FALSE}
ALL$Sample.type.2 <- substr(ALL$ID2, 13, 15) #Keep only sample type numbers.
table(ALL$Sample.type.2) #6949 primary; 0 normal; 0 recurrent.
```



## PART 5: CALCULATE THE BALT SCORE OF EACH SAMPLE

#Step 1: Create a new variable that is be the Balt score.
```{r message = FALSE, warning = FALSE, results='hide'}
ALL$balt <- sqrt((max(ALL$ALT_EJ_repair)-ALL$ALT_EJ_repair)^2+
                   (min(ALL$Upregulated_TGF_beta)-ALL$Upregulated_TGF_beta)^2) - 
            sqrt((min(ALL$ALT_EJ_repair)-ALL$ALT_EJ_repair)^2+
                   (max(ALL$Upregulated_TGF_beta)-ALL$Upregulated_TGF_beta)^2)
ALL$balt <- ALL$balt * -1 
```


## PART 6: ANALYSIS OF THE CORRELATION BETWEEN TGFB AND ALTEJ SSGSEA SCORES 

#Step 1: Create a scatterplot of TGFB versus ALTEJ ssGSEA scores, coloring the cancer type.
```{r message = FALSE, warning = FALSE, results='hide'}
ggplot(ALL, aes(x=ALL$Upregulated_TGF_beta, y=ALL$ALT_EJ_repair, col=ALL$`TCGA PanCanAtlas Cancer Type Acronym`)) + 
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(x = "TGFβ ssGSEA score", y = "Alt-EJ ssGSEA score", color = "Cancer type") + 
  theme(legend.title=element_text(size=8)) 
```

#Step 2: Calculate Pearson correlation coefficient.
```{r message = FALSE, warning = FALSE}
cor.test(ALL$Upregulated_TGF_beta, ALL$ALT_EJ_repair, method = "pearson")
```

#Step 3: Calculate the Pearson correlation coefficient BY CANCER TYPE.
```{r message = FALSE, warning = FALSE, results='hide'}
library(broom)
TumCor <- ALL %>% 
  group_by(`TCGA PanCanAtlas Cancer Type Acronym`) %>% 
  do(tidy(cor.test(.$Upregulated_TGF_beta, .$ALT_EJ_repair, method = "pearson"))) #All anticorrelated. All significantly except PAAD and KIRC.
```



## PART 7: ANALYSIS OF THE CORRELATION OF "AGE VS BALT SCORE"  

#Step 1: Calculate Pearson correlation coefficient of age vs BAlt score.
```{r message = FALSE, warning = FALSE}
cor.test(ALL$balt, ALL$`Diagnosis Age`, method = "pearson") #PCC=-0.0396
```

#Step 2: Create a scatterplot of Age vs BAlt score.
```{r message = FALSE, warning = FALSE, results='hide'}
ggplot(ALL, aes(x=ALL$balt, y=ALL$`Diagnosis Age`)) + geom_point(col="grey30") +
  labs(x = "βAlt score", y = "Age (years)") + 
  theme(legend.title=element_text(size=8))  +
  geom_smooth(method = "lm", fill = NA, size=2)
```



## PART 8: SELECT PATIENTS WHO PROBABLY RECEIVED GENOTOXIC TREATMENT 

#Step 1-A: Create a subdataframe with patients treated with RT.
```{r message = FALSE, warning = FALSE, results='hide'}
table(ALL$`Radiation Therapy`)
ALL_RT <- ALL[which(ALL$`Radiation Therapy`=="Yes"),] #6.949-->1.737 patients. 
```
```{r message = FALSE, warning = FALSE}
dim(ALL_RT)
```

#Step 2-A: Create a new variable that represents the tertile that each patient belongs to, according to the value of the BAlt score.
```{r message = FALSE, warning = FALSE, results='hide'}
ALL_RT <- ALL_RT %>% mutate(tertile = ntile(balt,3)) 
table(ALL_RT$`TCGA PanCanAtlas Cancer Type Acronym`, ALL_RT$tertile)
ALL_RT$Tertile <- ifelse(ALL_RT$tertile == 2, NA, ALL_RT$tertile)
```

#Step 1-B: Create a subdataframe with patients who, based on their cancer type and stage, their standard of care treatment includes RT and/or genotoxic ChT.
```{r message = FALSE, warning = FALSE, results='hide'}
table(ALL$`Cancer Type`, ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code`)
table(ALL$`TCGA PanCanAtlas Cancer Type Acronym`, ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code`)
ALL_RTChT <- ALL[which(ALL$`Radiation Therapy` == "Yes" |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "BLCA" & ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IV") |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "COAD" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE I" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" ) |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "ESCA" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE I" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IB" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IIA") |
                            ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "GBM" |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "HNSC" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE I" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IB" ) |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "LUAD" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE I" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IB" ) |
                            (ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "LUSC" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE I" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IB" ) |
                            ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "OV" & !ALL$`Neoplasm Disease Stage American Joint Committee on Cancer Code` == "STAGE IA" |
                            ALL$`TCGA PanCanAtlas Cancer Type Acronym` == "PAAD" |
                            ALL$`Cancer Type`== "Non-Seminomatous Germ Cell Tumor"
),] #6.949-->3.577 patients.
```
```{r message = FALSE, warning = FALSE}
dim(ALL_RTChT)
```
```{r message = FALSE, warning = FALSE, results='hide'}
##Inclusion criteria. RT yes or:
##BLCA: Stage 4.
##BRCA: RT.
##COAD: Stage not 1.
##ESCA (Esophageal): Stage not 1 nor 2a.
##GBM: All.
##HNSC: Stage not 1.
##KIRC (Renal): RT.
##LIHC (Liver): RT.
##LUAD and LUSC: Stage not 1.
##OV: Stage not 1a.
##PAAD (Pancreas): All.
##PRAD (Prostate): RT.
##SKCM: RT. 
##TGCT seminoma: RT. 
##TGCT non-seminoma: All.
##THCA (Thyroid): RT.
##UCEC (Endometrial): RT.
```
```{r message = FALSE, warning = FALSE, results='hide'}
##Reasoning behind the inclusion criteria:
##BLCA: Stage 4 usually includes genotoxic ChT with platin agents and earlier stages are frequently treated with surgery alone or followed by BCG or intracavitary non-genotoxic ChT. 
##BRCA: Usually treated with surgery + RT +/- HT +/- Trastuzumab +/- ChT. ChT frequently consists of anthracyclines and taxanes, so many times it does not include genotoxic drugs. 
##COAD: Stage I may be treated with surgery alone. In other stages ChT usually includes genotoxic Platin agents (FOLFOX/FOLFIRI). 
##ESCA: Stages T1-2N0M0 (stages <IIB) may be treated with surgery alone. Otherwise RT, ChT or both are usually used. ChT usually includes genotoxic platin agents. 
##GBM: Standard treatment is usually surgery + RT + ChT. ChT usually consists of genotoxic Temozolamide. 
##HNSC: Stage I may be treated with surgery alone. Otherwise treatment usually includes RT, ChT or both (most frequently). ChT usually includes genotoxic platin agents (TPF for induction, CDDP in concomitance with RT, EXTREME in palliative cases).
##KIRC: RT is rarely indicated. Many times treated with surgery alone and the systemic treatment, when given, rarely includes genotoxic drugs. 
##LIHC: Many treated with surgery alone. Systemic treatment rarely includes genotoxic drugs, as the most frequent agent is sorafenib. 
##LUAD and LUSC: Stage I may be treated with surgery alone. Otherwise treatment usually includes RT, ChT or both (most frequently). ChT usually includes genotoxic platin agents.
##OV: Stage Ia may be treated with surgery alone and stage>Ia usually receive ChT. ChT usually includes genotoxic platin agents.
##PAAD: Rarely treated with surgery alone. RT and/or ChT are frequently added. ChT most usually includes pirimidine analogues (such as gemcitabine, capecitabine or 5FU) or FOLFIRINOX.
##PRAD: Usually treated with HT combined with surgery, RT or both. ChT is given only in some stage IV patients and many times it does not include genotoxic drugs, as it usually consists of taxane drugs. 
##SKCM: Many treated with surgery alone. Systemic treatment rarely includes genotoxic drugs. 
##TGCT: Treated with surgery. In seminomas, RT is usually added (and/or ChT if stage >I) and, in non-seminomas, ChT. ChT usually includes platin agents (BEP). 
##THCA: Most of them are treated with surgery alone. Few exceptions are treated with both RT and ChT, such as anaplastic carcinomas. 
##UCEC: Stage I may be treated with surgery alone. Otherwise RT is usually given. ChT is used mainly in stage IV and even so not always, as if possible HT is preferred.
```

#Step 2-B: Create a new variable that represents the tertile that each patient belongs to, according to the value of the BAlt score.
```{r message = FALSE, warning = FALSE, results='hide'}
ALL_RTChT <- ALL_RTChT %>% mutate(tertile = ntile(balt,3)) 
table(ALL_RTChT$`TCGA PanCanAtlas Cancer Type Acronym`, ALL_RTChT$tertile)
ALL_RTChT$Tertile <- ifelse(ALL_RTChT$tertile == 2, NA, ALL_RTChT$tertile)
```


## PART 9: OVERALL SURVIVAL CURVES OF THE BALT SCORE GROUPS 

#Step 1: Plot and compare the OS curves between Balt score top and bottom tertiles.
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT
my.surv.object <- Surv(time=ALL_RT$`Overall Survival (Months)`, event=ALL_RT$OSstatus)
my.fit<-survfit(my.surv.object~ALL_RT$Tertile)
my.fit
ggsurvplot(my.fit, data=ALL_RT, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "royalblue4"))
```
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT and/or probably genotoxic ChT
my.surv.object <- Surv(time=ALL_RTChT$`Overall Survival (Months)`, event=ALL_RTChT$OSstatus)
my.fit<-survfit(my.surv.object~ALL_RTChT$Tertile)
my.fit
ggsurvplot(my.fit, data=ALL_RTChT, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "royalblue4"))
```



## PART 10: CALCULATE THE HAZARD RATIOS OF THE SURVIVAL CURVES FROM PART 9 

#Step 1: In the variable tertile, make "high TGFB and low ALTEJ" the reference level so that the Hazard Ratios are expressed as 0.x instead of 1.x. 
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT
ALL_RT$Tertile <- ifelse(ALL_RT$Tertile == "1", "high TGFB and low ALTEJ", 
                            ifelse(ALL_RT$Tertile == "3", "low TGFB and high ALTEJ", NA))
ALL_RT$Tertile <- factor(ALL_RT$Tertile, levels = c("high TGFB and low ALTEJ", "low TGFB and high ALTEJ"))
table(ALL_RT$Tertile)
```
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT and/or genotoxic ChT
ALL_RTChT$Tertile <- ifelse(ALL_RTChT$Tertile == "1", "high TGFB and low ALTEJ", 
                               ifelse(ALL_RTChT$Tertile == "3", "low TGFB and high ALTEJ", NA))
ALL_RTChT$Tertile <- factor(ALL_RTChT$Tertile, levels = c("high TGFB and low ALTEJ", "low TGFB and high ALTEJ"))
table(ALL_RTChT$Tertile)
```

#Step 2: Calculate the OS hazard ratio of BAlt tertile 1 versus 3.
```{r message = FALSE, warning = FALSE}
##In patients treated with RT
my.surv.object <- Surv(time=ALL_RT$`Overall Survival (Months)`, event=ALL_RT$OSstatus)
cox<-coxph(my.surv.object ~  ALL_RT$Tertile)
summary(cox) #HR=0.56
```
```{r message = FALSE, warning = FALSE}
##In patients treated with RT and/or genotoxic ChT
my.surv.object <- Surv(time=ALL_RTChT$`Overall Survival (Months)`, event=ALL_RTChT$OSstatus)
cox<-coxph(my.surv.object ~  ALL_RTChT$Tertile)
summary(cox) #HR=0.60
```


## PART 11: SCATTERPLOTS SHOWING THE TWO BALT SCORE GROUPS COMPARED IN PART 9 

#Step 1: Create a scatterplot of TGFB versus ALTEJ ssGSEA scores, coloring the BAlt tertiles.
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT
ggplot(ALL_RT, aes(x=ALL_RT$Upregulated_TGF_beta, y=ALL_RT$ALT_EJ_repair, col=ALL_RT$Tertile)) + 
  geom_point() +
  scale_color_manual(values=c("orange", "royalblue4"), 
                     labels = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
                     na.translate=TRUE, na.value="grey") + 
  labs(x = "TGFβ ssGSEA score", y = "Alt-EJ ssGSEA score", color = "βAlt score groups") + 
  theme(legend.title=element_text(size=8)) 
```
```{r message = FALSE, warning = FALSE}
cor.test(ALL_RT$Upregulated_TGF_beta, ALL_RT$ALT_EJ_repair, method = "pearson") #PCC=-0.23
```

##In patients treated with RT and/or genotoxic ChT
```{r message = FALSE, warning = FALSE, results='hide'}
ggplot(ALL_RTChT, aes(x=ALL_RTChT$Upregulated_TGF_beta, y=ALL_RTChT$ALT_EJ_repair, col=ALL_RTChT$Tertile)) + 
  geom_point() +
  scale_color_manual(values=c("orange", "royalblue4"), 
                     labels = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
                     na.translate=TRUE, na.value="grey") + 
  labs(x = "TGFβ ssGSEA score", y = "Alt-EJ ssGSEA score", color = "βAlt score groups") + 
  theme(legend.title=element_text(size=8)) 
```
```{r message = FALSE, warning = FALSE}
cor.test(ALL_RTChT$Upregulated_TGF_beta, ALL_RTChT$ALT_EJ_repair, method = "pearson") #PCC=-0.16
```


## PART 12: MULTIVARIATE COX OF THE ASSOCIATION OF THE BALT SCORE WITH OS, ADJUSTED FOR AGE AND STAGE 

#Step 1: Group stages into 4 groups (1, 2, 3 and 4) 
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT
table(ALL_RT$`Neoplasm Disease Stage American Joint Committee on Cancer Code`)
ALL_RT$stage <- ALL_RT$`Neoplasm Disease Stage American Joint Committee on Cancer Code`
ALL_RT$stage <- ifelse(ALL_RT$stage == "STAGE I", "1",
                          ifelse(ALL_RT$stage == "STAGE IA", "1",
                          ifelse(ALL_RT$stage == "STAGE IB", "1",
                          ifelse(ALL_RT$stage == "STAGE II", "2",
                          ifelse(ALL_RT$stage == "STAGE IIA", "2",
                          ifelse(ALL_RT$stage == "STAGE IIB", "2",
                          ifelse(ALL_RT$stage == "STAGE IIC", "2",
                          ifelse(ALL_RT$stage == "STAGE III", "3",
                          ifelse(ALL_RT$stage == "STAGE IIIA", "3",
                          ifelse(ALL_RT$stage == "STAGE IIIB", "3",
                          ifelse(ALL_RT$stage == "STAGE IIIC", "3",
                          ifelse(ALL_RT$stage == "STAGE IS", "1",
                          ifelse(ALL_RT$stage == "STAGE IV", "4",
                          ifelse(ALL_RT$stage == "STAGE IVA", "4",
                          ifelse(ALL_RT$stage == "STAGE IVB", "4",
                          ifelse(ALL_RT$stage == "STAGE IVC", "4", NA))))))))))))))))
table(ALL_RT$stage)
```
```{r message = FALSE, warning = FALSE, results='hide'}
##In patients treated with RT and/or genotoxic ChT
table(ALL_RTChT$`Neoplasm Disease Stage American Joint Committee on Cancer Code`)
ALL_RTChT$stage <- ALL_RTChT$`Neoplasm Disease Stage American Joint Committee on Cancer Code`
ALL_RTChT$stage <- ifelse(ALL_RTChT$stage == "STAGE I", "1",
                            ifelse(ALL_RTChT$stage == "STAGE IA", "1",
                            ifelse(ALL_RTChT$stage == "STAGE IB", "1",
                            ifelse(ALL_RTChT$stage == "STAGE II", "2",
                            ifelse(ALL_RTChT$stage == "STAGE IIA", "2",
                            ifelse(ALL_RTChT$stage == "STAGE IIB", "2",
                            ifelse(ALL_RTChT$stage == "STAGE IIC", "2",
                            ifelse(ALL_RTChT$stage == "STAGE III", "3",
                            ifelse(ALL_RTChT$stage == "STAGE IIIA", "3",
                            ifelse(ALL_RTChT$stage == "STAGE IIIB", "3",
                            ifelse(ALL_RTChT$stage == "STAGE IIIC", "3",
                            ifelse(ALL_RTChT$stage == "STAGE IS", "1",
                            ifelse(ALL_RTChT$stage == "STAGE IV", "4",
                            ifelse(ALL_RTChT$stage == "STAGE IVA", "4",
                            ifelse(ALL_RTChT$stage == "STAGE IVB", "4",
                            ifelse(ALL_RTChT$stage == "STAGE IVC", "4", NA))))))))))))))))
table(ALL_RTChT$stage)
```

#Step 2: Calculate survival multivariate Cox regression adjusted for age and stage of the BAlt score continuous variable. 
```{r message = FALSE, warning = FALSE}
##In patients treated with RT
my.surv.object <- Surv(time=ALL_RT$`Overall Survival (Months)`, event=ALL_RT$OSstatus)
cox<-coxph(my.surv.object ~  ALL_RT$balt + ALL_RT$`Diagnosis Age` + ALL_RT$stage)
summary(cox) #P=0.003
```
##In patients treated with RT and/or genotoxic ChT
```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL_RTChT$`Overall Survival (Months)`, event=ALL_RTChT$OSstatus)
cox<-coxph(my.surv.object ~  ALL_RTChT$balt + ALL_RTChT$`Diagnosis Age` + ALL_RTChT$stage)
summary(cox) #P<0.001
```


