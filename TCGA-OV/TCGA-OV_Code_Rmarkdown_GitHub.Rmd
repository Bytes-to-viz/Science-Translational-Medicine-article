---
title: "TCGA-OV"
Author: "Bytes-to-viz"
Date: "August 2020"
output:
  html_document: default
  word_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Code from the Science Translational Medicine article about predicting patients' prognosis based on their transcriptomic phenotype

This R Markdown document is part of a series containing the code that was used to perform some of the analyses from the article "Loss of TGFβ signaling increases alternative end-joining DNA repair that sensitizes to genotoxic therapies across cancer types", published on Science Translational Medicine on 2021 (link: https://stm.sciencemag.org/content/13/580/eabc4465).

Specifically, this is the code that was used to analyze the ovarian cancer (OV) dataset from The Cancer Genome Atlas (TCGA). Many of the results from this analysis can be seen in Figure 4 of the article.


## Part 0: Prepare the environment

Step 1: Load (+/- install) the necessary packages.

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
```


## Part 1: Import the input files

Step 1: Import the file with gene expression information.

```{r message = FALSE, warning = FALSE}
GDC_ov_genes <- read_excel("Input/GDC_ov_genes.xlsx") 
#This file was downloaded from GDC using the R package “TCGAbiolinks” in January 2020.
#Contains expression info measured by Agilent Microarray (AgilentG4502A_07_3). 
#Dimensions: 562 samples x 16210 genes. Some samples are from the same patients and some are from normal tissue. 
```

Step 2: Import the files with clinical information. 

```{r message = FALSE, warning = FALSE}
ov_tcga_pan_can_atlas_2018_clinical_data <- read.delim("Input/ov_tcga_pan_can_atlas_2018_clinical_data.tsv") 
#This file was downloaded on April 2 2020 from cBioPortal from the project "Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas)".
#This file contains all the clinical information that we will need, except the "Figo stage".
#Dimensions: 585 samples x 96 variables.

```

```{r message = FALSE, warning = FALSE}
GDC_ov_clinical <- read_excel("Input/GDC_ov_clinical.xlsx") 
#This file was downloaded from GDC using the R package “TCGAbiolinks” in January 2020.
#We will use this file to have "Figo stage" information from the patients.
#Dimensions: 562 samples x 63 variables.
```

Step 3: Import the TGFB and ALTEJ genelists. 

```{r message = FALSE, warning = FALSE}
TGFBgeneset <- read.table("Input/TGFBUPgeneset.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
ALTEJgeneset <- read.table("Input/ALTEJgeneset.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
```

Step 4: Turn the genelistS into lists and create a list containing both of them.

```{r message = FALSE, warning = FALSE}
TGFBlist <- as.list(TGFBgeneset[,1:1]) 
ALTEJlist <- as.list(ALTEJgeneset[,1:1]) 
Bothgenelists <- list(TGFBUPgeneset=TGFBlist, ALTEJgeneset=ALTEJlist)
```


# Part 2: Calculate the TGFβ and alt-EJ ssGSEA scores of each sample

Step 1: Turn the dataframe with gene expression into a matrix.

```{r message = FALSE, warning = FALSE, results='hide'}
rownames(GDC_ov_genes) <- GDC_ov_genes$GeneID
GDC_ov_genes$GeneID <- NULL
GDC_ov_genes <- as.matrix(GDC_ov_genes)
```

Step 2: Convert gene expression values into Z-scores ((x - mean(Xcolumn)) / sd(Xcolumn)).

```{r message = FALSE, warning = FALSE, results='hide'}
genesT <- t(GDC_ov_genes) #Transpose the matrix.
genesT<-scale(genesT, center = TRUE, scale = TRUE) #Calculate the Z-scores. IMPORTANT: To calculate the Z-score of each gene, each gene has to be a column.
```

Step 3: Calculate the ssGSEA scores of TGFB and ALTEJ signatures in each sample.

```{r message = FALSE, warning = FALSE, results='hide'}
genes<-t(genesT)
ssgsea <- gsva(genes, Bothgenelists, method=c("ssgsea"))
```


Step 4: Traspose the matrix with the ssGSEA results and turn it into a dataframe.

```{r message = FALSE, warning = FALSE, results='hide'}
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
```


# Part 3: Calculate the βAlt score of each sample 

Step 1: Create a new variable that is be the Balt score.

```{r message = FALSE, warning = FALSE, results='hide'}
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                    (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                    (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
```

# Part 4: Merge all the information in one single dataframe 

Step 1: Merge the "ssgsea" and "ov_tcga_pan_can_atlas_2018_clinical_data"dataframes.

```{r message = FALSE, warning = FALSE, results='hide'}
ssgsea$ID1 <- rownames(ssgsea) #Create a variable that are sample IDs.
ssgsea$ID2 <- substr(ssgsea$ID1, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(ssgsea, ov_tcga_pan_can_atlas_2018_clinical_data, 
             by.x = "ID2", by.y = "Sample.ID",
             all.x=TRUE, all.y=FALSE) #There are 535 patients with all the data and 562 patients with gene expression but not necessarilly clinical data.
dim(ALL) #562 samples.

```


Step 2: Add the normalized (z-score transformed) TGFβ1 expression to the "ALL" dataframe.

```{r message = FALSE, warning = FALSE, results='hide'}
TGFB1 <- as.data.frame(genesT[,c("TGFB1")])
TGFB1$TGFB1 <- TGFB1$`genesT[, c("TGFB1")]`
TGFB1$`genesT[, c("TGFB1")]` <- NULL
TGFB1$ID1 <- rownames(TGFB1) 
TGFB1$ID2 <- substr(TGFB1$ID1, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(TGFB1, ALL, 
             by.x = "ID2", by.y = "ID2",
             all.x=FALSE, all.y=TRUE) #562p + 562p --> 564p. 
dim(ALL) #564 samples.
```

Step 3: Add the "Figo stage" to the "ALL" dataframe.

```{r message = FALSE, warning = FALSE, results='hide'}
Figo <- as.data.frame(GDC_ov_clinical[, c("sample", "figo_stage")])
Figo$ID2 <- substr(Figo$sample, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(Figo, ALL, 
             by.x = "ID2", by.y = "ID2",
             all.x=FALSE, all.y=TRUE) #564p + 562p --> 568p. 
dim(ALL) #568 samples.
```


# Part 5: Eliminate duplicated, recurrent and normal tissue samples 

Step 1: Create a variable that shows the sample type. 


```{r message = FALSE, warning = FALSE, results='hide'}
ALL$Sample.type.2 <- substr(ALL$sample, 13, 15) #Keep only sample type numbers.
table(ALL$Sample.type.2) #548 primary; 3 normal; 17 recurrent.
```

Step 2: Remove samples from normal tissue.

```{r message = FALSE, warning = FALSE, results='hide'}
ALL <- ALL[-which(ALL$Sample.type.2 == "-11"), ] #568->565p. 
```

Step 3: Remove samples from recurrent tumors.

```{r message = FALSE, warning = FALSE, results='hide'}
ALL <- ALL[-which(ALL$Sample.type.2 == "-02"), ] #565->548p.
```

Step 4: Remove duplicated patients. 

```{r message = FALSE, warning = FALSE, results='hide'}
ALL$Patient.ID.2 <- substr(ALL$sample, 0, 12) #Keep only the first 12 characters of patient IDs.
Duplicated <- ALL[duplicated(ALL$Patient.ID.2), ] #7 duplicated patients. All 7 samples are from the same 1 patient.
ALL <- distinct(ALL, ALL$Patient.ID.2, .keep_all=TRUE) #548->541p. 
rownames(ALL) <- ALL$Patient.ID.2
```

Step 5: Create a new variable that represents the tertile that each patient belongs to, according to the value of the βAlt score.


```{r message = FALSE, warning = FALSE, results='hide'}
ALL$tertile <- cut(ALL$balt, 
                   quantile(ALL$balt, c(0, 1/3, 2/3, 1), na.rm=TRUE, include.lowest=TRUE), right=FALSE,
                   labels = c("low", "middle", "high"))
ALL$Tertile <- ifelse(ALL$tertile == "low", "high TGFB and low ALTEJ", 
                      ifelse(ALL$tertile == "high", "low TGFB and high ALTEJ", NA))
ALL$Tertile <- as.character(ALL$Tertile)
```


# Part 6: Analysis of the correlation between TGFβ and alt-EJ ssGSEA scores 

Step 1: Create a scatterplot of TGFβ versus ALTEJ ssGSEA scores, coloring the βAlt tertiles.


```{r message = FALSE, warning = FALSE, results='hide'}
ggplot(ALL, aes(x=ALL$TGFBUPgeneset, y=ALL$ALTEJgeneset, col=ALL$Tertile)) + geom_point() +
  scale_color_manual(values=c("orange", "blue"), 
                     labels = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
                     na.translate=TRUE, na.value="grey") + 
  labs(x = "TGFβ ssGSEA score", y = "Alt-EJ ssGSEA score", color = "βAlt score groups") + 
  theme(legend.title=element_text(size=8)) 
```

Step 2: Calculate Pearson correlation coefficient.

```{r message = FALSE, warning = FALSE}
cor.test(ALL$TGFBUPgeneset, ALL$ALTEJgeneset, method = "pearson") #PCC=-0.32
```


# Part 7: Analysis of the correlation of "TGFβ1 vs TGFβ signature" and "age vs βAlt score"

Step 1: Calculate Pearson correlation coefficient of TGFβ1 vs TGFβ ssGSEA.


```{r message = FALSE, warning = FALSE}
cor.test(ALL$TGFBUPgeneset, ALL$TGFB1, method = "pearson") #PCC=0.312
```


Step 2: Create a scatterplot of TGFβ1 versus TGFβ ssGSEA score.

```{r message = FALSE, warning = FALSE}
ggplot(ALL, aes(x=ALL$TGFBUPgeneset, y=ALL$TGFB1)) + geom_point(col="grey30") +
  labs(x = "TGFβ ssGSEA score", y = "TGFβ1") + 
  theme(legend.title=element_text(size=8)) 
```


Step 3: Calculate Pearson correlation coefficient of age vs βAlt score.

```{r message = FALSE, warning = FALSE}
cor.test(ALL$balt, ALL$Diagnosis.Age, method = "pearson") #Non-significant.
```



# Part 8: Analysis of the correlation between the βAlt score and genomic alterations 

Step 1: Plot boxplots showing the mutational load of the BAlt score top and bottom tertiles.

```{r message = FALSE, warning = FALSE}
ALL$Fraction.Genome.Altered <- as.numeric(ALL$Fraction.Genome.Altered)
boxplot(ALL$Fraction.Genome.Altered ~ ALL$Tertile,
        xlab="βAlt score groups", ylab="Fraction of the genome altered", 
        names = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
        col=c("orange", "blue"), 
        outpch=NA, cex.axis=0.8) 
stripchart(ALL$Fraction.Genome.Altered ~ ALL$Tertile, 
           vertical=TRUE, add=TRUE, method="jitter", col=c("dark orange", "blue"), 
           pch=21, bg="bisque") #Add dots to it representing each sample.

```


Step 2: Do a U Mann Whitney Test to compare Genomic alterations between the top and bottom βAlt tertiles.
```{r message = FALSE, warning = FALSE}
wilcox.test(ALL$Fraction.Genome.Altered ~ ALL$Tertile, data=ALL) 
```


# Part 9: Overall survival and progression free survival curves of the βAlt score groups 

Step 1: Prepare the OS variables so that they are usable.
```{r message = FALSE, warning = FALSE,  results='hide'}
table(ALL$Overall.Survival.Status)
ALL$OSstatus <- ifelse(ALL$Overall.Survival.Status == "DECEASED", 1, 
                       ifelse(ALL$Overall.Survival.Status == "LIVING", 0, NA)) 
class(ALL$OSstatus) #[1] "numeric"
class(ALL$Overall.Survival..Months.) #[1] "numeric"
```


Step 2: Prepare the PFS variables so that they are usable.
```{r message = FALSE, warning = FALSE,  results='hide'}
table(ALL$Progression.Free.Status) 
ALL$PFSstatus <- ifelse(ALL$Progression.Free.Status == "PROGRESSION", 1, 
                        ifelse(ALL$Progression.Free.Status == "CENSORED", 0, NA))
class(ALL$PFSstatus) #[1] "numeric"
class(ALL$Progress.Free.Survival..Months.) #[1] "numeric"
```


Step 3: Plot and compare the OS curves between βAlt score top and bottom tertiles.
```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
my.fit<-survfit(my.surv.object~ALL$Tertile)
my.fit
ggsurvplot(my.fit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "blue"))
```


Step 4: Plot and compare the PFS curves between βAlt score top and bottom tertiles.
```{r message = FALSE, warning = FALSE,  results='hide'}
my.DFsurv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
my.DFfit<-survfit(my.DFsurv.object~ALL$Tertile)
my.DFfit
ggsurvplot(my.DFfit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Progression Free Fraction", xlim=c(0, 60),   
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="PFS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "blue"))
```


# Part 10: Calculate the hazard ratios of the OS and PFS from part 9

Step 1: In the variable tertile, make "high TGFβ and low ALTEJ" the reference level so that the Hazard Ratios are 0.x instead of 1.x. 
```{r message = FALSE, warning = FALSE,  results='hide'}
table(ALL$Tertile)
ALL$Tertile <- factor(ALL$Tertile, levels = c("high TGFB and low ALTEJ","low TGFB and high ALTEJ"))
table(ALL$Tertile)
```


Step 2: Calculate the OS hazard ratio of βAlt tertile 1 versus 3.
```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
cox<-coxph(my.surv.object ~  ALL$Tertile)
summary(cox) #HR=0.67.
```

Step 3: Calculate the PFS hazard ratio of βAlt tertile 1 versus 3.
```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
cox<-coxph(my.surv.object ~  ALL$Tertile)
summary(cox) #HR=0.68.
```


# Part 11: Multivariate Cox of the association of the βAlt score with OS and PFS, adjusted for age and stage 

Step 1: Group Figo stages into 3 groups. 

```{r message = FALSE, warning = FALSE,  results='hide'}
ALL$Figo <- ifelse(ALL$figo_stage == "Stage IA", "1-2",
                   ifelse(ALL$figo_stage == "Stage IB", "1-2",
                   ifelse(ALL$figo_stage == "Stage IC", "1-2",
                   ifelse(ALL$figo_stage == "Stage IIA", "1-2",
                   ifelse(ALL$figo_stage == "Stage IIB", "1-2",
                   ifelse(ALL$figo_stage == "Stage IIC", "1-2",
                   ifelse(ALL$figo_stage == "Stage IIIA", "3",
                   ifelse(ALL$figo_stage == "Stage IIIB", "3",
                   ifelse(ALL$figo_stage == "Stage IIIC", "3",
                   ifelse(ALL$figo_stage == "Stage IV", "4", NA))))))))))
table(ALL$Figo)
```


Step 2: Calculate βAlt score (as a continuous variable) OS multivariate Cox adjusted for age and stage. 

```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
cox<-coxph(my.surv.object ~  ALL$balt + ALL$Diagnosis.Age + ALL$Figo)
summary(cox) #p=0.035.
```


Step 3: Calculate βAlt score (as a continuous variable) PFS multivariate Cox adjusted for age and stage. 

```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
cox<-coxph(my.surv.object ~  ALL$balt + ALL$Diagnosis.Age + ALL$Figo)
summary(cox) #p=0.002.
```


# Part 12: TGFβ1 association with overall survival and progresion free survival 

Step 1: Split TGFβ1 into tertiles.

```{r message = FALSE, warning = FALSE,  results='hide'}
ALL$TGFB1tertile <- cut(ALL$TGFB1, 
                        quantile(ALL$TGFB1, c(0, 1/3, 2/3, 1), na.rm=TRUE, include.lowest=TRUE), 
                        labels = c("low", "middle", "high"))
ALL$TGFB1Tertile <- ifelse(ALL$TGFB1tertile == "low", "low TGFB1", 
                           ifelse(ALL$TGFB1tertile == "high", "high TGFB1", NA))
ALL$TGFB1Tertile <- as.character(ALL$TGFB1Tertile)
```

Step 2: Plot and compare the OS curves between TGFβ1 top and bottom tertiles.

```{r message = FALSE, warning = FALSE}
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
my.fit<-survfit(my.surv.object~ALL$TGFB1Tertile)
my.fit
ggsurvplot(my.fit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ1", "low TGFβ1"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their TGFβ1 expression", 
           font.main=13, palette=c("orange", "blue"))
```


Step 3: Plot and compare the PFS curves between TGFβ1 top and bottom tertiles.

```{r message = FALSE, warning = FALSE}
my.DFsurv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
my.DFfit<-survfit(my.DFsurv.object~ALL$TGFB1Tertile)
my.DFfit
ggsurvplot(my.DFfit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Progression Free Fraction", xlim=c(0, 60),   
           legend.labs=c("high TGFβ1", "low TGFβ1"),
           legend.title=" ", title="PFS in tertiles 1 versus 3 of patients according to their TGFβ1 expression", 
           font.main=13, palette=c("orange", "blue")) 
```








