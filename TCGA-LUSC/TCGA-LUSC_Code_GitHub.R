#R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: TCGA-LUSC ANALYSIS



# PART 0: PREPARE THE ENVIRONMENT --------------------------------

#Step 1: Load the necessary packages. 
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



# PART 1: IMPORT THE FILES FOR THE ANALYSIS -----------------------------------------

#Step 1: Import the file with gene expression information.
Xena_lusc_genes <- read_excel("Input/Xena_lusc_genes.xlsx")
##This file was downloaded from "UCSC Xena cancer browser" using R package "UCSCXenaTools" in August 2020.
##Contains expression info measured by IlluminaHiSeq_RNASeqV2. Unit: pan-cancer normalized log2(norm_count+1).
##Some patients are duplicated and some samples are from normal tissue. 
##Dimensions: 553 samples x 20530 genes.

#Step 2: Import the files with clinical information. 
lusc_tcga_pan_can_atlas_2018_clinical_data <- read.delim("Input/cbioportal_lusc_tcga_pan_can_atlas_2018_clinical_data.tsv")
##This file was downloaded on on January 2020 from cBioPortal from the project "LUSC TCGA PanCancer Atlas".
##This file contains the main clinical information, such as overall survival and progression free survival.
##Dimensions: 487 samples x 96 variables.
Xena_lusc_clinical <- read_excel("Input/Xena_lusc_clinical.xlsx")
##This file was downloaded from "UCSC Xena cancer browser" using R package "UCSCXenaTools".
##This file contains some additional clinical information.
##Dimensions: 626 samples x 118 variables.

#Step 3: Import the TGFB and ALTEJ genelists. 
TGFBgeneset <- read.table("Input/TGFBUPgeneset.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
ALTEJgeneset <- read.table("Input/ALTEJgeneset.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

#Step 4: Turn the genelistS into lists and create a list containing both of them.
TGFBlist <- as.list(TGFBgeneset[,1:1]) 
ALTEJlist <- as.list(ALTEJgeneset[,1:1]) 
Bothgenelists <- list(TGFBUPgeneset=TGFBlist, ALTEJgeneset=ALTEJlist)



# PART 2: CALCULATE THE TGFB AND ALTEJ SSGSEA SCORES OF EACH SAMPLE -----------------------------------------

#Step 1: Turn the datafrane with gene expression into a matrix.
rownames(Xena_lusc_genes) <- Xena_lusc_genes$sample
Xena_lusc_genes$sample <- NULL
Xena_lusc_genes <- as.matrix(Xena_lusc_genes)

#Step 2: Convert gene expression values into Z-scores ((x - mean(Xcolumn)) / sd(Xcolumn)).
genesT <- t(Xena_lusc_genes) #Transpose the matrix.
#Calculate the Z-scores. IMPORTANT: To calculate the Z-score of each gene, each gene has to be a column.
genesT<-scale(genesT, center = TRUE, scale = TRUE)

#Step 3: Calculate the ssGSEA scores of TGFB and ALTEJ signatures in each sample.
genes<-t(genesT)
ssgsea <- gsva(genes, Bothgenelists, method=c("ssgsea"))

#Step 4: Traspose the matrix with the ssGSEA results and turn it into a dataframe.
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)



# PART 3: PUT ALL THE INFORMATION THAT WE WILL NEED FOR FURTHER ANALYSES IN ONE SINGLE DATAFRAME -----------------------------------

#Step 1: Merge the "ssgsea" and "Xena_lusc_clinical" dataframes.
ssgsea$ID1 <- rownames(ssgsea) #Create a variable that are sample IDs.
ssgsea$ID2 <- substr(ssgsea$ID1, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(ssgsea, Xena_lusc_clinical, 
             by.x = "ID2", by.y = "sampleID",
             all.x=TRUE, all.y=FALSE) #553 patients with gene expression, all of them with clinical data.

#Step 2: Add normalized (z-score transformed) TGFB1 expression to the "ALL" dataframe.
TGFB1 <- as.data.frame(genesT[,c("TGFB1")])
TGFB1$TGFB1 <- TGFB1$`genesT[, c("TGFB1")]`
TGFB1$`genesT[, c("TGFB1")]` <- NULL
TGFB1$ID1 <- rownames(TGFB1) 
TGFB1$ID2 <- substr(TGFB1$ID1, 0, 15) #Keep only the first 15 characters of sample IDs.
ALL <- merge(TGFB1, ALL, 
             by.x = "ID2", by.y = "ID2",
             all.x=FALSE, all.y=TRUE) #553p + 553p --> 553 patients. 

#Step 3: Add "lusc_tcga_pan_can_atlas_2018_clinical_data" to the "ALL" dataframe.
ALL <- merge(ALL, lusc_tcga_pan_can_atlas_2018_clinical_data,
             by.x = "ID2", by.y = "Sample.ID",
             all.x=TRUE, all.y=FALSE) #553p + 487p --> 553 patients. 484 patients were in common.
rownames(ALL) <- ALL$ID2



# PART 4: ELIMINATE DUPLICATED, RECURRENT AND NORMAL TISSUE SAMPLES -----------------------------------

#Step 1: Remove duplicated patients. 
Duplicated <- ALL[duplicated(ALL$`_PATIENT`), ] #51 duplicated patients. All 51 samples are from normal tissue.
ALL <- distinct(ALL, ALL$`_PATIENT`, .keep_all=TRUE) #553->502p. 

#Step 2: Create a variable that shows the sample type. 
ALL$Sample.type.2 <- substr(ALL$ID2, 13, 15) #Keep only sample type numbers.
table(ALL$Sample.type.2) #502 primary; 0 normal; 0 recurrent.



# PART 5: CALCULATE THE βALT SCORE OF EACH SAMPLE -----------------------------------

#Step 1: Create a new variable that is be the Balt score.
ALL$balt <- sqrt((max(ALL$ALTEJgeneset)-ALL$ALTEJgeneset)^2+
                 (min(ALL$TGFBUPgeneset)-ALL$TGFBUPgeneset)^2) - 
            sqrt((min(ALL$ALTEJgeneset)-ALL$ALTEJgeneset)^2+
                 (max(ALL$TGFBUPgeneset)-ALL$TGFBUPgeneset)^2)
ALL$balt <- ALL$balt * -1 

#Step 2: Create a new variable that represents the tertile that each patient belongs to, according to the value of the BAlt score.
ALL$tertile <- cut(ALL$balt, 
                   quantile(ALL$balt, c(0, 1/3, 2/3, 1)), na.rm=TRUE, include.lowest=TRUE, 
                   labels = c("low", "middle", "high"))
ALL$Tertile <- ifelse(ALL$tertile == "low", "high TGFB and low ALTEJ", 
                      ifelse(ALL$tertile == "high", "low TGFB and high ALTEJ", NA))
ALL$Tertile <- as.character(ALL$Tertile)



# PART 6: ANALYSIS OF THE CORRELATION BETWEEN TGFβ AND ALTEJ SSGSEA SCORES -----------------------

#Step 1: Create a scatterplot of TGFB versus ALTEJ ssGSEA scores, coloring the BAlt tertiles.
ggplot(ALL, aes(x=ALL$TGFBUPgeneset, y=ALL$ALTEJgeneset, col=ALL$Tertile)) + geom_point() +
  scale_color_manual(values=c("orange", "blue"), 
                     labels = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
                     na.translate=TRUE, na.value="grey") + 
  labs(x = "TGFβ ssGSEA score", y = "Alt-EJ ssGSEA score", color = "βAlt score groups") + 
  theme(legend.title=element_text(size=8)) 

#Step 2: Calculate Pearson correlation coefficient.
cor.test(ALL$TGFBUPgeneset, ALL$ALTEJgeneset, method = "pearson") #PCC=-0.43



# PART 7: ANALYSIS OF THE CORRELATION OF "TGFβ1 VS TGFβ SIGNATURE" AND "AGE VS βALT SCORE"  -----------------------

#Step 1: Calculate Pearson correlation coefficient of TGFB1 vs TGFB ssGSEA.
cor.test(ALL$TGFBUPgeneset, ALL$TGFB1, method = "pearson") #PCC=0.51

#Step 2: Create a scatterplot of TGFB1 versus TGFB ssGSEA score.
ggplot(ALL, aes(x=ALL$TGFBUPgeneset, y=ALL$TGFB1)) + geom_point(col="grey30") +
  labs(x = "TGFβ ssGSEA score", y = "TGFβ1") + 
  theme(legend.title=element_text(size=8)) 

#Step 3: Calculate Pearson correlation coefficient of age vs BAlt score.
cor.test(ALL$balt, ALL$Diagnosis.Age, method = "pearson") #PCC=-0.14

#Step 4: Create a scatterplot of Age vs BAlt score.
ggplot(ALL, aes(x=ALL$balt, y=ALL$Diagnosis.Age)) + geom_point(col="grey30") +
  labs(x = "βAlt score", y = "Age (years)") + 
  theme(legend.title=element_text(size=8))  +
  ggtitle("TCGA-LUSC")



# PART 8: ANALYSIS OF THE CORRELATION BETWEEN THE βALT SCORE AND GENOMIC ALTERATIONS -----------------------------

#Step 1: Plot boxplots showing the mutational load of the BAlt score top and bottom tertiles.
boxplot(ALL$Fraction.Genome.Altered ~ ALL$Tertile, 
        xlab="βAlt score groups", ylab="Fraction of the genome altered", 
        names = c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"), 
        col=c("orange", "blue"), 
        outpch=NA, cex.axis=0.8) + 
  stripchart(ALL$Fraction.Genome.Altered ~ ALL$Tertile, 
             vertical=TRUE, add=TRUE, method="jitter", col=c("dark orange", "blue"), 
             pch=21, bg="bisque") #Add dots to it representing each sample.

#Step 2: Do a U Mann Whitney Test to compare Genomic alterations between the top and bottom BAlt tertiles.
wilcox.test(ALL$Fraction.Genome.Altered ~ ALL$Tertile, data=ALL) 



# PART 9: IDENTIFICATION OF PATIENTS WHO PROBABLY RECEIVED RT AND/OR GENOTOXIC CHT -----------------------------

#Step 1: Remove stage I patients as they may had only been treated with surgery.
table(ALL$pathologic_stage) #3+90+151 are stage 1.
ALL <- ALL[-which(ALL$pathologic_stage == "Stage I" | ALL$pathologic_stage == "Stage IA" | ALL$pathologic_stage == "Stage IB"), ] #502->258 patients.

#Step 2: Eliminate patients in whom it is explicitly indicated that no RT nor ChT was given. 
frequencies <- table(ALL$Radiation.Therapy, ALL$radiation_therapy, ALL$additional_radiation_therapy, ALL$additional_pharmaceutical_therapy)
frequencies<-melt(frequencies)
names(frequencies) <- c("bioportal_RT", "Xena_RT", "Xena_additional_RT", "Xena_additional_pharma", "Freq") #27 did not receive neither ChT nor RT, so will be eliminated. 
ALL <- ALL[-which(ALL$radiation_therapy == "NO" & ALL$Radiation.Therapy == "No" & 
                  ALL$additional_pharmaceutical_therapy == "NO" & ALL$additional_radiation_therapy == "NO"), ] #258->231p.



# PART 10: OVERALL SURVIVAL AND PROGRESSION PREE SURVIVAL CURVES OF THE βALT SCORE GROUPS -----------------------------------------------

#Step 1: Edit the OS variables so that they are usable. These variables come from the cbioportal dataset.
table(ALL$Overall.Survival.Status) 
ALL$OSstatus <- ifelse(ALL$Overall.Survival.Status == "DECEASED", 1, 
                       ifelse(ALL$Overall.Survival.Status == "LIVING", 0, NA)) 
class(ALL$OSstatus) #[1] "numeric"
class(ALL$Overall.Survival..Months.) #[1] "numeric"

#Step 2: Edit the PFS variables so that they are usable. These variables come from the cbioportal dataset.
table(ALL$Progression.Free.Status) 
ALL$PFSstatus <- ifelse(ALL$Progression.Free.Status == "PROGRESSION", 1, 
                        ifelse(ALL$Progression.Free.Status == "CENSORED", 0, NA))
class(ALL$PFSstatus) #[1] "numeric"
class(ALL$Progress.Free.Survival..Months.) #[1] "numeric"

#Step 3: Plot and compare the OS curves between Balt score top and bottom tertiles.
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
my.fit<-survfit(my.surv.object~ALL$Tertile)
my.fit
ggsurvplot(my.fit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "blue"))

#Step 4: Plot and compare the PFS curves between Balt score top and bottom tertiles.
my.DFsurv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
my.DFfit<-survfit(my.DFsurv.object~ALL$Tertile)
my.DFfit
ggsurvplot(my.DFfit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Progression Free Fraction", xlim=c(0, 60),   
           legend.labs=c("high TGFβ & low Alt-EJ", "low TGFβ & high Alt-EJ"),
           legend.title=" ", title="PFS in tertiles 1 versus 3 of patients according to their βAlt score", 
           font.main=13, palette=c("orange", "blue"))



#PART 11: CALCULATE THE HAZARD RATIOS OF THE OS AND PFS CURVES FROM PART 10 -------------------------------------

#Step 1: In the variable tertile, make "high TGFB and low ALTEJ" the reference level so that the Hazard Ratios are 0.x instead of 1.x. 
table(ALL$Tertile)
ALL$Tertile <- factor(ALL$Tertile, levels = c("high TGFB and low ALTEJ", "low TGFB and high ALTEJ"))
table(ALL$Tertile)

#Step 2: Calculate the OS hazard ratio of BAlt tertile 1 versus 3.
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
cox<-coxph(my.surv.object ~  ALL$Tertile)
summary(cox) #HR=0.60.

#Step 3: Calculate the PFS hazard ratio of BAlt tertile 1 versus 3.
my.surv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
cox<-coxph(my.surv.object ~  ALL$Tertile)
summary(cox) #HR=0.51.



#PART 12: MULTIVARIATE COX OF THE ASSOCIATION OF THE βALT SCORE WITH OS AND PFS, ADJUSTED FOR AGE AND STAGE -------------------------------------

#Step 1: Group Figo stages into 3 groups (2, 3 and 4) 
table(ALL$pathologic_stage)
ALL$stage <- ifelse(ALL$pathologic_stage == "Stage II", "2",
                   ifelse(ALL$pathologic_stage == "Stage IIA", "2",
                   ifelse(ALL$pathologic_stage == "Stage IIB", "2",
                   ifelse(ALL$pathologic_stage == "Stage III", "3",
                   ifelse(ALL$pathologic_stage == "Stage IIIA", "3",
                   ifelse(ALL$pathologic_stage == "Stage IIIB", "3",
                   ifelse(ALL$pathologic_stage == "Stage IV", "4", NA)))))))
table(ALL$stage)

#Step 2: Calculate BAlt score(continuous variable) OS multivariate Cox adjusted for age and stage. 
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
cox<-coxph(my.surv.object ~  ALL$balt + ALL$Diagnosis.Age + ALL$stage)
summary(cox) #p=0.035.

#Step 3: Calculate BAlt score (continuous variable) PFS multivariate Cox adjusted for age and stage. 
my.surv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
cox<-coxph(my.surv.object ~  ALL$balt + ALL$Diagnosis.Age + ALL$stage)
summary(cox) #p=0.01.



# PART 13: TGFβ1 ASSOCIATION WITH OVERALL SURVIVAL AND PROGRESSION PREE SURVIVAL -----------------------------------------------

#Step 1: Split TGFB1 into tertiles.
ALL$TGFB1tertile <- cut(ALL$TGFB1, 
                        quantile(ALL$TGFB1, c(0, 1/3, 2/3, 1), na.rm=TRUE, include.lowest=TRUE), 
                        labels = c("low", "middle", "high"))
ALL$TGFB1Tertile <- ifelse(ALL$TGFB1tertile == "low", "low TGFB1", 
                           ifelse(ALL$TGFB1tertile == "high", "high TGFB1", NA))
ALL$TGFB1Tertile <- as.character(ALL$TGFB1Tertile)

#Step 2: Plot and compare the OS curves between TGFB1 top and bottom tertiles.
my.surv.object <- Surv(time=ALL$Overall.Survival..Months., event=ALL$OSstatus)
my.fit<-survfit(my.surv.object~ALL$TGFB1Tertile)
my.fit
ggsurvplot(my.fit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),  
           legend.labs=c("high TGFβ1", "low TGFβ1"),
           legend.title=" ", title="OS in tertiles 1 versus 3 of patients according to their TGFβ1 expression", 
           font.main=13, palette=c("orange", "blue"))

#Step 3: Plot and compare the PFS curves between TGFB1 top and bottom tertiles.
my.DFsurv.object <- Surv(time=ALL$Progress.Free.Survival..Months., event=ALL$PFSstatus)
my.DFfit<-survfit(my.DFsurv.object~ALL$TGFB1Tertile)
my.DFfit
ggsurvplot(my.DFfit, data=ALL, pval = TRUE, conf.int = TRUE, break.time.by = 12,
           xlab="Time (months)", ylab="Progression Free Fraction", xlim=c(0, 60),   
           legend.labs=c("high TGFβ1", "low TGFβ1"),
           legend.title=" ", title="PFS in tertiles 1 versus 3 of patients according to their TGFβ1 expression", 
           font.main=13, palette=c("orange", "blue")) 




