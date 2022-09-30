#installing required packages
BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("affy")
BiocManager::install("limma")
install.packages("ggplot2")
install.packages("pheatmap")

#loading the installed packages into R
library(GEOquery)
library(oligo)
library(affy)
library(limma)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#downloading the chosen GSE Supp files (Raw data)
getGEOSuppFiles('GSE50737')
#Depending on the network strength the download might fail a couple of times hence another alternative might be to directly download the supplementary (.tar) file from the GEO Accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

#setting the working directory where the file has been downloaded
setwd("C:/Users/Sasmitha/Documents/BL5631_Sasmitha Deepak Singh_CA1/GSE50737")

#extracting the .cel files
untar("GSE50737_RAW.tar", exdir = 'Data/')

#setting the working directory where the .cel files have been extracted
setwd("C:/Users/Sasmitha/Documents/BL5631_Sasmitha Deepak Singh_CA1/GSE50737/Data")

#listing and reading the .cel files from the set working directory
GSE50737_celdata <- read.celfiles(list.celfiles())

image(GSE50737_celdata[,1])

#Background correction, normalizing and calculation of expression
GSE50737_Eset <- oligo::rma(GSE50737_celdata)

#checking for histogram evidence of rma
hist(GSE50737_celdata, lwd=2, xlab='log intensity', which = 'all', main = "CEL file densities before RMA")
hist(GSE50737_Eset, main = "CEL file densities after RMA")

#cross-checking metadata and phenodata with the series on GEO to ensure that the order matches
GEO_gse50737 <- getGEO("GSE50737")
GEO_gse50737Eset <- GEO_gse50737[[1]]
pData(GEO_gse50737Eset)
pData(GSE50737_Eset)

#Naming the study groups after cross-checking
StudyGroups <- c("BenzenePoisoning", "BenzenePoisoning", "BenzenePoisoning", "BenzenePoisoning", "BenzeneExposed", "BenzeneExposed", "BenzeneExposed", "Control", "Control", "Control")

#Generating a model matrix with specified design
Design <- model.matrix(~0+factor(StudyGroups))
Design
colnames(Design) <- c("BenzeneExposed", "BenzenePoisoning", "Control")
Design

#Fitting the model to the design
fit <- lmFit(GSE50737_Eset, Design)
fitted.ebayes <- eBayes(fit)
options(digits=3) #Decimal setting
topTable(fitted.ebayes)

#In order to save a copy of the differential expression values
write.table(GOI, "Differential_Exp.txt", sep="\t") 

#Creating contrast matrix
contrast_matrix <- makeContrasts(BenzenePoisoning-Control, BenzeneExposed-Control, BenzenePoisoning-BenzeneExposed, levels=Design)

#Fitting the contrast matrix to make all pair-wise comparisons
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
results
summary(decideTests(fit2), lfc=1)

#pheatmap option1
GOI <- topTable(fit2, coef=1, adjust="BH")
GOI
EsetOfInterest <- GSE50737_Eset[rownames(GOI),]
pheatmap(exprs(EsetOfInterest))

#Volcano plot
GOI2 <- topTable(fit2,coef = 2, number =100, adjust.method = "none", p.value=0.05, lfc=0.5)
GOI2
volcanoplot(fit2, coef=2, main=sprintf("Top 100 features that pass our cutoffs",nrow(GOI2)), highlight=10) #As many as needed can be highlighted. Only 10 highlighted for clarity purposes.
points(GOI2[['logFC']],-log10(GOI2[['P.Value']]),col='red')

#ggplot for the ten samples in the experiment
ggplot(GOI, aes(GOI[['logFC']],-log10(GOI[['P.Value']]), colour = StudyGroups)) + 
  geom_point()


#Boxplot of the differential expression of our genes of interest (GOI) in each sample
par(mar=c(20,5,2,2))
boxplot(exprs(EsetOfInterest), main="Boxplots of selected differentially expressed genes",
        ylab="Expression",
        col="Blue",
        border="Black", las=2)
