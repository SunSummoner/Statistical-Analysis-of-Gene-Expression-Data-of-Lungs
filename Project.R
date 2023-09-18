library(GEOquery)
gds <- getGEO("GSM254626")
#for normal lungs
gsm <- getGEO(filename=system.file("extdata/GSM254626.txt.gz",package="GEOquery"))
gds <- getGEO("GSM254626")
head(Meta(gds))

Table(gds)[1:5,]
Columns(gsm)
Columns(gds)[,1:3]
gse <- getGEO("GSE10072", GSEMatrix = FALSE)
names(GSMList(gse))
GSMList(gse)[[1]]

gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL96'},GSMList(gse))
length(gsmlist)
probesets <- Table(GPLList(gse)[[1]])$ID
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
    + {tab <- Table(x)
    + mymatch <- match(probesets,tab$ID_REF)
    + return(tab$VALUE[mymatch])
    + }))

data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]

install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

my_id <- "GSE10072"
gse <- getGEO(my_id)

gse <- gse[[1]]
gse

pData(gse)
fData(gse) ## print the gene annotation
exprs(gse)##The expression data

summary(exprs(gse))


par("mar")
par(mar=c(1,1,1,1))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

library(dplyr)

sampleInfo <- pData(gse)
sampleInfo

sampleInfo

library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 
rownames(sampleInfo)
colnames(corMatrix)
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)   
library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()
