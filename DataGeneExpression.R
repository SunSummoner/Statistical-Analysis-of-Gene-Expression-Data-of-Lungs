install.packages("tidyverse")
getwd()
library("GEOquery")
library(SummarizedExperiment)

list.files("GEO")

tempdir()
#Temporary directory for GEO files

gse <- getGEO(GEO = "GSE103512")

gse

class(gse)
length(gse)

gse103512 <- as(gse$GSE103512_series_matrix.txt.gz , "SummarizedExperiment")

class(gse103512)
gse103512

isS4(gse103512)
typeof(gse103512)
attributes(gse103512)[["class"]]

class(gse103512)
is(gse103512, "SummarizedExperiment")

getClass("SummarizedExperiment")

slotNames(gse103512)

metadata(gse103512)

class(metadata(gse103512))

names(metadata(gse103512))


metadata(gse103512)$formula <- exprs ~ cancer.type.ch1 + normal.ch1
#Adding design matrix formula for differential expression analysis

names(metadata(gse103512))

metadata(gse103512)[["experimentData"]]

class(metadata(gse103512)[["experimentData"]])

isS4(metadata(gse103512)[["experimentData"]])

    

miame <- metadata(gse103512)[["experimentData"]]


miame

expinfo(miame)

#expinfo(miame)[c("Title", "url")]

abstract(miame)

pubMedIds(miame)

names(otherInfo(miame))

otherInfo(miame)[c("relation", "overall_design")]

assays(gse103512)

class(assay(gse103512, 'exprs'))

dim(assay(gse103512, 'exprs'))


tibble::as_tibble(assay(gse103512, "exprs"))



my_id <- "GSE10072"

gse10072 <- getGEO(my_id)

my_id <- "GSE10072"

gse <- gse[[1]]

gse

pData(gse10072)

fData(gse10072)

exprs(gse10072)

gse10072 <- as(gse10072$GSE10072_series_matrix.txt.gz)

gse10072 <- as(gse10072$GSE10072_series_matrix.txt.gz, "ExpressionSet")

summary(exprs(gse10072))




