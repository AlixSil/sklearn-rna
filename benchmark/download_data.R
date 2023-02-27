library(TCGAbiolinks)
library(dplyr)

query.clinical <- GDCquery(
    project = "TCGA-ACC",
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
)
GDCdownload(query.clinical)
clinicaldat <- GDCprepare(query)

a = clinicaldat@colData[, c("patient", "days_to_last_follow_up", "days_to_death", "vital_status")]



query.exp <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
)
GDCdownload(query.exp)
expdat <- GDCprepare(
    query = query.exp,
    summarizedExperiment = FALSE,
)

f = paste("unstranded", rownames(a), sep="_")
expdat[1:.N, ..f]

#TCGA-OR-A5LD-01A-11R-A29S-07

