library(TCGAbiolinks)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
project =  args[1]


clinicaldat <- GDCquery_clinic(project, type="clinical", save.csv=FALSE)

clinicaldat.formatted = clinicaldat[, c("bcr_patient_barcode", "days_to_last_follow_up", "days_to_death", "vital_status")]


query.exp <- GDCquery(
    project = project, 
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

GDCdownload(query.exp)
expdat <- GDCprepare(
    query = query.exp,
    summarizedExperiment = FALSE,
)

#Keep only raw data

raw_unstranded_samples <- expdat %>% names() %>% 
    grep("^unstranded", ., value = TRUE)


samples_handling = tibble(column_name = raw_unstranded_samples)

samples_handling$patient <- samples_handling$column_name %>% 
    sub('unstranded_([^-]+-[^-]+-[^-]+).*', "\\1", .) 

samples_handling$duplicated <- duplicated(samples_handling$patient)

#TODO handle cases of duplication


samples_handling <- samples_handling %>% filter(patient %in% clinicaldat.formatted$bcr_patient_barcode)
clinicaldat.formatted <- clinicaldat.formatted %>% filter(bcr_patient_barcode %in% samples_handling$patient)

expdat.filter <- expdat %>% select(.,all_of(c("gene_id", "gene_name", "gene_type", samples_handling$column_name)))

renaming_vector <- samples_handling$column_name
names(renaming_vector) <- samples_handling$patient

expdat.filter <- expdat.filter %>% rename(., all_of(renaming_vector)) %>% filter(row_number() < n()-3)


output_dir = paste("Formatted_data/", project, sep="")
dir.create(output_dir, showWarnings=FALSE)

write.csv(expdat.filter, file = paste(output_dir, "raw_counts.csv", sep="/"), row.names = FALSE)
write.csv(clinicaldat.formatted, file = paste(output_dir, "clinical.csv", sep="/"), row.names = FALSE)