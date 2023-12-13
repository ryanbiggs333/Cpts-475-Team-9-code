library(tidyverse)
library(dplyr)

# four columns: label, Donor_ID, tumor_gene_uni, tissue
# suggestions and notes from mentor:
#  -Tissue is the cancer type
#  -tumor_gene_uni: which gene has caused the cancer in that tissue
#  -Look for patterns, for example, PTPRT could be exclusive to breast

gene_tissue_data <- read.csv("icgc_468_gene_tissue.csv", header = TRUE)

gene_tissue_raw_data <- read.csv("icgc_data_468.csv", header = TRUE)

gene_sequence <- read.csv("donor_art_seq_468.csv", header = TRUE)

# the mentor suggested to go to the Project_Code column and 
#extract the first word using the - as the separator
#we also extract the second word for potential use in the future

gene_tissue_raw_data <- gene_tissue_raw_data %>%
    separate(Project_Code, into = c("Project_Type", "Cancer_Type"))


donor_sequence_data <- read.csv("donor_art_seq_468.csv", header = TRUE)

# put that word into a new column in the donor_art_seq_468.csv file, 
# using the Donor_id column from the second csv file and to match it to 
# the sample id column from the donor_art_seq_468.csv file.
# repeat same process for second word
donor_sequence_data <- donor_sequence_data %>%
  left_join(select(gene_tissue_raw_data, Donor_ID, Project_Type), by = c("sample.id" = "Donor_ID"))
donor_sequence_data <- donor_sequence_data %>%
  left_join(select(gene_tissue_raw_data, Donor_ID, Cancer_Type), by = c("sample.id" = "Donor_ID"))

view(donor_sequence_data)
write.csv(donor_sequence_data, file = "v2.csv", row.names = FALSE)


selected <- select(gene_tissue_data, Donor_ID, tumor_gene_uni)

# Use left_join to merge d1 with the selected columns of d2
# Specify the join condition by matching 'sample_id' from d1 with 'donor_id' from d2
merged_data <- gene_sequence %>%
  left_join(selected, by = c("sample id" = "Donor_ID"))

# -------------------------------------------------------------------------------------------------------------------------------------------
# *Let's investigate how often a gene is associated with a tissue*
gene_frequency <- gene_tissue_data %>% 
  separate_rows(tumor_gene_uni, sep = ", ") %>% 
  count(tumor_gene_uni, tissue)
print(gene_frequency)

# this shows which genes are predominant in causing certain cancers
gene_frequency %>%
  arrange(desc(n)) %>%
  top_n(10)
# Notes:
#   -TP53 gene is frequently associated with multiple cancer types:
#    (Pancreatic, Esophageal, and Liver)
#   -KRAS, GRIN2A, PTPRT, CTNNB1, BRAF appear to have convincing associations
#    with specific cancer types
#
# Potential next step: Investigate each of these genes individually, find a 
# relation with cancer biology, maybe research on the web.

# now let's try investigating which genes are predominant in causing certain cancers
# for a specific tissues
gene_frequency %>%
  filter(tissue == "Panc") %>%
  arrange(desc(n))
# Notes:
#   -KRAS gene is the most prevalent
#   -TP53 is the next gene most showed up and next, PTPRT(26) and CDKN2A(20) show up almost equally
#   -"A hallmark of pancreatic cancer is the prevalence of oncogenic mutation in the KRAS gene"
#     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8380752/#:~:text=A%20hallmark%20of%20pancreatic%20cancer,major%20target%20for%20therapeutic%20intervention.

gene_frequency %>%
  filter(tissue == "Skin") %>%
  arrange(desc(n))
# Notes:
#   -GRIN2A and PTPRT gene is the most prevalent
#   -Followed by BRAF, with EGFR and NRAS coming in next with almost the same frequency
#   -"...protein tyrosine phosphatase receptor type T (PTPRT), which is a well‐known phosphatase and mutates frequently in melanoma..."
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8817076/

gene_frequency %>%
  filter(tissue == "Eso") %>%
  arrange(desc(n))
# Notes:
#   -TP53 and PTPRT gene is the most prevalent
#   -"In esophagogastric cancer, the combination mutation frequency of PTPRT and TP53 in metastatic 
#    cancer is significantly higher than that in primary cancer"
#    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9250438/

gene_frequency %>%
  filter(tissue == "Liver") %>%
  arrange(desc(n))
# Notes:
#   -CTNNB1 and TP53 gene is the most prevalent
#   -couldn't find an article showing these two genes have a connection with liver cancer :(

# General notes about these four tissues:
#   -TP53 seems to be connected to many tissues
#   -certain genes seem to be more tissue-specific, KRAS for pan and CTNNB1 for liver

# Moving along, let's investigate rare occurrences, genes that appear once across all tissues
gene_frequency %>%
  filter(n == 1) %>%
  print(n=74)

# Notes:
#   - KEAP1 is unique to Myeloid and GRIN2A is unique to Bone

# Code below is for identifying genes that appear across multiple tissue types
gene_frequency %>%
  group_by(tumor_gene_uni) %>%
  filter(n() > 1) %>%
  arrange(desc(n())) %>%
  print(n=205)


