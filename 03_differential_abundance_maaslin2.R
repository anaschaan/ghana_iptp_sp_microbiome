############################################################
# 03_differential_abundance_maaslin2.R
# Differential abundance analysis with MaAsLin2
############################################################

library(phyloseq)
library(dplyr)
library(tidyr)
library(Maaslin2)

# Load phyloseq object

ps <- readRDS("phyloseq_Ghana_V4_iptp.rds")
# Prepare metadata

metadata <- sample_data(ps) %>%
  data.frame() %>%
  tibble::rownames_to_column("sample_id")

# Keep only samples with IPTp-SP category and parity information
metadata <- metadata %>%
  filter(
    !is.na(sp.category3),
    sp.category3 != "NA",
    !is.na(n.birth.category),
    n.birth.category != "NA"
  )

# Set reference levels
metadata$sp.category3 <- factor(
  metadata$sp.category3,
  levels = c("Three or more", "Less than three")
)

metadata$n.birth.category <- factor(metadata$n.birth.category)

#Prep abundance table
prepare_maaslin_input <- function(ps_object, tax_level = NULL, sample_ids) {
  
  ps_sub <- prune_samples(sample_ids, ps_object)
  
  if (!is.null(tax_level)) {
    ps_sub <- tax_glom(ps_sub, taxrank = tax_level)
  }
  
  abundance <- as.data.frame(otu_table(ps_sub))
  
  if (taxa_are_rows(ps_sub)) {
    abundance <- t(abundance)
  }
  
  abundance <- as.data.frame(abundance)
  
  return(abundance)
}

# Run MaAsLin2

run_maaslin <- function(ps_object, metadata, output_name, tax_level = NULL) {
  
  abundance <- prepare_maaslin_input(
    ps_object = ps_object,
    tax_level = tax_level,
    sample_ids = metadata$sample_id
  )
  
  # Match sample order between abundance and metadata
  metadata_matched <- metadata %>%
    filter(sample_id %in% rownames(abundance)) %>%
    arrange(match(sample_id, rownames(abundance)))
  
  abundance <- abundance[metadata_matched$sample_id, ]
  
  fit_data <- metadata_matched %>%
    tibble::column_to_rownames("sample_id")
  
  Maaslin2(
    input_data = abundance,
    input_metadata = fit_data,
    output = file.path("results/maaslin2", output_name),
    fixed_effects = c("sp.category3", "n.birth.category"),
    reference = c("sp.category3,Three or more"),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE,
    min_prevalence = 0.1
  )
}

# Run ASV-level differential abundance

run_maaslin(
  ps_object = ps,
  metadata = metadata,
  output_name = "asv_iptp_adjusted_parity",
  tax_level = NULL
)

