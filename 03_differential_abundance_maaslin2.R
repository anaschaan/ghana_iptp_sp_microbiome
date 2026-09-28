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

# Prep abundance table
# Stool samples with an IPTp-SP category. The prevalence filter is applied after any sample exclusion (so it matches each model's own sample set).
prepare_maaslin_input <- function(ps_object, tax_level = NULL,
                                  exclude_partial = FALSE) {
  
  ps_sub <- ps_object %>%
    subset_samples(sample_type == "stool") %>%
    subset_samples(!is.na(sp.category3))
  
  # Sensitivity analysis: drop the single participant with one dose (ipt == 1)
  if (exclude_partial) {
    ps_sub <- subset_samples(ps_sub, !is.na(ipt) & ipt != 1)
  }
  
  if (!is.null(tax_level)) {
    ps_sub <- tax_glom(ps_sub, taxrank = tax_level)
  }
  
  ps_sub <- filter_taxa(ps_sub, function(x) sum(x > 0) > 0.1 * length(x), TRUE)
  
  abundance <- as.data.frame(otu_table(ps_sub))
  if (taxa_are_rows(ps_sub)) {
    abundance <- as.data.frame(t(abundance))
  }
  
  metadata <- data.frame(sample_data(ps_sub))
  
  # Setting reference levels
  metadata$sp.category3 <- factor(
    metadata$sp.category3,
    levels = c("Less than three", "Three or more")
  )
  metadata$parity_collapsed <- factor(
    metadata$parity_collapsed,
    levels = c("Nulliparous", "Parous")
  )
  metadata$parity <- factor(
    metadata$parity,
    levels = c("Nulliparous", "Primiparous", "Multiparous")
  )
  
  list(abundance = abundance, metadata = metadata)
}

# Run MaAsLin2

run_maaslin <- function(ps_object, output_name,
                        fixed_effects = "sp.category3",
                        reference = "sp.category3,Less than three",
                        tax_level = NULL,
                        exclude_partial = FALSE) {
  
  input <- prepare_maaslin_input(ps_object, tax_level, exclude_partial)
  
  # Keep only samples with complete data for the model covariates
  keep <- complete.cases(input$metadata[, fixed_effects, drop = FALSE])
  abundance <- input$abundance[keep, ]
  metadata <- input$metadata[keep, ]
  
  Maaslin2(
    input_data = abundance,
    input_metadata = metadata,
    output = file.path("results/maaslin2", output_name),
    fixed_effects = fixed_effects,
    reference = reference,
    normalization = "CLR",
    transform = "NONE",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE
  )
}

# ASV-level differential abundance, gut

# Unadjusted
run_maaslin(ps, "asv_iptp_unadjusted")

# Adjusted for parity, 2 levels
run_maaslin(
  ps, "asv_iptp_adjusted_parity_2lvl",
  fixed_effects = c("sp.category3", "parity_collapsed"),
  reference = c("sp.category3,Less than three",
                "parity_collapsed,Nulliparous")
)

# Adjusted for parity, 3 levels
run_maaslin(
  ps, "asv_iptp_adjusted_parity_3lvl",
  fixed_effects = c("sp.category3", "parity"),
  reference = c("sp.category3,Less than three",
                "parity,Nulliparous")
)

# Sensitivity analysis: unadjusted, excluding the participant with one dose
run_maaslin(ps, "asv_iptp_unadjusted_excluding_partial_dose",
            exclude_partial = TRUE)
