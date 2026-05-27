############################################################
# 01_load_and_prepare_phyloseq.R
# Load ASV table, taxonomy, tree, and metadata into phyloseq
############################################################

# Load libraries
library(phyloseq)
library(dplyr)

# -----------------------------
# Define input files
# -----------------------------

asv_file <- "annotation/table.tsv"
taxonomy_file <- "taxonomy-formatted-final.tsv"
tree_file <- "annotation/seq_rooted_tree_unzip/data/tree.nwk"
metadata_file <- "metadata-ghana.csv"

# -----------------------------
# Read input data
# -----------------------------

asv_table <- read.table(
  asv_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

taxonomy_table <- read.table(
  taxonomy_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

tree <- read_tree(tree_file)

metadata <- read.csv(
  metadata_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# -----------------------------
# Create phyloseq object
# -----------------------------

ASV <- otu_table(as.matrix(asv_table), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(taxonomy_table))
META <- sample_data(metadata)

ps <- phyloseq(ASV, TAX, META, tree)

# -----------------------------
# Keep only Ghana IPTp samples
# -----------------------------

ps <- subset_samples(ps, batch == "ghana")

# -----------------------------
# Add IPTp-SP uptake category
# -----------------------------
# Less than three doses: ipt <= 2
# Three or more doses: ipt >= 3

sample_data(ps)$sp.category3 <- ifelse(
  is.na(sample_data(ps)$ipt),
  NA,
  ifelse(sample_data(ps)$ipt <= 2, "Less than three", "Three or more")
)

# -----------------------------
# Save processed object
# -----------------------------

saveRDS(ps, "phyloseq_Ghana_V4_iptp.rds")
