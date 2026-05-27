############################################################
# 02_diversity_and_composition.R
# Rarefaction, alpha diversity, beta diversity, and
# taxonomic composition analyses
############################################################

library(phyloseq)
library(picante)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# -----------------------------
# Load processed phyloseq object
# -----------------------------

ps <- readRDS("phyloseq_Ghana_V4_iptp.rds")

# -----------------------------
# Rarefaction
# -----------------------------

set.seed(1)
ps_rare <- rarefy_even_depth(
  ps,
  sample.size = 12000,
  rngseed = 1,
  replace = FALSE
)

saveRDS(ps_rare, "results/phyloseq_Ghana_V4_iptp_rarefied_12000.rds")

# Sample overview

sample_overview <- sample_data(ps) %>%
  data.frame() %>%
  count(sample_type)

paired_samples <- sample_data(ps) %>%
  data.frame() %>%
  group_by(individual) %>%
  summarise(
    sample_types = paste(sort(unique(sample_type)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(grepl("stool", sample_types) & grepl("vaginal", sample_types))

write.csv(sample_overview, "results/sample_overview.csv", row.names = FALSE)
write.csv(paired_samples, "results/paired_samples.csv", row.names = FALSE)

# Alpha diversity

alpha_div <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)

otu_mat <- as(otu_table(ps_rare), "matrix")

if (!taxa_are_rows(ps_rare)) {
  otu_mat <- t(otu_mat)
}

faith_pd <- pd(
  samp = t(otu_mat),
  tree = phy_tree(ps_rare),
  include.root = TRUE
)

alpha_div$sample_id <- rownames(alpha_div)
faith_pd$sample_id <- rownames(faith_pd)

metadata <- sample_data(ps_rare) %>%
  data.frame() %>%
  tibble::rownames_to_column("sample_id")

alpha_div_meta <- alpha_div %>%
  left_join(faith_pd, by = "sample_id") %>%
  left_join(metadata, by = "sample_id")

write.csv(alpha_div_meta, "results/alpha_diversity_metrics.csv", row.names = FALSE)

# Beta diversity

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
saveRDS(ps_rel, "results/phyloseq_Ghana_V4_iptp_relative_abundance.rds")

# Unweighted UniFrac at ASV level

unifrac_dist <- phyloseq::distance(ps_rel, method = "unifrac")
unifrac_pcoa <- ordinate(ps_rel, method = "PCoA", distance = unifrac_dist)

# Sample type comparison

metadata_rel <- sample_data(ps_rel) %>% data.frame()

adonis_sample_type <- adonis2(
  unifrac_dist ~ sample_type,
  data = metadata_rel,
  permutations = 999
)

beta_disp_sample_type <- betadisper(
  unifrac_dist,
  metadata_rel$sample_type
)

beta_disp_sample_type_test <- anova(beta_disp_sample_type)

write.csv(
  as.data.frame(adonis_sample_type),
  "results/permanova_sample_type_unifrac.csv"
)

write.csv(
  as.data.frame(beta_disp_sample_type_test),
  "results/betadisper_sample_type_unifrac.csv"
)

p_pcoa_sample_type <- plot_ordination(ps_rel, unifrac_pcoa, color = "sample_type") +
  geom_point(size = 3) +
  theme_minimal() +
  labs(color = "Sample type")

ggsave(
  "figures/pcoa_unifrac_sample_type.pdf",
  p_pcoa_sample_type,
  width = 5,
  height = 4
)

# IPTp-SP comparisons separately for stool and vaginal samples

run_iptp_beta <- function(ps_object, body_site, output_prefix) {
  
  ps_sub <- subset_samples(
    ps_object,
    sample_type == body_site & !is.na(sp.category3) & sp.category3 != "NA"
  )
  
  dist <- phyloseq::distance(ps_sub, method = "unifrac")
  ord <- ordinate(ps_sub, method = "PCoA", distance = dist)
  meta <- sample_data(ps_sub) %>% data.frame()
  
  adonis_res <- adonis2(
    dist ~ sp.category3,
    data = meta,
    permutations = 999
  )
  
  bd <- betadisper(dist, meta$sp.category3)
  bd_test <- anova(bd)
  
  write.csv(
    as.data.frame(adonis_res),
    paste0("results/permanova_", output_prefix, "_iptp_unifrac.csv")
  )
  
  write.csv(
    as.data.frame(bd_test),
    paste0("results/betadisper_", output_prefix, "_iptp_unifrac.csv")
  )
  
  p_pcoa <- plot_ordination(ps_sub, ord, color = "sp.category3") +
    geom_point(size = 3) +
    stat_ellipse() +
    theme_minimal() +
    labs(color = "IPTp-SP uptake")
  
  ggsave(
    paste0("figures/pcoa_unifrac_", output_prefix, "_iptp.pdf"),
    p_pcoa,
    width = 5,
    height = 4
  )
}

run_iptp_beta(ps_rel, body_site = "stool", output_prefix = "stool")
run_iptp_beta(ps_rel, body_site = "vaginal", output_prefix = "vaginal")

# -----------------------------
# Taxonomic composition
# -----------------------------

plot_taxonomic_composition <- function(ps_object, tax_level, top_n = 15, output_file) {
  
  ps_glom <- tax_glom(ps_object, taxrank = tax_level)
  
  tax_df <- as.data.frame(tax_table(ps_glom)) %>%
    tibble::rownames_to_column("taxon_id")
  
  otu_df <- as.data.frame(otu_table(ps_glom))
  
  if (taxa_are_rows(ps_glom)) {
    otu_df <- otu_df %>%
      tibble::rownames_to_column("taxon_id") %>%
      pivot_longer(
        cols = -taxon_id,
        names_to = "sample_id",
        values_to = "abundance"
      )
  } else {
    otu_df <- otu_df %>%
      tibble::rownames_to_column("sample_id") %>%
      pivot_longer(
        cols = -sample_id,
        names_to = "taxon_id",
        values_to = "abundance"
      )
  }
  
  meta_df <- sample_data(ps_glom) %>%
    data.frame() %>%
    tibble::rownames_to_column("sample_id")
  
  composition_df <- otu_df %>%
    left_join(tax_df, by = "taxon_id") %>%
    left_join(meta_df, by = "sample_id") %>%
    mutate(
      taxon = ifelse(is.na(.data[[tax_level]]) | .data[[tax_level]] == "", "Unclassified", .data[[tax_level]])
    )
  
  top_taxa <- composition_df %>%
    group_by(taxon) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop") %>%
    slice_max(total_abundance, n = top_n) %>%
    pull(taxon)
  
  composition_df <- composition_df %>%
    mutate(taxon = ifelse(taxon %in% top_taxa, taxon, "Other")) %>%
    group_by(sample_id, sample_type, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  write.csv(
    composition_df,
    paste0("results/composition_", tolower(tax_level), ".csv"),
    row.names = FALSE
  )
  
  p <- ggplot(composition_df, aes(x = sample_id, y = abundance, fill = taxon)) +
    geom_col(width = 0.9) +
    facet_wrap(~ sample_type, scales = "free_x") +
    theme_minimal() +
    labs(
      x = "Sample",
      y = "Relative abundance",
      fill = tax_level
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ggsave(output_file, p, width = 12, height = 6)
}

plot_taxonomic_composition(
  ps_rel,
  tax_level = "Genus",
  top_n = 15,
  output_file = "figures/genus_level_composition.pdf"
)

plot_taxonomic_composition(
  ps_rel,
  tax_level = "Family",
  top_n = 15,
  output_file = "figures/family_level_composition.pdf"
)

plot_taxonomic_composition(
  ps_rel,
  tax_level = "Phylum",
  top_n = 15,
  output_file = "figures/phylum_level_composition.pdf"
)
