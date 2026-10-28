# ---- Epigenetic trajectory clustering (degPatterns) ----
# Clusters EpiAgeSig peaks into trajectories across 10-year age bins, then
# hierarchically clusters the resulting trajectory centroids into broader
# regulatory phases.

library(DESeq2)
library(DEGreport)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

# ---- Build DESeq2 dataset from EpiAgeSig peak counts ----

countData <- readRDS("/path/to/epiagesig_counts.rds")
colData   <- readRDS("/path/to/sample_metadata.rds")

colData$age_group <- cut(colData$age,
                         breaks = c(-Inf, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90),
                         labels = c("fetal", "0-9", "10-19", "20-29", "30-39",
                                    "40-49", "50-59", "60-69", "70-79", "80-89"),
                         right = FALSE)
colData$age_group <- as.factor(colData$age_group)

stopifnot(all(colnames(countData) == rownames(colData)))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 0 + age_group)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# variance-stabilised signal, used as input to degPatterns
rld <- vst(dds, blind = FALSE)
ma  <- assay(rld)

# ---- Cluster peaks into trajectories across age bins ----

patterns <- degPatterns(
  ma,
  metadata = as.data.frame(colData(dds)),
  minc = 50,
  summarize = "merge",
  time = "age_group",
  reduce = FALSE,
  nClusters = 10,
  skipDendrogram = FALSE
)

# ---- Hierarchically cluster trajectory centroids into regulatory phases ----

age_order <- c("fetal", "0-9", "10-19", "20-29", "30-39",
               "40-49", "50-59", "60-69", "70-79", "80-89")

df <- as.data.frame(patterns$normalized) %>%
  mutate(genes = as.character(genes), cluster = as.character(cluster),
         age_group = as.character(age_group), value = as.numeric(value))

trajectory_profiles <- df %>%
  mutate(age_group = factor(age_group, levels = age_order)) %>%
  group_by(cluster, age_group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

traj_mat <- trajectory_profiles %>%
  arrange(cluster, age_group) %>%
  pivot_wider(names_from = age_group, values_from = mean_value) %>%
  column_to_rownames("cluster") %>%
  as.matrix()
traj_mat <- traj_mat[complete.cases(traj_mat), ]

traj_mat_scaled <- t(scale(t(traj_mat)))
traj_dist   <- dist(traj_mat_scaled, method = "euclidean")
traj_hclust <- stats::hclust(traj_dist, method = "ward.D2")

# three broader regulatory phases
trajectory_groups <- cutree(traj_hclust, k = 3)
