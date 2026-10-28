# ---- AD-associated regulatory changes versus ageing ----
# Compares AD-versus-control differential H3K27ac peaks with EpiAgeSig
# peaks, mapping both onto a shared consensus peak set, then classifies each
# overlapping peak into one of four quadrants by the direction of change in
# AD and with age.

library(data.table)
library(GenomicRanges)
library(dplyr)

classified_cpgs <- fread("/path/to/epiagesig_classified_peaks.txt")  # CpG, class (LI/LD/NLI/NLD), gam_slope, Chr, Start, End
ad_deg          <- fread("/path/to/ad_vs_control_differential_peaks.txt")  # genes, direction (UP/DOWN), logFC, Chr, Start, End

# ---- Build a shared consensus peak set from both inputs ----

consensus_peaks_gr <- GRanges(
  seqnames = c(classified_cpgs$Chr, ad_deg$Chr),
  ranges = IRanges(start = c(classified_cpgs$Start, ad_deg$Start),
                    end   = c(classified_cpgs$End,   ad_deg$End))
)
consensus_gr <- GenomicRanges::reduce(consensus_peaks_gr, min.gapwidth = 31)
consensus_gr$ConsensusID <- paste0("ConsensusPeak_", seq_along(consensus_gr))

# ---- Define the four directional peak sets ----

classified_cpgs <- classified_cpgs[, .(PeakID = CpG, Directions = class, Chr, Start, End)]
ad_deg          <- ad_deg[, .(PeakID = genes, Directions = direction, Chr, Start, End)]

peak_sets <- list(
  EpiAgeSig_UP   = classified_cpgs[Directions %in% c("LI", "NLI")],
  EpiAgeSig_DOWN = classified_cpgs[Directions %in% c("LD", "NLD")],
  AD_UP          = ad_deg[Directions == "UP"],
  AD_DOWN        = ad_deg[Directions == "DOWN"]
)

# ---- Map each set onto the consensus peak set ----

overlap_results <- lapply(names(peak_sets), function(set_name) {
  x <- peak_sets[[set_name]]
  query_gr <- GRanges(seqnames = x$Chr, ranges = IRanges(x$Start, x$End), peak_id = x$PeakID)
  hits <- findOverlaps(query = query_gr, subject = consensus_gr)
  data.table(
    peak_set = set_name,
    peak_id = mcols(query_gr)$peak_id[queryHits(hits)],
    consensus_peak_id = consensus_gr$ConsensusID[subjectHits(hits)]
  )
})
names(overlap_results) <- names(peak_sets)

# ---- Attach effect-size scores: GAM endpoint slope for ageing, logFC for AD ----

ageing_scores <- classified_cpgs_full[, .(PeakID, score = gam_slope)]  # from the linear/non-linear classification step
ad_scores     <- ad_deg_full[, .(PeakID, score = logFC)]

values_to_use <- list(EpiAgeSig_UP = ageing_scores, EpiAgeSig_DOWN = ageing_scores,
                       AD_UP = ad_scores, AD_DOWN = ad_scores)

combined_results <- rbindlist(lapply(names(overlap_results), function(nm) {
  merged <- merge(overlap_results[[nm]], values_to_use[[nm]], by.x = "peak_id", by.y = "PeakID", all.x = TRUE)
  data.table(consensus_peak_id = merged$consensus_peak_id, score = merged$score, name = nm)
}))

# only peaks appearing in exactly one EpiAgeSig set AND one AD set
combined_results <- combined_results[, if (uniqueN(name) == 2) .SD, by = consensus_peak_id]

# where a peak has multiple matches within a comparison, keep the strongest
combined_results[, comparison := fifelse(name %in% c("EpiAgeSig_UP", "EpiAgeSig_DOWN"), "Ageing", "AD")]
strongest <- combined_results[, .SD[which.max(abs(score))], by = .(consensus_peak_id, comparison)]

plot_df <- dcast(strongest, consensus_peak_id ~ comparison, value.var = "score")

# ---- Quadrant classification ----

plot_df[, quadrant := fifelse(Ageing > 0 & AD > 0, "Up in ageing / Up in AD",
                      fifelse(Ageing > 0 & AD < 0, "Up in ageing / Down in AD",
                      fifelse(Ageing < 0 & AD > 0, "Down in ageing / Up in AD",
                                                   "Down in ageing / Down in AD")))]
