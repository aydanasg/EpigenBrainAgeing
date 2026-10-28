# ---- Differential H3K27ac analysis pipeline ----
# Starting point: a DGEList object (edgeR) with raw peak counts and sample
# metadata (at minimum a `group` column for developmental stage, plus any
# covariates to include in the design)

library(edgeR)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
library(gprofiler2)

dgList <- readRDS("path/to/dgList.rds")

# ---- Filtering and normalisation ----
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep, ]

dgList <- calcNormFactors(dgList, method = "TMM")

# ---- Design matrix and dispersion ----
# covariates (e.g. tissue type, sex) are included here so TMM/dispersion
# estimation accounts for them, even though the final contrasts below only
# test the stage (group) term
design <- model.matrix(~0 + group + covariate1 + covariate2, data = dgList$samples)
colnames(design) <- gsub("group|covariate1|covariate2", "", colnames(design))
dgList <- estimateDisp(dgList, design, robust = TRUE)

# ---- Model fit ----
# fit uses a stage-only design: each stage's contrast is tested without the
# covariate terms, since covariates were already accounted for at the
# dispersion-estimation step above
design_fit <- model.matrix(~0 + group, data = dgList$samples)
colnames(design_fit) <- gsub("group", "", colnames(design_fit))
fit <- glmQLFit(y = dgList, design = design_fit, robust = TRUE)

# ---- One-vs-average-of-rest contrasts per stage ----
contrasts <- makeContrasts(
  stage1VSall = stage1 - (stage2 + stage3 + stage4)/3,
  stage2VSall = stage2 - (stage1 + stage3 + stage4)/3,
  stage3VSall = stage3 - (stage1 + stage2 + stage4)/3,
  stage4VSall = stage4 - (stage1 + stage2 + stage3)/3,
  levels = design_fit
)

# peak coordinates (GeneID, Chr, Start, End) for the consensus peak set
consensus_peaks <- fread("path/to/consensus_peaks.saf", sep = "\t")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb); genes_gr$gene_id <- names(genes_gr)
prom_gr <- promoters(genes_gr, upstream = 2000, downstream = 500)
sym <- mapIds(org.Hs.eg.db, keys = genes_gr$gene_id, column = "SYMBOL", keytype = "ENTREZID")
prom_gr$symbol <- sym[prom_gr$gene_id]

# promoter-capture Hi-C loops, used as tier-b evidence below
pcHiC <- fread("path/to/pcHiC_loops.bed")
loops <- pcHiC[, c("sourceChrom","sourceStart","sourceEnd","targetChrom","targetStart","targetEnd")]
a1 <- GRanges(loops$sourceChrom, IRanges(loops$sourceStart, loops$sourceEnd))
a2 <- GRanges(loops$targetChrom, IRanges(loops$targetStart, loops$targetEnd))
oa1 <- findOverlaps(a1, prom_gr); oa2 <- findOverlaps(a2, prom_gr)
lg1 <- data.table(loop = queryHits(oa1), gene = prom_gr$symbol[subjectHits(oa1)])
lg2 <- data.table(loop = queryHits(oa2), gene = prom_gr$symbol[subjectHits(oa2)])

stages <- colnames(contrasts)
fdr_cutoff <- 0.05
top_n_peaks <- 200   # peaks per direction feeding pathway analysis

for (q in stages) {

  # ---- differential peaks for this stage ----
  qlf <- glmQLFTest(glmfit = fit, contrast = contrasts[, q])
  top_deg <- as.data.frame(topTags(qlf, n = Inf, sort.by = "PValue"))
  top_deg <- top_deg[top_deg$FDR < fdr_cutoff, ]
  top_deg <- merge(top_deg, consensus_peaks, by.x = "genes", by.y = "GeneID")
  if (nrow(top_deg) == 0) next

  # ---- up/down defined by the sign of logFC ----
  top_deg$direction <- ifelse(top_deg$logFC > 0, "up", "down")

  # ---- three-tier peak-to-gene annotation ----
  peaks_gr <- GRanges(top_deg$Chr, IRanges(top_deg$Start, top_deg$End), PeakID = top_deg$genes)

  pa <- annotatePeak(peaks_gr, tssRegion = c(-2000, 500), TxDb = txdb, level = "gene", annoDb = "org.Hs.eg.db")
  ag <- as.GRanges(pa)
  prom <- ag[grepl("Promoter", ag$annotation)]      # tier a: promoter overlap
  distal <- ag[!grepl("Promoter", ag$annotation)]

  a_df <- data.table(PeakID = prom$PeakID, gene = prom$SYMBOL, tier = "a")

  op1 <- findOverlaps(distal, a1); op2 <- findOverlaps(distal, a2)   # tier b: pcHi-C loop
  p1 <- data.table(PeakID = distal$PeakID[queryHits(op1)], loop = subjectHits(op1))
  p2 <- data.table(PeakID = distal$PeakID[queryHits(op2)], loop = subjectHits(op2))
  b_df <- unique(rbind(merge(p1, lg2, by = "loop")[, .(PeakID, gene)],
                       merge(p2, lg1, by = "loop")[, .(PeakID, gene)]))
  b_df <- b_df[!is.na(gene)]; b_df[, tier := "b"]

  unresolved <- distal[!distal$PeakID %in% b_df$PeakID]
  op3 <- findOverlaps(unresolved, genes_gr, maxgap = 10000)         # tier c: nearest gene <=10kb
  c_df <- data.table(PeakID = unresolved$PeakID[queryHits(op3)],
                     gene = sym[genes_gr$gene_id[subjectHits(op3)]], tier = "c")
  c_df <- c_df[!is.na(gene)]

  gene_annotation <- rbind(a_df, b_df, c_df)[!is.na(gene) & gene != ""]
  gene_annotation <- merge(gene_annotation, top_deg[, c("genes","direction","FDR")],
                            by.x = "PeakID", by.y = "genes")

  # ---- top N genes per direction (by FDR), fed into pathway enrichment ----
  sel_genes <- function(dir) {
    dt <- unique(gene_annotation[gene_annotation$direction == dir, c("PeakID","FDR")])
    dt <- dt[order(dt$FDR), ]
    top_peaks <- head(dt$PeakID, top_n_peaks)
    unique(gene_annotation$gene[gene_annotation$PeakID %in% top_peaks])
  }
  up_genes   <- sel_genes("up")
  down_genes <- sel_genes("down")

  # ---- pathway enrichment (GO:BP, KEGG, Reactome) ----
  run_pathways <- function(genes) {
    if (length(genes) < 1) return(NULL)
    gost(query = genes, organism = "hsapiens", ordered_query = TRUE,
         evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
         sources = c("GO:BP","KEGG","REAC"))$result
  }
  pathways_up   <- run_pathways(up_genes)
  pathways_down <- run_pathways(down_genes)
}
