# remotes::install_github("aydanasg/hmagmaR")
library(data.table)
library(hmagmaR)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# ---- Build a stage-specific H-MAGMA annotation file ----

stage <- c("stage1", "stage2", "stage3", "stage4")

# required reference files
snpgeneexon     <- fread("/path/to/required_files/exonic_promoter_snps_chipSeeker.txt")
snps            <- fread("/path/to/required_files/g1000_snps.txt")
annotated_genes <- fread("/path/to/required_files/annotated_genes_TxDb.Hsapiens.UCSC.hg19.knownGene.org.Hs.eg.db.txt")

# promoter-capture Hi-C loops (pooled across stages given limited per-stage
# sample size)
hic <- fread("/path/to/pcHiC_loops.bed")
hic <- hic[, c("sourceChrom", "sourceStart", "sourceEnd", "targetChrom", "targetStart", "targetEnd")]
hic <- hic[!grepl("chrX|chrY", hic$sourceChrom) & !grepl("chrX|chrY", hic$targetChrom)]
colnames(hic) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")

for (s in seq_along(stage)) {

  # stage-specific active H3K27ac regulatory regions
  regulatoryRegions <- fread(paste0("/path/to/bedfiles/", stage[s], "_consensus.bed"))

  hmagmaR::AnnotationFileHmagma(
    hic = hic,
    regulatoryRegions = regulatoryRegions,
    snps = snps,
    annotated_genes = annotated_genes,
    snpgeneexon = snpgeneexon,
    AnnotationFile = paste0("/path/to/annotation_files/olig2_ageing.", stage[s])
  )
}

# ---- Gene-level H-MAGMA association per disease and stage ----

disease <- c("AD_Jansen2019", "AD_Kunkle2019", "AD_Bellenguez2022", "PD_GP22025",
             "SCZ_Trubetskoy2022", "SVD_Sargurupremraj2022", "PVS_Duperron2023",
             "Covid_Castineira2023", "MS_Andlauer2016", "ALS_Rheenen2021",
             "Longevity_Deelan2019_90th")

for (d in seq_along(disease)) {
  for (s in seq_along(stage)) {

    hmagmaR::GeneLevelAnalysis_hmagma(
      magma = "/path/to/magma/magma",
      g1000 = "/path/to/g1000_reference/g1000_eur",
      gwas = paste0("/path/to/munged_gwas/", disease[d], "_munged.tsv"),
      AnnotationFile = paste0("/path/to/annotation_files/olig2_ageing.", stage[s], ".transcript.annot"),
      output = paste0("/path/to/magma_output/olig2_ageing.", disease[d], ".", stage[s])
    )
  }
}
