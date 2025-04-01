# Getting the normalized RNASeq counts from Assignment 1
destfile = "precomputed/upreg_genes_gost.rds"
if (!file.exists(destfile)) {
  download.file(
    url = "https://github.com/bcb420-2025/Koji_Wong/raw/refs/heads/main/A3_KojiWong/precomputed/upreg_genes_gost.rds",
    destfile = destfile
  )
}
# Load in the normalized count matrix from our rds object file
upreg_genes_gost <- readRDS(destfile)
head(upreg_genes_gost)