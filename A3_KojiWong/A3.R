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


destfile = "precomputed/qlf_output_hits.rds"
if (!file.exists(destfile)) {
  download.file(
    url = "https://github.com/bcb420-2025/Koji_Wong/raw/refs/heads/main/A3_KojiWong/precomputed/qlf_output_hits.rds",
    destfile = destfile
  )
}
# Load in the normalized count matrix from our rds object file
qlf_output_hits <- readRDS(destfile)
analysis_table <- qlf_output_hits$table
head(analysis_table)
analysis_table$rank <- -log(
  analysis_table$PValue, base = 10
) * sign(analysis_table$logFC)

sorted_pvalue_analysis_table <- analysis_table[order(analysis_table$PValue), ]
head(sorted_pvalue_analysis_table)
sorted_rank_analysis_table <- analysis_table[order(analysis_table$rank), ]
head(sorted_rank_analysis_table)

analysis_table$hgnc_symbol <- rownames(analysis_table)

ranked_genes_path = "./data/GFPPositive_GFPNegative.rnk"
rank_list <- as.matrix(analysis_table[,c("hgnc_symbol", "rank")])
write.table(
  rank_list,
  file = ranked_genes_path,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
)

analysis_name <- "GFPNegative_vs_GFPPositive"
output_dir <- "./data"
gsea_jar <- file.path(getwd(), "GSEA_4.4.0/gsea-cli.sh")

if (!file.exists("./gsea_output.txt")) {
  dest_gmt_file = "./data/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_symbol.gmt"
  command <- paste("",gsea_jar,  
                   "GSEAPreRanked -gmx", dest_gmt_file,
                   "-rnk" , ranked_genes_path,
                   "-collapse false -nperm 1000 -scoring_scheme weighted",
                   "-rpt_label ", analysis_name,
                   "  -plot_top_x 20 -rnd_seed 12345  -set_max 200",
                   " -set_min 15 -zip_report false ",
                   " -out" ,output_dir,
                   " > gsea_output.txt",sep=" ")
  system(command, intern = TRUE)
}