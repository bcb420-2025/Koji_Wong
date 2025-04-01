library(data.table)
library(knitr)

#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  
           install.packages("RCurl")}, 
         finally = library("RCurl"))

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

ranked_path <- "./data/GFPPositive_GFPNegative.rnk"
analysis_table$hgnc_symbol <- rownames(analysis_table)
ranks <- as.matrix(analysis_table[,c("hgnc_symbol", "rank")])

write.table(
  ranks,
  file = ranked_path,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
)

# GSEA parameters
working_dir <- file.path(getwd(), "data")
output_dir <- file.path(getwd(), "generated")
analysis_name <- "GFPPositive_vs_GFPNegative"
rnk_file = "GFPPositive_GFPNegative.rnk"
rank_list <- as.matrix(analysis_table[,c("hgnc_symbol", "rank")])
run_gsea <- FALSE
gsea_run_id <- "1743482506447"
dest_gmt_file = ""

write.table(
  rank_list,
  file = rnk_file,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t",
)

analysis_name <- "GFPNegative_vs_GFPPositive"
output_dir <- "./data"

gsea_jar <- file.path(getwd(), "GSEA/gsea-cli.sh")
dest_gmt_file = ""

## download the .gmt file
if(dest_gmt_file == ""){
  gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"
  
  #list all the files on the server
  filenames = getURL(gmt_url)
  tc = textConnection(filenames)
  contents = readLines(tc)
  close(tc)
  
  #get the gmt that has all the pathways and does not include terms 
  # inferred from electronic annotations(IEA)
  #start with gmt file that has pathways only and GO Biological Process only.
  rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_noPFOCR_no_GO_iea.*.)(.gmt)(?=\">)",
                contents, perl = TRUE)
  gmt_file = unlist(regmatches(contents, rx))
  
  dest_gmt_file <- file.path(output_dir,gmt_file )
  
  #check if this gmt file already exists
  if(!file.exists(dest_gmt_file)){
    download.file(
      paste(gmt_url,gmt_file,sep=""),
      destfile=dest_gmt_file
    )
  }
}

if(run_gsea){
  command <- paste("",gsea_jar,  
                   "GSEAPreRanked -gmx", dest_gmt_file, 
                   "-rnk" ,file.path(working_dir,rnk_file), 
                   "-collapse false -nperm 1000 -scoring_scheme weighted", 
                   "-rpt_label ",analysis_name,
                   "  -plot_top_x 20 -rnd_seed 12345  -set_max 200",  
                   " -set_min 15 -zip_report false ",
                   " -out" ,output_dir, 
                   " > gsea_output.txt",sep=" ")
  system(command)
}

result_dir <- paste("data/", analysis_name, ".GseaPreranked.", gsea_run_id, sep="")
na_pos_file <- file.path(getwd(), 
                         result_dir, 
                         paste("gsea_report_for_na_pos_", 
                               gsea_run_id, 
                               ".tsv",
                               sep = ""))
na_neg_file <- file.path(getwd(), 
                         result_dir, 
                         paste("gsea_report_for_na_neg_", 
                               gsea_run_id, 
                               ".tsv",
                               sep = ""))
positive_ranked <- as.data.frame(fread(na_pos_file))
negative_ranked <- as.data.frame(fread(na_neg_file))

kable(positive_ranked[, 1:5])
