if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("Pheatmap")


set.seed(123)
mat <- matrix(sample(1:100, 50, replace=TRUE), nrow=5, ncol=10)
colnames(mat) <- c("cond1", "cond2", "cond3", "cond4", "cond5", "ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5")
rownames(mat) <- paste0("gene", 1:5)
fold_change <- rowMeans(mat[, 1:5]) / rowMeans(mat[, 6:10])
fold_change

