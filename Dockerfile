# arm64 for macbook M chip
FROM risserlin/bcb420-base-image:winter2025-arm64

# Install R packages
RUN R -e "install.packages(c('DESeq2', 'pheatmap', 'enrichplot'), repos='http://cran.rstudio.com/')"

EXPOSE 8787
