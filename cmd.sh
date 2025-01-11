docker build --tag 'rstudio:latest' .
docker run -e PASSWORD=changeit -v "$(pwd)":/home/rstudio/projects -p 8787:8787 \
'rstudio:latest'