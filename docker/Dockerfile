FROM satijalab/seurat:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('harmony','data.table','cowplot'), repos='https://cloud.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('SingleCellExperiment','batchelor'), ask=FALSE, update=FALSE)"
RUN R -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('satijalab/seurat-wrappers')"

WORKDIR /workspace
