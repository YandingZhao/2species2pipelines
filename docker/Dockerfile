FROM satijalab/seurat:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    ca-certificates \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/rpy && \
    /opt/rpy/bin/pip install --no-cache-dir \
    anndata==0.10.9 \
    h5py==3.16.0 \
    numpy==1.26.4 \
    scipy==1.14.1

ENV RETICULATE_PYTHON=/opt/rpy/bin/python

RUN R -e "install.packages(c('harmony','data.table','cowplot'), repos='https://cloud.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('SingleCellExperiment','batchelor','zellkonverter'), ask=FALSE, update=FALSE)"
RUN R -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('satijalab/seurat-wrappers')"

WORKDIR /workspace
