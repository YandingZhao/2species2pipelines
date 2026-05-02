FROM satijalab/seurat:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    ca-certificates \
    g++ \
    make \
    python3-dev \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/rpy && \
    /opt/rpy/bin/pip install --no-cache-dir --default-timeout=300 --retries=10 \
    anndata==0.10.9 \
    bbknn==1.6.0 \
    h5py==3.16.0 \
    jax==0.4.28 \
    jaxlib==0.4.28 \
    leidenalg==0.10.2 \
    numpy==1.26.4 \
    pandas==2.2.3 \
    scanorama==1.7.4 \
    scanpy==1.10.3 \
    scikit-learn==1.5.2 \
    scipy==1.13.1 \
    scvi-tools==1.1.6

ENV RETICULATE_PYTHON=/opt/rpy/bin/python
ENV PATH="/opt/rpy/bin:${PATH}"

RUN R -e "install.packages(c('harmony','data.table','cowplot'), repos='https://cloud.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('SingleCellExperiment','batchelor'), ask=FALSE, update=FALSE)"
RUN R -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('satijalab/seurat-wrappers')"
RUN R -e "if (!requireNamespace('orthogene', quietly=TRUE)) remotes::install_github('neurogenomics/orthogene')"

WORKDIR /workspace
