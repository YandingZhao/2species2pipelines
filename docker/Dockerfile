FROM satijalab/seurat:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('harmony','data.table','cowplot'), repos='https://cloud.r-project.org')"

WORKDIR /workspace
