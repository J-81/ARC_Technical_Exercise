FROM mambaorg/micromamba:1-lunar

# Do not pick up python packages from $HOME
ENV PYTHONNUSERSITE=1

# Install tools available via Anaconda using micromamba
COPY assets/conda.yml assets/conda.yml
RUN micromamba env create -f assets/conda.yml -y

# Add mOTU database version 3.0.1 to image
# TODO Hold off until storage requirement for building are availabile (>30 GB)
# COPY db_mOTU /opt/conda/envs/main/lib/python3.9/site-packages/motus/db_mOTU
# USER root
# RUN chmod -R a+rX /opt/conda/envs/main/lib/python3.9/site-packages/motus/db_mOTU
# USER mambauser

# Add mOTU database at available time of build
RUN eval "$(micromamba shell hook --shell bash)" \
    && micromamba activate main \
    && motus downloadDB

# Install BBMap tools
COPY assets/BBMap_39.01.tar.gz assets/BBMap_39.01.tar.gz
RUN tar -xzf assets/BBMap_39.01.tar.gz -C /opt/conda/envs/main/bin/

# Copy in nextflow workflow code and data
COPY modules modules
COPY main.nf main.nf
COPY bin bin

USER root
RUN apt-get update && apt-get install unzip

