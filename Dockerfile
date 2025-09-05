FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV VEP_VERSION=109.3

# Set working directory
WORKDIR /app

# Update and install absolutely essential packages only
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Install mamba
RUN conda install -c conda-forge mamba -y

# Create conda environment with specified packages
RUN mamba create -n emumadz -c conda-forge -c bioconda -c defaults \
    python \
    pandas==2.3.0 \
    pysam \
    bcftools==1.19 \
    samtools==1.19.2 \
    snpeff==5.2 \
    perl \
    perl-dbi \
    perl-dbd-mysql \
    perl-xml-parser \
    perl-archive-zip \
    perl-archive-extract \
    perl-compress-raw-zlib \
    perl-file-copy-recursive \
    perl-module-build \
    perl-bioperl \
    perl-bio-db-hts \
    perl-set-intervaltree \
    perl-json \
    perl-perlio-gzip \    
    && mamba clean -all

# Activate the environment
SHELL ["conda", "run", "-n", "emumadz", "/bin/bash", "-c"]

# Install GATK 4.5.0.0 manually
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    && unzip gatk-4.5.0.0.zip \
    && mv gatk-4.5.0.0 /opt/gatk \
    && rm gatk-4.5.0.0.zip

# Download and install VEP
RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/${VEP_VERSION}.tar.gz && \
    tar -xzf ${VEP_VERSION}.tar.gz && \
    cd ensembl-vep-release-${VEP_VERSION} && \
    perl INSTALL.pl --NO_UPDATE --NO_HTSLIB --NO_TEST && \
    rm -rf /app/${VEP_VERSION}.tar.gz

# Make conda environment the default and add tools to PATH
RUN echo "conda activate emumadz" >> ~/.bashrc
ENV CONDA_DEFAULT_ENV=emumadz
ENV PATH="/opt/gatk:/app/ensembl-vep-release-${VEP_VERSION}:/opt/conda/envs/emumadz/bin:$PATH"

# Verify installations
# Test PATH for each tool
RUN python --version 2>&1 || echo "Python version failed"
RUN bcftools --version 2>&1 || echo "bcftools version failed" 
RUN samtools --version 2>&1 || echo "samtools version failed"
RUN gatk --version 2>&1 || echo "GATK version failed"
RUN snpEff -version 2>&1 || echo "snpEff version failed"
RUN vep --help 2>&1 | head -5 || echo "VEP help failed"
RUN python -c "import pysam; print(f'pysam version: {pysam.__version__}')"
RUN python -c "import pandas; print(f'pandas version: {pandas.__version__}')"

# Set default command
CMD ["/bin/bash"]