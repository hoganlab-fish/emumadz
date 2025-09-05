FROM condaforge/mambaforge:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV VEP_VERSION=109.3

# Set working directory
WORKDIR /app

# Create conda environment with specified packages
RUN mamba create -n emumadz -c bioconda -c conda-forge -c defaults \
    python \
    pandas \
    htslib \
    make \
    wget \
    curl \
    unzip \    
    gcc_linux-64 \
    gxx_linux-64 \    
    ensembl-utils \
    pysam>=0.16.0 \
    bcftools>=1.19 \
    samtools>=1.19.2 \
    snpeff=5.2 \
    perl>=5.22 \
    perl-dbi>=1.643 \
    perl-archive-zip>=1.6.8 \
    perl-dbd-mysql>=4.050 \
    perl-set-intervaltree>=0.12 \
    perl-json>=4.10 \
    perl-perlio-gzip>=0.20 \
    perl-bio-bigfile>=1.07 \
    perl-list-moreutils>=0.430 \
    htslib \
    mysql-connector-c \
    && mamba clean -all

# Activate the environment
SHELL ["conda", "run", "-n", "emumadz", "/bin/bash", "-c"]

# Install cpanm first
RUN curl -L https://cpanmin.us | perl - App::cpanminus

# Install specific Perl modules with required versions
RUN cpanm \
    DBD::mysql \
    Set::IntervalTree \
    JSON \
    PerlIO::gzip \
    Bio::DB::BigFile

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
    cd .. && \
    rm -rf ${VEP_VERSION}.tar.gz

# Make conda environment the default and add tools to PATH
RUN echo "conda activate emumadz" >> ~/.bashrc
ENV CONDA_DEFAULT_ENV=emumadz
ENV PATH="/opt/gatk:/app/ensembl-vep-release-${VEP_VERSION}:/opt/conda/envs/emumadz/bin:$PATH"

# Verify installations
RUN python --version
RUN bcftools --version
RUN samtools --version
RUN gatk --version
RUN snpEff -version
RUN vep --help | head -5
RUN python -c "import pandas; print(f'pandas version: {pandas.__version__}')"
RUN python -c "import pysam; print(f'pysam version: {pysam.__version__}')"

# Set default command
CMD ["/bin/bash"]