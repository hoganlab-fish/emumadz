FROM condaforge/mambaforge:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV VEP_VERSION=109.3

# Set working directory
WORKDIR /app

# Create conda environment with specified packages
# RUN mamba create -n emumadz -c bioconda -c conda-forge -c defaults \
#     python \
#     pandas \
#     htslib \
#     make \
#     wget \
#     curl \
#     unzip \    
#     gcc_linux-64 \
#     gxx_linux-64 \    
#     ensembl-utils \
#     pysam>=0.16.0 \
#     bcftools>=1.19 \
#     samtools>=1.19.2 \
#     snpeff=5.2 \
#     perl>=5.22 \
#     perl-dbi>=1.643 \
#     perl-archive-zip>=1.6.8 \
#     perl-dbd-mysql>=4.050 \
#     perl-set-intervaltree>=0.12 \
#     perl-json>=4.10 \
#     perl-perlio-gzip>=0.20 \
#     perl-bio-bigfile>=1.07 \
#     perl-list-moreutils>=0.430 \
#     htslib \
#     mysql-connector-c \
#     && mamba clean -all

RUN mamba create -n emumadz -c bioconda -c conda-forge -c defaults \
    _libgcc_mutex \
    _openmp_mutex \
    _r-mutex \
    alabaster \
    alsa-lib \
    amqp \
    appdirs \
    argcomplete \
    argtable2 \
    babel \
    backports \
    backports.functools_lru_cache \
    bcrypt \
    billiard \
    binutils_impl_linux-64 \
    binutils_linux-64 \
    blast \
    blast-legacy \
    boost \
    boost-cpp \
    brotli-python \
    bwidget \
    bzip2 \
    c-ares \
    ca-certificates \
    cached-property \
    cached_property \
    cd-hit \
    celery \
    certifi \
    cffi \
    charset-normalizer \
    click \
    click-didyoumean \
    click-plugins \
    click-repl \
    clustalo \
    clustalw \
    colorama \
    cryptography \
    curl \
    cyrus-sasl \
    docker-py \
    docutils \
    expat \
    font-ttf-dejavu-sans-mono \
    font-ttf-inconsolata \
    font-ttf-source-code-pro \
    font-ttf-ubuntu \
    fontconfig \
    fonts-conda-ecosystem \
    fonts-conda-forge \
    freetype \
    fribidi \
    gawk \
    gcc_impl_linux-64 \
    gcc_linux-64 \
    gfortran_impl_linux-64 \
    giflib \
    glib \
    glib-tools \
    gmp \
    gnutls \
    graphite2 \
    gunicorn \
    gxx_impl_linux-64 \
    harfbuzz \
    hmmer \
    htslib \
    icu \
    idna \
    imagesize \
    importlib-metadata \
    itsdangerous \
    jinja2 \
    jq \
    kernel-headers_linux-64 \
    keyutils \
    kombu \
    krb5 \
    lcms2 \
    ld_impl_linux-64 \
    lerc \
    libabseil \
    libasprintf \
    libblas \
    libcblas \
    libcups \
    libcurl \
    libdb \
    libdeflate \
    libedit \
    libev \
    libevent \
    libexpat \
    libffi \
    libgcc \
    libgcc-devel_linux-64 \
    libgcc-ng \
    libgettextpo \
    libgfortran \
    libgfortran-ng \
    libgfortran5 \
    libglib \
    libgomp \
    libiconv \
    libjpeg-turbo \
    liblapack \
    liblzma \
    liblzma-devel \
    libnghttp2 \
    libnsl \
    libntlm \
    libopenblas \
    libpng \
    libprotobuf \
    libsanitizer \
    libsodium \
    libsqlite \
    libssh2 \
    libstdcxx \
    libstdcxx-devel_linux-64 \
    libstdcxx-ng \
    libtiff \
    libuuid \
    libuv \
    libwebp-base \
    libxcb \
    libxcrypt \
    libxml2 \
    libzlib \
    lz4-c \
    mafft \
    make \
    markdown-it-py \
    markupsafe \
    mdit-py-plugins \
    mdurl \
    mpfr \
    muscle \
    mysql \
    mysql-client \
    mysql-common \
    mysql-connector-c \
    mysql-devel \
    mysql-libs \
    mysql-server \
    myst-parser \
    ncurses \
    nettle \
    numpy \
    oniguruma \
    openjdk \
    openssl \
    packaging \
    paml \
    pandas \
    pango \
    parallel \
    paramiko \
    patchelf \
    pcre \
    pcre2 \
    perl \
    perl-algorithm-diff \
    perl-archive-tar \
    perl-archive-zip \
    perl-b-cow \
    perl-base \
    perl-bio-asn1-entrezgene \
    perl-bio-bigfile \
    perl-bio-coordinate \
    perl-bio-featureio \
    perl-bio-samtools \
    perl-bio-searchio-hmmer \
    perl-bio-tools-phylo-paml \
    perl-bio-tools-run-alignment-clustalw \
    perl-bio-tools-run-alignment-tcoffee \
    perl-bioperl \
    perl-bioperl-core \
    perl-bioperl-run \
    perl-business-isbn \
    perl-business-isbn-data \
    perl-capture-tiny \
    perl-carp \
    perl-class-data-inheritable \
    perl-clone \
    perl-common-sense \
    perl-compress-raw-bzip2 \
    perl-compress-raw-zlib \
    perl-constant \
    perl-data-dump \
    perl-data-dumper \
    perl-db_file \
    perl-dbd-mysql \
    perl-dbi \
    perl-devel-checklib \
    perl-devel-stacktrace \
    perl-digest-hmac \
    perl-digest-md5 \
    perl-encode \
    perl-encode-locale \
    perl-exception-class \
    perl-exporter \
    perl-exporter-tiny \
    perl-extutils-cppguess \
    perl-extutils-makemaker \
    perl-file-listing \
    perl-file-path \
    perl-file-slurper \
    perl-file-sort \
    perl-file-spec \
    perl-file-temp \
    perl-getopt-long \
    perl-html-parser \
    perl-html-tagset \
    perl-http-cookies \
    perl-http-daemon \
    perl-http-date \
    perl-http-message \
    perl-http-negotiate \
    perl-inc-latest \
    perl-io-compress \
    perl-io-html \
    perl-io-socket-ssl \
    perl-io-string \
    perl-io-tty \
    perl-io-zlib \
    perl-ipc-run \
    perl-json \
    perl-json-xs \
    perl-libwww-perl \
    perl-libxml-perl \
    perl-list-moreutils \
    perl-list-moreutils-xs \
    perl-lwp-mediatypes \
    perl-lwp-protocol-https \
    perl-mime-base64 \
    perl-module-build \
    perl-module-load \
    perl-mozilla-ca \
    perl-net-http \
    perl-net-ssleay \
    perl-ntlm \
    perl-parent \
    perl-pathtools \
    perl-perlio-gzip \
    perl-scalar-list-utils \
    perl-set-intervaltree \
    perl-socket \
    perl-storable \
    perl-sub-uplevel \
    perl-test-deep \
    perl-test-differences \
    perl-test-exception \
    perl-test-fatal \
    perl-test-most \
    perl-test-needs \
    perl-test-requiresinternet \
    perl-test-warn \
    perl-test-warnings \
    perl-text-diff \
    perl-time-hires \
    perl-time-local \
    perl-timedate \
    perl-tree-dag_node \
    perl-try-tiny \
    perl-types-serialiser \
    perl-uri \
    perl-url-encode \
    perl-www-robotrules \
    perl-xml-dom \
    perl-xml-dom-xpath \
    perl-xml-parser \
    perl-xml-regexp \
    perl-xml-xpathengine \
    pip \
    pixman \
    poa \
    prompt-toolkit \
    prompt_toolkit \
    provean \
    pthread-stubs \
    pybiolib \
    pycparser \
    pycryptodome \
    pygments \
    pyjwt \
    pynacl \
    pysam \
    pysocks \
    python \
    python-dateutil \
    python_abi \
    pytz \
    pyvcf \
    pywin32-on-windows \
    pyyaml \
    readline \
    requests \
    rich \
    samtools \
    sed \
    setuptools \
    six \
    snowballstemmer \
    snpeff \
    snpsift \
    sqlite \
    sysroot_linux-64 \
    t-coffee \
    tabix \
    tk \
    tktable \
    toml \
    tomlkit \
    tqdm \
    tree \
    typing-extensions \
    typing_extensions \
    tzdata \
    urllib3 \
    viennarna \
    vine \
    wcwidth \
    websocket-client \
    werkzeug \
    wheel \
    xmltodict \
    xorg-fixesproto \
    xorg-inputproto \
    xorg-kbproto \
    xorg-libice \
    xorg-libsm \
    xorg-libx11 \
    xorg-libxau \
    xorg-libxdmcp \
    xorg-libxext \
    xorg-libxfixes \
    xorg-libxi \
    xorg-libxrender \
    xorg-libxt \
    xorg-libxtst \
    xorg-recordproto \
    xorg-renderproto \
    xorg-xextproto \
    xorg-xproto \
    xz \
    xz-gpl-tools \
    xz-tools \
    yaml \
    yq \
    zipp \
    zlib \
    zstd \
    && mamba clean -all
    
# Activate the environment
SHELL ["conda", "run", "-n", "emumadz", "/bin/bash", "-c"]

# Install cpanm first
RUN curl -L https://cpanmin.us | perl - App::cpanminus

# Install specific Perl modules with required versions
# RUN cpanm --mirror http://cpan.cpan.org \
#     Digest::SHA1 \
#     Bio::DB::GenBank \
#     Bio::DB::GenPept \
#     Cache::FileCache \
#     Bio::EnsEMBL::Registry \
#     Bio::Perl \
#     Set::IntervalTree \
#     JSON \
#     PerlIO::gzip \
#     Bio::DB::BigFile

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
