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
    _libgcc_mutex=0.1=conda_forge \
    _openmp_mutex=4.5=2_gnu \
    _r-mutex=1.0.1=anacondar_1 \
    alabaster=0.7.13=pyhd8ed1ab_0 \
    alsa-lib=1.2.6.1=h7f98852_0 \
    amqp=5.2.0=pyhd8ed1ab_1 \
    appdirs=1.4.4=pyh9f0ad1d_0 \
    argcomplete=3.0.8=pyhd8ed1ab_0 \
    argtable2=2.13=hd590300_1004 \
    babel=2.14.0=pyhd8ed1ab_0 \
    backports=1.0=pyhd8ed1ab_4 \
    backports.functools_lru_cache=2.0.0=pyhd8ed1ab_0 \
    bcrypt=3.2.2=py37h540881e_0 \
    billiard=3.6.4.0=py37h540881e_2 \
    binutils_impl_linux-64=2.43=h4bf12b8_4 \
    binutils_linux-64=2.43=h4852527_4 \
    blast=2.5.0=hc0b0e79_3 \
    blast-legacy=2.2.26=h9ee0642_4 \
    boost=1.78.0=py37h48bf904_0 \
    boost-cpp=1.78.0=h2c5509c_4 \
    brotli-python=1.0.9=py37hd23a5d3_7 \
    bwidget=1.10.1=ha770c72_0 \
    bzip2=1.0.8=h4bc722e_7 \
    c-ares=1.34.5=hb9d3cd8_0 \
    ca-certificates=2025.4.26=hbd8a1cb_0 \
    cached-property=1.5.2=hd8ed1ab_1 \
    cached_property=1.5.2=pyha770c72_1 \
    cd-hit=4.8.1=h43eeafb_11 \
    celery=5.2.7=pyhd8ed1ab_0 \
    certifi=2024.8.30=pyhd8ed1ab_0 \
    cffi=1.15.1=py37h43b0acd_1 \
    charset-normalizer=3.4.0=pyhd8ed1ab_0 \
    click=8.1.3=py37h89c1867_0 \
    click-didyoumean=0.3.1=pyhd8ed1ab_0 \
    click-plugins=1.1.1=py_0 \
    click-repl=0.3.0=pyhd8ed1ab_0 \
    clustalo=1.2.4=h503566f_9 \
    clustalw=2.1=h9948957_12 \
    colorama=0.4.6=pyhd8ed1ab_0 \
    cryptography=38.0.2=py37h5994e8b_1 \
    curl=8.14.1=h14366f7_1 \
    cyrus-sasl=2.1.28=hd9c7081_0 \
    docker-py=7.0.0=pyhd8ed1ab_0 \
    docutils=0.19=py37h89c1867_0 \
    expat=2.7.0=h5888daf_0 \
    font-ttf-dejavu-sans-mono=2.37=hab24e00_0 \
    font-ttf-inconsolata=3.000=h77eed37_0 \
    font-ttf-source-code-pro=2.038=h77eed37_0 \
    font-ttf-ubuntu=0.83=h77eed37_3 \
    fontconfig=2.14.2=h14ed4e7_0 \
    fonts-conda-ecosystem=1=0 \
    fonts-conda-forge=1=0 \
    freetype=2.12.1=h267a509_2 \
    fribidi=1.0.10=h36c2ea0_0 \
    gawk=5.3.1=hcd3d067_0 \
    gcc_impl_linux-64=13.4.0=h69c5793_4 \
    gcc_linux-64=13.4.0=h621f4e2_11 \
    gfortran_impl_linux-64=13.4.0=h847f9e2_4 \
    giflib=5.2.2=hd590300_0 \
    glib=2.80.2=hf974151_0 \
    glib-tools=2.80.2=hb6ce0ca_0 \
    gmp=6.3.0=hac33072_2 \
    gnutls=3.6.13=h79a8f9a_0 \
    graphite2=1.3.14=h295c915_1 \
    gunicorn=20.1.0=py37h89c1867_2 \
    gxx_impl_linux-64=13.4.0=haf17267_4 \
    harfbuzz=9.0.0=hfac3d4d_0 \
    hmmer=3.4=h503566f_3 \
    htslib=1.21=h5efdd21_0 \
    icu=73.2=h59595ed_0 \
    idna=3.10=pyhd8ed1ab_0 \
    imagesize=1.4.1=pyhd8ed1ab_0 \
    importlib-metadata=4.11.4=py37h89c1867_0 \
    itsdangerous=2.1.2=pyhd8ed1ab_0 \
    jinja2=3.1.4=pyhd8ed1ab_0 \
    jq=1.7.1=hd590300_0 \
    kernel-headers_linux-64=3.10.0=he073ed8_18 \
    keyutils=1.6.1=h166bdaf_0 \
    kombu=5.2.4=py37h89c1867_1 \
    krb5=1.21.3=h659f571_0 \
    lcms2=2.16=hb7c19ff_0 \
    ld_impl_linux-64=2.43=h712a8e2_4 \
    lerc=4.0.0=h0aef613_1 \
    libabseil=20240116.2=cxx17_he02047a_1 \
    libasprintf=0.25.1=h3f43e3d_1 \
    libblas=3.9.0=31_h59b9bed_openblas \
    libcblas=3.9.0=31_he106b2a_openblas \
    libcups=2.3.3=h4637d8d_4 \
    libcurl=8.14.1=hc1efc7f_1 \
    libdb=6.2.32=h9c3ff4c_0 \
    libdeflate=1.20=hd590300_0 \
    libedit=3.1.20250104=pl5321h7949ede_0 \
    libev=4.33=hd590300_2 \
    libevent=2.1.12=hf998b51_1 \
    libexpat=2.7.0=h5888daf_0 \
    libffi=3.4.6=h2dba641_1 \
    libgcc=14.2.0=h767d61c_2 \
    libgcc-devel_linux-64=13.4.0=hba01cd7_104 \
    libgcc-ng=14.2.0=h69a702a_2 \
    libgettextpo=0.25.1=h3f43e3d_1 \
    libgfortran=14.2.0=h69a702a_2 \
    libgfortran-ng=14.2.0=h69a702a_2 \
    libgfortran5=14.2.0=hf1ad2bd_2 \
    libglib=2.80.2=hf974151_0 \
    libgomp=14.2.0=h767d61c_2 \
    libiconv=1.18=h4ce23a2_1 \
    libjpeg-turbo=3.1.0=hb9d3cd8_0 \
    liblapack=3.9.0=31_h7ac8fdf_openblas \
    liblzma=5.8.1=hb9d3cd8_0 \
    liblzma-devel=5.8.1=hb9d3cd8_0 \
    libnghttp2=1.58.0=h47da74e_1 \
    libnsl=2.0.1=hd590300_0 \
    libntlm=1.8=hb9d3cd8_0 \
    libopenblas=0.3.29=pthreads_h94d23a6_0 \
    libpng=1.6.43=h2797004_0 \
    libprotobuf=4.25.3=h08a7969_0 \
    libsanitizer=13.4.0=h14bf0c3_4 \
    libsodium=1.0.18=h36c2ea0_1 \
    libsqlite=3.46.0=hde9e2c9_0 \
    libssh2=1.11.1=h251f7ec_0 \
    libstdcxx=14.2.0=h8f9b012_2 \
    libstdcxx-devel_linux-64=13.4.0=hba01cd7_104 \
    libstdcxx-ng=14.2.0=h4852527_2 \
    libtiff=4.6.0=h1dd3fc0_3 \
    libuuid=2.38.1=h0b41bf4_0 \
    libuv=1.51.0=hb9d3cd8_0 \
    libwebp-base=1.5.0=h851e524_0 \
    libxcb=1.15=h0b41bf4_0 \
    libxcrypt=4.4.36=hd590300_1 \
    libxml2=2.12.7=hc051c1a_1 \
    libzlib=1.2.13=h4ab18f5_6 \
    lz4-c=1.9.4=hcb278e6_0 \
    mafft=7.526=h4bc722e_0 \
    make=4.4.1=hb9d3cd8_2 \
    markdown-it-py=2.2.0=pyhd8ed1ab_0 \
    markupsafe=2.1.1=py37h540881e_1 \
    mdit-py-plugins=0.3.5=pyhd8ed1ab_0 \
    mdurl=0.1.2=pyhd8ed1ab_0 \
    mpfr=4.2.1=h90cbb55_3 \
    muscle=5.3=h9948957_3 \
    mysql=8.3.0=h27aab58_4 \
    mysql-client=8.3.0=h545f5f4_4 \
    mysql-common=8.3.0=hf1915f5_4 \
    mysql-connector-c=6.1.11=h659d440_1008 \
    mysql-devel=8.3.0=hf1915f5_4 \
    mysql-libs=8.3.0=hca2cd23_4 \
    mysql-server=8.3.0=h27b9c54_4 \
    myst-parser=1.0.0=pyhd8ed1ab_0 \
    ncurses=6.5=h2d0b736_3 \
    nettle=3.4.1=h1bed415_1002 \
    numpy=1.21.6=py37h976b520_0 \
    oniguruma=6.9.10=hb9d3cd8_0 \
    openjdk=11.0.1=h516909a_1016 \
    openssl=3.5.2=h26f9b46_0 \
    packaging=23.2=pyhd8ed1ab_0 \
    paml=4.10.9=h7b50bb2_0 \
    pandas=1.3.5=py37he8f5f7f_0 \
    pango=1.54.0=h4c5309f_1 \
    parallel=20250622=ha770c72_0 \
    paramiko=3.5.0=pyhd8ed1ab_0 \
    patchelf=0.17.2=h58526e2_0 \
    pcre=8.45=h9c3ff4c_0 \
    pcre2=10.43=hcad00b1_0 \
    perl=5.32.1=7_hd590300_perl5 \
    perl-algorithm-diff=1.201=pl5321hd8ed1ab_0 \
    perl-archive-tar=2.40=pl5321hdfd78af_0 \
    perl-archive-zip=1.68=pl5321hdfd78af_0 \
    perl-b-cow=0.007=pl5321hb9d3cd8_1 \
    perl-base=2.23=pl5321hd8ed1ab_0 \
    perl-bio-asn1-entrezgene=1.73=pl5321hdfd78af_3 \
    perl-bio-bigfile=1.07=pl5321h41f7678_5 \
    perl-bio-coordinate=1.007001=pl5321hdfd78af_3 \
    perl-bio-featureio=1.6.905=pl5321hdfd78af_4 \
    perl-bio-samtools=1.43=pl5321he4a0461_4 \
    perl-bio-searchio-hmmer=1.7.3=pl5321hdfd78af_0 \
    perl-bio-tools-phylo-paml=1.7.3=pl5321hdfd78af_3 \
    perl-bio-tools-run-alignment-clustalw=1.7.4=pl5321hdfd78af_3 \
    perl-bio-tools-run-alignment-tcoffee=1.7.4=pl5321hdfd78af_5 \
    perl-bioperl=1.7.8=hdfd78af_1 \
    perl-bioperl-core=1.7.8=pl5321hdfd78af_1 \
    perl-bioperl-run=1.007003=pl5321hdfd78af_0 \
    perl-business-isbn=3.007=pl5321hd8ed1ab_0 \
    perl-business-isbn-data=20210112.006=pl5321hd8ed1ab_0 \
    perl-capture-tiny=0.48=pl5321ha770c72_1 \
    perl-carp=1.50=pl5321hd8ed1ab_0 \
    perl-class-data-inheritable=0.09=pl5321ha770c72_0 \
    perl-clone=0.46=pl5321hb9d3cd8_1 \
    perl-common-sense=3.75=pl5321hd8ed1ab_0 \
    perl-compress-raw-bzip2=2.201=pl5321hbf60520_1 \
    perl-compress-raw-zlib=2.202=pl5321hadc24fc_0 \
    perl-constant=1.33=pl5321hd8ed1ab_0 \
    perl-data-dump=1.25=pl5321h7b50bb2_2 \
    perl-data-dumper=2.183=pl5321hb9d3cd8_1 \
    perl-db_file=1.858=pl5321hb9d3cd8_0 \
    perl-dbd-mysql=4.050=pl5321h9948957_3 \
    perl-dbi=1.643=pl5321hb9d3cd8_1 \
    perl-devel-checklib=1.16=pl5321h7b50bb2_1 \
    perl-devel-stacktrace=2.04=pl5321h296ab09_0 \
    perl-digest-hmac=1.04=pl5321hdfd78af_0 \
    perl-digest-md5=2.59=pl5321hb9d3cd8_3 \
    perl-encode=3.21=pl5321hb9d3cd8_1 \
    perl-encode-locale=1.05=pl5321hdfd78af_7 \
    perl-exception-class=1.45=pl5321ha770c72_0 \
    perl-exporter=5.74=pl5321hd8ed1ab_0 \
    perl-exporter-tiny=1.002002=pl5321hd8ed1ab_0 \
    perl-extutils-cppguess=0.26=pl5321h9948957_3 \
    perl-extutils-makemaker=7.70=pl5321hd8ed1ab_0 \
    perl-file-listing=6.16=pl5321hdfd78af_0 \
    perl-file-path=2.18=pl5321hd8ed1ab_0 \
    perl-file-slurper=0.014=pl5321hdfd78af_0 \
    perl-file-sort=1.01=pl5321hdfd78af_3 \
    perl-file-spec=3.48_01=pl5321hdfd78af_2 \
    perl-file-temp=0.2304=pl5321hd8ed1ab_0 \
    perl-getopt-long=2.58=pl5321hdfd78af_0 \
    perl-html-parser=3.83=pl5321h9948957_1 \
    perl-html-tagset=3.24=pl5321hdfd78af_0 \
    perl-http-cookies=6.11=pl5321hdfd78af_0 \
    perl-http-daemon=6.16=pl5321hdfd78af_0 \
    perl-http-date=6.06=pl5321hdfd78af_0 \
    perl-http-message=7.00=pl5321hdfd78af_0 \
    perl-http-negotiate=6.01=pl5321hdfd78af_4 \
    perl-inc-latest=0.500=pl5321ha770c72_0 \
    perl-io-compress=2.201=pl5321h503566f_5 \
    perl-io-html=1.004=pl5321hdfd78af_0 \
    perl-io-socket-ssl=2.075=pl5321hd8ed1ab_0 \
    perl-io-string=1.08=pl5321hdfd78af_4 \
    perl-io-tty=1.20=pl5321hb9d3cd8_3 \
    perl-io-zlib=1.14=pl5321hdfd78af_0 \
    perl-ipc-run=20200505.0=pl5321hdfd78af_0 \
    perl-json=4.10=pl5321hdfd78af_1 \
    perl-json-xs=4.03=pl5321h9948957_4 \
    perl-libwww-perl=6.68=pl5321hdfd78af_0 \
    perl-libxml-perl=0.08=pl5321hdfd78af_3 \
    perl-list-moreutils=0.430=pl5321hdfd78af_0 \
    perl-list-moreutils-xs=0.430=pl5321h7b50bb2_5 \
    perl-lwp-mediatypes=6.04=pl5321hdfd78af_1 \
    perl-lwp-protocol-https=6.10=pl5321hdfd78af_0 \
    perl-mime-base64=3.16=pl5321hb9d3cd8_3 \
    perl-module-build=0.4234=pl5321ha770c72_1 \
    perl-module-load=0.34=pl5321hdfd78af_0 \
    perl-mozilla-ca=20211001=pl5321hdfd78af_0 \
    perl-net-http=6.22=pl5321hdfd78af_0 \
    perl-net-ssleay=1.92=pl5321hf14f497_1 \
    perl-ntlm=1.09=pl5321hdfd78af_5 \
    perl-parent=0.243=pl5321hd8ed1ab_0 \
    perl-pathtools=3.75=pl5321hb9d3cd8_2 \
    perl-perlio-gzip=0.20=pl5321he4a0461_6 \
    perl-scalar-list-utils=1.69=pl5321hb9d3cd8_0 \
    perl-set-intervaltree=0.12=pl5321h503566f_5 \
    perl-socket=2.027=pl5321h5c03b87_6 \
    perl-storable=3.15=pl5321hb9d3cd8_2 \
    perl-sub-uplevel=0.2800=pl5321hb9d3cd8_0 \
    perl-test-deep=1.130=pl5321hd8ed1ab_0 \
    perl-test-differences=0.72=pl5321ha770c72_0 \
    perl-test-exception=0.43=pl5321hd8ed1ab_0 \
    perl-test-fatal=0.016=pl5321hdfd78af_0 \
    perl-test-most=0.38=pl5321hdfd78af_0 \
    perl-test-needs=0.002009=pl5321hd8ed1ab_0 \
    perl-test-requiresinternet=0.05=pl5321hdfd78af_1 \
    perl-test-warn=0.37=pl5321hd8ed1ab_0 \
    perl-test-warnings=0.031=pl5321ha770c72_0 \
    perl-text-diff=1.45=pl5321hd8ed1ab_0 \
    perl-time-hires=1.9764=pl5321h7b50bb2_6 \
    perl-time-local=1.35=pl5321hdfd78af_0 \
    perl-timedate=2.33=pl5321hdfd78af_2 \
    perl-tree-dag_node=1.35=pl5321hdfd78af_0 \
    perl-try-tiny=0.32=pl5321ha770c72_1 \
    perl-types-serialiser=1.01=pl5321hdfd78af_0 \
    perl-uri=5.17=pl5321ha770c72_0 \
    perl-url-encode=0.03=pl5321h9ee0642_1 \
    perl-www-robotrules=6.02=pl5321hdfd78af_4 \
    perl-xml-dom=1.46=pl5321hdfd78af_1 \
    perl-xml-dom-xpath=0.14=pl5321hdfd78af_2 \
    perl-xml-parser=2.44_01=pl5321hf2e2c51_1004 \
    perl-xml-regexp=0.04=pl5321hdfd78af_3 \
    perl-xml-xpathengine=0.14=pl5321hdfd78af_3 \
    pip=24.0=pyhd8ed1ab_0 \
    pixman=0.46.0=h29eaf8c_0 \
    poa=2.0=h7b50bb2_6 \
    prompt-toolkit=3.0.48=pyha770c72_0 \
    prompt_toolkit=3.0.48=hd8ed1ab_1 \
    provean=1.1.5=h503566f_3 \
    pthread-stubs=0.4=hb9d3cd8_1002 \
    pybiolib=1.2.735=pyhdfd78af_0 \
    pycparser=2.21=pyhd8ed1ab_0 \
    pycryptodome=3.15.0=py37h294ce10_0 \
    pygments=2.17.2=pyhd8ed1ab_0 \
    pyjwt=2.8.0=pyhd8ed1ab_1 \
    pynacl=1.5.0=py37h540881e_1 \
    pysam=0.16.0=py37ha9a96c6_0 \
    pysocks=1.7.1=py37h89c1867_5 \
    python=3.7.12=hf930737_100_cpython \
    python-dateutil=2.9.0=pyhd8ed1ab_0 \
    python_abi=3.7=4_cp37m \
    pytz=2024.2=pyhd8ed1ab_0 \
    pyvcf=0.6.8=py37hc8dfbb8_1002 \
    pywin32-on-windows=0.1.0=pyh1179c8e_3 \
    pyyaml=6.0=py37h540881e_4 \
    readline=8.2=h8c095d6_2 \
    requests=2.29.0=pyhd8ed1ab_0 \
    rich=13.8.1=pyhd8ed1ab_0 \
    samtools=1.21=h50ea8bc_0 \
    sed=4.9=h6688a6e_0 \
    setuptools=59.8.0=py37h89c1867_1 \
    six=1.16.0=pyh6c4a22f_0 \
    snowballstemmer=2.2.0=pyhd8ed1ab_0 \
    snpeff=5.2=hdfd78af_1 \
    snpsift=5.2=hdfd78af_0 \
    sqlite=3.46.0=h6d4b2fc_0 \
    sysroot_linux-64=2.17=h0157908_18 \
    t-coffee=12.00.7fb08c2=h26a2512_0 \
    tabix=1.11=hdfd78af_0 \
    tk=8.6.13=noxft_h4845f30_101 \
    tktable=2.10=h8bc8fbc_6 \
    toml=0.10.2=pyhd8ed1ab_0 \
    tomlkit=0.12.5=pyha770c72_0 \
    tqdm=4.67.1=pyhd8ed1ab_0 \
    tree=2.2.1=hb9d3cd8_0 \
    typing-extensions=4.7.1=hd8ed1ab_0 \
    typing_extensions=4.7.1=pyha770c72_0 \
    tzdata=2025b=h78e105d_0 \
    urllib3=1.26.19=pyhd8ed1ab_0 \
    viennarna=2.4.7=py27_2 \
    vine=5.1.0=pyhd8ed1ab_0 \
    wcwidth=0.2.10=pyhd8ed1ab_0 \
    websocket-client=1.6.1=pyhd8ed1ab_0 \
    werkzeug=2.2.3=pyhd8ed1ab_0 \
    wheel=0.42.0=pyhd8ed1ab_0 \
    xmltodict=0.14.2=pyhd8ed1ab_0 \
    xorg-fixesproto=5.0=hb9d3cd8_1003 \
    xorg-inputproto=2.3.2=hb9d3cd8_1003 \
    xorg-kbproto=1.0.7=hb9d3cd8_1003 \
    xorg-libice=1.1.2=hb9d3cd8_0 \
    xorg-libsm=1.2.6=he73a12e_0 \
    xorg-libx11=1.8.9=h8ee46fc_0 \
    xorg-libxau=1.0.12=hb9d3cd8_0 \
    xorg-libxdmcp=1.1.5=hb9d3cd8_0 \
    xorg-libxext=1.3.4=h0b41bf4_2 \
    xorg-libxfixes=5.0.3=h7f98852_1004 \
    xorg-libxi=1.7.10=h4bc722e_1 \
    xorg-libxrender=0.9.11=hd590300_0 \
    xorg-libxt=1.3.0=hd590300_1 \
    xorg-libxtst=1.2.5=h4bc722e_0 \
    xorg-recordproto=1.14.2=hb9d3cd8_1003 \
    xorg-renderproto=0.11.1=hb9d3cd8_1003 \
    xorg-xextproto=7.3.0=hb9d3cd8_1004 \
    xorg-xproto=7.0.31=hb9d3cd8_1008 \
    xz=5.2.6=h166bdaf_0 \
    xz-gpl-tools=5.8.1=hbcc6ac9_0 \
    xz-tools=5.8.1=hb9d3cd8_0 \
    yaml=0.2.5=h7f98852_2 \
    yq=3.4.3=pyhd8ed1ab_0 \
    zipp=3.15.0=pyhd8ed1ab_0 \
    zlib=1.2.13=h4ab18f5_6 \
    zstd=1.5.6=ha6fb4c9_0 \
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
