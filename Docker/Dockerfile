FROM ubuntu:20.04
MAINTAINER bhaas@broadinstitute.org

ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev


RUN pip install numpy


RUN curl -L https://cpanmin.us | perl - App::cpanminus

## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

ENV DATA /usr/local/data
RUN mkdir $DATA


## perl lib installations

RUN cpanm install PerlIO::gzip
RUN cpanm install Set::IntervalTree
RUN cpanm install DB_File
RUN cpanm install URI::Escape
RUN cpanm install Carp::Assert
RUN cpanm install JSON::XS.pm


######################
## Tool installations:
######################



## Bowtie2
WORKDIR $SRC
ENV BOWTIE2_VERSION 2.4.4
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip/download -O bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-${BOWTIE2_VERSION}-linux-x86_64

## Samtools
WORKDIR $SRC
ENV SAMTOOLS_VERSION 1.9
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install


## Seqtk
WORKDIR $SRC
RUN git clone https://github.com/lh3/seqtk.git && \
             cd seqtk &&  make && cp seqtk $BIN/


########
## Minimap2

ENV MM2_ver=2.28
WORKDIR $SRC
RUN curl -L https://github.com/lh3/minimap2/releases/download/v${MM2_ver}/minimap2-${MM2_ver}_x64-linux.tar.bz2 | tar -jxvf - && \
    mv ./minimap2-${MM2_ver}_x64-linux/minimap2 $BIN/


##############
## STAR

ENV STAR_VERSION=2.7.11a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar xvf htslib-1.17.tar.bz2 && \
    cd htslib-1.17 && \
    ./configure && make && make install


RUN ln -s /usr/bin/python3 /usr/bin/python

RUN pip install pysam

RUN apt-get install -y git && apt-get clean

WORKDIR $SRC


# K8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
    cp k8-0.2.4/k8-`uname -s` $BIN/k8
    
RUN pip install pandas

RUN pip install requests igv-reports==1.11.0

ENV ctat_LR_fusion_version 1.0.2-postdev1
ENV CTAT_LR_FUSION_CO 38d6b8c4b6db259f3e6d6a083801bc323b677526

RUN git clone --recursive https://github.com/TrinityCTAT/CTAT-LR-fusion.git  && \
	cd CTAT-LR-fusion && \
    git checkout ${CTAT_LR_FUSION_CO} && \
    git submodule init && git submodule update && \
    git submodule foreach --recursive git submodule init && \
    git submodule foreach --recursive git submodule update && \
    make


RUN mv $SRC/CTAT-LR-fusion/* $BIN/

