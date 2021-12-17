## with 'DOCKER_BUILDKIT=1 docker build --secret id=mysecret,src=mysecret.txt ./', where mysecret.txt is the 

# Start from the NVIDIA CUDA Docker image built on top of Ubuntu 20
FROM nvidia/cuda:11.4.1-runtime-ubuntu20.04

# basic setup -- update the instance
# ################################################
RUN apt-get update
RUN apt-get install -y apt-transport-https
RUN apt install -y python3.8
RUN apt install -y python3-pip python-is-python3

RUN mkdir -p /usr/share/man/man1

# add basic tools to the instance
# ################################################
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
RUN apt-get install -y wget unzip net-tools emboss gnupg2 bash make gcc g++ zlib1g-dev libboost-all-dev curl libhdf5-serial-dev libssl-dev git build-essential bzip2 g++ zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev libcurl4-gnutls-dev libssl-dev curl make cmake wget python3-all-dev
# libcurl3 libhdf5-cpp-11 libssl1.0.0

RUN mkdir /ont/
RUN python3 --version 1>&2 
RUN pip install nanofilt
RUN pip install pycoqc

# Add the oxford Nanopore repository to the container
# ################################################

RUN apt-get update
RUN apt-get install wget lsb-release
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt focal-stable non-free" |  tee /etc/apt/sources.list.d/nanoporetech.sources.list
RUN apt-get update
RUN apt update
RUN apt install -y ont-guppy

# Add Rerio models
# ################################################
RUN git clone https://github.com/nanoporetech/rerio /ont/rerio
RUN /ont/rerio/download_model.py /ont/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001

# Install porechop
# ################################################
RUN git clone https://github.com/rrwick/Porechop.git /ont/Porechop 

COPY adapters.py /ont/Porechop/porechop/

RUN cd /ont/Porechop && \
    python3 setup.py install && \
    porechop -h

# Install Canu
# ################################################
RUN cd /ont/ && \
    wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz && \
    tar xf canu-2.2.Linux-amd64.tar.xz
ENV PATH="/ont/canu-2.2/bin/:${PATH}"

# Install Racon
RUN git clone --recursive https://github.com/lbcb-sci/racon.git /ont/racon && \
    cd /ont/racon && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make install


# install miniasm and minimap2 
# ################################################
RUN cd /ont/ && git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
RUN cd /ont/ && git clone https://github.com/lh3/miniasm  && (cd miniasm  && make)
ENV PATH="/ont/minimap2/:${PATH}"
ENV PATH="/ont/miniasm/:${PATH}"

# install nextpolish
# ################################################
RUN apt install -y libncurses5-dev libncursesw5-dev libbz2-dev libclang-dev liblzma-dev
RUN pip install psutil
RUN wget --quiet https://github.com/Nextomics/NextPolish/releases/download/v1.4.0/NextPolish.tgz -O /ont/NextPolish.tgz && \
    cd /ont/ && \
    tar -vxzf NextPolish.tgz && cd NextPolish && make
ENV PATH="/ont/NextPolish/:${PATH}"

# Install fast5_subset
# ################################################
RUN pip install ont-fast5-api

# java 1.8 is needed by canu
# ################################################
RUN apt update
RUN apt install -y openjdk-8-jdk

# Medaka for sequence correction - lots of deps
# here that are a bit of a nightmare to work out
# ################################################
RUN pip install medaka

RUN wget -O /ont/bcftools-1.11.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
    cd /ont/ && tar xf bcftools-1.11.tar.bz2 && cd bcftools-1.11 && make && make install
RUN wget -O /ont/htslib-1.11.tar.bz2 https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    cd /ont/ && tar xf htslib-1.11.tar.bz2 && cd htslib-1.11 && make && make install

RUN wget -O /ont/samtools-1.11.tar.bz2 https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    cd /ont/ && tar xf samtools-1.11.tar.bz2 && cd samtools-1.11 && make && make install

RUN pip install numpy==1.19.5

RUN pip install pydivsufsort

RUN pip install megalodon

RUN mkdir /og_methyl/

COPY ont-guppy_4.4.2_linux64.tar.gz /og_methyl/

RUN cd /og_methyl/ && tar -xf ont-guppy_4.4.2_linux64.tar.gz 

RUN git clone https://github.com/novoalab/modPhred /og_methyl/modPhred && \
    pip install matplotlib pysam pandas seaborn && \
    pip install pyguppyclient==0.0.9

RUN pip install ont-pyguppy-client-lib==5.0.11

# cp the github pipeline files 
# ################################################
RUN mkdir /plasmidseq/
RUN mkdir /plasmid_test89

RUN --mount=type=secret,id=mysecret git clone https://femiliani:$(cat /run/secrets/mysecret)@github.com/mckennalab/PlasmidSeq.git /plasmidseq/


## setup the annotation tool
RUN pip install typing-extensions==3.10.0.0
RUN wget https://github.com/barricklab/pLannotate/archive/refs/tags/v1.0.3.tar.gz -O /plasmidseq/v1.0.3.tar.gz && \
    cd /plasmidseq/ && tar -xzf v1.0.3.tar.gz && mv pLannotate-1.0.3 pLannotate && cd pLannotate && python setup.py install
RUN wget https://github.com/barricklab/pLannotate/releases/download/v1.0.3/BLAST_dbs.tar.gz -O /plasmidseq/pLannotate/BLAST_dbs.tar.gz && \
    cd /plasmidseq/pLannotate/ && tar -xzf BLAST_dbs.tar.gz

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz -O /plasmidseq/ncbi-blast-2.12.0+-x64-linux.tar.gz && \
    cd /plasmidseq/ && tar xf ncbi-blast-2.12.0+-x64-linux.tar.gz
    
ENV PATH="/plasmidseq/ncbi-blast-2.12.0+/bin/:${PATH}"

RUN apt install infernal

RUN mkdir /plasmidseq/diamond

RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.0.12/diamond-linux64.tar.gz -O /plasmidseq/diamond/diamond-linux64.tar.gz

RUN cd /plasmidseq/diamond/ && tar xf diamond-linux64.tar.gz 

ENV PATH="/plasmidseq/diamond/:${PATH}"

# add Flye to the docker container
RUN wget https://github.com/fenderglass/Flye/archive/refs/tags/2.9.tar.gz -O /plasmidseq/Flye-2.9.tar.gz && cd /plasmidseq/ && tar xf Flye-2.9.tar.gz && cd /plasmidseq/Flye-2.9 && make
ENV PATH="/plasmidseq/Flye-2.9/bin/:${PATH}"

RUN pip install pysam
RUN pip install paralleltask