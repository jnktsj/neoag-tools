# Dockerfile for neoag-tools

FROM ubuntu:22.04

LABEL maintainer="Junko Tsuji <jtsuji@broadinstitute.org>"

WORKDIR /root/

RUN apt-get update && \
    # apt-get upgrade -y && \
    apt-get install -y \
      # build-essential \
      # software-properties-common \
      # libz-dev \
      # libncurses-dev \
      # liblzma-dev \
      # libbz2-dev \
      # openjdk-8-jdk \
      # libcurl4-openssl-dev \
      python3 \
      python3-pip \
      # wget \
      # curl \
      # git \
      # unzip && \
      && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

# # Clone neoantigen phasing and translation tools
# RUN git clone https://github.com/jnktsj/neoag-tools.git && \
#     mv neoag-tools main && pip3 install -r /root/main/requirements.txt

COPY . /root/neoag-tools
RUN mv neoag-tools main && pip3 install -r /root/main/requirements.txt

# Deleting unneeded caches
RUN rm -rf /var/lib/apt/lists/*

# Copy Dockerfile to workdir
COPY Dockerfile /root/
