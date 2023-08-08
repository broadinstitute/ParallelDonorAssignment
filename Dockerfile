FROM python:3.8
# GCLOUD

WORKDIR /build

# From https://github.com/getzlab/CGA_standard_image/blob/master/build.sh
RUN mkdir /gcsdk && \
  wget -O gcs.tgz https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-318.0.0-linux-x86_64.tar.gz && \
  tar xzf gcs.tgz -C /gcsdk && \
  /gcsdk/google-cloud-sdk/install.sh --usage-reporting false --path-update true --quiet && \
  ln -s /gcsdk/google-cloud-sdk/bin/* /usr/bin

# SAMTOOLS
RUN apt-get install -y libssl-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libncurses5-dev zlib1g-dev
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xf htslib-1.16.tar.bz2 && rm htslib-1.16.tar.bz2 && cd htslib-1.16 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean
RUN wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 && \
    tar -xf samtools-1.16.tar.bz2 && rm samtools-1.16.tar.bz2 && cd samtools-1.16 && \
    ./configure --with-htslib=system && make && make install && make clean
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    tar -xf bcftools-1.16.tar.bz2 && rm bcftools-1.16.tar.bz2 && cd bcftools-1.16 && \
    ./configure --with-htslib=system && make && make install && make clean

RUN apt install zlib1g

WORKDIR /app
ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY donor_assigment ./donor_assigment
