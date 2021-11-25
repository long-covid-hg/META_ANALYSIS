FROM ubuntu:18.04

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install mercurial --yes
RUN apt-get update && apt-get install -qqy python3 python3-venv python3-pip && pip3 install gitpython numpy
RUN apt-get install  gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install bzip2 curl unzip gawk git  make gcc g++ libz-dev libbz2-dev liblzma-dev  littler libssl-dev --yes


#QCTOOL
RUN wget https://code.enkre.net/qctool/zip/v2.0.7/qctool.tgz && \
unzip qctool.tgz && \
cd qctool && ./waf-1.5.18 configure && ./waf-1.5.18 && cd ..
#move binaries
RUN cp qctool/build/release/qctool* /usr/local/bin/qctool
 
#BGEN
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && \
tar -xvzf bgen.tgz && mv bgen.tgz bgen && \
cd bgen && ./waf configure && ./waf && cd ..
# move binaries
RUN cp /bgen/build/apps/bgenix /usr/local/bin/bgenix && cp /bgen/build/apps/cat-bgen /usr/local/bin/cat-bgen

# SAMTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && \
tar -xvjf htslib-1.10.2.tar.bz2 && cd htslib-1.10.2 && ./configure && make && make install && cd .. && rm  htslib-1.*tar*

RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && tar -xvjf bcftools-1.10.2.tar.bz2 && cd bcftools-1.10.2 && \
./configure && make && make install && cd ..  && rm  bcftools-*bz2 && export BCFTOOLS_PLUGINS=~/bcftools-1.10.2/plugins/

#PLINK 1.9
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_* -d /usr/local/bin && rm plink*zip

#PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip


ADD Scripts ../Scripts
