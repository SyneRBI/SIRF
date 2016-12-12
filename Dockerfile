# Copyright Casper da Costa-Luis Dec 2016
# Mozilla Public License 2.0 (https://mozilla.org/MPL/2.0/)

FROM ubuntu:16.04
MAINTAINER Casper da Costa-Luis <imaging@caspersci.uk.to>

ENV DEBIAN_FRONTEND noninteractive
# Set locale, suppress warnings
RUN localedef -v -c -i en_GB -f UTF-8 en_GB.UTF-8 || :

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y apt-utils

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y \
  bash-completion      \
  build-essential      \
  cmake                \
  make                 \
  automake autoconf    \
  sudo                 \
  man                  \
  git

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y \
  net-tools            \
  openssh-server       \
  vim

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y \
  python python-dev python-pip
RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y \
  gcc g++
RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y \
  libboost-dev

RUN sed -i -r \
  "s/^.*PermitRootLogin prohibit-password/PermitRootLogin without-password/" \
  /etc/ssh/sshd_config

ARG mainUser=caspersci

ADD passwd.txt /src/
RUN useradd --create-home --shell /bin/bash --groups sudo $mainUser \
  -p $(cat /src/passwd.txt | openssl passwd -1 -stdin)
RUN sudo rm -rf /src/

RUN echo "$mainUser ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers.d/$mainUser

USER $mainUser
ENV HOME /home/$mainUser
WORKDIR $HOME

RUN sudo -H pip install virtualenv
RUN virtualenv py2
ARG PIPINST="py2/bin/python -m pip install -U"
RUN $PIPINST pip
RUN $PIPINST cython
RUN $PIPINST docopt tqdm docutils
RUN $PIPINST numpy
RUN $PIPINST scipy
RUN $PIPINST matplotlib
RUN mkdir .local


RUN git clone -b matlab https://github.com/CCPPETMR/swig
RUN sudo apt-get update && sudo apt-get dist-upgrade -y && sudo apt-get install -y \
  libpcre3-dev libpcre++-dev byacc
RUN . py2/bin/activate && cd swig && ./autogen.sh && \
  ./configure --enable-cpp11-testing --prefix=$HOME/.local && \
  make -j8 && make install

RUN git clone https://github.com/CCPPETMR/STIR
RUN . py2/bin/activate && cd STIR && mkdir build && cd build && \
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DBUILD_DOCUMENTATION=0 \
    -DBUILD_SWIG_MATLAB=0 \
    -DBUILD_SWIG_OCTAVE=0 \
    -DBUILD_SWIG_PYTHON=1 \
    -DCMAKE_BUILD_TYPE=Release \
    -DDISABLE_AVW=1 \
    -DDISABLE_CERN_ROOT_SUPPORT=1 \
    -DDISABLE_ITK=1 \
    -DDISABLE_LLN_MATRIX=1 \
    -DDISABLE_RDF=1 \
    -DSTIR_MPI=0 \
    -DSTIR_OPENMP=1 && \
  make -j8 && make install
# TODO: git clone -b casper https://github.com/CCPPETMR/gadgetron

# RUN git clone --recursive https://github.com/CCPPETMR/SIRF
ADD . $HOME/SIRF
RUN sudo chown -R $mainUser SIRF
RUN . py2/bin/activate && cd SIRF && \
  mkdir dockerbuild && cd dockerbuild && \
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local \
  -DSTIR_LIBRARY_DIR=$HOME/.local/lib && \
  make -j8 && make install
RUN cp SIRF/.bash_aliases $HOME/

CMD cd $HOME && /bin/bash
