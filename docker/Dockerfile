FROM registry.access.redhat.com/ubi7:7.9
LABEL author="William Jenkins <wjenkins@ucsd.edu>"
RUN yum update -y
RUN yum install -y gcc-gfortran gdb git make vim
RUN git clone https://github.com/gerstoft/saga.git
ENV SAGADIR=/saga
ENV HOSTTYPE=x86_64
ENV FORTRAN=gfortran
COPY setup.sh $SAGADIR
WORKDIR $SAGADIR
RUN bash ./setup.sh
WORKDIR $SAGADIR/src
RUN make clean && make all
WORKDIR $SAGADIR