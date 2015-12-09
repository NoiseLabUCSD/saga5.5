#!/bin/sh
# Create machine dependent object and binary directories for OASES3D/SCATT
#
  if [ ! -d ./obj/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./obj/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating ${SAGADIR}/src/obj/${HOSTTYPE}-${FORTRAN} <<<'
  fi

  if [ ! -d ../bin/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ../bin/${HOSTTYPE}-${FORTRAN}
    cp -f  ../bin/script/* ../bin/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/bin/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./snap/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./snap/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/snap/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./snaprd/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./snaprd/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/snaprd/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./oases/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./oases/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/oases/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./tpem/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./tpem/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/tpem/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./popp/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./popp/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/popp/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./prosim/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./prosim/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/prosim/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./cprosim/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./cprosim/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/cprosim/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./gama/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./gama/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/cprosim/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./orca90/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./orca90/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/cprosim/${HOSTTYPE}-${FORTRAN} <<<'
  fi
  if [ ! -d ./misc/${HOSTTYPE}-${FORTRAN} ]; then
    mkdir ./misc/${HOSTTYPE}-${FORTRAN}
    echo '>>> Creating  ${SAGADIR}/misc/${HOSTTYPE}-${FORTRAN} <<<'
  fi
