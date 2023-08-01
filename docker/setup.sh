#!/bin/bash
# 
# William Jenkins
# Scripps Institution of Oceanography
# 21 November 2022
# 
# Script works in concert with Dockerfile to add SAGA configuration to .bashrc.

( echo "" ; echo "# >> Begin SAGA Configuration <<" ) >> ~/.bashrc
echo "export FORTRAN=$FORTRAN" >> ~/.bashrc
echo "export HOSTTYPE=$HOSTTYPE" >> ~/.bashrc
echo "export PATH=$SAGADIR/bin:$SAGADIR/bin/$HOSTTYPE-$FORTRAN:$PATH" >> ~/.bashrc
echo "# >> End SAGA Configuration <<" >> ~/.bashrc
