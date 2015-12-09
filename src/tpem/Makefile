# -----------------------------------------------------------------------------
# TPEM library for SAGA
# -----------------------------------------------------------------------------
#
PROG   = empty
#FC     = ${FCC}
#DFLAGS = 
#FFLAGS = $(FFLAGS_warn)
#CC     = cc
#CFLAGS = -C -O
MV= mv

RM     = rm -f
RANLIB = ranlib

$(OBJ)%.o : %.f  		
	$(FC) $(FFLAGS) -I $(MOD)  -c $*.f  $(OBJ_opt) $@
$(OBJ)%.o : %.f90	
	$(FC) $(FFLAG) -I $(MOD) -module $(MOD) -c $*.f90  $(OBJ_opt) $@

TARGET = ..
DFLAGS = 
OBJ    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)tpem.a 
MOD    = ../obj/${HOSTTYPE}-${FORTRAN}/
# -------------------------------------------------------------------------

all:	$(XLIB) 
# $(PROG)
MODS = 
OBJO = 	$(OBJ)tpemsubs.o $(OBJ)tpemsubs2.o $(OBJ)sinfft.o $(OBJ)seofgen.o $(OBJ)refparam.o $(OBJ)trfhead.o

PRJS =  $(OBJO:.o=.prj)

$(XLIB):$(MODS) $(OBJO)
	ar rvu $(XLIB) $(OBJO)
	$(RANLIB) $(XLIB)
#	$(MV) *.o $(LIB)

$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) $(PROG).o $(XLIB) 
	strip $(PROG)

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)

dircmp:
	@echo 'Tracing == ' `pwd`
	@trace

### Dependencies:
$(OBJ)sinfft.o:	sinfft.f    fftsiz.inc
	$(FC) $(FFLAGS)  -c sinfft.f
	@ $(MV) ${@F} $@

$(OBJ)tpem.o:		tpem.f      tpem.inc
	$(FC) $(FFLAGS)  -c tpem.f
	@ $(MV) ${@F} $@

$(OBJ)tpemsubs.o:	tpemsubs.f  tpem.inc
	$(FC) $(FFLAGS)  -c tpemsubs.f
	@ $(MV) ${@F} $@

$(OBJ)tpemsubs2.o:	tpemsubs2.f tpem.inc
	$(FC) $(FFLAGS)  -c tpemsubs2.f
	@ $(MV) ${@F} $@

$(OBJ)seofgen.o:	seofgen.f
	$(FC) $(FFLAGS)  -c 	seofgen.f
	@ $(MV) ${@F} $@

$(OBJ)refparam.o:	refparam.f tpem.h tpem.inc
	$(FC) $(FFLAGS)  -c refparam.f
	@ $(MV) ${@F} $@

$(OBJ)trfhead.o:	trfhead.f 
	$(FC) $(FFLAGS)  -c trfhead.f
	@ $(MV) ${@F} $@

