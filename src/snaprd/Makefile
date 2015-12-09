# -----------------------------------------------------------------------------
# SNAPRD library for SAGA
# -----------------------------------------------------------------------------
#
PROG   =  empty
#FC     = ${FCC}
#DFLAGS = 
#FFLAGS = -O $(DFLAGS)
#CC     = cc
#CFLAGS = -C -O

RM     = rm -f
#MV     = mv
RANLIB = ranlib
#INFO   = 

# Suffix rules
# ------------
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.f.o:;		$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
.f.prj:;	ftnchek -project -noextern -library $<

TARGET = ..
DFLAGS = 
OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)snaprd.a 
MV = mv
# -------------------------------------------------------------------------

all:	$(XLIB)
# $(PROG)

OBJO =	$(OBJ)cohrnt.o $(OBJ)exchng.o $(OBJ)file.o   $(OBJ)modfun.o\
	$(OBJ)range.o  $(OBJ)refl1.o  $(OBJ)revise.o $(OBJ)rfform.o\
	$(OBJ)ssnap2.o $(OBJ)ssnap3.o $(OBJ)tlran.o  $(OBJ)tlranpg.o\
	$(OBJ)snaprd.o

PRJS =  $(OBJO:.o=.prj)

$(XLIB):$(OBJO)
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


$(OBJ)cohrnt.o: cohrnt.f
	$(FC)  $(FFLAGS)  -c cohrnt.f
	@ $(MV) ${@F} $@


$(OBJ)exchng.o: exchng.f
	$(FC)  $(FFLAGS)  -c exchng.f
	@ $(MV) ${@F} $@


$(OBJ)file.o: file.f
	$(FC)  $(FFLAGS)  -c file.f
	@ $(MV) ${@F} $@

$(OBJ)modfun.o: modfun.f
	$(FC)  $(FFLAGS)  -c modfun.f
	@ $(MV) ${@F} $@

$(OBJ)range.o: range.f
	$(FC)  $(FFLAGS)  -c range.f
	@ $(MV) ${@F} $@

$(OBJ)refl1.o: refl1.f
	$(FC)  $(FFLAGS)  -c refl1.f
	@ $(MV) ${@F} $@


$(OBJ)revise.o: revise.f
	$(FC)  $(FFLAGS)  -c revise.f
	@ $(MV) ${@F} $@

$(OBJ)rfform.o: rfform.f
	$(FC)  $(FFLAGS)  -c rfform.f
	@ $(MV) ${@F} $@

$(OBJ)ssnap2.o: ssnap2.f
	$(FC)  $(FFLAGS)  -c ssnap2.f
	@ $(MV) ${@F} $@

$(OBJ)ssnap3.o: ssnap3.f
	$(FC)  $(FFLAGS)  -c ssnap3.f
	@ $(MV) ${@F} $@

$(OBJ)tlran.o: tlran.f
	$(FC)  $(FFLAGS)  -c tlran.f
	@ $(MV) ${@F} $@

$(OBJ)tlranpg.o: tlranpg.f
	$(FC)  $(FFLAGS)  -c tlranpg.f
	@ $(MV) ${@F} $@

$(OBJ)snaprd.o: snaprd.f
	$(FC)  $(FFLAGS)  -c snaprd.f
	@ $(MV) ${@F} $@
