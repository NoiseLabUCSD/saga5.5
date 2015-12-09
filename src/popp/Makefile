# -----------------------------------------------------------------------------
# POPP library for SAGA
# -----------------------------------------------------------------------------
#
PROG   = empty
#FC     = ${FCC}
#DFLAGS = 
#FFLAGS = -O $(DFLAGS)
#CC     = cc
#CFLAGS = -C -O

MV     = mv
RM     = rm -f
RANLIB = ranlib
#INFO   = 
#INFO  = mailback

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
OBJ    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)libpopp.a 
# -------------------------------------------------------------------------

all:	$(XLIB) 
#$(PROG)

OBJO =	$(OBJ)attdep.o	$(OBJ)attenu.o   $(OBJ)bdmds.o   $(OBJ)datack.o  \
	$(OBJ)dopar.o    $(OBJ)error.o    $(OBJ)intin.o   $(OBJ)intout.o \
	$(OBJ)intpts.o   $(OBJ)itpdat.o   $(OBJ)kirkfb.o  $(OBJ)kirkff.o \
	$(OBJ)layer.o    $(OBJ)linear.o   $(OBJ)normal.o  $(OBJ)omodes.o \
	$(OBJ)output.o   $(OBJ)qtrial.o  $(OBJ)turnpt.o   $(OBJ)xmatch.o

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

$(OBJ)attdep.o:	attdep.f      
	$(FC) $(FFLAGS)  -c attdep.f
	@ $(MV) ${@F} $@

$(OBJ)attenu.o:	attenu.f      
	$(FC) $(FFLAGS)  -c attenu.f
	@ $(MV) ${@F} $@

$(OBJ)bdmds.o:	bdmds.f      
	$(FC) $(FFLAGS)  -c bdmds.f
	@ $(MV) ${@F} $@

$(OBJ).o:	.f      
	$(FC) $(FFLAGS)  -c .f
	@ $(MV) ${@F} $@

$(OBJ)datack.o:	datack.f      
	$(FC) $(FFLAGS)  -c datack.f
	@ $(MV) ${@F} $@

$(OBJ)dopar.o:	dopar.f      
	$(FC) $(FFLAGS)  -c dopar.f
	@ $(MV) ${@F} $@

$(OBJ)error.o:	error.f      
	$(FC) $(FFLAGS)  -c error.f
	@ $(MV) ${@F} $@

$(OBJ)intin.o:	intin.f      
	$(FC) $(FFLAGS)  -c intin.f
	@ $(MV) ${@F} $@

$(OBJ)intout.o:	intout.f      
	$(FC) $(FFLAGS)  -c intout.f
	@ $(MV) ${@F} $@

$(OBJ)intpts.o:	intpts.f      
	$(FC) $(FFLAGS)  -c intpts.f
	@ $(MV) ${@F} $@

$(OBJ)itpdat.o:	itpdat.f      
	$(FC) $(FFLAGS)  -c itpdat.f
	@ $(MV) ${@F} $@

$(OBJ)kirkfb.o:	kirkfb.f      
	$(FC) $(FFLAGS)  -c kirkfb.f
	@ $(MV) ${@F} $@

$(OBJ)kirkff.o:	kirkff.f      
	$(FC) $(FFLAGS)  -c kirkff.f
	@ $(MV) ${@F} $@

$(OBJ)layer.o:	layer.f      
	$(FC) $(FFLAGS)  -c layer.f
	@ $(MV) ${@F} $@

$(OBJ)linear.o:	linear.f      
	$(FC) $(FFLAGS)  -c linear.f
	@ $(MV) ${@F} $@

$(OBJ)normal.o:	normal.f      
	$(FC) $(FFLAGS)  -c normal.f
	@ $(MV) ${@F} $@
$(OBJ)omodes.o:	omodes.f      
	$(FC) $(FFLAGS)  -c omodes.f
	@ $(MV) ${@F} $@

$(OBJ)output.o:	output.f      
	$(FC) $(FFLAGS)  -c output.f
	@ $(MV) ${@F} $@
$(OBJ)qtrial.o:	qtrial.f      
	$(FC) $(FFLAGS)  -c qtrial.f
	@ $(MV) ${@F} $@

$(OBJ)turnpt.o:	turnpt.f      
	$(FC) $(FFLAGS)  -c turnpt.f
	@ $(MV) ${@F} $@

$(OBJ)xmatch.o:	xmatch.f      
	$(FC) $(FFLAGS)  -c xmatch.f
	@ $(MV) ${@F} $@

