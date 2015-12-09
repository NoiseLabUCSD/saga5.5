# -----------------------------------------------------------------------------
# SLAVESNAP library for SAGA
# -----------------------------------------------------------------------------
#
PROG   = slave 
FC     = ${FCC}
DFLAGS = 
FFLAGS = -O $(DFLAGS)
CC     = cc
CFLAGS = -C -O

RM     = rm -f
MV     = mv
RANLIB = ranlib
INFO   = 

# Suffix rules
# ------------
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.f.o:;		$(FC) $(FFLAGS) -c $*.f $@
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
.f.prj:;	ftnchek -project -noextern -library $<

TARGET = ..
DFLAGS = 
OBJ    = $(TARGET)/obj/
LIB    = $(TARGET)/obj/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)snap.a 
# -------------------------------------------------------------------------

all:	$(PROG)

PRJS =  $(OBJO:.o=.prj)

$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) slavesnap.o $(XLIB) 
	strip $(PROG)

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)

dircmp:
	@echo 'Tracing == ' `pwd`
	@trace

### Dependencies:



