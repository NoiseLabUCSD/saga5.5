# -----------------------------------------------------------------------------
# RAM library for SAGA
# -----------------------------------------------------------------------------

PROG   = raminit
#FC     = ${FCC}
#DFLAGS = 
#FFLAGS = -O $(DFLAGS)
#CC     = cc
#CFLAGS = -C -O
MV= mv
RM     = rm -f
#MV     = mv
#RANLIB = ranlib
#INFO   = 

# Suffix rules
# ------------
FFLAGS= -C
FFLAGS= 
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.f.o:;		$(FC) $(FFLAGS) -c  $*.f  $(OBJ_opt)  $@
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F -o $@
.f.prj:;	ftnchek -project -noextern -library $<

TARGET = ..
DFLAGS = 
OBJ    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)ramgeo.a 
BIN    = 

# -------------------------------------------------------------------------

all:	$(XLIB) 

#$(PROG)
OBJO =	$(OBJ)ramgeo1.5.o 

$(XLIB):$(OBJO)
	ar rvu $(XLIB) $(OBJO)
	$(RANLIB) $(XLIB)

$(PROG):$(PROG).o $(XLIB)
	$(FC)  -o   $(PROG) $(OBJO) $(XLIB) 

clean: 
	$(RM) *.a *~ *.*~ core *.o *.prj .f $(PROG)

### Dependencies:
$(OBJ)ram1.5.o:	ram1.5.f ram.h
	$(FC) $(FFLAGS)  -c ram1.5.f
	@ $(MV) ${@F} $@
$(OBJ)ramgeo1.5.o:	ramgeo1.5.f ram.h
	$(FC) $(FFLAGS)  -c ramgeo1.5.f
	@ $(MV) ${@F} $@
$(OBJ)main_ram.o:	main_ram.f ram.h
	$(FC) $(FFLAGS)  -c ram_main.f
	@ $(MV) ${@F} $@
$(OBJ)raminit.o:	raminit.f ram.h comramgeo.h
	$(FC) -c raminit.f
	@ $(MV) ${@F} $@




