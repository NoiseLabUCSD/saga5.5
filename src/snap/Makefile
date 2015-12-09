# -----------------------------------------------------------------------------
# SNAP library for SAGA
# -----------------------------------------------------------------------------
#
PROG   = empty 
RM     = rm -f
#MV     = mv
RANLIB = ranlib
#INFO   = 

# Suffix rules
# ------------
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.f.o:;		$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@; mv $*.o $(OBJ)
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
.f.prj:;	ftnchek -project -noextern -library $<

TARGET = ..
DFLAGS = 
OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)snap.a 
MV = mv
# -------------------------------------------------------------------------

all:	$(XLIB) 
#$(PROG)

OBJO =	$(OBJ)accrcy.o    $(OBJ)brent.o   $(OBJ)charac.o   $(OBJ)convrg.o \
	$(OBJ)czs.o    	  $(OBJ)eigvec.o  $(OBJ)field.o    $(OBJ)init.o  \
	$(OBJ)isobr.o     $(OBJ)isoint.o  $(OBJ)lagrange.o $(OBJ)modfun.o \
	$(OBJ)newton.o    $(OBJ)porter.o  $(OBJ)pressr.o   $(OBJ)refl1.o  \
	$(OBJ)sturm.o     $(OBJ)velocity.o $(OBJ)snapsub.o $(OBJ)getprofiles.o

PRJS =  $(OBJO:.o=.prj)

$(XLIB):$(OBJO)
	ar rvu $(XLIB) $(OBJO)
	$(RANLIB) $(XLIB)

$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) $(PROG).o $(XLIB) 
	strip $(PROG)

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)

dircmp:
	@echo 'Tracing == ' `pwd`
	@trace

### Dependencies:

$(OBJ)accrcy.o: accrcy.f
	$(FC)  $(FFLAGS)  -c accrcy.f
	@ $(MV) ${@F} $@

$(OBJ)brent.o: brent.f
	$(FC)  $(FFLAGS)  -c brent.f
	@ $(MV) ${@F} $@

$(OBJ)charac.o: charac.f
	$(FC)  $(FFLAGS)  -c charac.f
	@ $(MV) ${@F} $@

$(OBJ)convrg.o: convrg.f
	$(FC)  $(FFLAGS)  -c convrg.f
	@ $(MV) ${@F} $@

$(OBJ)czs.o: czs.f
	$(FC)  $(FFLAGS)  -c czs.f
	@ $(MV) ${@F} $@

$(OBJ)eigvec.o: eigvec.f
	$(FC)  $(FFLAGS)  -c eigvec.f
	@ $(MV) ${@F} $@

$(OBJ)field.o: field.f
	$(FC)  $(FFLAGS)  -c field.f
	@ $(MV) ${@F} $@

$(OBJ)init.o: init.f
	$(FC)  $(FFLAGS)  -c init.f
	@ $(MV) ${@F} $@

$(OBJ)isobr.o: isobr.f
	$(FC)  $(FFLAGS)  -c isobr.f
	@ $(MV) ${@F} $@

$(OBJ)isoint.o: isoint.f
	$(FC)  $(FFLAGS)  -c isoint.f
	@ $(MV) ${@F} $@

$(OBJ)lagrange.o: lagrange.f
	$(FC)  $(FFLAGS)  -c lagrange.f
	@ $(MV) ${@F} $@

$(OBJ)modfun.o: modfun.f
	$(FC)  $(FFLAGS)  -c modfun.f
	@ $(MV) ${@F} $@

$(OBJ)newton.o: newton.f
	$(FC)  $(FFLAGS)  -c newton.f
	@ $(MV) ${@F} $@

$(OBJ)porter.o: porter.f
	$(FC)  $(FFLAGS)  -c porter.f
	@ $(MV) ${@F} $@

$(OBJ)pressr.o: pressr.f
	$(FC)  $(FFLAGS)  -c pressr.f
	@ $(MV) ${@F} $@

$(OBJ)refl1.o: refl1.f
	$(FC)  $(FFLAGS)  -c refl1.f
	@ $(MV) ${@F} $@

$(OBJ)sturm.o: sturm.f
	$(FC)  $(FFLAGS)  -c sturm.f
	@ $(MV) ${@F} $@

$(OBJ)velocity.o: velocity.f
	$(FC)  $(FFLAGS)  -c velocity.f
	@ $(MV) ${@F} $@

$(OBJ)snapsub.o: snapsub.f
	$(FC)  $(FFLAGS)  -c snapsub.f
	@ $(MV) ${@F} $@

$(OBJ)getprofiles.o: getprofiles.f90
	$(FC90)  $(FFLAGS90)  -c getprofiles.f90 
	@ $(MV) ${@F} $@

