# ----------------------------------------------------------------------------
# OASES library for SAGA
# ----------------------------------------------------------------------------
#
PROG   = nothing
#FC     = ${FC}
#DFLAGS = 
#FFLAGS = -O $(DFLAGS)
#CC     = cc
#CFLAGS = -C -O
#FFLAGS_oasd=
#FFLAGS_csvdc=
RM     = rm -f
MV     = mv -f
RANLIB = ranlib
#INFO   = 

# Suffix rules
# ------------
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.f.o:;		$(FC) $(FFLAGS) -c $*.f $(OBJ_opt) $@
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
.f.prj:;	ftnchek -project -noextern -library $<

LIB    = ../obj/${HOSTTYPE}-${FORTRAN}/
BIN    = ../bin/
LIBA   = $(LIB)apmath.a 
LIBO   = $(LIB)oaslib.a
XLIB   = $(LIBA) $(LIBO)
OBJDIR=${HOSTTYPE}-${FORTRAN}/
# -------------------------------------------------------------------------

all:	$(XLIB)	


OBJM =	$(OBJDIR)apmath.o  $(OBJDIR)linpack.o $(OBJDIR)csvdc.o $(OBJDIR)oasqvms11.o  $(OBJDIR)oashun11.o  $(OBJDIR)cbgemrva.o $(OBJDIR)indexx.o $(OBJDIR)recipes.o $(OBJDIR)ran2.o

OBJO = 	$(OBJDIR)oasavms12.o  $(OBJDIR)oasbvms10.o $(OBJDIR)oascvms10.o  $(OBJDIR)oasdvms11.o  $(OBJDIR)oasevms12.o  \
	$(OBJDIR)oasfvms11.o $(OBJDIR)oasgvms11.o  $(OBJDIR)oasivms12.o  $(OBJDIR)oasjvms11.o  $(OBJDIR)oaskvms11.o  \
	$(OBJDIR)oasmvms13.o  $(OBJDIR)oasnvms13.o

PRJS =  $(OBJO:.o=.prj)

$(LIBA):$(OBJM)
	ar rvu $(LIBA) $(OBJM)
	$(RANLIB) $(LIBA)

$(LIBO):$(OBJO)
	ar rvu $(LIBO) $(OBJO)
	$(RANLIB) $(LIBO)
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
$(OBJDIR)oasavms12.o: oasavms12.f compar.f comnp.f comnla.f comnrd.f 
	$(FC) $(FFLAGS) -c oasavms12.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasbvms10.o: oasbvms10.f compar.f comnp.f comnla.f comnrd.f 
	$(FC) $(FFLAGS) -c oasbvms10.f
	@ $(MV) ${@F} $@
$(OBJDIR)oascvms10.o: oascvms10.f compar.f comnp.f comnla.f comnrd.f 
	$(FC) $(FFLAGS) -c oascvms10.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasdvms11.o: oasdvms11.f compar.f comnp.f comnla.f comnrd.f comti.f complo.f
	$(FC) $(FFLAGS_oasd) -c oasdvms11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasevms12.o: oasevms12.f compar.f comnp.f comnla.f comnrd.f
	$(FC) $(FFLAGS) -c oasevms12.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasfvms11.o: oasfvms11.f compar.f comnp.f comnla.f comnrd.f comfip.f complo.f
	$(FC) $(FFLAGS) -c oasfvms11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasgvms11.o: oasgvms11.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS) -c oasgvms11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oashun11.o: oashun11.f
	$(FC) $(WFLAGS) $(FFLAGS) $(FFLAG1) -c oashun11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasivms12.o: oasivms12.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS) -c oasivms12.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasjvms11.o: oasjvms11.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS) -c oasjvms11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oaskvms11.o: oaskvms11.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS_oasd) -c oaskvms11.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasmvms13.o: oasmvms13.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS) -c oasmvms13.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasnvms13.o: oasnvms13.f compar.f comnp.f comnla.f comnrd.f complo.f
	$(FC) $(FFLAGS) -c oasnvms13.f
	@ $(MV) ${@F} $@

$(OBJDIR)apmath.o: apmath.f
	$(FC)  $(FFLAGS)  -c apmath.f
	$(MV) ${@F} $@

$(OBJDIR)recipes.o: recipes.f
	$(FC)  $(FFLAGS_oasd)  -c recipes.f
#	$(FC)  $(FFLAGS)  -c recipes.f
	$(MV) ${@F} $@

$(OBJDIR)linpack.o: linpack.f
	$(FC) $(FFLAGS)  -c linpack.f
	@ $(MV) ${@F} $@
$(OBJDIR)csvdc.o: csvdc.f
	$(FC) $(FFLAGS_csvdc)  -c csvdc.f
	@ $(MV) ${@F} $@

$(OBJDIR)oasqvms11.o: oasqvms11.f
	$(FC) $(FFLAGS)  -c oasqvms11.f
	@ $(MV) ${@F} $@
$(OBJDIR)cbgemrva.o: cbgemrva.f
	$(FC) $(FFLAGS)  -c cbgemrva.f
	@ $(MV) ${@F} $@
$(OBJDIR)indexx.o: indexx.f
	$(FC) $(FFLAGS)  -c indexx.f
	@ $(MV) ${@F} $@
$(OBJDIR)ran2.o: ran2.f
	$(FC) $(FFLAGS)  -c ran2.f
	@ $(MV) ${@F} $@


