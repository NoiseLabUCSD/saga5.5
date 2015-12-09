# -----------------------------------------------------------------------------
# GAMA library for SAGA
# -----------------------------------------------------------------------------
#.SUFFIXES :
#.SUFFIXES : .f .c .o .prj
#.f.o:;		$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@; mv $*.o $(OBJ)
#.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
#.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
#.f.prj:;	ftnchek -project -noextern -library $<
PROG  = gama
TARGET = ..
DFLAGS = 
OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)libgama.a 
MV = mv
.SUFFIXES :
$(OBJ)%.o : %.f  		
	$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@
$(OBJ)%.o : %.f90	
	$(FC90) $(FFLAG90) -c $*.f90  $(OBJ_opt) $@


OBJS =   \
 $(OBJ)acouspl.o $(OBJ)addlay.o $(OBJ)airy.o $(OBJ)axlabel.o $(OBJ)bargr.o \
 $(OBJ)bd2calc.o $(OBJ)bdfill.o $(OBJ)bdfit.o $(OBJ)bdparsh.o $(OBJ)bdsort.o\
 $(OBJ)besl0.o $(OBJ)botlos.o $(OBJ)bsvpfit.o $(OBJ)cdcwine.o $(OBJ)chdisc.o\
 $(OBJ)chpeak.o $(OBJ)ckdup.o $(OBJ)cperr.o $(OBJ)cubext.o $(OBJ)cubroot.o\
 $(OBJ)dawrite.o $(OBJ)delta.o $(OBJ)discpts.o $(OBJ)divrat.o $(OBJ)dopshif.o\
 $(OBJ)dr2calc.o $(OBJ)dr2lay.o $(OBJ)drswich.o $(OBJ)dubloma.o $(OBJ)eigend.o\
 $(OBJ)eigini.o $(OBJ)eigwrit.o $(OBJ)expo.o $(OBJ)fft_tools.o $(OBJ)found.o\
 $(OBJ)geom.o $(OBJ)saga_getdat1.o $(OBJ)saga_getdat2.o $(OBJ)saga_getdat3.o\
 $(OBJ)getloss.o $(OBJ)saga_getsvp.o $(OBJ)horcaus.o $(OBJ)interp.o $(OBJ)laycrit.o\
 $(OBJ)layers.o $(OBJ)layparm.o $(OBJ)linefit.o $(OBJ)lnfit.o $(OBJ)loma.o $(OBJ)maxup.o\
 $(OBJ)mem_lim.o $(OBJ)mname.o $(OBJ)morelay.o $(OBJ)newrdep.o $(OBJ)newsdep.o $(OBJ)ocean.o\
 $(OBJ)pasc_make.o $(OBJ)pathfil.o $(OBJ)plcalls.o $(OBJ)plotpic.o $(OBJ)plplot.o\
 $(OBJ)polfit.o $(OBJ)pr.o $(OBJ)prhead.o $(OBJ)proflay.o $(OBJ)prwarn.o\
 $(OBJ)pwreal.o $(OBJ)pwrefco.o $(OBJ)quadfit.o $(OBJ)raycalc.o $(OBJ)rayplot.o\
 $(OBJ)rays.o $(OBJ)rcheck.o $(OBJ)rdrcalc.o $(OBJ)rdrlay.o $(OBJ)rdrsum.o\
 $(OBJ)readerr.o $(OBJ)record.o $(OBJ)rnglay.o $(OBJ)rtc_calc.o $(OBJ)rtc_multi.o\
 $(OBJ)rthplot.o $(OBJ)rtpres.o $(OBJ)rzlay.o $(OBJ)setbp.o $(OBJ)setop.o\
 $(OBJ)shadcal.o $(OBJ)shadow.o $(OBJ)smcaus.o $(OBJ)smfind.o $(OBJ)sonify.o\
 $(OBJ)srcstep.o $(OBJ)srlev.o $(OBJ)sscalc.o $(OBJ)star.o $(OBJ)sublay.o\
 $(OBJ)timatt.o $(OBJ)tmatlay.o $(OBJ)traduc.o $(OBJ)traxpl.o\
 $(OBJ)winefit.o $(OBJ)zwine.o $(OBJ)trfhead.o $(OBJ)svp_read_gama.o\
 $(OBJ)getdat2tilt.o $(OBJ)getdat3tilt.o $(OBJ)geomtilt.o\
 $(OBJ)forwinit.o $(OBJ)write_info.o $(OBJ)write_gama_svp.o $(OBJ)write_gama_opt.o


all:	$(XLIB) 

$(XLIB): $(OBJS)
	ar rvu $(XLIB)  $(OBJS)
	$(RANLIB) $(XLIB)


$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) $(PROG).o $(XLIB) 

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)

$(OBJ)sub.o:	${OBJO}	${COMS}
	$(RM) -f sub.f
	cat 	${OBJO}> sub.f
	$(FC)  $(FFLAGS_GAMA)  -c sub.f
	$(MV) sub.o $(OBJ)sub.o

