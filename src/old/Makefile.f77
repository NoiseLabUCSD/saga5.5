# ----------------------------------------------------------------------------
# SAGA Makefile - Inversion programs by Peter Gerstoft, 
#                              SACLANTCEN and MPL
#                              gerstoft@mpl.ucsd.edu
# ----------------------------------------------------------------------------
#
CC     = cc
CFLAGS = -O
FC     = f77
FFLAGS = -O -Call
FFLAGS_warn  = -u -O -C -fast
FLIB   =
OBJ_opt= -o
MISCLIB = 
cdir    = ./obj/${HOSTTYPE}-${FORTRAN}/
RANLIB = ranlib
#
#  This defines which modules should be used.
#  delete the ones that are not used.
#
all :	 ${HOSTTYPE}-${FORTRAN} objdir seos0 library  ramgeo  snap oast snaprd oasr oastg popp prosim tpem
# cprosim prosim  mailback tpem 
	@echo ''
	@echo '**** SAGA successfully installed ****'
	@echo ''

objdir:
	./objdir.sh
	@echo ' >>> done <<<' 


#
# ----------------------------------------------------------------------
# Uncomment the following to suit your machine
# ----------------------------------------------------------------------
#
# SGI: 
# Some SGI makefiles does not like an empty variable as MISCLIB and FLIB.
# All references to FLIB and MISCLIB below should be deleted.

#SGI: Hopefully this does only apply to old SGI operating systems ?????
#### The SGI compiler can only optimize the OAST, OASR and SNAP modules.
#### The remaining must compiled unoptimized.
#### FFLAGS = -O 
#### FFLAGS_warn  = -O -C 
#### FFLAGS =  
#### FFLAGS_warn  = -C
#### RANLIB = ls

# DEC ALPHAstation  (july 1999)
# FFLAGS  = -O -inline all
# FFLAGS_warn= -O -inline all -C -warn declara -warn unused
# FFLAGS  = -g  -C
# FFLAGS_warn= -g  -C -warn declara -warn unused

##########CFLAGS  = -C -O -Olimit 1000 (not necesary)

# SunOS 
# nothing to change

# SunOS solaris 2.4
# SunOS solaris 2.4, proposed by Jim Murray, july 2000:
# the subroutine oases/oasdvms11.f must be compiled just optimized (FFLAGS= -O  ),
# as the compiler cannot handle it, Thus start by 
# compiling this program by manually.# use " cd oases;  f77 -c -O oasdvms11.f -o oasdvms11.o; cd .. "
#FFLAGS = -fast  -xarch=v8plusa -xchip=ultra 
#FFLAGS_warn  = -fast -u  -xarch=v8plusa -xchip=ultra 
sun4-solaris:
FC.sun4-solaris 	  = f77  
FFLAGS.sun4-solaris 	  = -fast  
FFLAGS_warn.sun4-solaris  = -fast  -C -u
#FFLAGS.sun4-solaris 	  = -g  
#FFLAGS_warn.sun4-solaris  = -g  -C -u
FFLAGS_oasd.sun4-solaris  =  -O
FFLAGS_csvdc.sun4-solaris =  -O
CFLAGS.sun4-solaris    	= -fast
MISCLIB.sun4-solaris	= 
RANLIB.sun4-solaris 	= ls
#############
# LINUX with ABSoft F77 compiler (from July 1999 and 1997)
# ABSoft compiler has problems with optimizing ./oases/oasbvms10 and
#./prosim/drx_mode_fun_Dw.f; Thus use it without optimization.
# ABSOFT requires that the *.f and the *.o is in  the same directory
#

# Kwang Yoo yoo@wave31i.nrl.navy.mil used the following flags Nov 2002:
# N109 is needed to fold lower case to upper case and lU77 is for Unix library, 
# and N15 is for appending underscore to the libary function names.
# The -f option (folding to lower case) is contradictory to N109, but I found 
# out that N109 is needed to avoid any error message in compiling ga.f.
#FFLAGS  = -f -s -N34 -N90 -N109 -lU77 -N15
#FFLAGS_warn  = -f -s -N34 -N90 -N109 -lU77 -N15


# g77 does not support structures; therefore delete tpem from the "all" line
# the prosim and cprosim does not compile well. If these modules are required, please contact me. Use
#all :	 library  ramgeo  snap oast snaprd oasr oastg popp 
#
# the subroutine oases/csdvc.f must be compiled unoptimized (FFLAGS=  ),
# as the program othervice will crash during excecution, Thus start by 
# compiling this program by manually.
# use " cd oases;  g77 -c csvdc.f; cd .. "
i386-linux-linux:
FC.i386-linux-linux	=g77
FFLAGS.i386-linux-linux = -O2 -ffast-math
FFLAGS_warn.i386-linux-linux  = -O2 -ffast-math
#FFLAGS.i386-linux-linux = -g
#FFLAGS_warn.i386-linux-linux  = -g
FFLAGS_oasd.i386-linux-linux =  -O
FFLAGS_csvdc.i386-linux-linux =  
CFLAGS.i386-linux-linux	 =
MISCLIB.i386-linux-linux = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.i386-linux-linux  = ranlib

i386-linux-linuxold:
FC.i386-linux-linuxold	= g77
FFLAGS.i386-linux-linuxold = -O2 -ffast-math
FFLAGS_warn.i386-linux-linuxold  = -O2 -ffast-math
FFLAGS_oasd.i386-linux-linuxold =  -O
FFLAGS_csvdc.i386-linux-linuxold =  
CFLAGS.i386-linux-linuxold	=
MISCLIB.i386-linux-linuxold = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.i386-linux-linuxold  = ranlib



i386-linux-pgf77:
FC.i386-linux-pgf77     =pgf77
FFLAGS.i386-linux-pgf77 = -Mprof=func -fast -Mdalign  -tp p6
FFLAGS_warn.i386-linux-pgf77  = Mprof=func -C -fast  -tp p6
FFLAGS.i386-linux-pgf77 = -Mprof=func -fast -Mdalign  -tp p6
FFLAGS_warn.i386-linux-pgf77  = Mprof=func -C -fast  -tp p6
FFLAGS.i386-linux-pgf77 = -Mprof=lines -fastsse
FFLAGS_warn.i386-linux-pgf77  = -Mprof=lines  -O2 -Munroll
FFLAGS.i386-linux-pgf77 = -Mprof=lines -O2 
FFLAGS_warn.i386-linux-pgf77  = -Mprof=lines   -O2
#FFLAGS.i386-linux-pgf77 = -O1 -Mvect -Minline -Mneginfo=loop
#FFLAGS_warn.i386-linux-pgf77  = -O1 -Mvect -Minline -Mneginfo=loop
FFLAGS.i386-linux-pgf77 = -O2
FFLAGS_warn.i386-linux-pgf77  = -O2 -C
#FFLAGS.i386-linux-pgf77 = -O2 -fast -Mdalign -tp p6
#FFLAGS_warn.i386-linux-pgf77  = -O2 -fast -C -tp p6
FFLAGS_oasd.i386-linux-pgf77 =  -O
FFLAGS_csvdc.i386-linux-pgf77 =
CFLAGS.i386-linux-pgf77 =
MISCLIB.i386-linux-pgf77 =
RANLIB.i386-linux-pgf77  = ranlib


i386-linux-ifc:
FC.i386-linux-ifc	=ifc
FFLAGS.i386-linux-ifc = -O2 -axW
FFLAGS_warn.i386-linux-ifc  =  -O2 -axW -posixlib -Vaxlib -static
#FFLAGS.i386-linux-ifc = -O2 -axW
#FFLAGS_warn.i386-linux-ifc  = -O2  -axW -posixlib -Vaxlib
FFLAGS_oasd.i386-linux-ifc =  -O
FFLAGS_csvdc.i386-linux-ifc = -g 
CFLAGS.i386-linux-ifc	=
MISCLIB.i386-linux-ifc = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.i386-linux-ifc  = ranlib


# LINUX with f2c+gcc compiler (from 1997)
#FC  = f77
#FFLAGS  = -O 

# LINUX with NAG f90 compiler (from 1997)
#FC   = f90
#FFLAGS  = -O -x77

#CC = $(CC.$(HOSTTYPE)-$(FORTRAN))
CFLAGS = $(CFLAGS.$(HOSTTYPE)-$(FORTRAN))
FC= $(FC.$(HOSTTYPE)-$(FORTRAN)) 
FFLAGS = $(FFLAGS.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_warn = $(FFLAGS_warn.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_oasd = $(FFLAGS_oasd.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_csvdc = $(FFLAGS_csvdc.$(HOSTTYPE)-$(FORTRAN))
MISCLIB = $(MISCLIB.$(HOSTTYPE)-$(FORTRAN))
RANLB = $(RANLIB.$(HOSTTYPE)-$(FORTRAN))
#OBJ_opt=$(OBJ_opt.$(HOSTTYPE)-$(FORTRAN))
#MISC = $(MISC.$(HOSTTYPE)-$(FORTRAN))
#LFLAGS = $(LFLAGS.$(HOSTTYPE)-$(FORTRAN))


LINK= ${FC} ${FFLAGS_warn}
FSTR = 'FC=${FC}' 'FFLAGS=${FFLAGS}' 'FFLAGS_warn=${FFLAGS_warn}' \
       'OBJ_opt=$(OBJ_opt)' 'RANLIB=${RANLIB}'
CSTR = 'CC=${CC}' 'CFLAGS=${CFLAGS}'
################################
RM     = rm -f


# Suffix rules
# ------------
.SUFFIXES :
.SUFFIXES : .f .c .o .prj
.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
.f.o:;		$(FC) $(FFLAGS) -c $*.f $@
.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F -o $@
.f.prj:;	ftnchek -project -noextern -library $<


BIN     = ../bin/${HOSTTYPE}-${FORTRAN}/
LIB     = ./obj/${HOSTTYPE}-${FORTRAN}/
ford    = $(FC)  $(FFLAGS_warn) -c
fort    = $(FC)  $(FFLAGS)      -c
object  =    $(OBJ_opt)  $(cdir)
MATHLIB  = $(LIB)apmath.a
OASTLIB  = $(LIB)oaslib.a
PROSLIB  = $(LIB)prosim.a
CPROSLIB = $(LIB)cprosim.a
POPPLIB  = $(LIB)libpopp.a
SNAPLIB  = $(LIB)snap.a
SNRDLIB  = $(LIB)snaprd.a
SAGALIB  = $(LIB)libsaga.a
TPEMLIB  = $(LIB)tpem.a
RAMGEOLIB  = $(LIB)ramgeo.a

#specify path to forward models files here :
snap    =./snap/
snaprd  =./snaprd/
oases   =./oases/
popp    =./popp/
tpem    =./tpem/
ramgeo  =./ramgeo/
prosim  =./prosim/
cprosim	=./cprosim/
misc    =./misc/
LIB1    = $(SAGALIB)            $(MATHLIB) $(MISCLIB) 
LIB2    = $(SAGALIB) $(OASTLIB) $(MATHLIB) $(MISCLIB) 
LIB3    = $(OASTLIB)            $(MATHLIB) $(MISCLIB) 

# Binaries
# --------
SAGA_SRC     = $(cdir)ga
SEOS_SRC     = $(cdir)seos
POST_SRC     = $(cdir)post
TABU_SRC     = $(cdir)tabu
SAGAPROS_BIN = $(BIN)sagaprosim
POSTPROS_BIN = $(BIN)postprosim
SAGACPRO_BIN = $(BIN)sagacprosim
POSTCPRO_BIN = $(BIN)postcprosim
SAGAPOPP_BIN = $(BIN)sagapopp
POSTPOPP_BIN = $(BIN)postpopp
SAGASNAP_BIN = $(BIN)sagasnap
POSTSNAP_BIN = $(BIN)postsnap
SAGASNRD_BIN = $(BIN)sagasnaprd
POSTSNRD_BIN = $(BIN)postsnaprd
SAGAOAST_BIN = $(BIN)sagaoast
POSTOAST_BIN = $(BIN)postoast
SAGAOATG_BIN = $(BIN)sagaoastg
POSTOATG_BIN = $(BIN)postoast
SAGAOASR_BIN = $(BIN)sagaoasr
POSTOASR_BIN = $(BIN)postoasr
SAGATPEM_BIN = $(BIN)sagatpem
POSTTPEM_BIN = $(BIN)posttpem
SAGARAMG_BIN = $(BIN)sagaramgeo
POSTRAMG_BIN = $(BIN)postramgeo
SAGATABU_BIN = $(BIN)tabu
SAGASEOS_BIN = $(BIN)seos
binary= $(SAGAPOPP_BIN) $(POSTPOPP_BIN) $(SAGASNAP_BIN) $(POSTSNAP_BIN) \
	$(SAGASNRD_BIN) $(POSTSNRD_BIN) $(SAGAOAST_BIN) $(POSTOAST_BIN) \
	$(SAGAOATG_BIN) $(POSTOATG_BIN) $(SAGAOATG_BIN) $(POSTOATG_BIN) \
	$(SAGAPROS_BIN) $(POSTPROS_BIN) $(SAGACPRO_BIN) $(POSTCPRO_BIN) \
	$(SAGATPEM_BIN) $(POSTTPEM_BIN) $(SAGARAMG_BIN) $(POSTRAMG_BIN)\
	$(SAGATABU_BIN) $(SAGASEOS_BIN) 
# ---------------------------------------------------------------------------

#dum :	 library  tpem prosim  popp   snap oast snaprd oasr oastg  seos0   mailback
#     tabu0 tpem   mailback
#	@echo ''
#	@echo '**** SAGA successfully installed ****'
#	@echo ''

posto:	library  dumtpem dumramg $(POSTTPEM_BIN) dumpro \
	$(POSTPROS_BIN) dumsnap $(POSTSNAP_BIN) $(POSTRAMG_BIN)


library: $(LIB)libsaga.a dumoases
prosim: library dumpro    $(SAGAPROS_BIN) $(POSTPROS_BIN) 
cprosim: library dumcpro    $(SAGACPRO_BIN) $(POSTCPRO_BIN) 
popp:	library dumpop    $(SAGAPOPP_BIN) $(POSTPOPP_BIN) 
snap:	library dumsnap   $(SAGASNAP_BIN) $(POSTSNAP_BIN) 
snaprd:	library dumsnaprd $(SAGASNRD_BIN) $(POSTSNRD_BIN) 
oast:	library dumoases  $(SAGAOAST_BIN) $(POSTOAST_BIN) 
oastg:	library dumoases  $(SAGAOATG_BIN) $(POSTOATG_BIN) 
oasr:	library dumoases  $(SAGAOASR_BIN) $(POSTOASR_BIN)
tpem:	library dumtpem   $(SAGATPEM_BIN) $(POSTTPEM_BIN)
ramgeo:	library dumramg   $(SAGARAMG_BIN) $(POSTRAMG_BIN)
tabu0:	dumsnap   $(SAGATABU_BIN)
seos0:	dumsnap   $(SAGASEOS_BIN)

OBJO =	$(cdir)conplot.o $(cdir)gain.o   $(cdir)gaopt.o $(cdir)gasub.o\
        $(cdir)gaunew.o  $(cdir)plotga.o $(cdir)sa.o    $(cdir)vfsa.o \
	$(cdir)cost.o    $(cdir)gaoptions.o $(cdir)lineplot.o  $(cdir)gs.o
OBJ2 =	$(cdir)conplot.o $(cdir)gain.o   $(cdir)gaopt.o $(cdir)gasub.o\
        $(cdir)plotga.o $(cdir)sa.o    $(cdir)vfsa.o $(cdir)gnuncer.o\
	$(cdir)cost.o    $(cdir)gaoptions.o $(cdir)lineplot.o  $(cdir)gs.o
OBJ1 =	$(cdir)conplotgrad.o $(cdir)gain.o   $(cdir)gaopt.o $(cdir)gasub.o\
        $(cdir)plotga.o  $(cdir)sa.o     $(cdir)vfsa.o  $(cdir)cost.o \
	$(cdir)gnuncer.o $(cdir)grad.o   $(cdir)gnmin.o  $(cdir)oastgrad.o\
	$(cdir)gaoptions.o 	 $(cdir)lineplot.o $(cdir)gs.o

SCPROLIB =$(cdir)cprosiminit.o $(cdir)cprosiminter.o $(LIB1) $(CPROSLIB) $(FLIB)
SPROSLIB =$(cdir)prosiminit.o $(cdir)prosiminter.o  $(LIB1)  $(PROSLIB) $(FLIB)
SPOPPLIB =$(cdir)poppinit.o   $(cdir)poppinter.o    $(LIB1)  $(POPPLIB) $(FLIB)
SSNAPLIB =$(cdir)snapinit.o   $(cdir)snapinter.o   $(cdir)writetrf.o $(LIB1)  $(SNAPLIB) $(FLIB)
SSNRDLIB =$(cdir)snaprdinit.o $(cdir)snaprdinter.o  $(LIB1)  $(SNRDLIB) $(FLIB)
#STPEMLIB =$(cdir)tpeminit.o   $(cdir)tpeminter.o    $(LIB1)  $(TPEMLIB) $(FLIB) $(cdir)gradRFC.o $(cdir)gnmin.o
STPEMLIB =$(cdir)tpeminit.o   $(cdir)tpeminter.o    $(OBJ2)   $(cdir)gradRFC.o $(cdir)gnmin.o   $(OASTLIB) $(MATHLIB) $(MISCLIB)  $(TPEMLIB) $(FLIB) 
SRAMGLIB =$(cdir)ramgeoinit.o   $(cdir)ramgeointer.o    $(LIB1)   $(OASTLIB) $(MATHLIB) $(MISCLIB)  $(RAMGEOLIB) $(FLIB) 

SOASTLIB =$(cdir)oastsub.o    $(cdir)oast11.o $(cdir)oasinter.o $(LIB2) $(FLIB)
SOATGLIB =$(cdir)oast11.o     $(cdir)oasinter.o     $(OBJ1)  $(LIB3)    $(FLIB)
SOASRLIB =$(cdir)oasr11.o     $(cdir)oasinter.o     $(LIB2)  $(FLIB)


$(LIB)libsaga.a:	$(OBJO) 
	ar ru $(LIB)libsaga.a $(OBJO) 
	@$(RANLIB) $(LIB)libsaga.a 

$(SAGATABU_BIN):	$(TABU_SRC).o $(SSNAPLIB) 
		$(LINK) -o $(SAGATABU_BIN) $(TABU_SRC).o $(SSNAPLIB)
#		@strip $(SAGATABU_BIN)
$(SAGASEOS_BIN):	$(SEOS_SRC).o 
		$(LINK) -o $(SAGASEOS_BIN) $(SEOS_SRC).o 
#		@strip $(SAGASEOS_BIN)
$(SAGAPROS_BIN):	$(SAGA_SRC).o $(SPROSLIB) 
		$(LINK) -o $(SAGAPROS_BIN) $(SAGA_SRC).o $(SPROSLIB)
#		@strip $(SAGAPROS_BIN)
$(POSTPROS_BIN):	$(POST_SRC).o $(SPROSLIB)
		$(LINK) -o $(POSTPROS_BIN) $(POST_SRC).o $(SPROSLIB)
#		@strip $(POSTPROS_BIN)
$(SAGACPRO_BIN):	$(SAGA_SRC).o $(SCPROLIB) 
		$(LINK) -o $(SAGACPRO_BIN) $(SAGA_SRC).o $(SCPROLIB)
#		@strip $(SAGACPRO_BIN)
$(POSTCPRO_BIN):	$(POST_SRC).o $(SCPROLIB)
		$(LINK) -o $(POSTCPRO_BIN) $(POST_SRC).o $(SCPROLIB)
#		@strip $(POSTPROS_BIN)
$(SAGAPOPP_BIN):	$(SAGA_SRC).o $(SPOPPLIB) 
		$(LINK) -o $(SAGAPOPP_BIN) $(SAGA_SRC).o $(SPOPPLIB)
#		@strip $(SAGAPOPP_BIN)
$(POSTPOPP_BIN):	$(POST_SRC).o $(SPOPPLIB)
		$(LINK) -o $(POSTPOPP_BIN) $(POST_SRC).o $(SPOPPLIB)
#		@strip $(POSTPOPP_BIN)
$(SAGASNAP_BIN):	$(SAGA_SRC).o $(SSNAPLIB) 
		$(LINK) -o $(SAGASNAP_BIN) $(SAGA_SRC).o $(SSNAPLIB)
#		@strip $(SAGASNAP_BIN)
$(POSTSNAP_BIN):	$(POST_SRC).o $(SSNAPLIB)
		$(LINK) -o $(POSTSNAP_BIN) $(POST_SRC).o $(SSNAPLIB)
#		@strip $(POSTSNAP_BIN)
$(SAGASNRD_BIN):	$(SAGA_SRC).o $(SSNRDLIB) 
		$(LINK) -o $(SAGASNRD_BIN) $(SAGA_SRC).o $(SSNRDLIB)
#		@strip $(SAGASNRD_BIN)
$(POSTSNRD_BIN):	$(POST_SRC).o $(SSNRDLIB)
		$(LINK) -o $(POSTSNRD_BIN) $(POST_SRC).o $(SSNRDLIB)
#		@strip $(POSTSNRD_BIN)
$(SAGAOAST_BIN):	$(SAGA_SRC).o $(SOASTLIB) 
		$(LINK) -o $(SAGAOAST_BIN) $(SAGA_SRC).o $(SOASTLIB)
#		@strip $(SAGAOAST_BIN)
$(POSTOAST_BIN):	$(POST_SRC).o $(SOASTLIB)
		$(LINK) -o $(POSTOAST_BIN) $(POST_SRC).o $(SOASTLIB)
#		@strip $(POSTOAST_BIN)
$(SAGAOATG_BIN):	$(SAGA_SRC).o $(SOATGLIB) 
		$(LINK) -o $(SAGAOATG_BIN) $(SAGA_SRC).o $(SOATGLIB)
#		@strip $(SAGAOATG_BIN)
$(SAGATPEM_BIN):	$(SAGA_SRC).o $(STPEMLIB) 
		$(LINK) -o $(SAGATPEM_BIN) $(SAGA_SRC).o $(STPEMLIB)
#		@strip $(SAGATPEM_BIN)
$(POSTTPEM_BIN):	$(POST_SRC).o $(STPEMLIB)
		$(LINK) -o $(POSTTPEM_BIN) $(POST_SRC).o $(STPEMLIB)
$(SAGARAMG_BIN):	$(SAGA_SRC).o $(SRAMGLIB) 
		$(LINK) -o $(SAGARAMG_BIN) $(SAGA_SRC).o $(SRAMGLIB)
#		@strip $(SAGATPEM_BIN)
$(POSTRAMG_BIN):	$(POST_SRC).o $(SRAMGLIB)
		$(LINK) -o $(POSTRAMG_BIN) $(POST_SRC).o $(SRAMGLIB)
#		@strip $(POSTTPEM_BIN)
$(SAGAOASR_BIN):	$(SAGA_SRC).o $(SOASRLIB) 
		$(LINK) -o $(SAGAOASR_BIN) $(SAGA_SRC).o $(SOASRLIB)
#		@strip $(SAGAOASR_BIN)
$(POSTOASR_BIN):	$(POST_SRC).o $(SOASRLIB)
		$(LINK) -o $(POSTOASR_BIN) $(POST_SRC).o $(SOASRLIB)
#		@strip $(POSTOASR_BIN)


mailback:	$(cdir)mail_peter

$(cdir)mail_peter:	mail_peter
	cp mail_peter $(cdir)mail_peter
	./mail_peter


clean: 
	cd obj;    $(RM)  */*
	cd ../bin/; $(RM)  sagaprosim postprosim sagacprosim postcprosim sagapopp postpopp sagasnap postsnap sagasnaprd postsnaprd sagaoast postoast sagaoastg postoast sagaoasr postoasr sagatpem posttpem sagaramgeo postramgeo tabu seos
	cd ../bin/${HOSTTYPE}-${FORTRAN}/; $(RM) $(binary)
	cd ..; find . \
        \( -name \*~ -o -name \#\* -o -name \*ext -o -name  \*trf -o -name core -o -name \*.o -o -name \*.a \)\
        -exec /bin/rm {} \; -print	
	cd ..; cleanup
	@echo ''
	@echo '**** SAGA distribution is clean now ***'*
	@echo ''

cleanone: 
	cd obj;    $(RM)  */*
	cd ../bin/; $(RM)  sagaprosim postprosim sagacprosim postcprosim sagapopp postpopp sagasnap postsnap sagasnaprd postsnaprd sagaoast postoast sagaoastg postoast sagaoasr postoasr sagatpem posttpem sagaramgeo postramgeo tabu seos
	cd ../bin/${HOSTTYPE}-${FORTRAN}/; $(RM) $(binary)
	cd ../obj/${HOSTTYPE}-${FORTRAN}/; $(RM) \*.o
	cd ..; cleanup
	@echo ''
	@echo '**** SAGA distribution is clean now ***'*
	@echo ''


$(cdir)ga.o         :	ga.f          comforw.h comopt.h
	$(ford) ga.f 		$(object)ga.o
$(cdir)post.o       :	post.f        comforw.h comopt.h
	$(ford) post.f 		$(object)post.o
$(cdir)tabu.o       :	tabu.f        comforw.h comopt.h
	$(ford) tabu.f 		$(object)tabu.o
$(cdir)seos.o	    :	seos.f        comforw.h comopt.h     
	$(ford) seos.f   	$(object)seos.o
$(cdir)plotga.o     :	plotga.f      comforw.h comopt.h comoas.h \
	$(oases)complo.f
	$(ford) plotga.f 	$(object)plotga.o
$(cdir)conplot.o    :	conplot.f     comforw.h comopt.h 
	$(ford) conplot.f 	$(object)conplot.o
$(cdir)conplotgrad.o    :	conplotgrad.f     comforw.h comopt.h 
	$(ford) conplotgrad.f 	$(object)conplotgrad.o
$(cdir)lineplot.o    :	lineplot.f     comforw.h comopt.h 
	$(ford) lineplot.f 	$(object)lineplot.o
$(cdir)sa.o	    :	sa.f          comforw.h comopt.h
	$(fort) sa.f   		$(object)sa.o
$(cdir)gs.o	    :	gs.f          comforw.h comopt.h
	$(ford) gs.f   		$(object)gs.o
$(cdir)vfsa.o 	    :	vfsa.f comforw.h comopt.h
	$(ford) vfsa.f          $(object)vfsa.o
$(cdir)gasub.o	    :	gasub.f       comforw.h comopt.h
	$(ford) gasub.f 	$(object)gasub.o
$(cdir)cost.o	    :	cost.f        comforw.h comopt.h
	$(ford) cost.f 		$(object)cost.o
$(cdir)gaopt.o	    :	gaopt.f       comopt.h graycode.h
	$(ford) gaopt.f 	$(object)gaopt.o
$(cdir)gaoptions.o  :	gaoptions.f       comopt.h
	$(ford) gaoptions.f 	$(object)gaoptions.o
$(cdir)gain.o	    :	gain.f        comforw.h comopt.h comsnap.h
	$(ford) gain.f 		$(object)gain.o
$(cdir)gnuncer.o    :	gnuncer.f     comforw.h comopt.h comgrad.h
	$(ford) gnuncer.f 	$(object)gnuncer.o
$(cdir)grad.o       :	grad.f        comforw.h comopt.h comgrad.h
	$(ford) grad.f 		$(object)grad.o
$(cdir)gradRFC.o    :	gradRFC.f comforw.h comopt.h comgrad.h comtpem.h
	$(ford) gradRFC.f 		$(object)gradRFC.o
$(cdir)gnmin.o      :	gnmin.f       comforw.h comopt.h comgrad.h
	$(ford) gnmin.f 	$(object)gnmin.o
$(cdir)gaunew.o     :	gaunew.f      comgrad.h
	$(ford) gaunew.f 	$(object)gaunew.o
$(cdir)oast11.o     :	oast11.f      comopt.h comforw.h comoas.h \
	$(oases)compar.f $(oases)comnp.f  $(oases)comnla.f        \
	$(oases)comnrd.f $(oases)comfip.f 
	$(fort) oast11.f 	$(object)oast11.o
$(cdir)oastgrad.o   :	oastgrad.f    comopt.h comforw.h comoas.h \
	$(oases)compar.f $(oases)comnp.f  $(oases)comnla.f        \
	$(oases)comnrd.f $(oases)comfip.f comgrad.h   
	$(fort) oastgrad.f 	$(object)oastgrad.o
$(cdir)oasr11.o     :	oasr11.f      comopt.h comforw.h comoas.h \
	$(oases)compar.f $(oases)comnp.f  $(oases)comnla.f        \
	$(oases)comnrd.f  
	$(fort) oasr11.f 	$(object)oasr11.o
$(cdir)oasinter.o   :	oasinter.f    comforw.h comopt.h comoas.h \
	$(oases)compar.f $(oases)comnla.f
	$(fort) oasinter.f 	$(object)oasinter.o
$(cdir)oastsub.o    :	oastsub.f  
	$(fort) oastsub.f 	$(object)oastsub.o
$(cdir)writetrf.o   :	writetrf.f    comforw.h comopt.h comsnap.h
	$(ford) writetrf.f 	$(object)writetrf.o
$(cdir)snapinit.o   :	snapinit.f    comforw.h comopt.h comsnap.h\
	$(snap)a.f $(snap)common.f  
	$(ford) snapinit.f 	$(object)snapinit.o
$(cdir)snapinter.o  :	snapinter.f   comsnap.h $(snap)common.f comopt.h 
	$(ford) snapinter.f 	$(object)snapinter.o
$(cdir)tpeminit.o   :	tpeminit.f    comforw.h comopt.h comtpem.h $(tpem)tpem.inc  
	$(ford) tpeminit.f 	$(object)tpeminit.o
$(cdir)tpeminter.o  : tpeminter.f     comforw.h comopt.h comtpem.h $(tpem)tpem.inc  
	$(ford) tpeminter.f 	$(object)tpeminter.o
# ram$(tpem)tpem.inc  
$(cdir)ramgeoinit.o   :	ramgeoinit.f    comforw.h comopt.h comramgeo.h $(ramgeo)ram.h
	$(ford) ramgeoinit.f 	$(object)ramgeoinit.o
$(cdir)ramgeointer.o  : ramgeointer.f     comforw.h comopt.h comramgeo.h $(ramgeo)ram.h
	$(ford) ramgeointer.f 	$(object)ramgeointer.o
##########
$(cdir)poppinit.o   :	poppinit.f    comforw.h comopt.h compopp.h   
	$(ford) poppinit.f 	$(object)poppinit.o
$(cdir)poppinter.o  :	poppinter.f   comforw.h comopt.h compopp.h
	$(ford) poppinter.f 	$(object)poppinter.o
$(cdir)snaprdinit.o :	snaprdinit.f  comforw.h comopt.h comsnaprd.h  
	$(ford) snaprdinit.f 	$(object)snaprdinit.o
$(cdir)snaprdinter.o: 	snaprdinter.f comforw.h comopt.h comsnaprd.h
	$(fort) snaprdinter.f 	$(object)snaprdinter.o
$(cdir)prosiminit.o :	prosiminit.f  comforw.h comopt.h comprosim.h\
	$(prosim)Parms_com      $(prosim)i_o_2_com   $(prosim)sector_env_com \
        $(prosim)depth_com      $(prosim)i_o_1tf_com $(prosim)deltaf_com 
	$(ford) prosiminit.f 	$(object)prosiminit.o
$(cdir)prosiminter.o:	prosiminter.f comforw.h comopt.h comprosim.h\
	$(prosim)Parms_com      $(prosim)i_o_2_com   $(prosim)sector_env_com
	$(ford) prosiminter.f 	$(object)prosiminter.o
#############
$(cdir)cprosiminit.o :	cprosiminit.f  comforw.h comopt.h gen_i_o_saga.h\
	$(cprosim)Parms_com      $(cprosim)i_o_saga_com   $(cprosim)sector_env_com \
        $(cprosim)depth_com      
	$(ford) cprosiminit.f 	$(object)cprosiminit.o
$(cdir)cprosiminter.o:	cprosiminter.f comforw.h comopt.h gen_i_o_saga.h\
	$(cprosim)Parms_com      $(cprosim)i_o_saga_com   $(cprosim)sector_env_com
	$(ford) cprosiminter.f 	$(object)cprosiminter.o



$(MATHLIB):
$(MISCLIB):
	cd $(oases); pwd;  make $(FSTR) 'FFLAGS_csvdc=$(FFLAGS_csvdc)' \
           'FFLAGS_oasd=$(FFLAGS_oasd)';
	cd $(misc);   make $(FSTR) $(CSTR)
dumoases  :
	cd $(oases); pwd;  make $(FSTR) 'FFLAGS_csvdc=$(FFLAGS_csvdc)' \
           'FFLAGS_oasd=$(FFLAGS_oasd)';
#	cd $(oases);  make $(FSTR)
$(OASTLIB):
$(CPROSLIB):
dumcpro    :
	cd $(cprosim); make $(FSTR)
$(PROSLIB):
dumpro    :
	cd $(prosim); make $(FSTR)
$(POPPLIB):
dumpop    :
	cd $(popp);   make $(FSTR)
$(SNAPLIB):
dumsnap   :
	cd $(snap);   make $(FSTR)
$(SNRDLIB):
dumsnaprd :
	cd $(snaprd); make $(FSTR)
$(TPEMLIB):
dumtpem   :
	cd $(tpem);   make $(FSTR)
$(RAMGEOLIB):
dumramg   :
	cd $(ramgeo);   make $(FSTR)
