# ----------------------------------------------------------------------------
# SAGA Makefile - Inversion programs by Peter Gerstoft, 
#                              MPL, Scripps Institute of Oceanography
#                              gerstoft@ucsd.edu
# ----------------------------------------------------------------------------
#
CC     = cc
CFLAGS = -O
FC     = f90
FFLAGS = -O -Call
FFLAGS_warn  = -u -O -C -fast
FLIB   =
OBJ_opt= -o
MISCLIB = 
cdir    = ./obj/${HOSTTYPE}-${FORTRAN}/
cdirorca= ./orca90//${HOSTTYPE}-${FORTRAN}/

RANLIB = ranlib
#
#  This defines which modules should be used.
#  delete the ones that are not used.
# old module: seos0 (!! oastg did not compile under x86_64-linux  ....) oastg
all :	library snap orca ramgeo oast snaprd oasr  popp prosim cprosim  gama tpem  resample 
	@echo ''
	@echo '**** SAGA successfully installed ****'
	@echo ''
#	@make -s $(cdir)mail_peter

objdir:
	./objdir.sh
	@echo ' >>> done <<<' 

# ----------------------------------------------------------------------
# Uncomment the following to suit your machine
# ----------------------------------------------------------------------

# intel fortran compiler release 10.0 ifort (August 2006)
# recipes.f must be compiled in 01 mode (it is done via FFLAGS_oasd.)
# intel fortran compiler release 10.0 ifort (August 2006)
# recipes.f must be compiled in 01 mode (it is done via FFLAGS_oasd.)
intel-mac-ifort:
FC.intel-mac-ifort	     = ifort 
FC90.intel-mac-ifort	     = ifort
FFLAGS.intel-mac-ifort      = -O 
FFLAGS90.intel-mac-ifort        = -O 
#FFLAGS_warn.intel-pc-ifort   =  -O0 -u  -posixlib -Vaxlib -CB
#FFLAGS_GAMA_warn.intel-pc-ifort = -O -r8 -align -w  -posixlib -Vaxlib 
FFLAGS_warn.intel-mac-ifort   =  -O0 -u    -Vaxlib -CB
FFLAGS_GAMA_warn.intel-mac-ifort = -O -r8 -align -w    -Vaxlib 
FFLAGS_GAMA.intel-mac-ifort   = -O2  -r8 -align  -w  
FFLAGS_oasd.intel-mac-ifort   =  -O1 
FFLAGS_csvdc.intel-mac-ifort  = -g 
CFLAGS.intel-mac-ifort        = -m64
MISCLIB.intel-mac-ifort        = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.intel-mac-ifort        = ranlib
MODULE.intel-mac-ifort      = -I$(MOD) -module $(MOD)
MODULE2.intel-mac-ifort     = -I ../$(MOD) -module ../$(MOD)

# intel fortran compiler release 10.0 ifort (August 2006)
# recipes.f must be compiled in 01 mode (it is done via FFLAGS_oasd.)
intel-pc-ifort:
FC.intel-pc-ifort	     = ifort 
FC90.intel-pc-ifort	     = ifort
FFLAGS.intel-pc-ifort      = -O 
FFLAGS90.intel-pc-ifort        = -O 
#FFLAGS_warn.intel-pc-ifort   =  -O0 -u  -posixlib -Vaxlib -CB
#FFLAGS_GAMA_warn.intel-pc-ifort = -O -r8 -align -w  -posixlib -Vaxlib 
FFLAGS_warn.intel-pc-ifort   =  -O0 -u    -Vaxlib -CB
FFLAGS_GAMA_warn.intel-pc-ifort = -O -r8 -align -w    -Vaxlib 
FFLAGS_GAMA.intel-pc-ifort   = -O2  -r8 -align  -w  
FFLAGS_oasd.intel-pc-ifort   =  -O1 
FFLAGS_csvdc.intel-pc-ifort  = -g 
CFLAGS.intel-pc-ifort        = -m64
MISCLIB.intel-pc-ifort        = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.intel-pc-ifort        = ranlib
MODULE.intel-pc-ifort      = -I$(MOD) -module $(MOD)
MODULE2.intel-pc-ifort     = -I ../$(MOD) -module ../$(MOD)

# -Qextname  allows function w/o underscore
# -qtune -qarch
# F:xlf_s defines that it schould follow my schema defined in the user file 
# /etc/opt/ibmcmp/xlf/8.1/xlf.cfg
#
# For the xlf compiler it is necessary to rename /snap/getprofiles.f90 to
# /snap/getprofiles.f90 (and change the name in the makefile)
powermac-xlf: 
FC.powermac-xlf	      = xlf  -F:xlf_s -qextname  -O 
FC90.powermac-xlf	      = xlf90 -qextname -O 
FFLAGS.powermac-xlf       = 
FFLAGS90.powermac-xlf       = -qsuffix=f=f90 
FFLAGS_warn.powermac-xlf  = 
FFLAGS_oasd.powermac-xlf  =  
FFLAGS_csvdc.powermac-xlf =  
CFLAGS.powermac-xlf	  =
MISCLIB.powermac-xlf      =  ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
#./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.powermac-xlf       = ranlib
MODULE.powermac-xlf       = -I$(MOD) -qmoddir=$(MOD)
MODULE2.powermac-xlf       = -I ../$(MOD) -qmoddir=../$(MOD)

# intel fortran compiler release 10.1 ifort (Aug 2008)
x86_64-linux-ifort:
FC.x86_64-linux-ifort	     = ifort -static
FC90.x86_64-linux-ifort	     = ifort
FFLAGS.x86_64-linux-ifort        = -O -axW  -w 
FFLAGS90.x86_64-linux-ifort        = -O -axW  -w 
FFLAGS_warn.x86_64-linux-ifort   =  -O -axW -u   -w -posixlib -Vaxlib -CB
FFLAGS_GAMA_warn.x86_64-linux-ifort = -O -r8 -align -w -axW -posixlib -Vaxlib 
FFLAGS_GAMA.x86_64-linux-ifort   =   -r8 -align  -axW -w 
FFLAGS_GAMA.x86_64-linux-ifort   = -O2  -r8 -align  -axW -w  
FFLAGS_oasd.x86_64-linux-ifort   =  -O0 -w95
FFLAGS_csvdc.x86_64-linux-ifort  = -g -w95
CFLAGS.x86_64-linux-ifort        =
MISCLIB.x86_64-linux-ifort       = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.x86_64-linux-ifort        = ranlib
MODULE.x86_64-linux-ifort      = -I$(MOD) -module $(MOD)
MODULE2.x86_64-linux-ifort       = -I ../$(MOD) -module ../$(MOD)

# GNU fortran compiler gfortran 4.8.4  (Oct 2015)
x86_64-linux-gfortran:
FC.x86_64-linux-gfortran	     = gfortran
FC90.x86_64-linux-gfortran	     = gfortran
FFLAGS.x86_64-linux-gfortran        = -O   
FFLAGS90.x86_64-linux-gfortran        = -O   
FFLAGS_warn.x86_64-linux-gfortran   =  -O -u  
FFLAGS_GAMA_warn.x86_64-linux-gfortran = -O  -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_GAMA.x86_64-linux-gfortran   =  -O   -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_GAMA.x86_64-linux-gfortran   = -O2  -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_oasd.x86_64-linux-gfortran   =  -O0 
FFLAGS_csvdc.x86_64-linux-gfortran  = -g 
CFLAGS.x86_64-linux-gfortran        =
MISCLIB.x86_64-linux-gfortran       = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.x86_64-linux-gfortran        = ranlib
MODULE.x86_64-linux-gfortran      = 
MODULE2.x86_64-linux-gfortran       = 

# GNU fortran compiler gfortran 5.1 on Imac  (Oct 2015)
x86_64-gfortran:
FC.x86_64-gfortran	     = gfortran 
FC90.x86_64-gfortran	     = gfortran
FFLAGS.x86_64-gfortran        = -O   
FFLAGS90.x86_64-gfortran        = -O   
FFLAGS_warn.x86_64-gfortran   =  -O   
FFLAGS_GAMA_warn.x86_64-gfortran = -O   -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_GAMA.x86_64-gfortran   =  -O -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_GAMA.x86_64-gfortran   = -O2    -fdefault-real-8  -Wall -fbounds-check -fcheck=all
FFLAGS_oasd.x86_64-gfortran   =  -O0 
FFLAGS_csvdc.x86_64-gfortran  = -g 
CFLAGS.x86_64-gfortran        =
MISCLIB.x86_64-gfortran       = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.x86_64-gfortran        = ranlib
MODULE.x86_64-gfortran      = 
MODULE2.x86_64-gfortran       = 

# intel fortran compiler release 9.0 ifort (Feb 2006)
i386-linux-ifort:
FC.i386-linux-ifort	     = ifort -static
FC90.i386-linux-ifort	     = ifort
FFLAGS.i386-linux-ifort        = -O -axW  -w 
FFLAGS90.i386-linux-ifort        = -O -axW  -w 
FFLAGS_warn.i386-linux-ifort   =  -O -axW -u   -w -posixlib -Vaxlib -CB
FFLAGS_GAMA_warn.i386-linux-ifort = -O -r8 -align -w -axW -posixlib -Vaxlib 
FFLAGS_GAMA.i386-linux-ifort   =   -r8 -align  -axW -w 
FFLAGS_GAMA.i386-linux-ifort   = -O2  -r8 -align  -axW -w  
FFLAGS_oasd.i386-linux-ifort   =  -O0 -w95
FFLAGS_csvdc.i386-linux-ifort  = -g -w95
CFLAGS.i386-linux-ifort        =
MISCLIB.i386-linux-ifort       = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.i386-linux-ifort        = ranlib
MODULE.i386-linux-ifort      = -I$(MOD) -module $(MOD)
MODULE2.i386-linux-ifort       = -I ../$(MOD) -module ../$(MOD)

# DEC ALPHAstation  (Oct 2005)
# FFLAGS  = -O -inline all
# FFLAGS_warn= -O -inline all -C -warn declara -warn unused
# FFLAGS  = -g  -C
# FFLAGS_warn= -g  -C -warn declara -warn unused
alpha-f90:
FC.alpha-f90 = f90
FC90.alpha-f90 = f90
FFLAGS.alpha-f90 = -O -inline all
FFLAGS90.alpha-f90 = -O -inline all
FFLAGS_warn.alpha-f90 = -O -inline all  -C -warn declara -warn unused
#FFLAGS.alpha-f90 = -g -C
#FFLAGS90.alpha-f90 = -g -C
#FFLAGS_warn.alpha-f90 = -g -C -warn declara -warn unused
FFLAGS_GAMA.alpha-f90      = -O  -r8 -align dcommons -inline all
FFLAGS_GAMA_warn.alpha-f90 = -O -r8 -align dcommons -C 
FFLAGS_oasd.alpha-f90 = -O
FFLAGS_csvdc.alpha-f90 = -O
MODULE.alpha-f90 = -module $(MOD) -I$(MOD)
MODULE2.alpha-f90 = -I../$(MOD) -module ../$(MOD)
MISCLIB.alpha-f90      = 
RANLIB.alpha-f90 = ranlib

# SunOS solaris 2.4, proposed by Jim Murray, july 2000:
# the subroutine oases/oasdvms11.f must be compiled just optimized (FFLAGS= -O  ),
# as the compiler cannot handle it, Thus start by 
# compiling this program by manually.# use " cd oases;  f77 -c -O oasdvms11.f -o oasdvms11.o; cd .. "
#FFLAGS = -fast  -xarch=v8plusa -xchip=ultra 
#FFLAGS_warn  = -fast -u  -xarch=v8plusa -xchip=ultra 
# oct2003: sun requires all module files also included as an inline *.o file
sun4-solaris:
FC.sun4-solaris 	  = f90  -Bstatic
FC90.sun4-solaris 	  = f90   -Bstatic
#FFLAGS.sun4-solaris 	  = -fast  
#FFLAGS90.sun4-solaris 	  = -fast  
#FFLAGS_warn.sun4-solaris  = -fast  -C -u
FFLAGS.sun4-solaris 	  =  
FFLAGS90.sun4-solaris 	  =  
FFLAGS_warn.sun4-solaris  = -C -u
#FFLAGS.sun4-solaris 	  = -g  
#FFLAGS_warn.sun4-solaris  = -g  -C -u
FFLAGS_oasd.sun4-solaris  =  -O
FFLAGS_csvdc.sun4-solaris =  -O
CFLAGS.sun4-solaris    	= -fast
MODULE.sun4-solaris = -M$(MOD) -I$(MOD)
MODULE2.sun4-solaris = -I../$(MOD) -M../$(MOD)
MODULE3.sun4-solaris = -M
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
FC.i386-linux-linux	      = g77 -ff90
FFLAGS.i386-linux-linux       = -O2 -ffast-math
FFLAGS_warn.i386-linux-linux  = -O2 -ffast-math -Wimplicit
#FFLAGS.i386-linux-linux      = -g
#FFLAGS_warn.i386-linux-linux = -g
FFLAGS_oasd.i386-linux-linux  =  -O
FFLAGS_csvdc.i386-linux-linux =  
CFLAGS.i386-linux-linux	      =
MISCLIB.i386-linux-linux      = ./obj/${HOSTTYPE}-${FORTRAN}/misc.a
RANLIB.i386-linux-linux       = ranlib

i386-linux-pgf77:
FC.i386-linux-pgf77     =pgf90
#FFLAGS.i386-linux-pgf77 = -Mprof=func -fast -Mdalign  -tp p6
#FFLAGS.i386-linux-pgf77 = -Mprof=func -fast -Mdalign  -tp p6
#FFLAGS.i386-linux-pgf77 = -Mprof=lines -fastsse
#FFLAGS.i386-linux-pgf77 = -Mprof=lines -O2 
#FFLAGS.i386-linux-pgf77 = -O1 -Mvect -Minline -Mneginfo=loop
FFLAGS_GAMA_warn.i386-linux-pgf77 =   -O2  -r8 -fast -Mdalign -tp p6
FFLAGS_GAMA.i386-linux-pgf77 = -O2  -r8 -fast -Mdalign -tp p6
FFLAGS.i386-linux-pgf77 = -O2
FFLAGS_warn.i386-linux-pgf77  = -O2 -C
FFLAGS_oasd.i386-linux-pgf77 =  -O
FFLAGS_csvdc.i386-linux-pgf77 =
CFLAGS.i386-linux-pgf77      =
MISCLIB.i386-linux-pgf77     =
RANLIB.i386-linux-pgf77      = ranlib


# g77 does not support structures; therefore delete tpem from the "all" line
# the prosim and cprosim does not compile well. If these modules are required, please contact me. Use
#all :	 library  ramgeo  snap oast snaprd oasr oastg popp 
#
# the subroutine oases/csdvc.f must be compiled unoptimized (FFLAGS=  ),
# as the program othervice will crash during excecution, Thus start by 
# compiling this program by manually.
# use " cd oases;  g77 -c csvdc.f; cd .. "





# LINUX with f2c+gcc compiler (from 1997)
#FC  = f77
#FFLAGS  = -O 

CFLAGS          = $(CFLAGS.$(HOSTTYPE)-$(FORTRAN))
FC              = $(FC.$(HOSTTYPE)-$(FORTRAN)) 
FC90            = $(FC90.$(HOSTTYPE)-$(FORTRAN)) 
FFLAGS          = $(FFLAGS.$(HOSTTYPE)-$(FORTRAN))
FFLAGS90        = $(FFLAGS90.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_warn     = $(FFLAGS_warn.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_oasd     = $(FFLAGS_oasd.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_csvdc    = $(FFLAGS_csvdc.$(HOSTTYPE)-$(FORTRAN))
MISCLIB         = $(MISCLIB.$(HOSTTYPE)-$(FORTRAN))
RANLIB           = $(RANLIB.$(HOSTTYPE)-$(FORTRAN))
MODULE          = $(MODULE.$(HOSTTYPE)-$(FORTRAN))
MODULE2         = $(MODULE2.$(HOSTTYPE)-$(FORTRAN))
MODULE3         = $(MODULE3.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_GAMA     = $(FFLAGS_GAMA.$(HOSTTYPE)-$(FORTRAN))
FFLAGS_GAMA_warn= $(FFLAGS_GAMA_warn.$(HOSTTYPE)-$(FORTRAN))

LINK   = ${FC}  ${FFLAGS_warn} 
FSTR   = 'FC=${FC}' 'FC90=${FC90}' 'FFLAGS=${FFLAGS}' 'FFLAGS90=${FFLAGS90}' 'FFLAGS_warn=${FFLAGS_warn}' 'MODULE=$(MODULE2)'   'OBJ_opt=$(OBJ_opt)' 'RANLIB=${RANLIB}' 
FSTRG   = 'FC=${FC}' 'FFLAGS=${FFLAGS_GAMA}' 'FFLAGS_warn=${FFLAGS_GAMA_warn}' \
	'OBJ_opt=$(OBJ_opt)' 'RANLIB=${RANLIB}' 
CSTR   = 'CC=${CC}' 'CFLAGS=${CFLAGS}'
RM     = rm -f

BIN       = ../bin/${HOSTTYPE}-${FORTRAN}/
LIB       = ./obj/${HOSTTYPE}-${FORTRAN}/
MOD       = ./obj/${HOSTTYPE}-${FORTRAN}/

ford      = $(FC)  $(FFLAGS_warn) -I$(MOD) $(MODULE)  -c
fort      = $(FC)  $(FFLAGS) -I$(MOD) $(MODULE)  -c
object    = $(OBJ_opt)  $(cdir)
MATHLIB   = $(LIB)apmath.a
OASTLIB   = $(LIB)oaslib.a
SAGALIB   = $(LIB)libsaga.a
LIB1      = $(SAGALIB)            $(MATHLIB) $(MISCLIB) 
LIB2      = $(SAGALIB) $(OASTLIB) $(MATHLIB) $(MISCLIB) 
LIB3      = $(OASTLIB)            $(MATHLIB) $(MISCLIB) 
INCLFIL   = 	comforw.h comopt.h  $(cdir)global.o

$(cdir)%.o : %.f $(INCLFIL)  		
	$(ford) -I$(MOD) -c $*.f  $(OBJ_opt) $@
$(cdir)%.o : %.f90
	$(FC90) $(FFLAGS90) $(MODULE)  -c $*.f90  $(OBJ_opt) $@

# --------

posto:	library  dumtpem dumramg $(BIN)posttpem  \
        dumsnap $(BIN)postsnap $(BIN)postramgeo
library: ${HOSTTYPE}-${FORTRAN} objdir $(LIB)libsaga.a dumoases  $(MISCLIB) 
prosim:  library dumpro    $(BIN)sagaprosim  $(BIN)postprosim
cprosim: library dumcpro   $(BIN)sagacprosim $(BIN)postcprosim
gama:    library dumgama   $(BIN)sagagama    $(BIN)postgama 
orca:    library dumorca   $(BIN)sagaorca    $(BIN)postorca 
popp:	 library dumpop    $(BIN)sagapopp    $(BIN)postpopp
snap:	 library dumsnap   $(BIN)sagasnap    $(BIN)postsnap
snaprd:	 library dumsnaprd $(BIN)sagasnaprd  $(BIN)postsnaprd
oast:	 library dumoases  $(BIN)sagaoast    $(BIN)postoast 
oastg:	 library dumoases  $(BIN)sagaoastg   $(BIN)postoast
oasr:	 library dumoases  $(BIN)sagaoasr    $(BIN)postoasr
tpem:	 library dumtpem   $(BIN)sagatpem    $(BIN)posttpem
ramgeo:	 library dumramg   $(BIN)sagaramgeo  $(BIN)postramgeo
tabu0:	 dumsnap   $(BIN)tabu
resample:	 dumsnap dumsnaprd dumramg	$(BIN)resamplsnap $(BIN)resamplsnaprd	$(BIN)resamplramgeo $(BIN)resamplorca
#resaramgeo:	 dumramg	$(BIN)resamplramgeo

seos0:	 dumsnap   $(BIN)seos

OBJO =	$(cdir)conplot.o $(cdir)gain.o $(cdir)gaopt.o $(cdir)gasub.o\
       $(cdir)gaunew.o  $(cdir)plotga.o $(cdir)sa.o $(cdir)vfsa.o \
	$(cdir)cost.o $(cdir)gaoptions.o $(cdir)lineplot.o  $(cdir)mh.o\
	$(cdir)powell.o	$(cdir)cov.o	 $(cdir)prior.o	$(cdir)writetrf.o
OBJ2 =	$(cdir)conplot.o $(cdir)gain.o   $(cdir)gaopt.o $(cdir)gasub.o\
        $(cdir)plotga.o $(cdir)sa.o    $(cdir)vfsa.o $(cdir)gnuncer.o\
	$(cdir)cost.o $(cdir)gaoptions.o $(cdir)lineplot.o  $(cdir)mh.o\
	$(cdir)powell.o $(cdir)cov.o	$(cdir)prior.o
OBJ1 =	$(cdir)conplotgrad.o $(cdir)gain.o   $(cdir)gaopt.o $(cdir)gasub.o\
        $(cdir)plotga.o  $(cdir)sa.o     $(cdir)vfsa.o  $(cdir)cost.o \
	$(cdir)gnuncer.o $(cdir)grad.o   $(cdir)gnmin.o  $(cdir)oastgrad.o\
	$(cdir)gaoptions.o $(cdir)lineplot.o $(cdir)mh.o\
	$(cdir)powell.o  $(cdir)cov.o	$(cdir)prior.o
SGAMALIB =$(cdir)gamainit.o $(cdir)gamainter.o $(LIB1) $(LIB)libgama.a $(FLIB)
SORCALIB =$(cdir)orcainit.o $(cdir)orcainter.o $(LIB)liborca.a $(LIB1) $(FLIB)
SCPROLIB =$(cdir)cprosiminit.o $(cdir)cprosiminter.o $(LIB1) $(LIB)cprosim.a $(FLIB)
SPROSLIB =$(cdir)prosiminit.o $(cdir)prosiminter.o $(LIB1) $(LIB)prosim.a $(FLIB)
SPOPPLIB =$(cdir)poppinit.o   $(cdir)poppinter.o $(LIB1) $(LIB)libpopp.a $(FLIB)
SSNAPLIB =$(cdir)snapinit.o   $(cdir)snapinter.o   \
 $(LIB1)  $(LIB)snap.a $(FLIB)
SSNRDLIB =$(cdir)snaprdinit.o $(cdir)snaprdinter.o $(LIB1) $(LIB)snaprd.a $(FLIB)
STPEMLIB =$(cdir)tpeminit.o   $(cdir)tpeminter.o    $(OBJ2) $(cdir)gradRFC.o \
      $(cdir)gnmin.o $(OASTLIB) $(MATHLIB) $(MISCLIB) $(LIB)tpem.a $(FLIB) 
SRAMGLIB =$(cdir)ramgeoinit.o   $(cdir)ramgeointer.o   $(cdir)writetrf.o $(LIB1)\
   $(OASTLIB) $(MATHLIB) $(MISCLIB)  $(LIB)ramgeo.a $(FLIB) 
SOASTLIB =$(cdir)oastsub.o    $(cdir)oast11.o $(cdir)oasinter.o $(LIB2) $(FLIB)
SOATGLIB =$(cdir)oast11.o     $(cdir)oasinter.o  $(OBJ1)  $(LIB3)    $(FLIB)
SOASRLIB =$(cdir)oasr11.o     $(cdir)oasinter.o  $(LIB2)  $(FLIB)

$(LIB)libsaga.a:	$(OBJO)
	ar ru $(LIB)libsaga.a $(OBJO)
	@$(RANLIB) $(LIB)libsaga.a 

$(BIN)tabu:	      $(cdir)tabu.o $(SSNAPLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)tabu.o $(SSNAPLIB)
$(BIN)resamplsnap:	$(cdir)resamp.o	$(SSNAPLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)resamp.o $(SSNAPLIB)
$(BIN)resamplsnaprd:	$(cdir)resamp.o	$(SSNRDLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)resamp.o $(SSNRDLIB) 
$(BIN)resamplramgeo:	$(cdir)resamp.o	$(SRAMGLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)resamp.o $(SRAMGLIB)
$(BIN)resamplorca:	$(cdir)resamp.o	$(SORCALIB) $(cdirorca)parms_com.o  $(cdirorca)gen_com.o $(cdirorca)i_o_com.o $(cdirorca)scairy_com.o $(cdirorca)fgs_com.o $(cdir)global.o
	$(LINK) -o $@ $(cdir)global.o $(cdir)resamp.o $(SORCALIB) $(cdirorca)parms_com.o  $(cdirorca)gen_com.o $(cdirorca)i_o_com.o $(cdirorca)scairy_com.o $(cdirorca)fgs_com.o
$(BIN)seos:	      $(cdir)seos.o
	$(LINK) -o $@ $(cdir)global.o $(cdir)seos.o
$(BIN)sagaprosim:     $(cdir)ga.o $(SPROSLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SPROSLIB)
$(BIN)postprosim:     $(cdir)post.o $(SPROSLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SPROSLIB)
$(BIN)sagacprosim:    $(cdir)ga.o $(SCPROLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SCPROLIB)
$(BIN)postcprosim:    $(cdir)post.o $(SCPROLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SCPROLIB)
$(BIN)sagagama:	      $(cdir)ga.o $(SGAMALIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SGAMALIB)
$(BIN)postgama:	      $(cdir)post.o $(SGAMALIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SGAMALIB)
$(BIN)sagaorca:	      $(cdir)ga.o $(SORCALIB) $(cdirorca)parms_com.o  $(cdirorca)gen_com.o $(cdirorca)i_o_com.o $(cdirorca)scairy_com.o $(cdirorca)fgs_com.o $(cdir)global.o
	$(LINK) -o $@ $(cdirorca)parms_com.o  $(cdirorca)gen_com.o $(cdirorca)i_o_com.o $(cdirorca)scairy_com.o $(cdirorca)fgs_com.o $(cdir)global.o $(cdir)ga.o $(SORCALIB) 
$(BIN)postorca:	$(cdir)post.o $(SORCALIB)
	$(LINK)  -o $@ $(cdirorca)parms_com.o  $(cdirorca)gen_com.o $(cdirorca)i_o_com.o $(cdirorca)scairy_com.o $(cdirorca)fgs_com.o $(cdir)liborca.a $(cdir)global.o $(cdir)post.o $(SORCALIB) 
$(BIN)sagapopp:	      $(cdir)ga.o $(SPOPPLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SPOPPLIB)
$(BIN)postpopp:	      $(cdir)post.o $(SPOPPLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SPOPPLIB)
$(BIN)sagasnap:	      $(cdir)ga.o $(SSNAPLIB) 
	$(LINK)  -o $@ $(cdir)global.o $(cdir)ga.o $(SSNAPLIB)
$(BIN)postsnap:	      $(cdir)post.o $(SSNAPLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SSNAPLIB)
$(BIN)sagasnaprd:     $(cdir)ga.o $(SSNRDLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SSNRDLIB)
$(BIN)postsnaprd:     $(cdir)post.o $(SSNRDLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SSNRDLIB)
$(BIN)sagaoast:	      $(cdir)ga.o $(SOASTLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SOASTLIB)
$(BIN)postoast:	      $(cdir)post.o $(SOASTLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SOASTLIB)
$(BIN)sagaoastg:      $(cdir)ga.o $(SOATGLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SOATGLIB)
$(BIN)sagatpem:	      $(cdir)ga.o $(STPEMLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(STPEMLIB)
$(BIN)posttpem:	      $(cdir)post.o $(STPEMLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(STPEMLIB)
$(BIN)sagaramgeo:     $(cdir)ga.o $(SRAMGLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SRAMGLIB)
$(BIN)postramgeo:     $(cdir)post.o $(SRAMGLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SRAMGLIB)
$(BIN)sagaoasr:	      $(cdir)ga.o $(SOASRLIB) 
	$(LINK) -o $@ $(cdir)global.o $(cdir)ga.o $(SOASRLIB)
$(BIN)postoasr:	      $(cdir)post.o $(SOASRLIB)
	$(LINK) -o $@ $(cdir)global.o $(cdir)post.o $(SOASRLIB)

#$(cdir)mail_peter:	mail_peter
#	cp mail_peter $(cdir)mail_peter
#	./mail_peter
clean: 
	/bin/rm -rf  obj/*/* ../bin/*/*
	find . \( -name \*~ -o -name \#\* -o -name \*ext -o -name  \*trf -o -name core -o -name \*.o -o -name \*.mod  -o -name  \*.MOD -o -name \*.a \)\
        -exec /bin/rm {} \; -print	
	cd ..
	@echo ''
	@echo '**** SAGA distribution is clean now ***'*
	@echo ''
cleanone: 
	/bin/rm  -rf obj/*/* ../bin/${HOSTTYPE}-${FORTRAN}/* \
                  */${HOSTTYPE}-${FORTRAN}/* ifc?????? */ifc??????
	cd ..
	@echo ''
	@echo '**** SAGA distribution is clean now ***'*
	@echo ''

$(cdir)plotga.o: plotga.f	$(INCLFIL)  comoas.h oases/complo.f 
	$(ford) -I$(MOD)  -c plotga.f  $(OBJ_opt) $@
$(cdir)sa.o	    : sa.f  $(INCLFIL) 
	$(fort) sa.f   		$(object)sa.o
$(cdir)gaopt.o:  gaopt.f comopt.h graycode.h $(cdir)global.o
	$(ford) -I$(MOD) -c  gaopt.f  $(OBJ_opt) $@
$(cdir)gaoptions.o  :gaoptions.f comopt.h $(cdir)global.o
	$(ford) -I$(MOD)  -c gaoptions.f  $(OBJ_opt) $@
$(cdir)gradRFC.o    : gradRFC.f	$(INCLFIL) comtpem.h 
	$(ford) -I$(MOD)  -c gradRFC.f  $(OBJ_opt) $@
$(cdir)gnmin.o      :gnmin.f	$(INCLFIL)  comgrad.h 
	$(ford) -I$(MOD)  -c gnmin.f  $(OBJ_opt) $@
$(cdir)gaunew.o     :gaunew.f	comgrad.h
	$(ford) -I$(MOD)  -c gaunew.f  $(OBJ_opt) $@
$(cdir)oast11.o :oast11.f   $(INCLFIL) comoas.h oases/compar.f  \
	oases/comnp.f oases/comnla.f oases/comnrd.f oases/comfip.f 
	$(fort) oast11.f $(object)oast11.o
$(cdir)oastgrad.o : oastgrad.f   $(INCLFIL) comoas.h \
	oases/compar.f oases/comnp.f oases/comnla.f \
	oases/comnrd.f oases/comfip.f comgrad.h   
	$(fort) oastgrad.f $(object)oastgrad.o
$(cdir)oasr11.o   : oasr11.f  $(INCLFIL) comoas.h  \
	oases/compar.f oases/comnp.f  oases/comnla.f oases/comnrd.f  
	$(fort) oasr11.f $(object)oasr11.o
$(cdir)oasinter.o  : oasinter.f   $(INCLFIL) comoas.h \
	oases/compar.f oases/comnla.f
	$(fort) oasinter.f $(object)oasinter.o
$(cdir)oastsub.o    :	oastsub.f  
	$(fort) oastsub.f $(object)oastsub.o
$(cdir)writetrf.o : writetrf.f  $(INCLFIL) comsnap.h 
	$(ford) -I$(MOD)  -c writetrf.f  $(OBJ_opt) $@
$(cdir)snapinit.o: snapinit.f $(INCLFIL) comsnap.h ./snap/a.f ./snap/common.f  
	$(ford) -I$(MOD)  -c snapinit.f  $(OBJ_opt) $@
$(cdir)snapinter.o: snapinter.f $(INCLFIL) ./snap/common.f 
	$(ford) -I$(MOD)  -c snapinter.f  $(OBJ_opt) $@
$(cdir)tpeminit.o:  tpeminit.f $(INCLFIL) comtpem.h  tpem/tpem.inc  
	$(ford) -I$(MOD)  -c tpeminit.f  $(OBJ_opt) $@
$(cdir)tpeminter.o:  tpeminter.f $(INCLFIL) comtpem.h tpem/tpem.inc  
	$(ford) -I$(MOD)  -c tpeminter.f  $(OBJ_opt) $@
$(cdir)ramgeoinit.o:  ramgeoinit.f $(INCLFIL) comramgeo.h ramgeo/ram.h
	$(ford) -I$(MOD)  -c ramgeoinit.f  $(OBJ_opt) $@
$(cdir)ramgeointer.o: ramgeointer.f $(INCLFIL) comramgeo.h ramgeo/ram.h
	$(ford) -I$(MOD)  -c ramgeointer.f  $(OBJ_opt) $@
$(cdir)poppinit.o    :	poppinit.f $(INCLFIL) compopp.h   
	$(ford) -I$(MOD)  -c poppinit.f  $(OBJ_opt) $@
$(cdir)poppinter.o   : poppinter.f $(INCLFIL)	 compopp.h
	$(ford) -I$(MOD)  -c poppinter.f  $(OBJ_opt) $@
$(cdir)snaprdinit.o  :	snaprdinit.f $(INCLFIL) comsnaprd.h  
	$(ford) -I$(MOD)  -c snaprdinit.f  $(OBJ_opt) $@
$(cdir)snaprdinter.o: snaprdinter.f $(INCLFIL) comsnaprd.h
	$(fort) snaprdinter.f 	$(object)snaprdinter.o
$(cdir)prosiminit.o : prosiminit.f $(INCLFIL) comprosim.h\
	prosim/Parms_com  prosim/i_o_2_com   prosim/sector_env_com \
        prosim/depth_com  prosim/i_o_1tf_com prosim/deltaf_com 
	$(ford) -I$(MOD)  -c prosiminit.f  $(OBJ_opt) $@
$(cdir)prosiminter.o: prosiminter.f	$(INCLFIL) comprosim.h\
	prosim/Parms_com  prosim/i_o_2_com   prosim/sector_env_com
	$(ford) -I$(MOD)  -c  prosiminter.f  $(OBJ_opt) $@
$(cdir)cprosiminit.o : cprosiminit.f	$(INCLFIL)  gen_i_o_saga.h\
	cprosim/Parms_com cprosim/i_o_saga_com cprosim/sector_env_com \
        cprosim/depth_com      
	$(ford) -I$(MOD) -c cprosiminit.f  $(OBJ_opt) $@
$(cdir)cprosiminter.o: cprosiminter.f	$(INCLFIL) gen_i_o_saga.h\
	cprosim/Parms_com cprosim/i_o_saga_com cprosim/sector_env_com
	$(ford) -I$(MOD)  -c cprosiminter.f  $(OBJ_opt) $@
$(cdir)orcainit.o: orcainit.f  	$(INCLFIL)
	$(ford) $(MODULE3)orca90   orcainit.f 	$(object)orcainit.o
$(cdir)orcainter.o: orcainter.f $(INCLFIL)
	$(ford) $(MODULE3)orca90 orcainter.f 	$(object)orcainter.o

dumoases  :
	cd oases;  make $(FSTR) 'FFLAGS_csvdc=$(FFLAGS_csvdc)' \
           'FFLAGS_oasd=$(FFLAGS_oasd)';
	cd misc;   make $(FSTR) $(CSTR)
dumcpro   :
	cd cprosim; make $(FSTR)
dumpro    :
	cd prosim; make $(FSTR)
dumpop    :
	cd popp;   make $(FSTR)
dumsnap   :
	cd snap;   make $(FSTR)
dumsnaprd :
	cd snaprd; make $(FSTR)
dumtpem   :
	cd tpem;   make $(FSTR)
dumramg   :
	cd ramgeo; make $(FSTR)
dumgama   :
	cd gama;   make  $(FSTRG)
dumorca   :
	cd orca90; make $(FSTR)
