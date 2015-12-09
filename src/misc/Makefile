# -----------------------------------------------------------------------------
# Misc library for SAGA - Linux support
# -----------------------------------------------------------------------------

OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIBPATH =  ../obj/${HOSTTYPE}-${FORTRAN}/
MV = mv

.SUFFIXES :
$(OBJ)%.o : %.f  		
	$(FC) $(FFLAGS90) -c $*.f  $(OBJ_opt) $@
$(OBJ)%.o : %.c  		
	$(CC) $(CFLAGS) -c $*.c  $(OBJ_opt) $@
$(OBJ)%.o : %.f90	
	$(FC) $(FFLAG90) -c $*.f90  $(OBJ_opt) $@

MLIB  = misc.a
MXLIB  = $(LIBPATH)$(MLIB)
# ---------------------------------------------------------------------
all:	$(MXLIB) 

OBJS =	$(OBJ)itime.o $(OBJ)etime.o $(OBJ)idate.o $(OBJ)getpid.o $(OBJ)lenstr2.o  $(OBJ)getcmds.o  \
	$(OBJ)ran.o	$(OBJ)fork.o  $(OBJ)wait.o \
	$(OBJ)cdabs.o $(OBJ)cdcos.o $(OBJ)cdexp.o $(OBJ)cdsin.o  $(OBJ)cdsqrt.o   $(OBJ)cdsrt.o   $(OBJ)cdtan.o
#$(OBJ)chmod.o
$(MXLIB):	$(OBJS)
		ar ru $(MXLIB) $(OBJS)

clean: 
	$(RM) *~ *.*~ core *.o *.prj *.a *.pl?
