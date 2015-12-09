# -----------------------------------------------------------------------------
# CPROSIM library for SAGA
# -----------------------------------------------------------------------------


# Suffix rules
# ------------
#.SUFFIXES :
#.SUFFIXES : .f .c .o .prj
#.f.o:;		$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@; mv $*.o $(OBJ)
#.c.o:;		$(CC) $(CFLAGS) -c $*.c -o $@
#.F.o:;		$(FC) $(FFLAGS) -D$(DFLAGS) -c $*.F $@
#.f.prj:;	ftnchek -project -noextern -library $<

TARGET = ..
DFLAGS = 
OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
XLIB   = $(LIB)cprosim.a 
PROG   = cprosim
MV = mv
RM     = rm -f
RANLIB = ranlib
.SUFFIXES :
$(OBJ)%.o : %.f  		
	$(FC) $(FFLAGS) -c $*.f  $(OBJ_opt) $@
$(OBJ)%.o : %.f90	
	$(FC90) $(FFLAG90) -c $*.f90  $(OBJ_opt) $@


OBJS =  $(OBJ)airy_hsp.o $(OBJ)bb_field.o $(OBJ)bb_init.o $(OBJ)bb_fmarch.o $(OBJ)bb_mloop.o \
        $(OBJ)bb_modes.o $(OBJ)blug.o $(OBJ)contour_find.o $(OBJ)cw_modes.o $(OBJ)duct_calcs.o \
        $(OBJ)eig_findm.o $(OBJ)em_calc.o $(OBJ)ep_calc.o $(OBJ)bb_brute.o \
        $(OBJ)fix_path.o $(OBJ)freq_init.o $(OBJ)h_space.o $(OBJ)k_plane.o $(OBJ)tenv_read_pro.o \
        $(OBJ)mode_branch.o $(OBJ)mfun_fill.o $(OBJ)mode_field.o $(OBJ)mode_find.o \
        $(OBJ)mode_fun.o $(OBJ)mode_ortho.o $(OBJ)opt_read_pro.o $(OBJ)out_writex.o $(OBJ)phz_calc.o \
        $(OBJ)r1r2_calc.o $(OBJ)rp_calc.o $(OBJ)mode_traj.o $(OBJ)rp_flay.o $(OBJ)rp_nomm.o \
        $(OBJ)rp_nomm0.o $(OBJ)rp_sfint.o $(OBJ)rp_slay.o $(OBJ)rp_slay0.o $(OBJ)mode_cmplx_field.o \
        $(OBJ)scairy2.o $(OBJ)sr_geom.o $(OBJ)svp_check.o $(OBJ)svp_read_pro.o \
       $(OBJ)tools.o $(OBJ)traj_mag.o $(OBJ)vbl_init.o \
        $(OBJ)xkh_init.o $(OBJ)zmx_init.o \
        $(OBJ)xi_cut_set.o \
        $(OBJ)bb_field_reuse.o $(OBJ)clean_array.o
COMS = 	Parms_com  depth_com  gen_cpro_com  i_o_com  i_o_saga_com\
	scairy_com debug_com  gen_com    gen_saga_com  i_o_cpro_com  lab_com\
	sector_env_com

all:	$(XLIB) 

$(XLIB): $(OBJS)
	ar rvu $(XLIB)  $(OBJS)
	$(RANLIB) $(XLIB)
#	$(MV) *.o $(LIB)

$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) $(PROG).o $(XLIB) 

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)


