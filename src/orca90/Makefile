################################
############# ORCA #############
################################

#FC=/usr/bin/f90
# F90 FFLAGS:
#FFLAGS = -C -V 
#FFLAGS = -O2 -inline all -fast
TARGET = .
OBJ    = ${HOSTTYPE}-${FORTRAN}/
LIBPATH =  ../obj/${HOSTTYPE}-${FORTRAN}/
LIB    = $(TARGET)/obj/${HOSTTYPE}-${FORTRAN}/
BIN    = $(TARGET)/bin/
MOD    = ../obj/${HOSTTYPE}-${FORTRAN}/

$(OBJ)%.o : %.f  		
	$(FC) $(FFLAGS) -I $(MOD)  -c $*.f  $(OBJ_opt) $@
$(OBJ)%.o: %.f90	
	$(FC) $(FFLAG) -I $(MOD) -c $*.f90  $(OBJ_opt) $@
	mv $*.mod $(MOD)$*.mod

MODS = ${OBJ}parms_com.o ${OBJ}gen_com.o ${OBJ}i_o_com.o ${OBJ}scairy_com.o ${OBJ}fgs_com.o

# bb_brute.o rx_bb.o vary_env.o 
#$(OBJ)hdf_write_gen_dummy.o  $(OBJ)hdf_write_end_dummy.o 

OBJS =   $(OBJ)sr_geom.o $(OBJ)airy_hsp.o  $(OBJ)bb_field.o $(OBJ)bb_fmarch.o $(OBJ)bb_init.o \
$(OBJ)blug.o $(OBJ)contour_find.o $(OBJ)cw_modes.o $(OBJ)duct_calcs.o \
$(OBJ)eig_findm.o $(OBJ)eig_list.o $(OBJ)em_calc.o $(OBJ)env_write.o $(OBJ)ep_calc.o \
$(OBJ)fix_path.o $(OBJ)freq_init.o $(OBJ)h_space.o \
$(OBJ)k_plane.o $(OBJ)krak_out.o $(OBJ)mode_branch.o $(OBJ)mfun_fill.o $(OBJ)mode_field.o \
$(OBJ)mode_find.o $(OBJ)mode_fun.o $(OBJ)mode_ortho.o $(OBJ)opt_read.o $(OBJ)phz_calc.o \
$(OBJ)pw_refco.o $(OBJ)r1r2_calc.o $(OBJ)rp_calc.o $(OBJ)mode_traj.o $(OBJ)rp_flay.o $(OBJ)rp_nomm.o \
$(OBJ)rp_nomm0.o $(OBJ)rp_sfint.o $(OBJ)rp_slay.o $(OBJ)rp_slay0.o \
$(OBJ)scairy2.o $(OBJ)sigma_fill.o  $(OBJ)svp_check.o $(OBJ)svp_read.o \
$(OBJ)tools.o $(OBJ)traj_mag.o $(OBJ)vbl_init.o \
$(OBJ)xkh_init.o $(OBJ)zmx_init.o  $(OBJ)rx_ek_calc.o $(OBJ)rx_mf_lay.o \
$(OBJ)rx_mode_fun.o $(OBJ)rx_mode_fun_Dw.o $(OBJ)rx_modes.o $(OBJ)xi_cut_set.o $(OBJ)orca_sub.o
    
#             
all: ${MODS} ${LIBPATH}liborca.a

${LIBPATH}liborca.a: 	${MODS}	${OBJS}
	ar rv ${LIBPATH}liborca.a ${OBJS}



















