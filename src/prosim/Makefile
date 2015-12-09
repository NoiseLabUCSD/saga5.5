#FFLAGS = -g
#FFLAGS = -g -C -V
#FFLAGS =  -O4   -fast   -inline all 
#FFLAGS = -g -Xlist
#FFLAGS = -g  -C 
#FFLAGS = -O2 -C  -inline all 
RM     = rm -f
MV     = mv
RANLIB = ranlib

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
XLIB   = $(LIB)prosim.a 
MV = mv
PROG =prosim

OBJO =	$(OBJ)bb_init.o $(OBJ)blug.o $(OBJ)opt_read.o  \
        $(OBJ)rx_mode_field_tilt.o $(OBJ)rx_mode_save.o \
	$(OBJ)sr_geom_tilt.o $(OBJ)svp_check.o $(OBJ)svp_read.o \
	$(OBJ)rx_bb_tool.o $(OBJ)rx_bb.o \
	$(OBJ)rx_ek_calc.o $(OBJ)rx_mf_lay.o $(OBJ)rx_mode_fun_Dw.o\
	$(OBJ)rx_modes.o \
	$(OBJ)tenv_read_tilt.o $(OBJ)svp_adj.o $(OBJ)bb_fft_out.o

COM = Parms_com deltaf_com depth_com freq2_com freq_com gen3_com\
gen1_com i_o_svp_com i_o_opt_com i_o_1b_com i_o_1c_com i_o_1tf_com \
gen2_com i_o_2_com modrng_com sector_env_com
#

$(XLIB):$(OBJO)
	ar rvu $(XLIB) $(OBJO)
	$(RANLIB) $(XLIB)
#	$(MV) *.o $(LIB)

$(PROG):$(PROG).o $(XLIB)
	$(FC) -o  $(PROG) $(PROG).o $(XLIB) 
	strip $(PROG)

clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)

dircmp:
	@echo 'Tracing == ' `pwd`
	@trace

#


clean: 
	$(RM) *~ *.*~ core *.o *.prj .f $(PROG)



$(OBJ)bb_init.o:	bb_init.f $(COM)
	$(FC)  $(FFLAGS)  -c bb_init.f
	@ $(MV) ${@F} $@


$(OBJ)blug.o: blug.f $(COM)
	$(FC)  $(FFLAGS)  -c blug.f
	@ $(MV) ${@F} $@



$(OBJ)opt_read.o: opt_read.f $(COM)
	$(FC)  $(FFLAGS)  -c opt_read.f
	@ $(MV) ${@F} $@


$(OBJ)rx_mode_field_tilt.o: rx_mode_field_tilt.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_mode_field_tilt.f
	@ $(MV) ${@F} $@


$(OBJ)rx_mode_save.o: rx_mode_save.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_mode_save.f
	@ $(MV) ${@F} $@


$(OBJ)sr_geom_tilt.o: sr_geom_tilt.f $(COM)
	$(FC)  $(FFLAGS)  -c sr_geom_tilt.f
	@ $(MV) ${@F} $@


$(OBJ)svp_check.o: svp_check.f $(COM)
	$(FC)  $(FFLAGS)  -c svp_check.f
	@ $(MV) ${@F} $@


$(OBJ)svp_read.o: svp_read.f $(COM)
	$(FC)  $(FFLAGS)  -c svp_read.f
	@ $(MV) ${@F} $@


$(OBJ)rx_bb_tool.o: rx_bb_tool.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_bb_tool.f
	@ $(MV) ${@F} $@


$(OBJ)rx_bb.o: rx_bb.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_bb.f
	@ $(MV) ${@F} $@


$(OBJ)rx_ek_calc.o: rx_ek_calc.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_ek_calc.f
	@ $(MV) ${@F} $@


$(OBJ)rx_mf_lay.o: rx_mf_lay.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_mf_lay.f
	@ $(MV) ${@F} $@


$(OBJ)rx_mode_fun_Dw.o: rx_mode_fun_Dw.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_mode_fun_Dw.f
	@ $(MV) ${@F} $@

$(OBJ)rx_modes.o: rx_modes.f $(COM)
	$(FC)  $(FFLAGS)  -c rx_modes.f
	@ $(MV) ${@F} $@


$(OBJ)tenv_read_tilt.o: tenv_read_tilt.f $(COM)
	$(FC)  $(FFLAGS)  -c tenv_read_tilt.f
	@ $(MV) ${@F} $@


$(OBJ)svp_adj.o: svp_adj.f $(COM)
	$(FC)  $(FFLAGS)  -c svp_adj.f
	@ $(MV) ${@F} $@


$(OBJ)bb_fft_out.o: bb_fft_out.f $(COM)
	$(FC)  $(FFLAGS)  -c bb_fft_out.f
	@ $(MV) ${@F} $@

