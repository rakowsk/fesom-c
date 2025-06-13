#Makefile

include Makefile.in

.SUFFIXES: .f90 .F90 .o

LIBS = $(GOTMLIB)  $(NETCDFLIB)
INCS = $(GOTMINCS) $(NETCDFINC)

OBJECTS= fv_var.o               \
      fv_utilit.o               \
      fv_vert_coord.o           \
      fv_sbc.o                  \
      fv_obc.o                  \
      fv_obc_2d.o               \
      fv_ic.o                   \
      fv_read_run.o             \
      fv_ini.o                  \
      fv_mesh_array.o           \
      fv_gradients.o            \
      fv_3D_vel.o               \
      fv_advection.o            \
      fv_average_dynamic.o      \
      fv_mixing.o               \
      fv_output.o               \
      fv_pressure_new.o         \
      fv_viscosity.o            \
      fv_tracer.o               \
      fv_tide_potential.o       \
      fv_ncoutput.o             \
      fv_sediment.o      \
      fv_sediment_adv.o      \
      fv_main.o  


INCS_INI = $(METISINC)
LIBS_INI = $(METISLIB)

OBJ_C_INI = metis_wrapper.o  
SRC_INI =  fv_var.f90 \
	   fv_utilit.f90 \
	   fv_read_init.f90 \
	   fv_mesh_array.f90 \
	   fv_distribute_mesh.f90 \
           fv_comm.f90 \
           fv_init.f90


#--------------------------------------------------

$(EXE): $(OBJECTS)  $(GOTMDIR)/lib/libturbulence.a
	 $(LINK) $(FCFLAGS) -o $(EXE) $(OBJECTS) $(LIBS)

$(EXE_INI): $(OBJ_C_INI) $(SRC_INI)
	 $(LINK) $(FCFLAGS) -o $(EXE_INI) $(SRC_INI) $(OBJ_C_INI) $(LIBS_INI)

GOTM_turbulence: 
	$(MAKE) -C $(GOTMDIR) all

all:
	$(MAKE) GOTM_turbulence
	$(MAKE) $(EXE)

clean :
	rm -f *~ *.o *.mod $(EXE) $(EXE_INI)

allclean :
	$(MAKE) clean
	$(MAKE) -C $(GOTMDIR) clean


run : $(EXE)         
	./$(EXE)

#--------------------------------------------------

.f90.o: $(MODULES)
	$(FC)  $(FCFLAGS) $(INCS) -c $*.f90

.c.o:
	$(CC) $(CCFLAGS) $(INCS_INI) -c $*.c
