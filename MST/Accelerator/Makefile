#*********************************************************************
#    makefile : generate library: libMatAcc.a
#
#    Book-keeping:      date                       note
#                     04-3-13
#*********************************************************************
VPATH = $(TOP_DIR)/Accelerator
INC_PATH = $(VPATH)

.SUFFIXES:	.o .f .F .f90 .F90 .c .cu .cpp

%.o: %.F90
	if [ "$(PREPROC)" != "$(FC)" ]; then \
           ( if [ $(XLC_I) != 1 ] ; then \
              $(PREPROC) $(PREPROC_FFLAGS) $< > $(ODIR)/$*.f90; \
           else \
              $(PREPROC) $(PREPROC_FFLAGS) $< ; mv $*.i $(ODIR)/$*.f90; \
           fi ); \
           cd $(ODIR); $(FC) $(FFLAGS) $*.f90; \
           (if [ $(No_Remove) != 1 ] ; then rm $*.f90; fi); \
        else \
	   $(FC) $(PREPROC_FFLAGS) $(FFLAGS) $< ;\
	fi

%.o: %.F
	if [ "$(PREPROC)" != "$(FC)" ]; then \
           ( if [ $(XLC_I) != 1 ] ; then \
              $(PREPROC) $(PREPROC_FFLAGS) $< > $(ODIR)/$*.f; \
           else \
              $(PREPROC) $(PREPROC_FFLAGS) $< ; mv $*.i $(ODIR)/$*.f; \
           fi ); \
           cd $(ODIR); $(FC) $(F77FLAGS) $*.f; \
           (if [ $(No_Remove) != 1 ] ; then rm $*.f; fi); \
        else \
	   $(FC) $(PREPROC_FFLAGS) $(F77FLAGS) $< ;\
	fi

%.o: %.c
ifdef ACCEL_CXX
	$(CC) $(CFLAGS) -I$(INC_PATH) $<
else
	touch $@
endif

%.o: %.cpp
ifdef ACCEL_CXX
	$(CXX) $(CFLAGS) -I$(INC_PATH) $<
else
	touch $@
endif

%.o: %.cu
ifdef ACCEL_CXX
	$(ACCEL_CXX) $(ACCEL_OPT) -I$(INC_PATH) $<
else
	touch $@
endif

OBJ0 = \
initialize_matinv_accel.o \
finalize_matinv_accel.o \
initialize_matinv_cuda.o \
finalize_matinv_cuda.o \
zblock_lu_accel.o \
zblock_lu_cuda_c.o \
zblock_lu_CPU.o \
inverse.o \
DeviceStorage.o \
interface_dstore.o

$(ACCLIB): $(OBJ0)
	$(ARCHV) $@ $(OBJ0)

clean:
	cd $(ODIR); rm -f $(OBJ0) $(ACCLIB)
