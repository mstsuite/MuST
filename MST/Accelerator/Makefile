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
acclib.o \
util_mat.o \
LSMS_Accel.o \
cusolver_KKR_c.o \
invertMatrixKKR_CUDA.o

$(ACCLIB): $(OBJ0)
	$(ARCHV) $@ $(OBJ0)

clean:
	cd $(ODIR); rm -f $(OBJ0) $(ACCLIB)
