#######################################################################
#
#  This makefile creates the Fortran example programs for the
#  linear equation routines in SuperLU_DIST.
#
#  Creation date:   July 29, 2003   version 2.0
#
#######################################################################
VPATH = $(ODIR)/../slulib

MODULE = superlupara.o superlu_mod.o SuperLUSolverModule.o

LOBS = IOroutines.o

WRAPPER = superlu_c2f_z_wrap.o

.SUFFIXES: .F90 .f90 .c .o

include $(SLUPATH)/make.inc
INCLUDEDIR = -I$(SLUPATH)/SRC -I$(MPI2INCLUDE_PATH)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -DNO_MPI2 -c $< $(VERBOSE)

	
.f90.o:
	$(FC) $(INCLUDEDIR) $(FFLAGSSLU) $< $(VERBOSE)

.F90.o:
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
	fi ;

$(SLULIB): $(MODULE) $(LOBS) $(WRAPPER)
	cd $(ODIR)
	$(ARCHV) $(SLULIB) $(MODULE) $(LOBS) $(WRAPPER)
	ranlib $(SLULIB)

clear:
	cd $(ODIR); rm -f $(MODULE) $(LOBS) $(SLULIB) $(WRAPPER) superlu*


