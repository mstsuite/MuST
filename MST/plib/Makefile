#*********************************************************************
#    makefile : generate library: mstlib.a
#
#    Book-keeping:      date                       note
#                     10-14-97          
#                     04-23-01          
#*********************************************************************
VPATH = $(TOP_DIR)/plib

MODULE = \
MPPModule.o \
GroupCommModule.o \
ParallelIOModule.o

.SUFFIXES:	.o .f .f90 .F90

.f.o:
	$(FC) $(FFLAGS) $<

.f90.o:
	$(FC) $(FFLAGS) $<

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

#=====================================================================
$(MPPLIB): $(MODULE)
	cd $(ODIR); $(ARCHV) $(MPPLIB) $(MODULE)
	ranlib $(MPPLIB)

clear: $(MODULE)
	cd $(ODIR); rm -f $(MODULE)
#=====================================================================
