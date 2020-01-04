#*********************************************************************
#    makefile : generate library: mstlib.a
#
#    Book-keeping:      date                       note
#                     10-14-97          
#                     04-23-01          
#*********************************************************************
VPATH = $(ODIR)/../iolib

OBJS = iotk.o

.SUFFIXES:	.o .f .f90 .F90

.c.o:
	$(CC) $(CFLAGS) -I$(VPATH) $<
	if test "$(DEBUG)" = "1"; then cp $(VPATH)/$*.c $(ODIR) ; fi

.f.o:
	$(FC) $(F77FLAGS) $<
	if test "$(DEBUG)" = "1"; then cp $(VPATH)/$*.f $(ODIR) ; fi

.f90.o:
	$(FC) $(FFLAGS) $<
	if test "$(DEBUG)" = "1"; then cp $(VPATH)/$*.f90 $(ODIR) ; fi

.F90.o:
	if [ "$(PREPROC)" != "$(FC)" ]; then \
	   ( if [ $(XLC_I) != 1 ] ; then \
	        $(PREPROC) $(PREPROC_FFLAGS) $< > $(ODIR)/$*.f90; \
	     else \
	        $(PREPROC) $(PREPROC_FFLAGS) $< ; mv $*.i $(ODIR)/$*.f90; \
	     fi ); \
	   cd $(ODIR); $(FC) $(FFLAGS) $*.f90; \
	   ( if [ $(No_Remove) != 1 ] ; then rm $*.f90; fi); \
	else \
	   $(FC) $(PREPROC_FFLAGS) $(FFLAGS) $< ;\
	fi ;
	   
#=====================================================================
$(IOLIB): $(OBJS)
	$(ARCHV_LSMS) $(ARCHV_OPT) $(IOLIB) $(OBJS); ranlib $(IOLIB)

clear: $(OBJS)
	cd $(ODIR); rm -f $(OBJS)
#=====================================================================
