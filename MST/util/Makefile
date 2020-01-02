VPATH = $(ODIR)/../util

.SUFFIXES:	.o .f .c .f90 .F90

.f.o:
	cd $(ODIR); $(FC) $(F77FLAGS) $<

.F.o:
	if [ "$(PREPROC)" != "$(FC)" ]; then \
	   ( if [ $(XLC_I) != 1 ] ; then \
	      $(PREPROC) $(PREPROC_FFLAGS) $< > $(ODIR)/$*.f90; \
	   else \
	      $(PREPROC) $(PREPROC_FFLAGS) $< ; mv $*.i $(ODIR)/$*.f90; \
	   fi ); \
	   cd $(ODIR); $(FC) $(F77FLAGS) $*.f90; \
	   (if [ $(No_Remove) != 1 ] ; then rm $*.f90; fi); \
	else \
	      $(FC) $(PREPROC_FFLAGS) $(F77FLAGS) $< ;\
	fi ;

.f90.o:
	cd $(ODIR); $(FC) $(FFLAGS) $<

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

.c.o:
	cd $(ODIR); $(CC) $(CFLAGS) $<

UTILITIES = $(ODIR)/measureVoronoi \
            $(ODIR)/charge \
            $(ODIR)/chgaver \
            $(ODIR)/newss \
            $(ODIR)/newa \
            $(ODIR)/slice \
            $(ODIR)/extrpot \
            $(ODIR)/murn \
            $(ODIR)/eosfit \
            $(ODIR)/genap \
            $(ODIR)/setEvec \
            $(ODIR)/cart2direct \
            $(ODIR)/bireft

UTILITIES_O = $(ODIR)/measureVoronoi.o \
              $(ODIR)/charge.o \
              $(ODIR)/chgaver.o \
              $(ODIR)/newss.o \
              $(ODIR)/newa.o \
              $(ODIR)/slice.o \
              $(ODIR)/extrpot.o \
              $(ODIR)/murn.o \
              $(ODIR)/eosfit.o \
              $(ODIR)/SampleModule.o \
              $(ODIR)/genap.o \
              $(ODIR)/setEvec.o \
              $(ODIR)/cart2direct.o \
              $(ODIR)/bireft.o

SRC_OBJS = $(ODIR)/PublicTypeDefinitionsModule.o \
           $(ODIR)/PublicParamDefinitionsModule.o $(ODIR)/DefaultParamModule.o \
           $(ODIR)/InputModule.o $(ODIR)/ScfDataModule.o $(ODIR)/SystemModule.o \
           $(ODIR)/ProcMappingModule.o \
           $(ODIR)/PolyhedraModule.o $(ODIR)/Atom2ProcModule.o \
           $(ODIR)/PotentialTypeModule.o $(ODIR)/readInputInOtherFormat.o \
           $(ODIR)/readPositionData.o $(ODIR)/ChemElementModule.o \
           $(ODIR)/readOldInfoTable.o $(ODIR)/rdin_old_lsms.o \
           $(ODIR)/rdin_old_infoevec.o $(ODIR)/rdin_old_infotable.o \
           $(ODIR)/readMomentDirectionData.o $(ODIR)/lattice.o \
           $(ODIR)/MatrixInverseModule.o

all: $(UTILITIES) gaps

gaps:
	cd  $(VPATH)/gaps; make "MPICC=$(MPICC)" "CFLAGS=$(CFLAGS)" "ODIR=$(ODIR)" -f $(VPATH)/gaps/makefile

clear:
	rm -f $(UTILITIES_O)

clean:
	rm -f $(UTILITIES); make "ODIR=$(ODIR)" clear -f $(VPATH)/Makefile 
	cd $(VPATH)/gaps; make clean "ODIR=$(ODIR)" "MPICC=$(MPICC)" "CFLAGS=$(CFLAGS)" -f makefile

$(ODIR)/setEvec: setEvec.o $(SRC_OBJS) $(MSTLIB) $(MPPLIB)
	$(LD) -o $(ODIR)/setEvec $(ODIR)/setEvec.o \
$(SRC_OBJS) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/charge: charge.o
	$(LD) -o $(ODIR)/charge $(ODIR)/charge.o \
$(SRC_OBJS) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/genap: SampleModule.o generateAtomPosition.o
	$(LD) -o $(ODIR)/genap $(ODIR)/generateAtomPosition.o $(ODIR)/SampleModule.o \
$(SRC_OBJS) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/chgaver: chgaver.o
	$(LD) -o $(ODIR)/chgaver $(ODIR)/chgaver.o \
$(SRC_OBJS) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/newss: newss.o
	$(FC) -o $(ODIR)/newss $(ODIR)/newss.o $(ODIR)/ChemElementModule.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/newa: newa.o
	$(FC) -o $(ODIR)/newa $(ODIR)/newa.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/slice: slice.o
	$(LD) -o $(ODIR)/slice $(ODIR)/slice.o \
$(ODIR)/ChemElementModule.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/extrpot: extractPotential.o
	$(LD) -o $(ODIR)/extrpot $(ODIR)/extractPotential.o \
$(ODIR)/ChemElementModule.o $(IOLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/murn: murn_new.o
	$(LD) -o $(ODIR)/murn $(ODIR)/murn_new.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/eosfit: eosfit.o
	$(LD) -o $(ODIR)/eosfit $(ODIR)/eosfit.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/bireft: bireft.o
	$(LD) -o $(ODIR)/bireft $(ODIR)/bireft.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/measureVoronoi: measureVoronoi.o $(ODIR)/readInputs.o $(ODIR)/PolyhedraModule.o
	$(LD) -o $(ODIR)/measureVoronoi measureVoronoi.o $(ODIR)/readInputs.o \
$(SRC_OBJS) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/readInputs.o: ../driver/readInputs.F90
	$(FC) $(FFLAGS) $<

$(ODIR)/cart2direct: cart2direct.o $(MSTLIB)
	$(LD) -o $(ODIR)/cart2direct $(ODIR)/cart2direct.o $(MSTLIB) $(ADDLIBS)
