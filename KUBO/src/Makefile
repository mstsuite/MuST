VPATH = $(TOP_DIR)/src

BasicModules = \
../../MST/bin/PhysParamModule.o \
../../MST/bin/ChemElementModule.o \
../../MST/bin/PublicTypeDefinitionsModule.o \
../../MST/bin/PublicParamDefinitionsModule.o \
../../MST/bin/AdaptIntegrationModule.o \
../../MST/bin/MatrixBlockInversionModule.o \
../../MST/bin/BinvMatrix.o \
../../MST/bin/WriteFunctionModule.o \
../../MST/bin/DefaultParamModule.o \
../../MST/bin/LatticeModule.o \
../../MST/bin/IBZRotationModule.o

BasicRoutines = \
../../MST/bin/readPositionData.o \
../../MST/bin/readMomentDirectionData.o \
../../MST/bin/readInputInOtherFormat.o \
../../MST/bin/readOldInfoTable.o \
../../MST/bin/rdin_old_lsms.o \
../../MST/bin/rdin_old_infotable.o \
../../MST/bin/rdin_old_infoevec.o \
../../MST/bin/getValueAtPosi.o \
../../MST/bin/lattice.o

AppModules = \
../../MST/bin/InputModule.o \
../../MST/bin/ScfDataModule.o \
../../MST/bin/PotentialTypeModule.o \
../../MST/bin/SystemModule.o \
../../MST/bin/BZoneModule.o \
../../MST/bin/MatsubaraModule.o \
../../MST/bin/ContourModule.o \
../../MST/bin/ProcMappingModule.o \
../../MST/bin/Atom2ProcModule.o \
../../MST/bin/PolyhedraModule.o \
../../MST/bin/SystemVolumeModule.o \
../../MST/bin/OutputModule.o \
../../MST/bin/VPI_main.o \
../../MST/bin/IsoparametricIntegrationModule.o \
../../MST/bin/NeighborModule.o \
../../MST/bin/SendRecvTmatModule.o \
../../MST/bin/AtomModule.o \
../../MST/bin/StepFunctionModule.o \
../../MST/bin/MadelungModule.o \
../../MST/bin/SystemSymmetryModule.o \
../../MST/bin/LdaCorrectionModule.o \
../../MST/bin/RadialGridModule.o \
../../MST/bin/PotentialModule.o \
../../MST/bin/StrConstModule.o \
../../MST/bin/RSpaceStrConstModule.o \
../../MST/bin/TriangleModule.o \
../../MST/bin/PolygonModule.o \
../../MST/bin/SurfElementsModule.o \
../../MST/bin/SpinRotationModule.o \
../../MST/bin/SSSolverModule.o \
../../MST/bin/DiracSolverModule.o \
../../MST/bin/RelScattererModule.o \
../../MST/bin/RelSSSolverModule.o \
../../MST/bin/SCPolesModule.o \
../../MST/bin/SMatrixPolesModule.o \
../../MST/bin/KreinModule.o \
../../MST/bin/CoreStatesModule.o \
../../MST/bin/Uniform3DGridModule.o \
../../MST/bin/ParallelFFTModule.o \
../../MST/bin/ConstrainLocalMomentModule.o \
../../MST/bin/ExchCorrFunctionalModule.o \
../../MST/bin/ValenceDensityModule.o \
../../MST/bin/ChargeDensityModule.o \
../../MST/bin/ChargeDistributionModule.o \
../../MST/bin/ChargeScreeningModule.o \
../../MST/bin/PotentialGenerationModule.o \
../../MST/bin/TotalEnergyModule.o \
../../MST/bin/ConvergenceCheckModule.o \
../../MST/bin/MixingModule.o \
../../MST/bin/TestPotentialModule.o \
../../MST/bin/BookKeepingModule.o \
../../MST/bin/ForceModule.o \

#HDF5Module.o

Modules_2 = \
../../MST/bin/ClusterMatrixModule.o \
../../MST/bin/MediumHostModule.o \
../../MST/bin/SROModule.o \
../../MST/bin/CrystalMatrixModule.o \
../../MST/bin/AccelerateCPAModule.o \
../../MST/bin/EmbeddedClusterModule.o \
../../MST/bin/CPAMediumModule.o \
../../MST/bin/BandStructureModule.o \
../../MST/bin/MSSolverModule.o \
../../MST/bin/RelMSSolverModule.o \
../../MST/bin/RelGreenFunctionModule.o \
../../MST/bin/GFMethodModule.o

Modules_Kubo = \
KuboDataModule.o \
CurrentMatrixModule.o \
ConductivityModule.o

Routines = \
../../MST/bin/setupLizNeighbor.o \
../../MST/bin/buildSendTable.o \
../../MST/bin/constructDataOnGrid.o \
../../MST/bin/printDataOnGrid.o \
../../MST/bin/printDataOnLine.o \
../../MST/bin/updatePotential.o \
../../MST/bin/setupMixingScheme.o \
../../MST/bin/setupMixRealArrayList.o \
../../MST/bin/setupMixCmplxArrayList.o \
../../MST/bin/updateMixRealValues.o \
../../MST/bin/updateMixCmplxValues.o \
../../MST/bin/potredg.o \
../../MST/bin/getpotg.o \
../../MST/bin/putpotg.o \
../../MST/bin/keep.o \
../../MST/bin/averageDOS.o \
../../MST/bin/printSystemMovie.o\
../../MST/bin/print_version.o \
../../MST/bin/driverSystemMovie.o \
../../MST/bin/force_openmp.o \
../../MST/bin/print_threads.o \
../../MST/bin/convertGijToRel.o \
../../MST/bin/gaspari_gyorffy_formula.o \
../../MST/bin/printMomentVsCoreSplit.o \
../../MST/bin/polyfermi.o \
../../MST/bin/calNicholsonPoles.o \
../../MST/bin/fetchVisualDomainParameters.o \
../../MST/bin/printDensityOnGrid.o \
../../MST/bin/printPotentialOnGrid.o \
../../MST/bin/printScfResults.o

.SUFFIXES:	.o .f .c .f90 .F90

.f.o:
	$(FC) $(F77FLAGS) $<
	if test "$(No_Remove)" = "1"; then cp -f $< $(ODIR); fi

.f90.o:
	$(FC) $(FFLAGS) $<
	if test "$(No_Remove)" = "1"; then cp -f $< $(ODIR); fi

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
	      $(FC) $(PREPROC_FFLAGS) $(FFLAGS) -I$(TOP_DIR)/../MST/bin/ $< ;\
	fi ;

.c.o:
	$(CC) $(CFLAGS) $<
	if test "$(No_Remove)" = "1"; then cp -f $< $(ODIR); fi

all: $(ODIR)/../../MST/bin/preproc print_version $(ODIR)/$(EXEC_NAME)

print_version:
	cd $(TOP_DIR) && ../MST/src/get_svn.sh ./src/print_version_include.h
	rm -f $(ODIR)/print_version.o

$(ODIR)/preproc: preproc.o
	$(CC) -o $(ODIR)/../../MST/bin/preproc preproc.o

$(ODIR)/$(EXEC_NAME): inc_fftw $(Modules_Kubo) kubo.o $(MSTLIB) $(MPPLIB)
	$(LD) -o $@ $(ODIR)/kubo.o $(Routines) \
$(BasicModules) $(BasicRoutines) $(AppModules) $(Modules_2) $(Modules_Kubo) $(IOLIB) $(MSTLIB) $(SLULIB) $(DSUPERLULIB) $(MPPLIB) $(ADDLIBS)
	rm -f fftw3f.h

inc_fftw:
	if test "$(Use_FFTW)" = "1"; then ln -fs $(FFTW_INC)/fftw3.f fftw3f.h; fi

clear:
	cd $(ODIR); rm -f $(BasicModules) $(BasicRoutines) $(AppModules) $(Routines) preproc.o kubo.o

clean:
	make clear; rm -f $(ODIR)/$(EXEC_NAME) $(ODIR)/preproc

printSystemMovie.o: $(Modules_2)

DefaultParamModule.o: ../../MST/src/DefaultParameters.h
InputModule.o: DefaultParamModule.o

keep.o: git_version.h
