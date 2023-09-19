VPATH = $(TOP_DIR)/src

BasicModules = \
$(MST_ODIR)/PhysParamModule.o \
$(MST_ODIR)/ChemElementModule.o \
$(MST_ODIR)/PublicTypeDefinitionsModule.o \
$(MST_ODIR)/PublicParamDefinitionsModule.o \
$(MST_ODIR)/AdaptIntegrationModule.o \
$(MST_ODIR)/MatrixBlockInversionModule.o \
$(MST_ODIR)/BinvMatrix.o \
$(MST_ODIR)/WriteFunctionModule.o \
$(MST_ODIR)/DefaultParamModule.o \
$(MST_ODIR)/CmdLineOptionModule.o \
$(MST_ODIR)/LatticeModule.o \
$(MST_ODIR)/IBZRotationModule.o \
$(MST_ODIR)/CheckPointModule.o

BasicRoutines = \
$(MST_ODIR)/readPositionData.o \
$(MST_ODIR)/readMomentDirectionData.o \
$(MST_ODIR)/readInputInOtherFormat.o \
$(MST_ODIR)/readOldInfoTable.o \
$(MST_ODIR)/rdin_old_lsms.o \
$(MST_ODIR)/rdin_old_infotable.o \
$(MST_ODIR)/rdin_old_infoevec.o \
$(MST_ODIR)/getValueAtPosi.o \
$(MST_ODIR)/lattice.o \
$(MST_ODIR)/averageAcrossProcesses.o \
$(MST_ODIR)/setupSemiCircleContour.o \
$(MST_ODIR)/getFermiDiracFunc.o

AppModules = \
$(MST_ODIR)/InputModule.o \
$(MST_ODIR)/ScfDataModule.o \
$(MST_ODIR)/PotentialTypeModule.o \
$(MST_ODIR)/ExchCorrFunctionalModule.o \
$(MST_ODIR)/SystemModule.o \
$(MST_ODIR)/BZoneModule.o \
$(MST_ODIR)/MatsubaraApproxModule.o \
$(MST_ODIR)/ContourModule.o \
$(MST_ODIR)/ProcMappingModule.o \
$(MST_ODIR)/Atom2ProcModule.o \
$(MST_ODIR)/PolyhedraModule.o \
$(MST_ODIR)/SystemVolumeModule.o \
$(MST_ODIR)/OutputModule.o \
$(MST_ODIR)/VPI_main.o \
$(MST_ODIR)/IsoparametricIntegrationModule.o \
$(MST_ODIR)/NeighborModule.o \
$(MST_ODIR)/SendRecvTmatModule.o \
$(MST_ODIR)/AtomModule.o \
$(MST_ODIR)/StepFunctionModule.o \
$(MST_ODIR)/MadelungModule.o \
$(MST_ODIR)/SystemSymmetryModule.o \
$(MST_ODIR)/LdaCorrectionModule.o \
$(MST_ODIR)/RadialGridModule.o \
$(MST_ODIR)/PotentialModule.o \
$(MST_ODIR)/StrConstModule.o \
$(MST_ODIR)/RSpaceStrConstModule.o \
$(MST_ODIR)/TriangleModule.o \
$(MST_ODIR)/PolygonModule.o \
$(MST_ODIR)/SurfElementsModule.o \
$(MST_ODIR)/SpinRotationModule.o \
$(MST_ODIR)/SSSolverModule.o \
$(MST_ODIR)/DiracSolverModule.o \
$(MST_ODIR)/RelScattererModule.o \
$(MST_ODIR)/RelSSSolverModule.o \
$(MST_ODIR)/SingleScatteringDOSModule.o \
$(MST_ODIR)/SCPolesModule.o \
$(MST_ODIR)/SMatrixPolesModule.o \
$(MST_ODIR)/SineMatrixZerosModule.o \
$(MST_ODIR)/KreinModule.o \
$(MST_ODIR)/CoreStatesModule.o \
$(MST_ODIR)/Uniform3DGridModule.o \
$(MST_ODIR)/ParallelFFTModule.o \
$(MST_ODIR)/ConstrainLocalMomentModule.o \
$(MST_ODIR)/ValenceDensityModule.o \
$(MST_ODIR)/ChargeDensityModule.o \
$(MST_ODIR)/ChargeDistributionModule.o \
$(MST_ODIR)/ChargeScreeningModule.o \
$(MST_ODIR)/PotentialGenerationModule.o \
$(MST_ODIR)/TotalEnergyModule.o \
$(MST_ODIR)/ConvergenceCheckModule.o \
$(MST_ODIR)/MixingModule.o \
$(MST_ODIR)/TestPotentialModule.o \
$(MST_ODIR)/BookKeepingModule.o \
$(MST_ODIR)/ForceModule.o 

#HDF5Module.o

Modules_2 = \
$(MST_ODIR)/ClusterMatrixModule.o \
$(MST_ODIR)/MediumHostModule.o \
$(MST_ODIR)/SROModule.o \
$(MST_ODIR)/CrystalMatrixModule.o \
$(MST_ODIR)/AccelerateCPAModule.o \
$(MST_ODIR)/EmbeddedClusterModule.o \
$(MST_ODIR)/CPAMediumModule.o \
$(MST_ODIR)/BandStructureModule.o \
$(MST_ODIR)/MSSolverModule.o \
$(MST_ODIR)/RelMSSolverModule.o \
$(MST_ODIR)/MultiScatteringDOSModule.o \
$(MST_ODIR)/RelGreenFunctionModule.o \
$(MST_ODIR)/OrbitalBasisModule.o \
$(MST_ODIR)/LocalGFModule.o \
$(MST_ODIR)/GFMethodModule.o

Modules_Kubo = \
KuboDataModule.o \
CurrentMatrixModule.o \
CPAConductivityModule.o \
LSMSConductivityModule.o \
ConductivityModule.o

Routines = \
$(MST_ODIR)/setupLizNeighbor.o \
$(MST_ODIR)/buildSendTable.o \
$(MST_ODIR)/constructDataOnGrid.o \
$(MST_ODIR)/printDataOnGrid.o \
$(MST_ODIR)/printDataOnLine.o \
$(MST_ODIR)/updatePotential.o \
$(MST_ODIR)/setupMixingScheme.o \
$(MST_ODIR)/setupMixRealArrayList.o \
$(MST_ODIR)/setupMixCmplxArrayList.o \
$(MST_ODIR)/updateMixRealValues.o \
$(MST_ODIR)/updateMixCmplxValues.o \
$(MST_ODIR)/potredg.o \
$(MST_ODIR)/getpotg.o \
$(MST_ODIR)/putpotg.o \
$(MST_ODIR)/keep.o \
$(MST_ODIR)/averageDOS.o \
$(MST_ODIR)/printSystemMovie.o\
$(MST_ODIR)/print_version.o \
$(MST_ODIR)/driverSystemMovie.o \
$(MST_ODIR)/force_openmp.o \
$(MST_ODIR)/print_threads.o \
$(MST_ODIR)/convertGijToRel.o \
$(MST_ODIR)/gaspari_gyorffy_formula.o \
$(MST_ODIR)/printMomentVsCoreSplit.o \
$(MST_ODIR)/polyfermi.o \
$(MST_ODIR)/calNicholsonPoles.o \
$(MST_ODIR)/fetchVisualDomainParameters.o \
$(MST_ODIR)/printDensityOnGrid.o \
$(MST_ODIR)/printPotentialOnGrid.o \
$(MST_ODIR)/printScfResults.o \
$(MST_ODIR)/adjustEnergy.o

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
	   cd $(ODIR); $(FC) $(FFLAGS) -I$(MST_ODIR) -I$(TOP_DIR)/../MST/src $*.f90; \
	   (if [ $(No_Remove) != 1 ] ; then rm $*.f90; fi); \
	else \
	      $(FC) $(PREPROC_FFLAGS) $(FFLAGS) -I$(MST_ODIR) -I$(TOP_DIR)/../MST/src $< ;\
	fi ;

.c.o:
	$(CC) $(CFLAGS) $<
	if test "$(No_Remove)" = "1"; then cp -f $< $(ODIR); fi

all: $(ODIR)/preproc print_version $(ODIR)/$(EXEC_NAME)

print_version:
	cd $(TOP_DIR) && ../MST/src/get_svn.sh ./src/print_version_include.h
	rm -f $(ODIR)/print_version.o

$(ODIR)/preproc: preproc.o
	$(CC) -o $(ODIR)/preproc preproc.o

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
