VPATH = $(ODIR)/../src

BasicModules = \
PhysParamModule.o \
ChemElementModule.o \
PublicTypeDefinitionsModule.o \
PublicParamDefinitionsModule.o \
AdaptIntegrationModule.o \
MatrixBlockInversionModule.o \
BinvMatrix.o \
WriteFunctionModule.o \
DefaultParamModule.o

BasicRoutines = \
readPositionData.o \
readMomentDirectionData.o \
readInputInOtherFormat.o \
readOldInfoTable.o \
rdin_old_lsms.o \
rdin_old_infotable.o \
rdin_old_infoevec.o \
getValueAtPosi.o \
lattice.o

AppModules = \
InputModule.o \
ScfDataModule.o \
PotentialTypeModule.o \
SystemModule.o \
BZoneModule.o \
MatsubaraModule.o \
ContourModule.o \
ProcMappingModule.o \
Atom2ProcModule.o \
PolyhedraModule.o \
SystemVolumeModule.o \
OutputModule.o \
VPI_main.o \
IsoparametricIntegrationModule.o \
NeighborModule.o \
SendRecvTmatModule.o \
AtomModule.o \
StepFunctionModule.o \
MadelungModule.o \
SystemSymmetryModule.o \
LdaCorrectionModule.o \
RadialGridModule.o \
PotentialModule.o \
CoreStatesModule.o \
StrConstModule.o \
RSpaceStrConstModule.o \
TriangleModule.o \
PolygonModule.o \
SurfElementsModule.o \
SpinRotationModule.o \
SSSolverModule.o \
SCPolesModule.o \
Uniform3DGridModule.o \
ParallelFFTModule.o \
KreinModule.o \
DiracSolverModule.o \
RelScattererModule.o \
ConstrainLocalMomentModule.o \
ExchCorrFunctionalModule.o \
ValenceDensityModule.o \
ChargeDensityModule.o \
ChargeDistributionModule.o \
PotentialGenerationModule.o \
TotalEnergyModule.o \
ConvergenceCheckModule.o \
MixingModule.o \
TestPotentialModule.o \
BookKeepingModule.o \
ForceModule.o \
RelSSSolverModule.o \

#HDF5Module.o

Modules_2 = \
LatticeModule.o \
IBZRotationModule.o \
ClusterMatrixModule.o \
MediumHostModule.o \
SROModule.o \
CrystalMatrixModule.o \
AccelerateCPAModule.o \
EmbeddedClusterModule.o \
CPAMediumModule.o \
BandStructureModule.o \
MSSolverModule.o \
RelMSSolverModule.o \
RelGreenFunctionModule.o \
GFMethodModule.o

Routines = \
setupLizNeighbor.o \
buildSendTable.o \
constructDataOnGrid.o \
printDataOnGrid.o \
printDataOnLine.o \
updatePotential.o \
setupMixingScheme.o \
setupMixRealArrayList.o \
setupMixCmplxArrayList.o \
updateMixRealValues.o \
updateMixCmplxValues.o \
potredg.o \
getpotg.o \
putpotg.o \
keep.o \
averageDOS.o \
printSystemMovie.o\
print_version.o \
driverSystemMovie.o \
force_openmp.o \
print_threads.o \
convertGijToRel.o \
gaspari_gyorffy_formula.o \
printMomentVsCoreSplit.o \
polyfermi.o \
calNicholsonPoles.o \
fetchVisualDomainParameters.o \
printDensityOnGrid.o \
printPotentialOnGrid.o \
printScfResults.o

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
	      $(FC) $(PREPROC_FFLAGS) $(FFLAGS) $< ;\
	fi ;

.c.o:
	$(CC) $(CFLAGS) $<
	if test "$(No_Remove)" = "1"; then cp -f $< $(ODIR); fi

all: $(ODIR)/preproc print_version $(ODIR)/$(EXEC_NAME)

print_version:
	cd .. && ./src/get_svn.sh ./src/print_version_include.h
	rm -f $(ODIR)/print_version.o

$(ODIR)/preproc: preproc.o
	$(CC) -o $(ODIR)/preproc preproc.o

$(ODIR)/$(EXEC_NAME): inc_fftw $(BasicModules) $(BasicRoutines) $(AppModules) $(Modules_2) $(Routines) mst2.o $(MSTLIB) $(MPPLIB)
	$(LD) -o $@ $(ODIR)/mst2.o $(Routines) \
$(BasicModules) $(BasicRoutines) $(AppModules) $(Modules_2) $(IOLIB) $(MSTLIB) $(SLULIB) $(DSUPERLULIB) $(MPPLIB) $(ADDLIBS)
	rm -f fftw3f.h

inc_fftw:
	if test "$(Use_FFTW)" = "1"; then ln -fs $(FFTW_INC)/fftw3.f fftw3f.h; fi

clear:
	cd $(ODIR); rm -f $(BasicModules) $(BasicRoutines) $(AppModules) $(Routines) preproc.o mst2.o

clean:
	make clear; rm -f $(ODIR)/$(EXEC_NAME) $(ODIR)/preproc

printSystemMovie.o: $(Modules_2)

DefaultParamModule.o: DefaultParameters.h
InputModule.o: DefaultParamModule.o

keep.o: git_version.h
