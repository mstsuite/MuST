VPATH = ../driver

.SUFFIXES:	.o .f .f90 .F90

.f.o:
	$(FC) $(F77FLAGS) $<

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

BASE_ROUTINES = readInputs.o startProcess.o finishProcess.o setupRadGridAndCell.o

SRC_OBJS =  \
$(ODIR)/PublicTypeDefinitionsModule.o \
$(ODIR)/PublicParamDefinitionsModule.o \
$(ODIR)/AtomModule.o \
$(ODIR)/Atom2ProcModule.o \
$(ODIR)/buildSendTable.o \
$(ODIR)/BZoneModule.o \
$(ODIR)/IBZRotationModule.o \
$(ODIR)/ChemElementModule.o \
$(ODIR)/MatsubaraModule.o \
$(ODIR)/ContourModule.o \
$(ODIR)/CoreStatesModule.o \
$(ODIR)/GauntFactorsModule.o \
$(ODIR)/getpotg.o \
$(ODIR)/DefaultParamModule.o  \
$(ODIR)/InputModule.o  \
$(ODIR)/IsoparametricIntegrationModule.o \
$(ODIR)/lattice.o  \
$(ODIR)/LdaCorrectionModule.o \
$(ODIR)/MadelungModule.o  \
$(ODIR)/LatticeModule.o  \
$(ODIR)/NeighborModule.o  \
$(ODIR)/OutputModule.o \
$(ODIR)/WriteFunctionModule.o \
$(ODIR)/putpotg.o \
$(ODIR)/potredg.o \
$(ODIR)/PolygonModule.o \
$(ODIR)/PotentialModule.o \
$(ODIR)/ProcMappingModule.o  \
$(ODIR)/PolyhedraModule.o  \
$(ODIR)/PotentialTypeModule.o  \
$(ODIR)/RadialGridModule.o  \
$(ODIR)/rdin_old_lsms.o  \
$(ODIR)/rdin_old_infoevec.o  \
$(ODIR)/rdin_old_infotable.o \
$(ODIR)/readInputInOtherFormat.o  \
$(ODIR)/readPositionData.o  \
$(ODIR)/readOldInfoTable.o  \
$(ODIR)/readMomentDirectionData.o  \
$(ODIR)/convertGijToRel.o \
$(ODIR)/ScfDataModule.o  \
$(ODIR)/setupLizNeighbor.o \
$(ODIR)/SpinRotationModule.o \
$(ODIR)/StepFunctionModule.o  \
$(ODIR)/SSSolverModule.o \
$(ODIR)/RelSSSolverModule.o \
$(ODIR)/StrConstModule.o  \
$(ODIR)/CrystalMatrixModule.o \
$(ODIR)/SurfElementsModule.o \
$(ODIR)/SystemModule.o  \
$(ODIR)/SystemSymmetryModule.o \
$(ODIR)/SystemVolumeModule.o \
$(ODIR)/TestPotentialModule.o \
$(ODIR)/TriangleModule.o \
$(ODIR)/calNicholsonPoles.o \
$(ODIR)/polyfermi.o \
$(ODIR)/VPI_main.o

TEST_CASES = $(ODIR)/testStrConst \
             $(ODIR)/tst_SurfElements \
             $(ODIR)/tst_ylm \
             $(ODIR)/tst_bessel \
             $(ODIR)/testIntegration \
             $(ODIR)/tst_string \
             $(ODIR)/testNeighbor \
             $(ODIR)/tst_write \
             $(ODIR)/tst_read \
             $(ODIR)/tst_findFit \
             $(ODIR)/testMadelung \
             $(ODIR)/madpot \
             $(ODIR)/ft \
             $(ODIR)/ft_s \
             $(ODIR)/ft_p \
             $(ODIR)/pav \
             $(ODIR)/testParallelFFT \
             $(ODIR)/testSineMatrixPole \
             $(ODIR)/testSineMatrixPole2 \
             $(ODIR)/findResonance \
             $(ODIR)/testPZGETRI \
             $(ODIR)/testStepFunc \
             $(ODIR)/testStepFunc2 \
             $(ODIR)/testProcMapping \
             $(ODIR)/testMPP \
             $(ODIR)/testPolyhedra \
             $(ODIR)/testGroupComm \
             $(ODIR)/testBZone \
             $(ODIR)/testQuadraticMatrix \
             $(ODIR)/testParallelIO \
             $(ODIR)/testSSSolver \
             $(ODIR)/testBlockInverse \
             $(ODIR)/testIBZRotation \
             $(ODIR)/testFEGreen \
             $(ODIR)/testWronskian \
             $(ODIR)/testFFTGrid \
             $(ODIR)/ssk \
             $(ODIR)/assk \
             $(ODIR)/testCoreStates \
             $(ODIR)/findBoundStates \
             $(ODIR)/testPolyFermi \
             $(ODIR)/testBandStructure \
             $(ODIR)/testExchCorr

TEST_CASES_O = $(ODIR)/testStrConst.o \
               $(ODIR)/tst_SurfElements.o \
               $(ODIR)/tst_ylm.o \
               $(ODIR)/tst_bessel.o \
               $(ODIR)/testIntegration.o \
               $(ODIR)/tst_string.o \
               $(ODIR)/testNeighbor.o \
               $(ODIR)/tst_write.o \
               $(ODIR)/tst_read.o \
               $(ODIR)/tst_findFit.o \
               $(ODIR)/testMadelung.o \
               $(ODIR)/madpot.o \
               $(ODIR)/pav.o \
               $(ODIR)/testSineMatrixPole.o \
               $(ODIR)/testSineMatrixPole2.o \
               $(ODIR)/findResonance.o \
               $(ODIR)/testPZGETRI.o \
               $(ODIR)/testStepFunc.o \
               $(ODIR)/testStepFunc2.o \
               $(ODIR)/testProcMapping.o \
               $(ODIR)/testMPP.o \
               $(ODIR)/testPolyhedra.o \
               $(ODIR)/testGroupComm.o \
               $(ODIR)/testBZone.o \
               $(ODIR)/testQuadraticMatrix.o \
               $(ODIR)/testParallelIO.o \
               $(ODIR)/testSSSolver.o \
               $(ODIR)/testBlockInverse.o \
               $(ODIR)/testFEGreen.o
#
# Unavailable or missing
#             $(ODIR)/tst_ScreenTau
#             $(ODIR)/findFactors \
#             $(ODIR)/testSetupProcDim \

all: $(BASE_ROUTINES) $(TEST_CASES)

clear:
	rm -f $(TEST_CASES_O)

clean:
	rm -f $(TEST_CASES); make "ODIR=$(ODIR)" clear -f $(VPATH)/Makefile

$(ODIR)/tst_findFit: tst_findFit.o
	$(LD) -o $(ODIR)/tst_findFit tst_findFit.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/testStrConst: testStrConst.o $(SRC_OBJS)
	$(LD) -o $(ODIR)/testStrConst testStrConst.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testBZone: testBZone.o
	$(LD) -o $(ODIR)/testBZone testBZone.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testQuadraticMatrix: testQuadraticMatrix.o
	$(LD) -o $(ODIR)/testQuadraticMatrix $(ODIR)/testQuadraticMatrix.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_SurfElements: tst_SurfElements.o
	$(LD) -o $(ODIR)/tst_SurfElements tst_SurfElements.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_ylm: tst_ylm.o
	$(LD) -o $(ODIR)/tst_ylm tst_ylm.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_bessel: tst_bessel.o
	$(LD) -o $(ODIR)/tst_bessel tst_bessel.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/testIntegration: testIntegration.o
	$(LD) -o $(ODIR)/testIntegration testIntegration.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_string: tst_string.o $(ODIR)/StringModule.o
	$(LD) -o $(ODIR)/tst_string tst_string.o $(ODIR)/StringModule.o $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_write: tst_write.o
	$(LD) -o $(ODIR)/tst_write tst_write.o $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/tst_read: tst_read.o
	$(LD) -o $(ODIR)/tst_read tst_read.o $(IOLIB) $(MSTLIB)  $(ADDLIBS)

$(ODIR)/testNeighbor: testNeighbor.o
	$(LD) -o $(ODIR)/testNeighbor testNeighbor.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB)  $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testMadelung: testMadelung.o $(ODIR)/MadelungModule.o
	$(LD) -o $(ODIR)/testMadelung testMadelung.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/madpot: madpot.o
	$(LD) -o $(ODIR)/madpot madpot.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/pav: pot_average.o
	$(LD) -o $(ODIR)/pav pot_average.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)


$(ODIR)/ft: inc_fftw fourier_transform.o
	$(LD) -o $(ODIR)/ft fourier_transform.o $(MPPLIB) $(MSTLIB) $(ADDLIBS)
	rm -f fftw3f.h

$(ODIR)/ft_s: inc_fftw ft_serial.o
	$(LD) -o $(ODIR)/ft_s ft_serial.o $(MSTLIB) $(ADDLIBS)
	rm -f fftw3f.h

$(ODIR)/ft_m: inc_fftw ft_mpi-fftw.o
	$(LD) -o $(ODIR)/ft_m ft_mpi-fftw.o $(MPPLIB) $(MSTLIB) $(ADDLIBS)
	rm -f fftw3f.h

$(ODIR)/ft_p: inc_fftw ft_p3dfft.o
	$(LD) -o $(ODIR)/ft_p ft_p3dfft.o $(MPPLIB) $(MSTLIB) $(ADDLIBS)
	rm -f fftw3f.h

inc_fftw:
	if test "$(Use_FFTW)" == "1"; then ln -fs $(FFTW_INC)/fftw3.f fftw3f.h; fi

$(ODIR)/testParallelFFT: inc_fftw testParallelFFT.o $(ODIR)/Uniform3DGridModule.o $(ODIR)/ParallelFFTModule.o
	$(LD) -o $(ODIR)/testParallelFFT testParallelFFT.o $(ODIR)/Uniform3DGridModule.o $(ODIR)/ParallelFFTModule.o \
$(SRC_OBJS)  $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)
	rm -f fftw3f.h

$(ODIR)/testSineMatrixPole: testSineMatrixPole.o calScatteringPoles.o
	$(LD) -o $(ODIR)/testSineMatrixPole testSineMatrixPole.o $(ODIR)/calScatteringPoles.o $(ODIR)/readInputs.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testSineMatrixPole2: testSineMatrixPole2.o calScatteringPoles.o
	$(LD) -o $(ODIR)/testSineMatrixPole2 testSineMatrixPole2.o $(ODIR)/readInputs.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/findResonance:  findResonance.o
	$(LD) -o $(ODIR)/findResonance findResonance.o $(ODIR)/SCPolesModule.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testSSSolver: $(ODIR)/SCPolesModule.o testSSSolver.o
	$(LD) -o $(ODIR)/testSSSolver $(ODIR)/testSSSolver.o $(ODIR)/SCPolesModule.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testWronskian: testWronskian.o
	$(LD) -o $(ODIR)/testWronskian $(ODIR)/testWronskian.o $(ODIR)/SCPolesModule.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testFFTGrid: testFFTGrid.o $(ODIR)/Uniform3DGridModule.o
	$(LD) -o $(ODIR)/testFFTGrid $(ODIR)/testFFTGrid.o Uniform3DGridModule.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testIBZRotation: testIBZRotation.o
	$(LD) -o $(ODIR)/testIBZRotation $(ODIR)/testIBZRotation.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/ssk: SingleSiteKrein.o 
	$(LD) -o $(ODIR)/ssk $(ODIR)/SingleSiteKrein.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/assk: AdaptIntegrationModule.o adaptiveSSK.o 
	$(LD) -o $(ODIR)/assk $(ODIR)/adaptiveSSK.o $(ODIR)/AdaptIntegrationModule.o \
$(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testBlockInverse: testBlockInverse.o
	$(LD) -o $(ODIR)/testBlockInverse $(ODIR)/testBlockInverse.o $(ODIR)/MatrixBlockInversionModule.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testParallelIO: testParallelIO.o
	$(LD) -o $(ODIR)/testParallelIO testParallelIO.o $(ODIR)/readInputs.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testPZGETRI: testPZGETRI.o 
	$(LD) -o $(ODIR)/testPZGETRI testPZGETRI.o $(ADDLIBS) $(MPPLIB) $(MSTLIB)

$(ODIR)/testProcMapping: testProcMapping.o $(ODIR)/ProcMappingModule.o
	$(LD) -o $(ODIR)/testProcMapping testProcMapping.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testStepFunc: testStepFunc.o $(ODIR)/StepFunctionModule.o
	$(LD) -o $(ODIR)/testStepFunc $(ODIR)/testStepFunc.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testStepFunc2: testStepFunc2.o
	$(LD) -o $(ODIR)/testStepFunc2 $(ODIR)/testStepFunc2.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testMPP: testMPP.o
	$(LD) -o $(ODIR)/testMPP $(ODIR)/testMPP.o $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testPolyhedra: tst_polyhedra.o
	$(LD) -o $(ODIR)/testPolyhedra tst_polyhedra.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testGroupComm: testGroupComm.o
	$(LD) -o $(ODIR)/testGroupComm testGroupComm.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

#$(ODIR)/tst_ScreenTau: tst_ScreenTau.o $(BASE_ROUTINES)
#	$(LD) -o $(ODIR)/tst_ScreenTau tst_ScreenTau.o $(BASE_ROUTINES) $(ODIR)/TauScreenKKRModule.o $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

#$(ODIR)/findFactors: findFactors.o
#	$(LD) -o $(ODIR)/findFactors findFactors.o $(MSTLIB)

#$(ODIR)/testSetupProcDim: testSetupProcDim.o $(BASE_ROUTINES)
#	$(LD) -o $(ODIR)/testSetupProcDim testSetupProcDim.o setupProcDim.o $(MSTLIB)

$(ODIR)/testFEGreen: testFEGreen.o $(SRC_OBJS)
	$(LD) -o $(ODIR)/testFEGreen testFEGreen.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testCoreStates: testCoreStates.o $(SRC_OBJS)
	$(LD) -o $(ODIR)/testCoreStates testCoreStates.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/findBoundStates: findBoundStates.o $(SRC_OBJS)
	$(LD) -o $(ODIR)/findBoundStates findBoundStates.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testPolyFermi: testPolyFermi.o $(SRC_OBJS)
	$(LD) -o $(ODIR)/testPolyFermi testPolyFermi.o $(BASE_ROUTINES) $(SRC_OBJS) $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testBandStructure: testBandStructure.o $(SRC_OBJS) $(ODIR)/BandStructureModule.o
	$(LD) -o $(ODIR)/testBandStructure testBandStructure.o $(ODIR)/BandStructureModule.o $(BASE_ROUTINES) $(SRC_OBJS) \
$(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)

$(ODIR)/testExchCorr: testExchCorr.o $(SRC_OBJS) $(ODIR)/ExchCorrFunctionalModule.o
	$(LD) -o $(ODIR)/testExchCorr testExchCorr.o $(BASE_ROUTINES) $(SRC_OBJS) $(ODIR)/ExchCorrFunctionalModule.o $(IOLIB) $(MPPLIB) $(MSTLIB) $(ADDLIBS)
