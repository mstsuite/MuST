#*********************************************************************
#    makefile : generate library: mstlib.a
#
#    Book-keeping:      date                       note
#                     10-14-97          
#                     04-23-01          
#*********************************************************************
VPATH = $(TOP_DIR)/lib

MODULE = \
KindParamModule.o \
ErrorHandlerModule.o \
MathParamModule.o \
CmplxFunctionModule.o \
IntegerFactorsModule.o \
LegendreModule.o \
SphericalHarmonicsModule.o \
GauntFactorsModule.o \
TimerModule.o \
DerivativeModule.o \
IntegrationModule.o \
MatrixDeterminantModule.o \
MatrixInverseModule.o \
Matrix3dModule.o \
MatrixModule.o \
InterpolationModule.o \
SortModule.o \
BesselModule.o \
WriteMatrixModule.o \
StringModule.o \
VectorModule.o \
DataServiceCenter.o \
LUdcmpModule.o \
RelativityToolsModule.o \
QuadraticMatrixModule.o

LOBS = \
copyRealArray.o \
erfc.o \
gauleg.o \
gaussq.o \
genGaussianPoints.o \
fit.o \
newder.o \
cnewder.o \
hunt.o \
interp.o \
ord3v.o \
rtflsp.o \
rtsafe.o \
rtsec.o \
vcross.o \
zbrent.o \
zriddr.o \
chebev.o \
ChebyshevModule.o \
beschb.o \
bessjy.o \
sphbes.o \
csphjy.o \
lagpolfit_r.o \
lagpolfit_c.o \
ylag.o \
zeroout.o \
copyCharArray2String.o \
copyString2CharArray.o \
nocaseCompare.o \
isNumber.o \
isRealNumber.o \
isInteger.o \
getTokenPosition.o \
getNumTokens.o \
getToken.o \
trim_string.o \
cylag.o \
zaxpby.o \
zmar1.o \
zrandn.o \
zucpx.o \
zuqmx.o \
zutfx.o \
FlushFile.o \
ran2.o \
Fermi_function.o \
Gaussian_function.o \
ErrorStep_function.o \
findCubicFit.o \
findQuadraticFit.o \
fourn.o \
rlft3.o \
ratintr.o \
GaussianElim.o \
zgetf2_nopivot.o \
zgetrf_nopivot.o \
pzgetf2_nopivot.o \
pzgetrf_nopivot.o \
pz_LogTrace.o \
PrimeFactorsModule.o \
computeProdExpan.o \
computeRealProdExpan.o \
computeCmplxProdExpan.o \
BroydenLib.o \
AndersonLib.o \
newtc.o \
decompose3Index.o \
breakLine.o \
arrayTools.o

ifdef ESSL_WORKAROUND
  LOBS += essl_workaround.o
endif

.SUFFIXES:	.o .f .F .f90 .F90 
#.SUFFIXES:	.o .F .F90 

.f.o:
	$(FC) $(F77FLAGS) $<
	if [ $(No_Remove) != 0 ] ; then cp $(VPATH)/$*.f $(ODIR) ; fi

.F.o:
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
	fi ;

.f90.o:
	$(FC) $(FFLAGS) $<
	if [ $(No_Remove) != 0 ] ; then cp $(VPATH)/$*.f90 $(ODIR) ; fi

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
$(MSTLIB): $(MODULE) $(LOBS)
	$(ARCHV) $(MSTLIB) $(MODULE) $(LOBS)
	ranlib $(MSTLIB)

IntegrationModule.o: IntegrationModule.F90
	$(FC) $(PREPROC_FFLAGS) $(FFLAGS_ALT) $<

clear: $(MODULE) $(LOBS)
	cd $(ODIR); rm -f $(MODULE) $(LOBS)
#=====================================================================
