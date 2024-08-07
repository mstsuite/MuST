target_sources(mst_main
               PUBLIC
               mst2.F90
               )

target_sources(mst_lib
              PUBLIC
              AccelerateCPAModule.F90
              AdaptIntegrationModule.F90
              AccelerateCPAModule.F90
              AnaliticContinuationModule.F90
              Atom2ProcModule.F90
              AtomModule.F90
              BZoneModule.F90
              BandStructureModule.F90
              BinvMatrix.F90
              BookKeepingModule.F90
              CheckPointModule.F90
              CPAMediumModule.F90
              ChargeDensityModule.F90
              ChargeDistributionModule.F90
              ChargeScreeningModule.F90
              ChemElementModule.F90
              ClusterMatrixModule.F90
              ConstrainLocalMomentModule.F90
              ContourModule.F90
              ConvergenceCheckModule.F90
              CoreStatesModule.F90
              CrystalMatrixModule.F90
              DefaultParamModule.F90
              DefaultParameters.h
              DiracSolverModule.F90
              EmbeddedClusterModule.F90
              ExchCorrFunctionalModule.F90
              ForceModule.F90
              GFMethodModule.F90
              HDF5Module.F90
              IBZRotationModule.F90
              InputModule.F90
              IsoparametricIntegrationModule.F90
              KreinModule.F90
              LatticeModule.F90
              LdaCorrectionModule.F90
              MSSolverModule.F90
              MadelungModule.F90
              MatrixBlockInversionModule.F90
              MatsubaraModule.F90
              MediumHostModule.F90
              #MediumMatrixModule.F90
              MixingModule.F90
              MultiScatteringDOSModule.F90
              NeighborModule.F90
              OutputModule.F90
              ParallelFFTModule.F90
              PhysParamModule.F90
              PolygonModule.F90
              PolyhedraModule.F90
              PotentialGenerationModule.F90
              PotentialModule.F90
              PotentialTypeModule.F90
              ProcMappingModule.F90
              PublicParamDefinitionsModule.F90
              PublicTypeDefinitionsModule.F90
              RSpaceStrConstModule.F90
              RadialGridModule.F90
              RelGreenFunctionModule.F90
              RelMSSolverModule.F90
              RelSSSolverModule.F90
              RelScattererModule.F90
              SCPolesModule.F90
              SMatrixPolesModule.F90
              SROModule.F90
              SSSolverModule.F90
              ScfDataModule.F90
              #ScreenStrConstModule.F90
              SendRecvTmatModule.F90
              SineMatrixZerosModule.F90
              SingleScatteringDOSModule.F90
              SpinRotationModule.F90
              StepFunctionModule.F90
              StrConstModule.F90
              SurfElementsModule.F90
              SystemModule.F90
              SystemSymmetryModule.F90
              SystemVolumeModule.F90
              TestPotentialModule.F90
              TotalEnergyModule.F90
              TriangleModule.F90
              Uniform3DGridModule.F90
              VPI_main.f
              ValenceDensityModule.F90
              WriteFunctionModule.F90
              adjustEnergy.F90
              averageAcrossProcesses.F90
              averageDOS.F90
              buildSendTable.F90
              #calBandStructure.F90
              calNicholsonPoles.F90
              computeTraceComplex.F90
              constructDataOnGrid.F90
              convertGijToRel.F90
              driverSystemMovie.F90
              fetchVisualDomainParameters.F90
              force_openmp.F90
              gaspari_gyorffy_formula.F90
              getFermiDiracFunc.F90
              getValueAtPosi.F90
              getpotg.F90
              imp_inp.h
              keep.F90
              lattice.F90
              polyfermi.F90
              potredg.F90
              #preproc.c
              printDataOnGrid.F90
              printDataOnLine.F90
              printDensityOnGrid.F90
              printMomentVsCoreSplit.F90
              printPotentialOnGrid.F90
              printScfResults.F90
              printSystemMovie.F90
              print_threads.F90
              print_version.F90
              putpotg.F90
              rdin_old_infoevec.F90
              rdin_old_infotable.F90
              rdin_old_lsms.F90
              readInputInOtherFormat.F90
              readMomentDirectionData.F90
              readOldInfoTable.F90
              readPositionData.F90
              setupLizNeighbor.F90
              setupMixCmplxArrayList.F90
              setupMixRealArrayList.F90
              setupMixingScheme.F90
              setupSemiCircleContour.F90
              updateMixCmplxValues.F90
              updateMixRealValues.F90
              updatePotential.F90
              )

target_include_directories(mst_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

execute_process(COMMAND git describe --tags
        OUTPUT_VARIABLE GIT_VERSION_TAG
        ERROR_QUIET)

# Check whether we got any revision (which isn't
# always the case, e.g. when someone downloaded a zip
# file from Github instead of a checkout)
if ("${GIT_VERSION_TAG}" STREQUAL "")
    set(GIT_VERSION_TAG "N/A")
else()
    string(STRIP "${GIT_VERSION_TAG}" GIT_VERSION_TAG)
endif()

string(TIMESTAMP TIME "%A, %B %d %Y %H:%M:%S")

cmake_host_system_information(RESULT FQDN QUERY FQDN)

execute_process(COMMAND whoami
        OUTPUT_VARIABLE USERNAME
        ERROR_QUIET)

if ("${USERNAME}" STREQUAL "")
    set(USERNAME "N/A")
else()
    string(STRIP "${USERNAME}" USERNAME)
endif()

configure_file (
        "${CMAKE_CURRENT_SOURCE_DIR}/git_version.h.in"
        "${CMAKE_CURRENT_BINARY_DIR}/git_version.h"
)

configure_file (
        "${CMAKE_CURRENT_SOURCE_DIR}/imp_inp.h"
        "${CMAKE_CURRENT_BINARY_DIR}/imp_inp.h"
        COPYONLY
)

configure_file (
        "${CMAKE_CURRENT_SOURCE_DIR}/print_version_include.h.in"
        "${CMAKE_CURRENT_BINARY_DIR}/print_version_include.h"
)

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/DefaultParameters.h"
        "${CMAKE_CURRENT_BINARY_DIR}/DefaultParameters.h"
        COPYONLY
)
