all: lsms wl-lsms rewl-lsms

ifdef USE_GPTL
  gptl_def=-DUSE_GPTL
else
  gptl_def=
endif

LSMS_OBJ = SystemParameters.o read_input.o PotentialIO.o buildLIZandCommLists.o \
           energyContourIntegration.o solveSingleScatterers.o \
           calculateDensities.o calculateChemPot.o checkConsistency.o \
           lsmsClass.o calculateEvec.o initializeAtom.o mixing.o \
           ReplicaExchangeWL.o AlloyBankIO.o rotateToGlobal.o \
           write_restart.o

clean:
	rm -f *.o *.a $(ODIR)/lsms \
                      $(ODIR)/lsms_inter \
                      $(ODIR)/wl-lsms \
                      $(ODIR)/rewl-lsms \
		      zblock_lu_driver

zblock_lu_driver: zblock_lu_cuda_driver.cpp libLSMS.a
	$(CXX) -fopenmp -UBUILDKKRMATRIX_GPU $(INC_PATH) -DNO_HDF5 -o zblock_lu_driver zblock_lu_cuda_driver.cpp libLSMS.a $(LIB_DIR) -L../Accelerator -lAccelerator $(ADD_LIBS)

lsms: lsms.cpp libLSMS.a
	$(CXX) $(INC_PATH) $(gptl_def) -o lsms lsms.cpp libLSMS.a $(LIBS) $(ADD_LIBS)
	mv lsms $(ODIR)

lsms_inter: lsms_inter.cpp
	$(CXX) $(INC_PATH) -o lsms_inter lsms_inter.cpp $(LIBS)
	mv lsms_inter $(ODIR)

wl-lsms: wl_lsms.cpp libLSMS.a
	$(CXX) $(INC_PATH) -o wl-lsms wl_lsms.cpp libLSMS.a $(LIBS) $(ADD_LIBS)
	mv wl-lsms $(ODIR)

rewl-lsms: rewl_lsms.cpp libLSMS.a
	$(CXX) $(INC_PATH) -o rewl-lsms rewl_lsms.cpp libLSMS.a $(LIBS) $(ADD_LIBS)
	mv rewl-lsms $(ODIR)

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

libLSMS.a: $(LSMS_OBJ)
	$(AR) -r libLSMS.a $(LSMS_OBJ) 

