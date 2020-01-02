OBJ = calculateChargesPotential.o interpolatePotential.o getqm_mt.o getvmt.o newexchg.o alpha2_c.o newpot_c.o epcorr.o rs.o getXCName.o libxcInterface.o calculateASAro3.o

all: libPotential.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libPotential.a : $(OBJ)
	$(AR) -r libPotential.a $(OBJ) 
	cp libPotential.a $(TOP_DIR)/lib
