OBJ = calculateTotalEnergy.o janake_c.o zeropt_c.o

all: libTotalEnergy.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libTotalEnergy.a : $(OBJ)
	$(AR) -r libTotalEnergy.a $(OBJ) 
	cp libTotalEnergy.a $(TOP_DIR)/lib
