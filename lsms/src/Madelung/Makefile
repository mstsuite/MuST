OBJ = cal_madelung_matrix.o getkncut.o getrscut.o getstruc.o interf.o interfsmr.o \
      lattice.o madewd.o madewdj.o madsum.o ord3v.o pqintg_c.o calculateMadelungMatrices.o \
      lmfacts.o bessj.o

all: libMadelung.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libMadelung.a : $(OBJ)
	$(AR) -r libMadelung.a $(OBJ) 
	cp libMadelung.a $(TOP_DIR)/lib
