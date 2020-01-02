OBJ = setup_vorpol.o setup_boundary.o filter_edge.o chkedge.o invm3.o chkbnd.o celbnd.o \
      polyhedron.o sort.o sortidx.o rcritpts.o volvor.o stepyll.o caltnode.o intpl0.o \
      calsig.o intphi.o inter_dip.o inter.o interstitial.o sigma.o \
      setupVorpol.o

all: libVORPOL.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libVORPOL.a : $(OBJ)
	$(AR) -r libVORPOL.a $(OBJ) 
	cp libVORPOL.a $(TOP_DIR)/lib
