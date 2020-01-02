OBJ = getcor_c.o calculateCoreStates.o corslv_c.o deepst_c.o semcst_c.o \
      inws_c.o outws_c.o inwhnk_c.o invals_c.o richnk_c.o
#      inws_c.o outws_c.o inwhnk_c.o invals_c.o richnk_c.o srcore.o se.o \
#      deepst_c_check.o

all: libCore.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libCore.a : $(OBJ)
	$(AR) -r libCore.a $(OBJ) 
	cp libCore.a $(TOP_DIR)/lib
