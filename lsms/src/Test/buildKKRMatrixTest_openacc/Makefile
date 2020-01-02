export XLF_LIB = /sw/summitdev/xl/20170606-beta/xlf/16.1.0/lib

export TOP_DIR = $(shell pwd)/../../..
export INC_PATH =
export LIBS := -L$(OLCF_CUDA_ROOT)/lib64 -lgfortran -L$(OLCF_ESSL_ROOT)/lib64 -lessl -L$(XLF_LIB) -lxlf90_r -lxlfmath -lxl -lxlopt -Wl,-rpath=$(XLF_LIB)

#include $(TOP_DIR)/architecture.h
#include arch-minsky-gnu
include arch-minsky-pgi

export INC_PATH += -I$(TOP_DIR)/include -I$(TOP_DIR)/src -I${TOP_DIR}/include
export MISC = $(TOP_DIR)/src/Misc

COEFICIENTS = cgaunt_c.o ifacts_c.o quadrature.o  fstop.o zeroout.o \
              gaunt.o plglmax.o associatedLegendreFunction.o

all: buildKKRMatrixTest-acc

clean:
	rm -f *.o buildKKRMatrixTest

cgaunt_c.o: $(MISC)/cgaunt_c.f
	$(F77) -c $(MISC)/cgaunt_c.f
ifacts_c.o: $(MISC)/ifacts_c.f
	$(F77) -c $(MISC)/ifacts_c.f
#gauleg.o: $(MISC)/gauleg.f
#	$(F77) -c $(MISC)/gauleg.f
quadrature.o: quadrature.cpp
	$(CXX) -c quadrature.cpp
plglmax.o: $(MISC)/plglmax.f
	$(F77) -c $(MISC)/plglmax.f
associatedLegendreFunction.o: associatedLegendreFunction.cpp associatedLegendreFunction.hpp
	$(CXX) -c associatedLegendreFunction.cpp
fstop.o: $(MISC)/fstop.f
	$(F77) -c $(MISC)/fstop.f
zeroout.o: $(MISC)/zeroout.f
	$(F77) -c $(MISC)/zeroout.f
gaunt.o: $(MISC)/gaunt.f
	$(F77) -c $(MISC)/gaunt.f

makegij_c.o: makegij_c.f
	$(F77) -c makegij_c.f
setgij.o: setgij.f
	$(F77) -c setgij.f

makegij_new.o: makegij_new.cpp
	$(CXX) $(INC_PATH) -c makegij_new.cpp
plglmax_new.o: plglmax_new.cpp
	$(CXX) $(INC_PATH) -c plglmax_new.cpp

quadrature.o: quadrature.cpp
	$(CXX) $(INC_PATH) -c quadrature.cpp

associatedLegendreFunction.o: associatedLegendreFunction.cpp associatedLegendreFunction.hpp
	$(CXX) -c associatedLegendreFunction.cpp

buildKKRMatrix_original.o: buildKKRMatrix_original.cpp
	$(CXX) $(INC_PATH) -c buildKKRMatrix_original.cpp
buildKKRMatrix_nrel_ns2.o: buildKKRMatrix_nrel_ns2.cpp
	$(CXX) $(INC_PATH) -c buildKKRMatrix_nrel_ns2.cpp
kkrmat_acc-mix.o: kkrmat_acc-mix.cpp
	#$(CXX) $(INC_PATH) -ta=tesla:cuda8.0,managed,fastmath,flushz -Minfo=accel -c kkrmat_acc.cpp
	$(CXX) $(INC_PATH) -c kkrmat_acc-mix.cpp

c_kkrmat.o: c_kkrmat.c
	#$(CXX) $(INC_PATH) -ta=tesla:cuda8.0,managed,fastmath,flushz -Minfo=accel -c c_kkrmat.cpp
	#$(CC) $(INC_PATH) -c99 -Minfo=accel -ta=tesla:cuda8.0,fastmath,flushz -c c_kkrmat.c
	$(CC) $(INC_PATH) -c99 -Minfo=accel -ta=tesla:cc60,cuda8.0,fastmath,flushz,safecache -c c_kkrmat.c

buildKKRMatrixTest-acc: buildKKRMatrixTest.cpp $(COEFICIENTS) \
                    buildKKRMatrix_original.o makegij_c.o setgij.o \
										buildKKRMatrix_nrel_ns2.o \
                    kkrmat_acc-mix.o c_kkrmat.o makegij_new.o plglmax_new.o 
	$(CXX) $(INC_PATH) -mp -ta=tesla:cc60,cuda8.0,flushz,fastmath,safecache -Minfo=accel \
	  -o buildKKRMatrixTest-acc buildKKRMatrixTest.cpp \
    buildKKRMatrix_original.o buildKKRMatrix_nrel_ns2.o kkrmat_acc-mix.o c_kkrmat.o \
    makegij_c.o setgij.o \
    $(COEFICIENTS) $(ADD_LIBS) $(LIBS)
	#$(CXX) $(INC_PATH) \
	#$(CXX) $(INC_PATH) -ta=tesla:cuda8.0,fastmath,flushz -Minfo=accel \

