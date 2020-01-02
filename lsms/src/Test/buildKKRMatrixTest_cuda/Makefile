
export TOP_DIR = $(shell pwd)/../../..
export INC_PATH =
export LIBS :=

include $(TOP_DIR)/architecture.h

export INC_PATH += -I $(TOP_DIR)/include -I $(TOP_DIR)/src
export MISC = $(TOP_DIR)/src/Misc

COEFICIENTS = cgaunt_c.o ifacts_c.o quadrature.o  fstop.o zeroout.o \
              gaunt.o plglmax.o associatedLegendreFunction.o

all: buildKKRMatrixTest

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

buildKKRMatrix_original.o: buildKKRMatrix_original.cpp
	$(CXX) $(INC_PATH) -c buildKKRMatrix_original.cpp
buildKKRMatrix_nrel_ns2.o: buildKKRMatrix_nrel_ns2.cpp
	$(CXX) $(INC_PATH) -c buildKKRMatrix_nrel_ns2.cpp

buildKKRMatrixTest: buildKKRMatrixTest.cpp $(COEFICIENTS) \
                    buildKKRMatrix_original.o makegij_c.o setgij.o \
                    buildKKRMatrix_nrel_ns2.o makegij_new.o plglmax_new.o 
	$(CXX) $(INC_PATH) -o buildKKRMatrixTest buildKKRMatrixTest.cpp \
               buildKKRMatrix_original.o buildKKRMatrix_nrel_ns2.o \
               makegij_c.o setgij.o \
               $(COEFICIENTS) $(ADD_LIBS) $(LIBS)

