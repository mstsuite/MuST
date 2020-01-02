all: testPotentialShift

clean:
	rm -f *.o *.a testPotentialShift $(TOP_DIR)/bin/testPotentialShift

testPotentialShift: testPotentialShift.cpp
	$(CXX) $(INC_PATH) -o testPotentialShift testPotentialShift.cpp ../Main/libLSMS.a $(LIBS) $(ADD_LIBS)
	cp testPotentialShift $(TOP_DIR)/bin

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<


