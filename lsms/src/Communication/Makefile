OBJ = distributeAtoms.o LSMSCommunication.o REWLCommunication.o

all: libCommunication.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<


libCommunication.a : $(OBJ)
	$(AR) -r libCommunication.a $(OBJ) 
	cp libCommunication.a $(TOP_DIR)/lib
