all: libRadialGrid.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

libRadialGrid.a: RadialGrid.o
	$(AR) -r libRadialGrid.a RadialGrid.o
	cp libRadialGrid.a $(TOP_DIR)/lib
