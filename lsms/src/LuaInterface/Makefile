# OBJ = LuaInterface.o RadialGrid_lua.o RadialPotential_lua.o LuaSupport.o
OBJ = LuaInterface.o LuaSupport.o

all: libLSMSLua.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<


libLSMSLua.a : $(OBJ)
	$(AR) -r libLSMSLua.a $(OBJ) 
	cp libLSMSLua.a $(TOP_DIR)/lib
