export INC_PATH += -I $(TOP_DIR)/lua/include -I $(TOP_DIR)/include -I $(TOP_DIR)/src
export LIBS += -L$(TOP_DIR)/lib -lLSMSLua -lCommunication \
               -lMultipleScattering -lSingleSite -lCore -lVORPOL -lAccelerator \
               -lMadelung -lPotential -lTotalEnergy -lMisc

all: lsms wl-lsms
# all: lsms wl-lsms rewl-lsms (rewl-lsms has problems with the seed generation)

clean:
	cd RadialGrid && $(MAKE) clean
	cd VORPOL && $(MAKE) clean
	cd Madelung && $(MAKE) clean
	cd Potential && $(MAKE) clean
	cd SingleSite && $(MAKE) clean
	cd MultipleScattering && $(MAKE) clean
	cd Core && $(MAKE) clean
	cd LuaInterface && $(MAKE) clean
	cd Communication && $(MAKE) clean
	cd Misc && $(MAKE) clean
	cd Main && $(MAKE) clean
	cd Accelerator && $(MAKE) clean
	cd TotalEnergy && $(MAKE) clean

Libraries:
	cd RadialGrid && $(MAKE)
	cd VORPOL && $(MAKE)
	cd Madelung && $(MAKE)
	cd Potential && $(MAKE)
	cd SingleSite && $(MAKE)
	cd MultipleScattering && $(MAKE)
	cd Core && $(MAKE)
	cd LuaInterface && $(MAKE)
	cd Communication && $(MAKE)
	cd Misc && $(MAKE)
	cd Accelerator && $(MAKE)
	cd TotalEnergy && $(MAKE)

lsms_inter: Libraries
	cd Main && $(MAKE) lsms_inter

lsms: Libraries
	cd Main && $(MAKE) lsms

wl-lsms: Libraries
	cd Main && $(MAKE) wl-lsms

rewl-lsms: Libraries
	cd Main && $(MAKE) rewl-lsms

zblock_lu_driver: 
	cd Accelerator && $(MAKE)
	cd Main && $(MAKE) zblock_lu_driver

test: Libraries lsms wl-lsms rewl-lsms
	cd Test && $(MAKE)
