export TOP_DIR = $(shell pwd)
export INC_PATH =
export LIBS := $(TOP_DIR)/mjson/mjson.a

ifndef EXTERN_LIB_PATH
EXTERN_LIB_PATH = $(TOP_DIR)/../external
endif

ifneq (,$(wildcard ./architecture.h))
include architecture.h 
endif

ifndef SUFFIX
SUFFIX=.
else
OBJ_DIR = $(SUFFIX)
endif

ifndef HDF5_PATH
#export ERROR_TARGETS=no_hdf5
export ADDITIONAL_TARGETS+=no_hdf5
HDF5_PATH=$(EXTERN_LIB_PATH)/hdf5/HDF5/$(SUFFIX)
endif

ifndef LUA_PATH
export ADDITIONAL_TARGETS+=no_lua
LUA_PATH=$(EXTERN_LIB_PATH)/lua/LUA/$(SUFFIX)
endif

export ADD_LIBS += -Wl,-rpath $(HDF5_PATH)/lib -L$(HDF5_PATH)/lib -lhdf5 -L$(LUA_PATH)/lib -llua
export INC_PATH +=  -I$(HDF5_PATH)/include -I$(LUA_PATH)/include

ifndef LIBXC_PATH
export ADDITIONAL_TARGETS+=no_libxc
LIBXC_PATH=$(EXTERN_LIB_PATH)/libxc/LibXC/$(SUFFIX)
endif

export ADD_LIBS += -L$(LIBXC_PATH)/lib -lxc
export INC_PATH += -I$(LIBXC_PATH)/include
export OPT_DEFINES += -DUSE_LIBXC

ifdef HAS_BACKTRACE
  export OPT_DEFINES += -DHAS_BACKTRACE
endif

export CC += $(CFLAGS) $(OPT_OPENMP) $(OPT_DEFINES)
export CXX += $(CXXFLAGS) $(OPT_OPENMP) $(OPT_DEFINES)
export F77 += $(FFLAGS) $(OPT_OPENMP)
export ODIR=$(TOP_DIR)/bin/$(SUFFIX)

all: odir $(ERROR_TARGETS) $(ADDITIONAL_TARGETS) libmjson LSMS
# all: liblua $(ADDITIONAL_TARGETS) libmjson LSMS
# all: liblua libjson $(ADDITIONAL_TARGETS) libmjson LSMS
# all: liblua LSMS Documentation

.PHONY: $(ERROR_TARGETS) $(ADDITIONAL_TARGETS) liblua libxc clean LSMS Documentation libmjson
# .PHONY: libxc clean LSMS Documentation liblua libjson libmjson Tools

odir:
	@if [ ! -d $(ODIR) ]; then echo "$(ODIR) folder does not exist and I am creating one ..."; mkdir $(ODIR); fi

clean:
	cd mjson && $(MAKE) clean
	cd src && $(MAKE) clean
	cd doc && $(MAKE) clean
	cd lib && $(MAKE) clean
	cd Tools && $(MAKE) clean
#	cd lua && $(MAKE) clean
#	cd libjson && $(MAKE) clean

distclean: clean
	rm -f architecture.h $(ODIR)/lsms $(ODIR)/wl-lsms
	@if [ "${OBJ_DIR}" ] && [ -d bin/$(OBJ_DIR) ]; then rm -rf bin/$(OBJ_DIR); fi

LSMS: libmjson $(ADDITIONAL_TARGETS)
# LSMS: liblua libmjson $(ADDITIONAL_TARGETS)
	cd src && $(MAKE)

# Tools: liblua 
#	cd Tools && $(MAKE)

Documentation:
	cd doc && $(MAKE)

# liblua:
#	cd lua; $(MAKE); $(MAKE) local

libmjson:
	cd mjson && $(MAKE)

# libjson:
#	cd libjson && $(MAKE)

zblock_lu_driver: liblua 
	cd src && $(MAKE) zblock_lu_driver


test: $(ADDITIONAL_TARGETS)
# test: liblua $(ADDITIONAL_TARGETS) libmjson
	cd src && $(MAKE) test

install:
	ln -fs $(ODIR)/lsms $(PREFIX)/bin/$(SUFFIX)/lsms; ln -fs $(ODIR)/wl-lsms $(PREFIX)/bin/$(SUFFIX)/wl-lsms

no_hdf5:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(TOP_DIR)/architecture.h" hdf5
#	@echo "Error: HDF5_PATH is not defined in arch file" \"$(ArchName)\"; exit 1

no_lua liblua:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(TOP_DIR)/architecture.h" lua

no_libxc libxc:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(TOP_DIR)/architecture.h" libxc
