export TOP_DIR = $(shell pwd)
export INC_PATH =
export LIBS =

ifneq (,$(wildcard ./architecture.h))
include architecture.h 
else
include $(ARCH_FILE)
endif

ifndef SUFFIX
export SUFFIX=.
endif

ifndef HDF5_PATH
export ADDITIONAL_TARGETS=hdf5
export HDF5_PATH=$(TOP_DIR)/hdf5/HDF5/$(SUFFIX)
endif

ifndef LUA_PATH
export ADDITIONAL_TARGETS+=lua
export LUA_PATH=$(TOP_DIR)/lua/LUA/$(SUFFIX)
endif

ifndef LIBXC_PATH
export ADDITIONAL_TARGETS+=libxc
export LIBXC_PATH=$(TOP_DIR)/libxc/LibXC/$(SUFFIX)
endif

ifndef FFTW_PATH
#export ERROR_TARGETS=no_fftw
export ADDITIONAL_TARGETS+=fftw
export FFTW_PATH=$(TOP_DIR)/fftw/FFTW/$(SUFFIX)
endif

ifndef P3DFFT_PATH
export ADDITIONAL_TARGETS+=p3dfft
export P3DFFT_PATH=$(TOP_DIR)/p3dfft/P3DFFT/$(SUFFIX)
export P3DFFT_CONFIG_FLAGS += --with-fftw=$(FFTW_PATH)
endif

all: $(ERROR_TARGETS) $(ADDITIONAL_TARGETS)

# .PHONY: hdf5 lua libxc fftw p3dfft clean distclean
# Since we are not using lsms package for the moment, we temporarily disabled compiling hdf5 and lua --2/18/2022
.PHONY: libxc fftw p3dfft clean distclean

hdf5:
	if [ ! -d hdf5/HDF5 ]; then \
DIR_NAME="$(shell tar -tzf hdf5/hdf5.tar.gz | head -1 | cut -f1 -d"/")" && cd hdf5 && tar xfz hdf5.tar.gz && \
ln -fs $$DIR_NAME HDF5; fi
	if [ ! -e hdf5/HDF5/$(SUFFIX)/lib/libhdf5.a ]; then \
cd hdf5/HDF5 && echo $$PWD && ./configure --prefix=$$PWD/$(SUFFIX) $(HDF5_CONFIG_FLAGS) && make clean && make && make -i install; fi

lua:
	if [ ! -d lua/LUA ]; then \
DIR_NAME="$(shell tar -tzf lua/lua.tar.gz | head -1 | cut -f1 -d"/")" && cd lua && rm -f LUA && tar xfz lua.tar.gz && \
ln -fs $$DIR_NAME LUA; fi
	if [ ! -e lua/LUA/$(SUFFIX)/lib/liblua.a ]; then \
cd lua/LUA && echo $$PWD && make clean; make INSTALL_TOP=$$PWD/$(SUFFIX) PLAT=generic CC="$(CC) -std=gnu99" ; \
make INSTALL_TOP=$$PWD/$(SUFFIX) PLAT=generic install; fi

libxc:
	if [ ! -d libxc/LibXC ]; then \
DIR_NAME="$(shell tar -tzf libxc/libxc.tar.gz | head -1 | cut -f1 -d"/")" && cd libxc && rm -f LibXC && tar xfz libxc.tar.gz && \
ln -fs $$DIR_NAME LibXC; fi
	if [ ! -e libxc/LibXC/$(SUFFIX)/lib/libxc.a ]; then \
cd libxc/LibXC && echo $$PWD && ./configure --prefix=$$PWD/$(SUFFIX) $(LIBXC_CONFIG_FLAGS) && make clean && make && make install; fi

fftw:
	if [ ! -d fftw/FFTW ]; then \
DIR_NAME="$(shell tar -tzf fftw/fftw.tar.gz | head -1 | cut -f1 -d"/")" && cd fftw && rm -f FFTW && tar xfz fftw.tar.gz && \
ln -fs $$DIR_NAME FFTW; fi
	if [ ! -e fftw/FFTW/$(SUFFIX)/lib/libfftw3.a ]; then \
cd fftw/FFTW && echo $$PWD && ./configure --prefix=$$PWD/$(SUFFIX) --enable-fortran $(FFTW_CONFIG_FLAGS) && make clean && \
make && make install; fi

p3dfft:
	if [ ! -d p3dfft/P3DFFT ]; then \
DIR_NAME="$(shell tar -tzf p3dfft/p3dfft.tar.gz | head -1 | cut -f1 -d"/")" && cd p3dfft && rm -f P3DFFT && tar xfz p3dfft.tar.gz && \
ln -fs $$DIR_NAME P3DFFT; fi
	if [ ! -e p3dfft/P3DFFT/$(SUFFIX)/lib/libp3dfft.a ]; then \
cd p3dfft/P3DFFT && echo $$PWD && ./configure --prefix=$$PWD/$(SUFFIX) --enable-fftw $(P3DFFT_CONFIG_FLAGS) && make clean; \
if [ "${MPIFH}" ]; then echo $$PWD && sed -i 's/mpifh/$(MPIFH)/g' ./sample/C/Makefile; fi && \
echo $$PWD && make && make -i install; fi

clean:
	if [ -d hdf5/HDF5 ]; then cd hdf5/HDF5; make clean; rm -rf $(SUFFIX)/lib; fi
	if [ -d lua/LUA ]; then cd lua/LUA; make clean; rm -rf $(SUFFIX)/lib; fi
	if [ -d libxc/LibXC ]; then cd libxc/LibXC; make clean; rm -rf $(SUFFIX)/lib; fi
	if [ -d p3dfft/P3DFFT ]; then cd p3dfft/P3DFFT; make clean; rm -rf $(SUFFIX)/lib; fi
	if [ -d fftw/FFTW ]; then cd fftw/FFTW; make clean; rm -rf $(SUFFIX)/lib; fi

distclean:
	@if [ -d hdf5/HDF5 ]; then rm -rf "$(shell readlink -f hdf5/HDF5)" hdf5/HDF5; fi
	@if [ -d lua/LUA ]; then rm -rf "$(shell readlink -f lua/LUA)" lua/LUA; fi
	@if [ -d libxc/LibXC ]; then rm -rf "$(shell readlink -f libxc/LibXC)" libxc/LibXC; fi
	@if [ -d p3dfft/P3DFFT ]; then rm -rf "$(shell readlink -f p3dfft/P3DFFT)" p3dfft/P3DFFT p3dfft/share; fi
	@if [ -d fftw/FFTW ]; then rm -rf "$(shell readlink -f fftw/FFTW)" fftw/FFTW; fi
	rm -f architecture.h

no_fftw:
	@echo "Error: FFTW_PATH is not defined in arch file" \"$(ArchName)\"; exit 1
