export TOP_DIR = $(shell pwd)
export INC_PATH =
export LIBS =

ifneq (,$(wildcard ./architecture.h))
include architecture.h 
else
include $(ARCH_FILE)
endif

ifndef HDF5_PATH
export ADDITIONAL_TARGETS=hdf5
export HDF5_PATH=$(TOP_DIR)/hdf5/HDF5
endif

ifndef LUA_PATH
export ADDITIONAL_TARGETS+=lua
export LUA_PATH=$(TOP_DIR)/lua/LUA
endif

ifndef LIBXC_PATH
export ADDITIONAL_TARGETS+=libxc
export LIBXC_PATH=$(TOP_DIR)/libxc/LibXC
endif

ifndef FFTW_PATH
#export ERROR_TARGETS=no_fftw
export ADDITIONAL_TARGETS+=fftw
export LIBXC_PATH=$(TOP_DIR)/fftw/FFTW
endif

ifndef P3DFFT_PATH
export ADDITIONAL_TARGETS+=p3dfft
export LIBXC_PATH=$(TOP_DIR)/p3dfft/P3DFFT
export P3DFFT_CONFIG_FLAGS += --with-fftw=$(TOP_DIR)/fftw/FFTW
endif

all: $(ERROR_TARGETS) $(ADDITIONAL_TARGETS)

.PHONY: hdf5 lua libxc fftw p3dfft clean distclean

hdf5:
	if [ ! -e hdf5/HDF5/lib/libhdf5.a ]; then \
DIR_NAME="$(shell tar -tzf hdf5/hdf5.tar.gz | head -1 | cut -f1 -d"/")" && cd hdf5 && rm -f HDF5 && tar xfz hdf5.tar.gz && \
ln -fs $$DIR_NAME HDF5; cd HDF5 && echo $$PWD && ./configure --prefix=$$PWD $(HDF5_CONFIG_FLAGS) && make && make -i install; fi

lua:
	if [ ! -e lua/LUA/lib/liblua.a ]; then \
DIR_NAME="$(shell tar -tzf lua/lua.tar.gz | head -1 | cut -f1 -d"/")" && cd lua && rm -f LUA && tar xfz lua.tar.gz && \
ln -fs $$DIR_NAME LUA && cd LUA && echo $$PWD && make INSTALL_TOP=$$PWD PLAT=generic CC="$(CC) -std=gnu99" ; make INSTALL_TOP=$$PWD PLAT=generic install; fi

libxc:
	if [ ! -e libxc/LibXC/lib/libxc.a ]; then \
DIR_NAME="$(shell tar -tzf libxc/libxc.tar.gz | head -1 | cut -f1 -d"/")" && cd libxc && rm -f LibXC && tar xfz libxc.tar.gz && \
ln -fs $$DIR_NAME LibXC; cd LibXC && echo $$PWD && ./configure --prefix=$$PWD $(LIBXC_CONFIG_FLAGS) && make && make install; fi

fftw:
	if [ ! -e fftw/FFTW/lib/libfftw3.a ]; then \
DIR_NAME="$(shell tar -tzf fftw/fftw.tar.gz | head -1 | cut -f1 -d"/")" && cd fftw && rm -f FFTW && tar xfz fftw.tar.gz && \
ln -fs $$DIR_NAME FFTW; cd FFTW && echo $$PWD && ./configure --prefix=$$PWD --enable-fortran $(FFTW_CONFIG_FLAGS) && make && make install; fi

p3dfft: fftw
	if [ ! -e p3dfft/P3DFFT/lib/libp3dfft.a ]; then \
DIR_NAME="$(shell tar -tzf p3dfft/p3dfft.tar.gz | head -1 | cut -f1 -d"/")" && cd p3dfft && rm -f P3DFFT && tar xfz p3dfft.tar.gz && \
ln -fs $$DIR_NAME P3DFFT; cd P3DFFT && echo $$PWD && ./configure --prefix=$$PWD --enable-fftw $(P3DFFT_CONFIG_FLAGS) && make && make -i install; fi

clean:
	if [ -d hdf5/HDF5 ]; then cd hdf5/HDF5; make clean; rm -rf lib; fi
	if [ -d lua/LUA ]; then cd lua/LUA; make clean; rm -rf lib; fi
	if [ -d libxc/LibXC ]; then cd libxc/LibXC; make clean; rm -rf lib; fi
	if [ -d p3dfft/P3DFFT ]; then cd p3dfft/P3DFFT; make clean; rm -rf lib; fi
	if [ -d fftw/FFTW ]; then cd fftw/FFTW; make clean; rm -rf lib; fi

distclean:
	@if [ -d hdf5/HDF5 ]; then rm -rf "$(shell readlink -f hdf5/HDF5)" hdf5/HDF5; fi
	@if [ -d lua/LUA ]; then rm -rf "$(shell readlink -f lua/LUA)" lua/LUA; fi
	@if [ -d libxc/LibXC ]; then rm -rf "$(shell readlink -f libxc/LibXC)" libxc/LibXC; fi
	@if [ -d p3dfft/P3DFFT ]; then rm -rf "$(shell readlink -f p3dfft/P3DFFT)" p3dfft/P3DFFT p3dfft/share; fi
	@if [ -d fftw/FFTW ]; then rm -rf "$(shell readlink -f fftw/FFTW)" fftw/FFTW; fi
	rm -f architecture.h

no_fftw:
	@echo "Error: FFTW_PATH is not defined in arch file" \"$(ArchName)\"; exit 1
