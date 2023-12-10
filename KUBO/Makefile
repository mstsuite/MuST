#!/bin/sh

#=====================================================================
# SystemName is the filename found under ./arch
#=====================================================================
SystemName = osx-gnu-openmpi
# SystemName = bolt_pgi
# SystemName = bolt_pgi_mpich

#=====================================================================
# Paths and internal libraries setup - This unlikely needs to be changed
#=====================================================================
export TOP_DIR = $(shell pwd)
OODIR=$(TOP_DIR)/bin

ifndef EXTERN_LIB_PATH
   EXTERN_LIB_PATH=$(TOP_DIR)/../external
endif

ifneq (,$(wildcard ./architecture.h))
   KUBO = 1
   include architecture.h
   ARCH_FILE = $(TOP_DIR)/architecture.h
else
   include ./arch/$(SystemName)
   ARCH_FILE = $(TOP_DIR)/arch/$(SystemName)
endif

#=====================================================================
# No need to change following settings, unless it is necessary
#=====================================================================

ifdef SUFFIX
   ifeq "$(SUFFIX)" "."
      ODIR = $(OODIR)
      MST_ODIR = $(TOP_DIR)/../MST/bin
   else
      ODIR = $(OODIR)/$(SUFFIX)
      MST_ODIR = $(TOP_DIR)/../MST/bin/$(SUFFIX)
   endif
else
   ifdef suffix
      SUFFIX = $(suffix)
      ifeq "$(suffix)" "."
         ODIR = $(OODIR)
         MST_ODIR = $(TOP_DIR)/../MST/bin
      else
         ODIR = $(OODIR)/$(suffix)
         MST_ODIR = $(TOP_DIR)/../MST/bin/$(suffix)
      endif
   else
      SUFFIX = "."
      ODIR = $(OODIR)
      MST_ODIR = $(TOP_DIR)/../MST/bin
   endif
endif

ifndef EXEC_NAME
   EXEC_NAME = kubo
endif

ADD_LIBS += $(LIBS)

MST_DIR_LINK = KUBO
MSTLIB = $(MST_ODIR)/libmst.a
MPPLIB = $(MST_ODIR)/libmpp.a
IOLIB  = $(MST_ODIR)/iolib.a
IODIR  = $(TOP_DIR)/../MST/iolib 
SLULIB = $(MST_ODIR)/slulib.a
ACCLIB = $(MST_ODIR)/libMatAcc.a

SuperLU=0
ifeq "$(SuperLU)" "0"
   SLULIB =
   DSUPERLULIB =
endif
 
ifeq "$(Acceleration)" "1"
   ifndef ACCEL_PATH
      ERROR_CODE = error
      ERROR_MESSAGE = "ACCEL_PATH is undefined..."
   endif
   ifndef ACCEL
      ERROR_CODE = error
      ERROR_MESSAGE = "ACCEL is undefined..."
   endif
   ifndef ACCEL_CXX
      ERROR_CODE = error
      ERROR_MESSAGE = "ACCEL_CXX is undefined..."
   endif
   ACCEL_OPT = -I$(ACCEL_PATH)/include -I. -DLSMS2 -D$(ACCEL)
   ADD_LIBS += $(ACCLIB)
   FPPFLAGS += -DACCEL -DCUDA
endif

ifdef ARCHV_LSMS
   ARCHV = $(ARCHV_LSMS) $(ARCHV_OPT)
endif

ifndef ARCHV
   ARCHV = ar -r
endif

ifndef PREPROC
   PREPROC = $(FC)
endif

ifndef XLC_I
   XLC_I = 0
endif

ifndef NotUse_LIBXC
   ifneq "$(NotUse_LIBXC)" "0"
      ifndef LIBXC_PATH
         ADDITIONAL_TARGETS = libxc
         LIBXC_PATH = $(EXTERN_LIB_PATH)/libxc/LibXC/$(SUFFIX)
      endif
      ifeq (,$(findstring -DLIBXC,$(FPPFLAGS)))
# ADD_LIBS  += -L$(LIBXC_PATH)/lib -lxcf90 -lxc
         ADD_LIBS  += $(LIBXC_PATH)/lib/libxcf90.a $(LIBXC_PATH)/lib/libxc.a
         FPPFLAGS += -DLIBXC5 -I$(LIBXC_PATH)/include
      endif
   endif
endif

ifndef Use_FFTW
   Use_FFTW = 1
endif

ifndef FFTW_PATH
   ADDITIONAL_TARGETS += fftw
   FFTW_PATH = $(EXTERN_LIB_PATH)/fftw/FFTW/$(SUFFIX)
endif
FFTW_INC = $(FFTW_PATH)/include

ifndef NotUse_P3DFFT
   ifndef P3DFFT_PATH
      ADDITIONAL_TARGETS += p3dfft
      P3DFFT_PATH = $(EXTERN_LIB_PATH)/p3dfft/P3DFFT/$(SUFFIX)
   endif
   FFT_FLAGS = -DFFTW -I$(FFTW_INC) -DP3DFFT -I$(P3DFFT_PATH)/include
   ADD_LIBS += -L$(P3DFFT_PATH)/lib -lp3dfft -L$(FFTW_PATH)/lib -lfftw3_mpi -lfftw3
else
   ifeq "$(NotUse_P3DFFT)" "0"
      ifndef P3DFFT_PATH
         ADDITIONAL_TARGETS += p3dfft
         P3DFFT_PATH = $(EXTERN_LIB_PATH)/p3dfft/P3DFFT/$(SUFFIX)
      endif
      FFT_FLAGS = -DFFTW -I$(FFTW_INC) -DP3DFFT -I$(P3DFFT_PATH)/include
      ADD_LIBS += -L$(P3DFFT_PATH)/lib -lp3dfft -L$(FFTW_PATH)/lib -lfftw3_mpi -lfftw3
   else
      FFT_FLAGS = -DFFTW -I$(FFTW_INC)
      ADD_LIBS += -L$(FFTW_PATH)/lib -lfftw3_mpi -lfftw3
   endif
endif

PREPROC_FFLAGS = $(FPPDEFS) $(FPPFLAGS) $(FFT_FLAGS)
PREPROC_CFLAGS = $(CPPDEFS) -DUnderscore
CFLAGS += $(PREPROC_CFLAGS) -c
FFLAGS += -c
ACCEL_OPT += -c
F77FLAGS = $(FFLAGS)

UNAME_S := $(shell uname -s)
ifdef XDR_INCLUDE
   CFLAGS += $(XDR_INCLUDE)
   ADD_LIBS += $(XDR_LIB)
else ifneq (,$(wildcard /usr/include/rpc/xdr.h))
   CFLAGS += -I/usr/include/rpc
else ifneq (,$(wildcard /usr/include/tirpc/rpc/xdr.h))
   CFLAGS += -I/usr/include/tirpc -I/usr/include/tirpc/rpc
   ADD_LIBS += -ltirpc
else ifneq (,$(wildcard /usr/include/gssrpc/xdr.h))
   CFLAGS += -I/usr/include/gssrpc
   ADD_LIBS += -lgssrpc
else ifneq ($(UNAME_S),Darwin)
   ERROR_CODE = error
   ERROR_MESSAGE = "In your system, /usr/include/rpc or /usr/include/tirpc is missing. Please install librpc or libtirpc."
endif

#=====================================================================

all: $(ERROR_CODE) $(ADDITIONAL_TARGETS) check-bin check-git-tag main

#=====================================================================
# Check bin directory to make sure it exists
#=====================================================================
check-bin:
	if [ ! -d $(OODIR) ]; then echo $(OODIR) "folder does not exist and I am creating one ..."; mkdir $(OODIR); fi
	if [ ! -d $(ODIR) ]; then echo $(ODIR) "folder does not exist and I am creating one ..."; mkdir $(ODIR); fi

check-git-tag:
	if [ ! -d ../.git ]; then echo '"Develop"' > $(ODIR)/git_version; else [[ `(git tag | sort -V | tail -n 1 | wc -m)` -eq 0 ]] && echo '"Develop"' > $(ODIR)/git_version || git tag | sort -V | tail -n 1 | sed 's/^/\"/' | sed 's/$$/\"/' > $(ODIR)/git_version; fi
	echo 'write(myunit,'\''("# ",t5,a,t30,": ",a)'\'') "MuST source code version",' > $(ODIR)/prepend
	cat $(ODIR)/prepend $(ODIR)/git_version | ( tr -d '\n'; echo; ) > src/git_version.h
#	rm -f $(ODIR)/git_version $(ODIR)/prepend

gen_link:
	cd ../ ; rm -f $(MST_DIR_LINK); ln -s $(TOP_DIR) $(MST_DIR_LINK)
	@echo
	@echo '---------------------------------------------------------------------------------------------'
	@echo '*** WARNING Generic link KUBO is pointing now to '$(TOP_DIR)' ***'
	@echo '---------------------------------------------------------------------------------------------'
	@echo

main:
#	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "PREPROC=$(PREPROC)" "ODIR=$(ODIR)" \
#"XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" -f $(TOP_DIR)/lib/Makefile
#	if test "$(Acceleration)" = "1"; then cd $(ODIR); \
#make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "PREPROC=$(PREPROC)" "CC=$(CC)" "CXX=$(CXX)" \
#"ACCEL_CXX=$(ACCEL_CXX)" "ACCEL_OPT=$(ACCEL_OPT)" "CFLAGS=$(CFLAGS)" \
#"XLC_I=$(XLC_I)" "ODIR=$(ODIR)" "ARCHV=$(ARCHV)" "ACCLIB=$(ACCLIB)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" -f $(TOP_DIR)/Accelerator/Makefile; fi
#	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "IOLIB=$(IOLIB)" \
#"PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
#"CC=$(CC)" "CFLAGS=$(CFLAGS)" "ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" -f $(TOP_DIR)/iolib/Makefile
#	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "MPPLIB=$(MPPLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" -f $(TOP_DIR)/plib/Makefile
#	if test "$(SuperLU)" = "1"; then cd $(ODIR); \
#make "FC=$(FC)" "FFLAGSSLU=$(FFLAGS)" "SLULIB=$(SLULIB)" "MPI2INCLUDE_PATH=$(MPI2INCLUDE_PATH)" "SLUPATH=$(SLUPATH)" "PREPROC=$(PREPROC)" \
#"PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "ODIR=$(ODIR)" "ARCHV=$(ARCHV)" "XLC_I=$(XLC_I)" -f $(TOP_DIR)/slulib/Makefile; fi
	cd $(ODIR); make "MST_ODIR=$(MST_ODIR)" "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" \
"PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "IOLIB=$(IOLIB)" "DSUPERLULIB=$(DSUPERLULIB)" "MPPLIB=$(MPPLIB)" "SLULIB=$(SLULIB)" "SLUPATH=$(SLUPATH)" \
"ODIR=$(ODIR)" "ADDLIBS=$(ADD_LIBS)" "XLC_I=$(XLC_I)" "LD=$(LD)" "Use_FFTW=$(Use_FFTW)" "FFTW_INC=$(FFTW_INC)" "EXEC_NAME=$(EXEC_NAME)" \
-f $(TOP_DIR)/src/Makefile
	if test -d $(HOME)/bin; then echo "copy the executable $(EXEC_NAME) to $(HOME)/bin"; else mkdir $(HOME)/bin; fi
	cp $(ODIR)/$(EXEC_NAME) $(HOME)/bin

clear:
	cd $(ODIR); make "ODIR=$(ODIR)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/src/Makefile

clear_src:
	cd $(ODIR); make "ODIR=$(ODIR)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/src/Makefile

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod $(ODIR)/*.a

distclean:
	rm -rf $(OODIR) architecture.h

error:
	@echo $(ERROR_MESSAGE); exit 1

libxc:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(ARCH_FILE)" libxc

fftw:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(ARCH_FILE)" fftw

p3dfft:
	@cd $(EXTERN_LIB_PATH) && make "ARCH_FILE=$(ARCH_FILE)" p3dfft

install:
	ln -fs $(ODIR)/kubo $(PREFIX)/bin/$(SUFFIX)/kubo
