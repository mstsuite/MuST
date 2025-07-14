#!/bin/sh

#=====================================================================
# SystemName is the filename found under ./arch
#=====================================================================
SystemName = bolt_intel
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
   MST = 1
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
   else
      ODIR = $(OODIR)/$(SUFFIX)
   endif
else
   ifdef suffix
      SUFFIX = $(suffix)
      ifeq "$(suffix)" "."
         ODIR = $(OODIR)
      else
         ODIR = $(OODIR)/$(suffix)
      endif
   else
      SUFFIX = "."
      ODIR = $(OODIR)
   endif
endif

ifndef EXEC_NAME
   EXEC_NAME = mst2
endif

ADD_LIBS += $(LIBS)

MST_DIR_LINK = MST-2
MSTLIB = $(ODIR)/libmst.a
MPPLIB = $(ODIR)/libmpp.a
IOLIB  = $(ODIR)/iolib.a
IODIR  = $(TOP_DIR)/iolib 
SLULIB = $(ODIR)/slulib.a
ACCLIB = $(ODIR)/libMatAcc.a

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
   ifndef LIBXC_PATH
      ADDITIONAL_TARGETS = libxc
      LIBXC_PATH = $(EXTERN_LIB_PATH)/libxc/LibXC/$(SUFFIX)
   endif
   ifeq (,$(findstring -DLIBXC,$(FPPFLAGS)))
#     if "-DLIBXC" is not part of the FPPFLAGS string
#     ADD_LIBS  += -L$(LIBXC_PATH)/lib -lxcf90 -lxc
      ifneq (,$(wildcard $(LIBXC_PATH)/lib/libxcf90.a))
         ADD_LIBS  += $(LIBXC_PATH)/lib/libxcf90.a $(LIBXC_PATH)/lib/libxc.a
         FPPFLAGS += -DLIBXCF90 -I$(LIBXC_PATH)/include
      else ifneq (,$(wildcard $(LIBXC_PATH)/lib/libxcf03.a))
         ADD_LIBS  += $(LIBXC_PATH)/lib/libxcf03.a $(LIBXC_PATH)/lib/libxc.a
         FPPFLAGS += -DLIBXCF03 -I$(LIBXC_PATH)/include
      else
         ERROR_CODE = error
         ERROR_MESSAGE = "LibXC library file, libxcf03.a or libxcf90.a, is missing. Please check the installation of LibXC on your system."
      endif
   endif
else
   ifeq "$(NotUse_LIBXC)" "0"
      ifndef LIBXC_PATH
         ADDITIONAL_TARGETS = libxc
         LIBXC_PATH = $(EXTERN_LIB_PATH)/libxc/LibXC/$(SUFFIX)
      endif
      ifeq (,$(findstring -DLIBXC,$(FPPFLAGS)))
#        if "-DLIBXC" is not part of the FPPFLAGS string
#        ADD_LIBS  += -L$(LIBXC_PATH)/lib -lxcf90 -lxc
         ifneq (,$(wildcard $(LIBXC_PATH)/lib/libxcf90.a))
            ADD_LIBS  += $(LIBXC_PATH)/lib/libxcf90.a $(LIBXC_PATH)/lib/libxc.a
            FPPFLAGS += -DLIBXCF90 -I$(LIBXC_PATH)/include
         else ifneq (,$(wildcard $(LIBXC_PATH)/lib/libxcf03.a))
            ADD_LIBS  += $(LIBXC_PATH)/lib/libxcf03.a $(LIBXC_PATH)/lib/libxc.a
            FPPFLAGS += -DLIBXCF03 -I$(LIBXC_PATH)/include
         else
            ERROR_CODE = error
            ERROR_MESSAGE = "LibXC library file, libxcf03.a or libxcf90.a, is missing. Please check the installation of LibXC on your system."
         endif
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

ifdef FFLAGS_ALT
   FFLAGS_ALT += -c
else
   FFLAGS_ALT = $(FFLAGS)
endif

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

all: $(ERROR_CODE) $(ADDITIONAL_TARGETS) check-bin check-git-tag main check util

#=====================================================================
# Check bin directory to make sure it exists
#=====================================================================
check-bin:
	if [ ! -d $(OODIR) ]; then echo $(OODIR) "folder does not exist and I am creating one ..."; mkdir $(OODIR); fi
	if [ ! -d $(ODIR) ]; then echo $(ODIR) "folder does not exist and I am creating one ..."; mkdir $(ODIR); fi

check-git-tag:
	if [ ! -d ../.git ]; then echo '"Develop"' > $(ODIR)/git_version; else [[ `(git tag | sort -V | tail -n 1 | wc -m)` -eq 0 ]] && echo '"Develop"' > $(ODIR)/git_version || git tag | sort -V | tail -n 1 | sed 's/^/\"/' | sed 's/$$/\"/' > $(ODIR)/git_version; fi
	echo 'write(myunit,'\''("# ",t5,a,t30,": ",a)'\'') "MuST source code version",' > $(ODIR)/prepend
	echo 'code_version =' `cat $(ODIR)/git_version` > $(ODIR)/MuST_version
	cat $(ODIR)/prepend $(ODIR)/git_version | ( tr -d '\n'; echo; ) > $(ODIR)/git_version.h
#	rm -f $(ODIR)/git_version $(ODIR)/prepend

gen_link:
	cd ../ ; rm -f $(MST_DIR_LINK); ln -s $(TOP_DIR) $(MST_DIR_LINK)
	@echo
	@echo '---------------------------------------------------------------------------------------------'
	@echo '*** WARNING Generic link MST-2 is pointing now to '$(TOP_DIR)' ***'
	@echo '---------------------------------------------------------------------------------------------'
	@echo

main:
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "PREPROC=$(PREPROC)" "ODIR=$(ODIR)" \
"XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "FFLAGS_ALT=$(FFLAGS_ALT)" -f $(TOP_DIR)/lib/Makefile
	if test "$(Acceleration)" = "1"; then cd $(ODIR); \
make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "PREPROC=$(PREPROC)" "CC=$(CC)" "CXX=$(CXX)" \
"ACCEL_CXX=$(ACCEL_CXX)" "ACCEL_OPT=$(ACCEL_OPT)" "CFLAGS=$(CFLAGS)" \
"XLC_I=$(XLC_I)" "ODIR=$(ODIR)" "ARCHV=$(ARCHV)" "ACCLIB=$(ACCLIB)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" -f $(TOP_DIR)/Accelerator/Makefile; fi
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "IOLIB=$(IOLIB)" \
"PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" -f $(TOP_DIR)/iolib/Makefile
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "MPPLIB=$(MPPLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV=$(ARCHV)" -f $(TOP_DIR)/plib/Makefile
	if test "$(SuperLU)" = "1"; then cd $(ODIR); \
make "FC=$(FC)" "FFLAGSSLU=$(FFLAGS)" "SLULIB=$(SLULIB)" "MPI2INCLUDE_PATH=$(MPI2INCLUDE_PATH)" "SLUPATH=$(SLUPATH)" "PREPROC=$(PREPROC)" \
"PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "ODIR=$(ODIR)" "ARCHV=$(ARCHV)" "XLC_I=$(XLC_I)" -f $(TOP_DIR)/slulib/Makefile; fi
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" \
"PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "FFLAGS_ALT=$(FFLAGS_ALT)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "IOLIB=$(IOLIB)" "DSUPERLULIB=$(DSUPERLULIB)" "MPPLIB=$(MPPLIB)" "SLULIB=$(SLULIB)" "SLUPATH=$(SLUPATH)" \
"ODIR=$(ODIR)" "ADDLIBS=$(ADD_LIBS)" "XLC_I=$(XLC_I)" "LD=$(LD)" "Use_FFTW=$(Use_FFTW)" "FFTW_INC=$(FFTW_INC)" "EXEC_NAME=$(EXEC_NAME)" \
-f $(TOP_DIR)/src/Makefile
	if test -d $(HOME)/bin; then echo "copy the executable $(EXEC_NAME) to $(HOME)/bin"; else mkdir $(HOME)/bin; fi
	cp $(ODIR)/$(EXEC_NAME) $(HOME)/bin

check: main
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "IOLIB=$(IOLIB)" "MPPLIB=$(MPPLIB)" \
"PREPROC=$(PREPROC)" "ODIR=$(ODIR)" "ADDLIBS=$(ADD_LIBS)" "XLC_I=$(XLC_I)" \
"PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "Use_FFTW=$(Use_FFTW)" "FFTW_INC=$(FFTW_INC)" "LD=$(LD)" -f $(TOP_DIR)/driver/Makefile

util: check
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "MPICC=$(MPICC)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" \
"PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "IOLIB=$(IOLIB)" "MPPLIB=$(MPPLIB)" "ODIR=$(ODIR)" "LD=$(LD)" "ADDLIBS=$(ADD_LIBS)" "XLC_I=$(XLC_I)" \
-f $(TOP_DIR)/util/Makefile

clear:
	cd $(ODIR); make "ODIR=$(ODIR)" "MSTLIB=$(MSTLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/lib/Makefile
	cd $(ODIR); make "ODIR=$(ODIR)" "IOLIB=$(IOLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/lib/Makefile
	cd $(ODIR); make "ODIR=$(ODIR)" "MPPLIB=$(MPPLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/plib/Makefile 
	cd $(ODIR); make "ODIR=$(ODIR)" "EXEC_NAME=$(EXEC_NAME)" clear -f $(TOP_DIR)/src/Makefile

clear_slulib:
	cd $(ODIR); make "ODIR=$(ODIR)" "SLULIB=$(SLULIB)" clear -f $(TOP_DIR)/slulib/Makefile

clear_util:
	cd $(ODIR); make "ODIR=$(ODIR)" clear -f $(TOP_DIR)/util/Makefile

clean_util:
	cd $(ODIR); make "ODIR=$(ODIR)" clean -f $(TOP_DIR)/util/Makefile

clear_check:
	cd $(ODIR); make "ODIR=$(ODIR)" clear -f $(TOP_DIR)/driver/Makefile

clean_check:
	cd $(ODIR); make "ODIR=$(ODIR)" clean -f $(TOP_DIR)/driver/Makefile

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
	ln -fs $(ODIR)/assk $(PREFIX)/bin/$(SUFFIX)/assk; \
ln -fs $(ODIR)/bireft $(PREFIX)/bin/$(SUFFIX)/bireft; \
ln -fs $(ODIR)/cart2direct $(PREFIX)/bin/$(SUFFIX)/cart2direct; \
ln -fs $(ODIR)/charge $(PREFIX)/bin/$(SUFFIX)/charge; \
ln -fs $(ODIR)/chgaver $(PREFIX)/bin/$(SUFFIX)/chgaver; \
ln -fs $(ODIR)/eosfit $(PREFIX)/bin/$(SUFFIX)/eosfit; \
ln -fs $(ODIR)/extrpot $(PREFIX)/bin/$(SUFFIX)/extrpot; \
ln -fs $(ODIR)/findBoundStates $(PREFIX)/bin/$(SUFFIX)/findBoundStates; \
ln -fs $(ODIR)/findResonance $(PREFIX)/bin/$(SUFFIX)/findResonance; \
ln -fs $(ODIR)/genap $(PREFIX)/bin/$(SUFFIX)/genap; \
ln -fs $(ODIR)/measureVoronoi $(PREFIX)/bin/$(SUFFIX)/measureVoronoi; \
ln -fs $(ODIR)/$(EXEC_NAME) $(PREFIX)/bin/$(SUFFIX)/$(EXEC_NAME); \
ln -fs $(ODIR)/murn $(PREFIX)/bin/$(SUFFIX)/murn; \
ln -fs $(ODIR)/newa $(PREFIX)/bin/$(SUFFIX)/newa; \
ln -fs $(ODIR)/newss $(PREFIX)/bin/$(SUFFIX)/newss; \
ln -fs $(ODIR)/setEvec $(PREFIX)/bin/$(SUFFIX)/setEvec; \
ln -fs $(ODIR)/slice $(PREFIX)/bin/$(SUFFIX)/slice; \
ln -fs $(ODIR)/ssk $(PREFIX)/bin/$(SUFFIX)/ssk
