# ==================================================================================================
# Makefile for installing MuST
#    The current directory (../MuST) has following sub-directories
#        architecture/ : contains architecture files for building MuST. The recommended name convention
#                        of an achitecture file is as follows: 
#                              system_name-compiler_name-MPI_name-GPU_status
#                        where system_name can be summit, bridges, frontier, etc.; GPU_status can be
#                        nogpu, gpu, amd, nvidia, etc.; compiler name can be intel, pgi, etc.; MPI_name
#                        can be mpich, openmpi, etc. If the compiler and MPI library are made by the same 
#                        ventor, e.g., Intel, compiler_name and MPI_name can be made into one character 
#                        string, e.g., intel.
#                        A templete "linux-gnu-openmpi-nogpu" is given as an example architecture file.
#        external/ : contains external libraries required or optionally required by MuST, 
#                   e.g., FFTW, Lua, P3DFFT, and LibXC.
#        lsms/ : contains LSMS and WL-LSMS codes targeted to extreme performance on petascale/exascale systems
#        MST/ : contains MST packages with wide physics capabilities, e.g. FP/MT KKR/LSMS/KKR-CPA, etc.
#
#    The steps to build MuST:
#        1. Create an architecture file and place it under architecture/.
#        2. make architecture_file_name
#        3. make install
# ==================================================================================================
#!/bin/sh

ArchName = $(MAKECMDGOALS)

MuST_PATH = $(shell pwd)

%:
	@if [ ! -e ./architecture/$(ArchName) ]; then echo "Architecture file" \"$(ArchName)\" "does not exist under ./architecture directory"; \
exit 1; fi
	@if [ ! -d bin ]; then @echo "bin folder does not exist and I am creating one ..."; mkdir bin; fi
	@cd external; ln -fs ../architecture/$(ArchName) architecture.h; make "EXTERNAL=1"
	@cd lsms; ln -fs ../architecture/$(ArchName) architecture.h; make "EXTERN_LIB_PATH=$(MuST_PATH)/external" "ArchName=$(ArchName)"
	@cd MST; ln -fs ../architecture/$(ArchName) architecture.h; make "MST=1" "EXTERN_LIB_PATH=$(MuST_PATH)/external" "ArchName=$(ArchName)"

install:
	cd lsms; make "PREFIX=$(MuST_PATH)" install
	cd MST; make "PREFIX=$(MuST_PATH)" install
	@echo
	@echo '----------------------------------------------------------------------------------------------'
	@echo '*** Generic links pointing to the generated executables are created under ./bin directory *** '
	@echo '----------------------------------------------------------------------------------------------'
	@echo

mst MST:
	@cd MST; make "MST=1" "EXTERN_LIB_PATH=$(MuST_PATH)/external"

lsms:
	@cd lsms; make "EXTERN_LIB_PATH=$(MuST_PATH)/external"

clean: clean-lsms clean-MST clean-external
	rm -f bin/*

clean-external:
	cd external; make "EXTERN_LIB_PATH=$(MuST_PATH)/external" clean

clean-lsms:
	cd lsms; make clean
	rm -f bin/lsms bin/wl-lsms

clean-MST:
	cd MST; make deepclean

distclean:
	@cd external; make distclean
	@cd lsms; make distclean
	@cd MST; make distclean
	@rm -f bin/*
