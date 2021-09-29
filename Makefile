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

ifndef SUFFIX
   ifdef suffix
      SUFFIX = $(suffix)
      SUFFIX_d = $(suffix)
   else
      SUFFIX = .
   endif
else
   SUFFIX_d = $(SUFFIX)
endif

%:
	@if [ ! -e ./architecture/$(ArchName) ]; then echo "Architecture file" \"$(ArchName)\" "does not exist under ./architecture directory"; \
exit 1; fi
	@if [ ! -d bin ]; then @echo "bin folder does not exist and I am creating one ..."; mkdir bin; fi
	@if [ "${SUFFIX_d}" ] && [ ! -d bin/$(SUFFIX_d) ]; then @echo "creating subdirectory $(SUFFIX_d) under bin ..."; mkdir bin/$(SUFFIX_d); fi
	@if [ "${SUFFIX_d}" ]; then echo $(SUFFIX_d) | tee bin/.SUFFIX; else echo "." | tee bin/.SUFFIX; fi
	@cd external; ln -fs ../architecture/$(ArchName) architecture.h; make "EXTERNAL=1" "SUFFIX=$(SUFFIX)"
	@cd lsms; ln -fs ../architecture/$(ArchName) architecture.h; make "EXTERN_LIB_PATH=$(MuST_PATH)/external" "ArchName=$(ArchName)" \
"SUFFIX=$(SUFFIX)"
	@cd MST; ln -fs ../architecture/$(ArchName) architecture.h; make "MST=1" "EXTERN_LIB_PATH=$(MuST_PATH)/external" "ArchName=$(ArchName)" \
"SUFFIX=$(SUFFIX)"
	@cd KUBO; ln -fs ../architecture/$(ArchName) architecture.h; make "KUBO=1" "EXTERN_LIB_PATH=$(MuST_PATH)/external" "ArchName=$(ArchName)" \
"SUFFIX=$(SUFFIX)"
	

install:
	$(eval SUFFIX := "$(shell tail -n 1 bin/.SUFFIX)")
	cd lsms; make "PREFIX=$(MuST_PATH)" "SUFFIX=$(SUFFIX)" install
	cd MST; make "PREFIX=$(MuST_PATH)" "SUFFIX=$(SUFFIX)" install
	cd KUBO; make "PREFIX=$(MuST_PATH)" "SUFFIX=$(SUFFIX)" install
	@echo
	@echo '----------------------------------------------------------------------------------------------'
	@echo '*** Generic links pointing to the generated executables are created under ./bin directory *** '
	@echo '----------------------------------------------------------------------------------------------'
	@echo

mst MST:
	@cd MST; make "MST=1" "EXTERN_LIB_PATH=$(MuST_PATH)/external" "SUFFIX=$(SUFFIX)"

lsms:
	@cd lsms; make "EXTERN_LIB_PATH=$(MuST_PATH)/external" "SUFFIX=$(SUFFIX)"

clean: clean-lsms clean-MST clean-external
	$(eval SUFFIX := "$(shell tail -n 1 bin/.SUFFIX)")
	@if [ "${SUFFIX}" ]; then rm -f bin/$(SUFFIX)/*; else rm -f bin/*; fi

clean-external:
	cd external; make "EXTERN_LIB_PATH=$(MuST_PATH)/external" "SUFFIX=$(SUFFIX)" clean

clean-lsms:
	cd lsms; make "SUFFIX=$(SUFFIX)" clean
	@if [ "${SUFFIX}" ]; then rm -f bin/$(SUFFIX)/lsms bin/$(SUFFIX)/wl-lsms; else rm -f bin/lsms bin/wl-lsms; fi

clean-MST:
	cd MST; make "SUFFIX=$(SUFFIX)" distclean

distclean:
	@cd external; make "SUFFIX=$(SUFFIX)" distclean
	@cd lsms; make "SUFFIX=$(SUFFIX)" distclean
	@cd MST; make "SUFFIX=$(SUFFIX)" distclean
	rm -rf bin
