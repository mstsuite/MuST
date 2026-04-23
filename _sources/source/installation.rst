************
Installation
************

Preferred method:
###################

1. Under architecture/ directory, create/modify an architecture file by following an existing example file;

2. In the top directory (MuST/), run the following commands to build executables
make architecture-file-name (e.g., make linux-intel-nogpu)
make install

Note --
make clean: delete the object, library, executable files under lsms and MST from installation
make distclean: delete the object, library, executable, and architecture.h files under lsms and MST from installation; also
                delete the executables under bin/.

Alternative method:
###################

The code MST (under MST/) and LSMS/WL-LSMS (under lsms/) can be built separately by running make under MST
and lsms. The executables can be found under MST/bin and lsms/bin, respectively. It requires to create
archietecture.h under MST and lsms using symbolic link. Steps are as follows:
* To build MST,
1. cd MST
2. set SystemName in Makefile (at line 6) to a proper name, or execute the following command:
   ln -s arch/architecture_file architecture.h
3. make
* To build LSMS/WL-LSMS,
1. cd lsms
2. ln -s arch/architecture_file architecture.h
3. make

Notes to the user of Fedora systems
###################################
MST may require using External Data Representation (XDR) library to store potential and charge density data.
Unfortunately, the latest Fedora Linux system does not place the library in conventional locations. Therefore,
before installing MuST or MST, please make sure that /usr/include/tirpc and /usr/include/tirpc/rpc exist. If not,
you need to ask your system administrator to istall libtirpc and librirpc-devel for you, or to run the following command
if you have the sys-admin privilige:
   sudo dnf install libtirpc libtirpc-devel
