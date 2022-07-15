## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

ACLOCAL_AMFLAGS = -I m4
SUBDIRS = pkgconfig scripts maple src testsuite cmake pylibxc

include_HEADERS = xc_version.h libxc.bib
EXTRA_DIST = ChangeLog.md PACKAGING setup.py config.h.cmake.in CMakeLists.txt CITATION

CLEANFILES = *~ *.bak

# These files are in the repository, and the command below should be run 
# to update them whenever references are changed.
# references.c references.h: get_references.pl ../libxc.bib
references:
	(cd $(srcdir)/src && ../scripts/get_references.pl ../libxc.bib)

rpm: dist
	@topdir=`rpmbuild --showrc | @GREP@ ': *_topdir' | cut -f 2`; \
	mv -f $(PACKAGE)-$(VERSION).tar.gz $$topdir/SOURCES; \
	cp build/$(PACKAGE).spec $$topdir/SPECS; \
	rpmbuild -ba $$topdir/SPECS/$(PACKAGE).spec; \
	echo "-------------------------------------------------------------"; \
        echo "  RPMs are built, and are located under:"; \
        echo "       $$topdir/RPMS/<arch>"; \
        echo "       $$topdir/SRPMS"; \
        echo "-------------------------------------------------------------"
#       scp $(RPM_DIR)/RPMS/*/libxc*.rpm marques@www.tddft.org:/server/www/tddft.org/programs/octopus/download/
