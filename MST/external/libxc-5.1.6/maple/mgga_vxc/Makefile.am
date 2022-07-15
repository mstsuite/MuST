## Copyright (C) 2017 M. Marques
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

mgga_order_4 = \
  mgga_x_tb09.mpl mgga_x_2d_prp10.mpl

dist_noinst_SCRIPTS = $(mgga_order_4)

maple2c:
	@for funct in $(mgga_order_4); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 4; \
	done

