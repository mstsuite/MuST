************
Installation
************

Preferred Method
################

1. Under the ``architecture/`` directory, create or modify an architecture file
   by following an existing example.

.. toctree::
   :maxdepth: 2

   architec

2. In the top-level directory (``MuST/``), run:

.. code-block:: bash

   make <architecture-file-name>
   make install

Example:

.. code-block:: bash

   make linux-intel-nogpu
   make install

.. note::

   - ``make clean``: Removes object, library, and executable files under
     ``lsms`` and ``MST``.
   - ``make distclean``: Removes object, library, executable, and
     ``architecture.h`` files under ``lsms`` and ``MST``, and also deletes
     executables under ``bin/``.

---

Alternative Method
##################

The ``MST`` (under ``MST/``) and ``LSMS/WL-LSMS`` (under ``lsms/``) components
can be built separately.

Executables will be located under ``MST/bin`` and ``lsms/bin``, respectively.
This method requires creating a symbolic link to ``architecture.h``.

**Build MST:**

1. Change directory:

   .. code-block:: bash

      cd MST

2. Set ``SystemName`` in the ``Makefile`` (line 6), or create a symbolic link:

   .. code-block:: bash

      ln -s arch/<architecture_file> architecture.h

3. Compile:

   .. code-block:: bash

      make

**Build LSMS / WL-LSMS:**

1. Change directory:

   .. code-block:: bash

      cd lsms

2. Create symbolic link:

   .. code-block:: bash

      ln -s arch/<architecture_file> architecture.h

3. Compile:

   .. code-block:: bash

      make

---

Notes for Fedora Systems
########################

``MST`` may require the External Data Representation (XDR) library to store
potential and charge density data.

On newer Fedora systems, this library may not be located in standard paths.
Ensure the following directories exist:

- ``/usr/include/tirpc``
- ``/usr/include/tirpc/rpc``

If they are missing, ask your system administrator to install the required
packages, or run the following command (with administrative privileges):

.. code-block:: bash

   sudo dnf install libtirpc libtirpc-devel
