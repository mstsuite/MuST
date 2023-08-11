"""
Find the LibXC shared object and imports it as core.
"""

import ctypes
import ctypes.util
import numpy as np
import os

# Attempt to load the compiled C code
core = None
__libxc_path = None

# First check the local folder
try:
    __libxc_path = os.path.abspath(os.path.dirname(__file__))
    core = np.ctypeslib.load_library("libxc", __libxc_path)
except OSError:
    # If no libxc is local, check LD_LIBRARY_PATHS's
    __libxc_path = ctypes.util.find_library("xc")

    # If we still havent found it, give up and throw an error
    if __libxc_path is None:
        raise ImportError(
            "LibXC Shared object not found, searched Python module local directory and library paths"
        )
    
    # Load the C object
    core = ctypes.CDLL(__libxc_path)


def get_core_path():
    """
    Returns the path of the loaded LibXC shared object.
    """

    return __libxc_path
