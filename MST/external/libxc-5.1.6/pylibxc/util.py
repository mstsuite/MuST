"""
Binds the LibXC utility functions.
"""

import ctypes
import numpy as np

from .core import core
from . import flags

### Set required ctypes bindings

core.xc_version.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
core.xc_version.restype = None

core.xc_version_string.restype = ctypes.c_char_p
core.xc_reference.restype = ctypes.c_char_p
core.xc_reference_doi.restype = ctypes.c_char_p

core.xc_functional_get_number.argtypes = (ctypes.c_char_p, )
core.xc_functional_get_number.restype = ctypes.c_int

core.xc_functional_get_name.argtypes = (ctypes.c_int, )
core.xc_functional_get_name.restype = ctypes.c_char_p

core.xc_family_from_id.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
core.xc_family_from_id.restype = ctypes.c_int

core.xc_available_functional_numbers.argtypes = (np.ctypeslib.ndpointer(dtype=np.intc, ndim=1, flags=("W", "C", "A")), )
core.xc_available_functional_numbers_by_name.argtypes = (np.ctypeslib.ndpointer(dtype=np.intc, ndim=1, flags=("W", "C", "A")), )

core.xc_available_functional_names.argtypes = (ctypes.POINTER(ctypes.c_char_p), )

### Build wrapper functions


def xc_version():
    """
    Returns the current LibXC version as semantic versioning tuple.

    Returns
    -------
    version : tuple
        The (major, minor, patch) version of the linked LibXC shared object.

    Examples
    --------
    >>> pylibxc.util.xc_version()
    (4, 0, 1)

    """
    major = ctypes.c_int()
    minor = ctypes.c_int()
    patch = ctypes.c_int()
    core.xc_version(major, minor, patch)
    return (major.value, minor.value, patch.value)


def xc_version_string():
    """
    Returns the current LibXC version as a string.

    Returns
    -------
    version : str
        The string representation of the current LibXC version.

    Examples
    --------
    >>> pylibxc.util.xc_version_string()
    "4.0.1"

    """
    return core.xc_version_string().decode("UTF-8")


def xc_reference():
    """
    Returns the reference for the current LibXC version as a string.

    Returns
    -------
    reference : str
        The string representation of the literature reference for LibXC.

    Examples
    --------
    >>> pylibxc.util.xc_reference()
    "S. Lehtola, C. Steigemann, M. J. Oliveira, and M. A. Marques, SoftwareX 7, 1 (2018)"

    """
    return core.xc_reference().decode("UTF-8")


def xc_reference_doi():
    """
    Returns the doi of the reference for the current LibXC version as a string.

    Returns
    -------
    doi : str
        The string representation of the doi of the literature reference for LibXC.

    Examples
    --------
    >>> pylibxc.util.xc_reference_doi()
    "10.1016/j.softx.2017.11.002"

    """
    return core.xc_reference_doi().decode("UTF-8")


def xc_functional_get_number(name):
    """
    Returns the functional ID from a given string.

    Parameters
    ----------
    name : str
        The functional name to find the ID of.

    Returns
    -------
    id : int
        The ID of the requested functional.

    Examples
    --------
    >>> pylibxc.util.xc_functional_get_number("XC_GGA_X_GAM")
    32
    """
    if not isinstance(name, str):
        raise TypeError("xc_functional_get_number: name must be a string. Got {}".format(name))
    return core.xc_functional_get_number(ctypes.c_char_p(name.encode()))


def xc_functional_get_name(func_id):
    """
    Returns the functional name from a ID.

    Parameters
    ----------
    func_id : int
        The LibXC functional ID

    Returns
    -------
    functional_name : str
        The functional_name of the requested functional.

    Examples
    --------
    >>> pylibxc.util.xc_functional_get_name(32)
    "gga_x_gam"
    """
    if not isinstance(func_id, (int, np.integer)):
        raise TypeError("xc_functional_get_name: func_id must be an int. Got {}".format(func_id))
    ret = core.xc_functional_get_name(func_id)
    if ret is not None:
        ret = ret.decode("UTF-8")
    return ret

def xc_family_from_id(func_id):
    """
    Returns the family class and family number (?).

    Parameters
    ----------
    func_id : int
        The LibXC functional ID

    Returns
    -------
    functional_family : int
        The family ID.
    functional_number : int
        The functional number within a family.

    Examples
    --------
    >>> pylibxc.util.xc_family_from_id(72)
    (4, 4)

    """
    if not isinstance(func_id, (int, np.integer)):
        raise TypeError("xc_family_from_id: func_id must be an int. Got {}".format(func_id))
    family = ctypes.c_int()
    number = ctypes.c_int()
    core.xc_family_from_id(func_id, ctypes.pointer(family), ctypes.pointer(number))

    return (family.value, number.value)


def xc_number_of_functionals():
    """
    Returns the totaly number of XC functionals available in LibXC.

    Returns
    -------
    number_of_functinals : int
        The total number of functionals available in LibXC.

    Examples
    --------
    >>> pylibxc.util.xc_number_of_functionals()
    447
    """

    return core.xc_number_of_functionals()


def xc_available_functional_numbers():
    """
    Returns a list of all available XC functional IDs


    Returns
    -------
    functional_ids : list of ints
        All available LibXC functional IDs.

    Examples
    --------
    >>> xc_func_list = pylibxc.util.xc_available_functional_numbers()
    np.array([1, 2, ..., 568, 569])
    """

    nfunc = xc_number_of_functionals()

    ret = np.zeros(nfunc, dtype=np.intc)
    core.xc_available_functional_numbers(ret)
    return ret


def xc_available_functional_names():
    """
    Returns a list of all available XC functional names

    Returns
    -------
    functional_names : list of strs
        All available LibXC functional names.

    Examples
    --------
    >>> xc_func_list = pylibxc.util.xc_available_functional_names()
    ['lda_x', 'lda_c_wigner', ..., 'hyb_mgga_x_revscan0', 'hyb_mgga_xc_b98']
    """

    # I give up trying to get char** working, someone else can pick it up.
    nfunc = xc_number_of_functionals()
    func_ids = np.zeros(nfunc, dtype=np.intc)
    core.xc_available_functional_numbers_by_name(func_ids)
    available_names = [xc_functional_get_name(x) for x in func_ids]
    return available_names
