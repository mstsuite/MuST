"""
Tests the utility LibXC functions.
"""

import pylibxc


def test_xc_version():

    assert tuple(map(int, pylibxc.version.__version__.split('.'))) == pylibxc.util.xc_version()


def test_xc_string():

    assert pylibxc.version.__version__ == pylibxc.util.xc_version_string()


def test_xc_functional_get_number():

    assert 32 == pylibxc.util.xc_functional_get_number("XC_GGA_X_GAM")
    assert 32 == pylibxc.util.xc_functional_get_number("gga_x_gam")
    assert 32 == pylibxc.util.xc_functional_get_number("gga_x_gam".upper())

    assert -1 == pylibxc.util.xc_functional_get_number("nothing")


def test_xc_functional_get_name():
    assert "gga_x_gam" == pylibxc.util.xc_functional_get_name(32)
    assert None is pylibxc.util.xc_functional_get_name(50000)


def test_xc_family_from_id():
    assert (4, 3) == pylibxc.util.xc_family_from_id(72)
    assert (0, 0) == pylibxc.util.xc_family_from_id(500000)


def test_xc_number_of_functionals():
    assert pylibxc.util.xc_number_of_functionals() > 400


def test_xc_available_functional_numbers():
    func_nums = pylibxc.util.xc_available_functional_numbers()

    assert len(func_nums) == pylibxc.util.xc_number_of_functionals()

    # Spot check some functionals
    assert (1 in func_nums)
    assert (76 in func_nums)
    assert (540 in func_nums)

    # Make sure the range on all values is reasonable
    assert all(x > 0 for x in func_nums)
    assert all(x < 2000 or x > 100000 for x in func_nums)


def test_xc_available_functional_names():
    func_names = pylibxc.util.xc_available_functional_names()

    assert len(func_names) == pylibxc.util.xc_number_of_functionals()

    # Spot check some functionals
    assert ("lda_c_vwn" in func_names)
    assert ("gga_c_pbe" in func_names)
    assert ("gga_x_lbm" in func_names)

    # Make sure the range on all values is reasonable
    assert all((("x" in x) or ("c" in x) or ("k" in x)) for x in func_names)
