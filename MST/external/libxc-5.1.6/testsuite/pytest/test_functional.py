"""
Tests the LibXCFunctional class.
"""

import pytest
import numpy as np
import ctypes
import pylibxc

compute_test_dim = 5

np.random.seed(0)


def _dict_array_comp(test, ref, keys):
    for key in keys:
        tmp = np.testing.assert_allclose(test[key], ref[key])
        # print(key, np.allclose(test[key], ref[key]), np.linalg.norm(test[key] - ref[key]))
    return True


_size_tuples = {"unpolarized": (1, 1, 1, 1), "polarized": (2, 3, 2, 2)}


def test_libxc_functional_build():

    pylibxc.LibXCFunctional(1, 1)
    pylibxc.LibXCFunctional(1, 2)

    pylibxc.LibXCFunctional("XC_LDA_C_VWN", "polarized")
    pylibxc.LibXCFunctional("lda_c_vwn", "unpolarized")

    # Check functional edge cases
    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("something", 1)

    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional(5000, 1)

    # Check spin edge cases
    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("lda_c_vwn", 10)

    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("lda_c_vwn", "something")


def test_libxc_functional_info():

    func = pylibxc.LibXCFunctional(1, 1)
    assert func.get_number() == 1
    assert func.get_kind() == 0
    assert func.get_name() == "Slater exchange"
    assert func.get_family() == 1
    # XC_FLAGS_3D + [HAVE_EXC | HAVE_VXC | HAVE_FXC | HAVE_KXC]
    assert func.get_flags() in [129, 131, 135, 143, 159]
    assert len(func.get_bibtex()) == 2
    assert len(func.get_references()) == 2
    assert len(func.get_doi()) == 2

    func = pylibxc.LibXCFunctional("XC_HYB_MGGA_XC_WB97M_V", 1)
    assert func.get_number() == 531
    assert func.get_kind() == 2
    assert func.get_name() == "wB97M-V exchange-correlation functional"
    assert func.get_family() == 64
    # XC_FLAGS_3D + XC_FLAGS_HYB_CAM + XC_FLAGS_VV10 + [HAVE_EXC | HAVE_VXC | HAVE_FXC | HAVE_KXC]
    assert func.get_flags() in [1409, 1411, 1415, 1423, 1439]
    assert len(func.get_bibtex()) == 1
    assert len(func.get_references()) == 1
    assert len(func.get_doi()) == 1

def test_ext_params():

    func = pylibxc.LibXCFunctional(1, 1)
    assert 0 == len(func.get_ext_param_descriptions())
    assert 0 == len(func.get_ext_param_default_values())
    func.set_dens_threshold(1.e-16)
    func.set_dens_threshold(5)
    with pytest.raises(ValueError):
        func.set_ext_params([])

    func = pylibxc.LibXCFunctional("XC_HYB_GGA_XC_HSE06", 1)
    assert 3 == len(func.get_ext_param_descriptions())
    assert 3 == len(func.get_ext_param_default_values())
    assert all("param" in x for x in func.get_ext_param_descriptions())
    func.set_dens_threshold(1.e-16)
    func.set_dens_threshold(5)

    # Segfaults, need to check it out
    func.set_ext_params([5, 3, 3])

    with pytest.raises(ValueError):
        func.set_ext_params([5, 3])

    with pytest.raises(ValueError):
        func.set_dens_threshold(-1)


@pytest.mark.parametrize("polar", [("unpolarized"), ("polarized")])
def test_lda_compute(polar):

    # Build input
    ndim = _size_tuples[polar]
    inp = {}
    inp["rho"] = np.random.random((compute_test_dim * ndim[0]))

    func = pylibxc.LibXCFunctional("lda_c_vwn", polar)

    # Check capabilities
    do_l = func._have_lxc
    do_k = func._have_kxc
    do_f = func._have_fxc
    do_v = func._have_vxc
    do_e = func._have_exc

    # Compute
    ret_full = func.compute(inp, do_exc=do_e, do_vxc=do_v, do_fxc=do_f, do_kxc=do_k, do_lxc=do_l)
    ret_ev = func.compute(inp, do_exc=do_e, do_vxc=do_v, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_e = func.compute(inp, do_exc=do_e, do_vxc=False, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_v = func.compute(inp, do_exc=False, do_vxc=do_v, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_f = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=do_f, do_kxc=False, do_lxc=False)
    ret_k = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=do_k, do_lxc=False)
    ret_l = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=False, do_lxc=do_l)

    # Test consistency
    if do_e:
        assert ret_full["zk"].size == compute_test_dim
    if do_v:
        assert ret_full["vrho"].size == compute_test_dim * ndim[0]

    if do_e:
        assert np.allclose(ret_full["zk"], ret_ev["zk"])
    if do_v:
        assert np.allclose(ret_full["vrho"], ret_ev["vrho"])

    if do_e:
        assert np.allclose(ret_full["zk"], ret_e["zk"])
    if do_v:
        assert np.allclose(ret_full["vrho"], ret_v["vrho"])
    if do_f:
        assert np.allclose(ret_full["v2rho2"], ret_f["v2rho2"])
    if do_k:
        assert np.allclose(ret_full["v3rho3"], ret_k["v3rho3"])
    if do_l:
        assert np.allclose(ret_full["v4rho4"], ret_l["v4rho4"])


@pytest.mark.parametrize("polar", [("unpolarized"), ("polarized")])
def test_gga_compute(polar):

    # Build input
    ndim = _size_tuples[polar]
    inp = {}
    inp["rho"] = np.random.random((compute_test_dim * ndim[0]))
    inp["sigma"] = np.random.random((compute_test_dim * ndim[1]))

    # Compute
    func = pylibxc.LibXCFunctional("gga_c_pbe", polar)

    # Check capabilities
    do_l = func._have_lxc
    do_k = func._have_kxc
    do_f = func._have_fxc
    do_v = func._have_vxc
    do_e = func._have_exc

    # Compute
    ret_full = func.compute(inp, do_exc=do_e, do_vxc=do_v, do_fxc=do_f, do_kxc=do_k, do_lxc=do_l)
    ret_ev = func.compute(inp, do_exc=do_e, do_vxc=do_v, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_e = func.compute(inp, do_exc=do_e, do_vxc=False, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_v = func.compute(inp, do_exc=False, do_vxc=do_v, do_fxc=False, do_kxc=False, do_lxc=False)
    ret_f = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=do_f, do_kxc=False, do_lxc=False)
    ret_k = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=do_k, do_lxc=False)
    ret_l = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=False, do_lxc=do_l)

    # Test consistency
    if do_e:
        assert ret_full["zk"].size == compute_test_dim
    if do_v:
        assert ret_full["vrho"].size == compute_test_dim * ndim[0]
        assert ret_full["vsigma"].size == compute_test_dim * ndim[1]

    if do_e and do_v:
        assert _dict_array_comp(ret_full, ret_ev, ["zk", "vrho", "vsigma"])

    if do_e:
        assert _dict_array_comp(ret_full, ret_e, ["zk"])
    if do_v:
        assert _dict_array_comp(ret_full, ret_v, ["vrho", "vsigma"])
    if do_f:
        assert _dict_array_comp(ret_full, ret_f, ["v2rho2", "v2rhosigma", "v2sigma2"])
    if do_k:
        assert _dict_array_comp(ret_full, ret_k, ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"])
    if do_l:
        assert _dict_array_comp(ret_full, ret_l, ["v4rho4", "v4rho3sigma", "v4rho2sigma2", "v4rhosigma3", "v4sigma4"])


@pytest.mark.parametrize("polar", [("unpolarized"), ("polarized")])
def test_mgga_compute(polar):

    # Build input
    ndim = _size_tuples[polar]

    inp = {}
    inp["rho"] = np.random.random((compute_test_dim * ndim[0]))
    inp["sigma"] = np.random.random((compute_test_dim * ndim[1]))
    inp["tau"] = np.random.random((compute_test_dim * ndim[3]))

    # Compute
    func = pylibxc.LibXCFunctional("mgga_c_tpss", polar)

    # Test consistency
    ret_ev = func.compute(inp, do_exc=True, do_vxc=True)
    ret_e = func.compute(inp, do_exc=True, do_vxc=False)
    ret_v = func.compute(inp, do_exc=False, do_vxc=True)

    assert ret_ev["zk"].size == compute_test_dim
    assert ret_ev["vrho"].size == compute_test_dim * ndim[0]
    assert ret_ev["vsigma"].size == compute_test_dim * ndim[1]

    assert _dict_array_comp(ret_ev, ret_e, ["zk"])
    assert _dict_array_comp(ret_ev, ret_v, ["vrho", "vsigma", "vtau"])


@pytest.mark.parametrize("polar", [("unpolarized"), ("polarized")])
def test_mgga_lapl_compute(polar):

    # Build input
    ndim = _size_tuples[polar]

    inp = {}
    inp["rho"] = np.random.random((compute_test_dim * ndim[0]))
    inp["sigma"] = np.random.random((compute_test_dim * ndim[1]))
    inp["tau"] = np.random.random((compute_test_dim * ndim[3]))
    inp["lapl"] = np.random.random((compute_test_dim * ndim[3]))

    # Compute
    func = pylibxc.LibXCFunctional("mgga_x_br89", polar)

    # Test consistency
    ret_ev = func.compute(inp, do_exc=True, do_vxc=True)
    ret_e = func.compute(inp, do_exc=True, do_vxc=False)
    ret_v = func.compute(inp, do_exc=False, do_vxc=True)

    assert ret_ev["zk"].size == compute_test_dim
    assert ret_ev["vrho"].size == compute_test_dim * ndim[0]
    assert ret_ev["vsigma"].size == compute_test_dim * ndim[1]

    assert _dict_array_comp(ret_ev, ret_e, ["zk"])
    assert _dict_array_comp(ret_ev, ret_v, ["vrho", "vsigma", "vtau"])

    # Test exception
    del inp["lapl"]
    with pytest.raises(KeyError):
        func.compute(inp, do_exc=True, do_vxc=True)


@pytest.mark.parametrize("polar", [("unpolarized"), ("polarized")])
def test_deriv_flags(polar):
    func = pylibxc.LibXCFunctional("mgga_c_tpss", polar)

    ndim = _size_tuples[polar]
    inp = {}
    inp["rho"] = np.random.random((compute_test_dim * ndim[0]))
    inp["sigma"] = np.random.random((compute_test_dim * ndim[1]))
    inp["tau"] = np.random.random((compute_test_dim * ndim[3]))
    inp["lapl"] = np.random.random((compute_test_dim * ndim[3]))

    # disabled as mgga_c_tpss now has fxc
    # please fix this better!
    #with pytest.raises(ValueError):
    #    func.compute(inp, do_fxc=True)

def test_hyb_getters():

    func = pylibxc.LibXCFunctional("hyb_gga_xc_b3lyp", "unpolarized")
    assert pytest.approx(0.2) == func.get_hyb_exx_coef()

    with pytest.raises(ValueError):
        func.get_cam_coef()
    with pytest.raises(ValueError):
        func.get_vv10_coef()


def test_cam_getters():

    func = pylibxc.LibXCFunctional("hyb_gga_xc_cam_b3lyp", "unpolarized")

    with pytest.raises(ValueError):
        assert pytest.approx(0.65) == func.get_hyb_exx_coef()

    omega, alpha, beta = func.get_cam_coef()
    assert pytest.approx(0.33) == omega
    assert pytest.approx(0.65) == alpha
    assert pytest.approx(-0.46) == beta

    with pytest.raises(ValueError):
        func.get_vv10_coef()


def test_vv10_getters():

    func = pylibxc.LibXCFunctional("gga_xc_vv10", "unpolarized")

    b, C = func.get_vv10_coef()
    assert pytest.approx(5.9) == b
    assert pytest.approx(0.0093) == C

    with pytest.raises(ValueError):
        func.get_cam_coef()

    with pytest.raises(ValueError):
        func.get_hyb_exx_coef()
