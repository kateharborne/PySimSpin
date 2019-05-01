import PySimSpin as PyS
import numpy as np
from math import pi as pi
import pytest

cosmo = PyS.LambdaCDM(H0 = PyS.constant.H0, Om0 = PyS.constant.OMEGA_M, Ode0 = PyS.constant.OMEGA_L)

def test_obs_1():
    sami = PyS.Telescope(fov=15, ap_shape="Circular", central_wvl=4800, lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04)
    obs = PyS.Observation(telescope=sami, z = 0.1, inc_deg=90, cosmo=cosmo)
    assert obs.sbin == 30

def test_obs_2():
    sami = PyS.Telescope(fov=15, ap_shape="Circular", central_wvl=4800, lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04)
    obs = PyS.Observation(telescope=sami, z = 0.1, inc_deg=90, cosmo=cosmo)
    assert obs.inc_rad == 0.5 * pi

def test_obs_3():
    with pytest.raises(PyS.ApShapeNotValidException):
        sami = PyS.Telescope(fov=15, ap_shape="Orange", central_wvl=4800, lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04)
        PyS.Observation(telescope=sami, z = 0.1, inc_deg=90, cosmo=cosmo)
    

def test_obs_4():
    sami = PyS.Telescope(fov=15, ap_shape="Circular", central_wvl=4800, lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04)
    simobs = PyS.ObsGalaxy(filename = "../../../data/SimSpin_example.hdf5", telescope=sami, z = 0.1, inc_deg=90, cosmo=cosmo)
    assert len(simobs.r_obs) == 25000 
