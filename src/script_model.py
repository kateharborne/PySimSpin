#!/usr/bin/env python
import PySimSpin as PyS

cosmo = PyS.LambdaCDM(H0 = PyS.constant.H0, Om0 = PyS.constant.OMEGA_M, Ode0 = PyS.constant.OMEGA_L)
sami = PyS.Telescope(fov=15, ap_shape="Circular", central_wvl=4800, lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04)
