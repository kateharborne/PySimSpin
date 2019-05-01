import PySimSpin as PyS
import numpy as np
from math import pi as pi
from math import isclose
import pytest

def test_sim_1(): # make sure error triggers when requesting particles that don't exist
    with pytest.raises(PyS.PtypeNotValidException):
        PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [1,2], m2l = [1,1])

def test_sim_2(): # make sure error triggers when no enough m2l ratio values are provided
    with pytest.raises(PyS.M2lNotValidException):
        PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [2,3], m2l = [2])

def test_sim_3(): # ensure that the correct amount of data is stored correctly (in test 3, 4 and 5)
    sim_3 = PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [2,3], m2l = [1,1])
    assert len(sim_3.x) == 25000 # when full particle data is requested

def test_sim_4():
    sim_4 = PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [2], m2l = [1])
    assert len(sim_4.x) == 10000 # when just the bulge is requested

def test_sim_5():
    sim_5 = PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [3], m2l = [1])
    assert len(sim_5.x) == 15000 # when just the disk is requested

def test_sim_6(): # checking that the ptype variable is generated correctly (has correct number of 3s)
    sim_6 = PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5", ptype = [2,3], m2l = [1,1])
    assert np.unique(sim_6.ptype, return_counts=True)[1][0] == 10000 

def test_sim_7(): # checking that the SSP data is stored correctly
    sim_7 = PyS.Simulation(filename = "../../../data/SimSpin_egSSP.hdf5", ssp="../../../data/SSP.hdf5")
    assert len(sim_7.ssp_wave) == 1221

def test_sim_8(): # check that the Lum property is NaN when SSP is provided
    sim_8 = PyS.Simulation(filename = "../../../data/SimSpin_egSSP.hdf5", ssp="../../../data/SSP.hdf5", ptype = [1,4])
    assert np.isnan(sim_8.Lum[190517])

def test_sim_9(): # ensure that defaults work as expected and make m2l == 1 for all ptypes
    sim_9 = PyS.Simulation(filename = "../../../data/SimSpin_example.hdf5")
    assert isclose(sim_9.Lum[0], sim_9.Mass[0] * 1e10, rel_tol=0.00001) # within 0.0001% 

def test_sim_10(): # checking that efaults work as expected and make m2l == NaN when SSP is present
    sim_10 = PyS.Simulation(filename = "../../../data/SimSpin_egSSP.hdf5", ssp="../../../data/SSP.hdf5")
    assert np.isnan(sim_10.Lum[310180])