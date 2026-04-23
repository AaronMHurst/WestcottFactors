import pytest
import unittest
import numpy as np
import pandas as pd

import westcott
gw = westcott.Westcott()

class KinematicsTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the `Kinematics` class of 
    the `westcott_gfactors.py` module."""

    bad_args = [2, 3.14, "x"]

    def test_disiplay_thermal_properties_returns_NoneType(self):
        self.assertIsNone(gw.display_thermal_properties())

    def test_disiplay_thermal_properties_raises_TypeError_passing_args(self):
        for arg in KinematicsTests.bad_args:
            with self.assertRaises(TypeError):
                gw.display_thermal_properties(arg)

    def test_disiplay_constants_returns_NoneType(self):
        self.assertIsNone(gw.display_constants())

    def test_disiplay_constants_raises_TypeError_passing_args(self):
        for arg in KinematicsTests.bad_args:
            with self.assertRaises(TypeError):
                gw.display_constants(arg)

    def test_vel_returns_float(self):
        for u in range(0,10000):
            v = gw.vel(u)
            self.assertIsInstance(v, float)

    def test_vel_raises_TypeError_passing_str(self):
        with self.assertRaises(TypeError):
            gw.vel("5.0")

    def test_phi_v_IdealGuide_returns_np_array(self):
        phi_ideal_293K = gw.phi_v_IdealGuide(293, np.linspace(1,10000,10000))
        self.assertIsInstance(phi_ideal_293K, np.ndarray)

    def test_phi_v_IdealGuide_raises_TypeError_passing_wrong_number_args(self):
        with self.assertRaises(TypeError):
            phi_ideal_293K = gw.phi_v_IdealGuide(293, np.linspace(1,10000,10000), 150)
        with self.assertRaises(TypeError):
            phi_ideal_293K = gw.phi_v_IdealGuide(293)

    def test_phi_v_IdealGuide_raises_TypeError_passing_wrong_type_args(self):
       
        with self.assertRaises(TypeError):
            phi_ideal_293K = gw.phi_v_IdealGuide(293, 800)
        with self.assertRaises(TypeError):
            phi_ideal_293K = gw.phi_v_IdealGuide("293", np.linspace(1,10000,10000))

    
                                              
