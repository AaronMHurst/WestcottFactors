import pytest
import unittest
import numpy as np
import pandas as pd

import westcott
gw = westcott.Westcott()

class gFactorsTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the `gFactors` class of the 
    `westcott_gfactors.py` module."""

    def test_gw_Maxwellian_returns_1_for_30Si_all_temperatures(self):
        expected_gW = 1.000
        endf_e, endf_cs = gw.sigma_ENDF('Si30')
        for T in range(20, 620, 20):
            gW = gw.gw_Maxwellian(T, endf_e, endf_cs)
            assert gW == pytest.approx(expected_gW, abs=0.001)

        # Check thermal T = 293 K
        gW = gw.gw_Maxwellian(293, endf_e, endf_cs)
        assert gW == pytest.approx(expected_gW, abs=0.001)

    def test_gw_IdealGuide_returns_1_for_30Si_all_temperatures(self):
        expected_gW = 1.000
        endf_e, endf_cs = gw.sigma_ENDF('Si30')
        for T in range(20, 620, 20):
            gW = gw.gw_IdealGuide(T, endf_e, endf_cs)
            assert gW == pytest.approx(expected_gW, abs=0.05)

        # Check thermal T = 293 K
        gW = gw.gw_IdealGuide(293, endf_e, endf_cs)
        assert gW == pytest.approx(expected_gW, abs=0.02)
        
            
            
        
