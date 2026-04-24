import pytest
import unittest
import numpy as np
import pandas as pd

import westcott
gw = westcott.Westcott()

class DataTests(unittest.TestCase):

    __doc__="""Unit tests for data files shipped with project."""

    # ENDF cross-section data
    def test_find_targets_returns_list_of_strings_containing_553_objects(self):
        t = gw.find_targets()
        self.assertIsInstance(t, list)
        assert len(t) == 554
        for s in t:
            self.assertIsInstance(s, str)

    # Resonance data
    def test_find_resonances_returns_list_of_strings_containing_458_objects(self):
        r = gw.find_resonances()
        self.assertIsInstance(r, list)
        assert len(r) == 458
        for s in r:
            self.assertIsInstance(s, str)

    def test_find_resonances_breit_wigner_returns_list_of_strings_containing_386_objects(self):
        bw = gw.find_resonances(res='BW')
        self.assertIsInstance(bw, list)
        assert len(bw) == 386
        for s in bw:
            self.assertIsInstance(s, str)

    def test_find_resonances_reich_moore_returns_list_of_strings_containing_72_objects(self):
        rm = gw.find_resonances(res='RM')
        self.assertIsInstance(rm, list)
        assert len(rm) == 72
        for s in rm:
            self.assertIsInstance(s, str)

    # Flux data
    def test_find_flux_returns_list_of_tuples_containing_4_objects(self):
        f = gw.find_flux()
        self.assertIsInstance(f, list)
        assert len(f) == 4
        for obj in f:
            self.assertIsInstance(obj, tuple)
            assert len(obj) == 2
            self.assertIsInstance(obj[0], int)
            self.assertIsInstance(obj[1], str)

            
        
    
