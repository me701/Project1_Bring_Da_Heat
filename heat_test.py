#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:28:53 2019

@author: brandon
"""

import scipy as sp
import numpy as np
from Heat1 import b_n, Heat_1D
import unittest

class TestHeat(unittest.TestCase):
    
    def setUp(self):
        self.L = np.pi
        self.T1 = 0
        self.T2 = 0
        self.n = 7
        self.fun = lambda x:100
        self.c = 1
        self.k = 10
        self.t = 1
    
    def test_Bn(self):
        bfound = sp.integrate.quad(b_n, 0, self.L, args=(self.L, self.T1, self.T2, self.n, self.fun))[0]
        bgiven = 200 / (self.n * np.pi) * (1 - np.cos(self.n * np.pi))
        self.assertAlmostEqual(bfound, bgiven, places=4)
        
    def test_DHeat(self):
        x = 2.0944 #20
        hfound = Heat_1D(self.c, self.L, self.t, self.T1, self.T2, self.fun)[0][0][20]
        i = 1
        expr = 0
        while i<self.k:
            expr = expr + (np.exp(-(i**2 * 0)) * np.sin(i * x) / i)
            i +=2
        hgiven = 400 / np.pi * expr
        self.assertAlmostEqual(hfound, hgiven, places=3)
        
if __name__ == '__main__':
    unittest.main()