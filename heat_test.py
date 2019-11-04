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
        self.k = 3
        self.t = 1
    
    def test_Bn(self):
        bfound = sp.integrate.quad(b_n, 0, self.L, args=(self.L, self.T1, self.T2, self.n, self.fun))[0]
        bgiven = 200 / (self.n * np.pi) * (1 - np.cos(self.n * np.pi))
        self.assertAlmostEqual(bfound, bgiven, places=4)
        
    def test_DHeat(self):
        xg = 2.27586 #20
        hfound1 = Heat_1D(self.c, self.L, self.t, self.T1, self.T2, self.fun)[]
        expr = np.exp(-(2 * self.k + 1)**2 * self.t / ((2 * self.k) + self.t )) * np.sin((2 * self.k) + 1) * xg
        hgiven = 400 / np.pi * expr
        print(hfound1)
        print(hgiven)
        self.assertAlmostEqual(hfound1, hgiven, places=4)
        
if __name__ == '__main__':
    unittest.main()