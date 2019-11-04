#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:28:53 2019

@author: brandon
"""

import scipy as sp
import numpy as np
from Heat1 import b_n, sumr, Heat_1D
import unittest

class heat_test(unittest.TestCase):
    
    def Bn(self):
        L = np.pi
        T1 = 100
        T2 = 0
        n = 1
        fun = lambda x: 100
        bfound = sp.integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun))[0]
        bgiven = 200 / (n * np.pi * (1 - np.cos(n * np.pi)))
        self.assertAlmostEqual(bfound, bgiven, places=4)
        
    def DHeat(self):
        c = 1
        L = np.pi
        t = 0
        T1 = 100
        T2 = 0
        x = L
        k = 3
        hfound = Heat_1D(c, L, t, T1, T2, fun = lambda x: 100)
        expr = np.exp(-(2 * k + 1)**2 * t / ((2 * k) + t )) * np.sin((2 * k) + 1) * x
        hgiven = 400 / np.pi * sumr(expr)
        self.assertAlmostEqual(hfound, hgiven, places=4)
        
if __name__ == '__main__':
    unittest.main()