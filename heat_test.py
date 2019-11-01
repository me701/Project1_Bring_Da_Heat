#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:28:53 2019

@author: brandon
"""

import sympy as sy
import numpy as np
from Heat1 import Heat_1D
import unittest

class heat_test(unittest.TestCase):
    
    def Bn(self):
        L = np.pi()
        x = L
        T1 = 100
        T2 = 0
        n = 1
        bfound = Heat_1D.b_n(x, L, T1, T2, n, fun = lambda x: 100 ) 
        bgiven = 200 / (n * np.pi()) * (1 - np.cos(n * np.pi()))
        return self.assertAlmostEqual(bfound, bgiven, places=4)
        
    def DHeat(self):
        c = 1
        L = np.pi()
        t = 0
        T1 = 100
        T2 = 0
        x = L
        k = sy.symbol('k')
        hfound = Heat_1D(c, L, t, T1, T2, fun = lambda x: 100)
        expr = np.exp(-(2 * k + 1)**2 * t / ((2 * k) + t )) * np.sin((2 * k) + 1) * x
        hgiven = 400 / np.pi() * sy.limit(expr, k, sy.oo)
        return self.assertAlmostEqual(hfound, hgiven, places=4)
        
if __name__ == '__main__':
    unittest.main()