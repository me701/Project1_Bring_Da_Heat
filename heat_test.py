#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:28:53 2019

@author: brandon
"""

import scipy as sp
import numpy as np
import Heat1
import unittest

class heat_test(unittest.TestCase):
    
    def Bn(self):
        c = 1
        L = np.pi()
        t = 0
        T1 = 100
        T2 = 0
        n = 1  "What is n value?"
        bfound = b_n(x, L, T1, T2, n, fun = lambda x: 100 ) 
        bgiven = 200 / (n * np.pi()) * (1 - np.cos(n * np.pi()))
        self.assertAlmostEqual(bfound, bgiven, places=4)
        
    def 1DHeat(self):
        c = 1
        L = np.pi()
        t = 0
        T1 = 100
        T2 = 0
        hfound = Heat_1D(c, L, t, T1, T2, fun = lambda x: 100)
#        hgiven = 400 / np.pi() sum np.exp(-(2 * k + 1)**2 * t / ((2 * k) + t )) * np.sin((2 * k) + 1) * x