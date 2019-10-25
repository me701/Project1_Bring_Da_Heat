#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Me701 Project 1
1D Heat Equation

Sean Cranford
T-Ying Lin
Brandon Yutzy
"""

import sympy as sy #run conda install sympy if you do not have enabled
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def sumr(a):
    if len(a) <= 2:
        return sum(a)
    else:
        return sumr(a[:len(a)//2]) + sumr(a[len(a)//2:])

def Heat(K, L, t, T1, T2):
    """
    This function gives the solution to the 1D heat equation. It outputs a graph of the solution at different upperbounds of the series, N.
    """
    
    pass