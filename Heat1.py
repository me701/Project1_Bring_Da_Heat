#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Me701 Project 1
1D Heat Equation

Sean Cranford
T-Ying Lin
Brandon Yutzy
"""

"""""""""""
  NOTE: If running in an IDE: go to tools, prefrences, iPython console, and then grahpics.
  Find graphics backend and switch to automatic in order to get best results with plots.
"""""""""""

import sympy as sy #run conda install sympy if you do not have enabled
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import pi


def sumr(a):
    """
    a recurrsive way to sum elements with little roundoff or truncation
    """
    if len(a) <= 2:
        return sum(a)
    else:
        return sumr(a[:len(a)//2]) + sumr(a[len(a)//2:])

def Heat_1D_Derivation(c, L, t, T1, T2, fun = lambda x: 100):
    """
    This function will contain all the sympy stuff to explain the 1D heat equation and will probably solve at a single value of x?
    """

def b_n(x, L, T1, T2, n, fun = lambda x: 100 ):
        """
        This function is used to calculate the fourier 
        coefficient b_n of the 1D heat equation:
            b_n = 2/L * integral((f(x) - u_1)*sin(n*pi*x/l)dx) evaluated from 0 to L
        
        Inputs:
            x: the position of interest and the variable of integration
                        
            L: the length of the member
        
            T1: boundary condition: u(0,t) = T1
        
            T2: boundary condition: u(L,t) = T2
                                   
            n: the iteration
            
            fun: boundary condition: u(x,0) = f(x
                
        b_n to be called with Heat_1D and ..... 
        """
        return (2/L)*(fun(x) - (x*(T1 - T2)/L + T1))*sp.sin(x*n*pi/L) 
    
def Heat_1D(c, L, t, T1, T2, fun = lambda x: 100):
    """
    This function gives the solution to the 1D heat equation:
        du(x,t)/dt = c^2 d^2(u(x,t))/dx^2
    The function only considers members with no internal heat sources. 
    
    Inputs:
        c: the diffusivity, k/(s*rho)
       
        L: the length of the member
        
        t: the time of interest
        
        T1: boundary condition: u(0,t) = T1
        
        T2: boundary condition: u(L,t) = T2
        
        fun: boundary condition: u(x,0) = f(x) 
        
            fun defaults to 100, this means that the bar is at a temperature of 100 when
            you start observing the system.
    """
    #plot the member, use annotations to mark the temperatures at ends
    #look into sympy seperation of variables 
    # f(x) is initial temp of bar - steady state solution
    
    
    if L <= 50: # this is actually probably too many points since evaluating and trying to run for x and n values
        point = 100*L
    else:
        point = 5000
    """sympy might actually make b_n solution quicker by computing integral once and then pluging in n for each value"""    
    X = sp.linspace(0,L,point)
    N1 = 0
    N2 = 0
    N3 = 0
    Sol1 = []
    Sol2 = []
    Sol3 = []
    
    for x in X:
        sn1 = [] #this is list that contains the values of the solution at each n up to N1
        sn2 = [] #this is list that contains the values of the solution at each n up to N2
        sn3 = [] #this is list that contains the values of the solution at each n up to N3
        
        for n in range(1,N1+1):
            lam1 = c*n*pi/L
            
            u1_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
            
            b1_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0]
    
            u1_2 = b1_n*sp.exp(-lam1**2*t)*sp.sin(n*pi*x/L)
            
            
            
            u1 = u1_1 + u1_2
            sn1.append(u1)
        
        for n in range(1,N2+1):
            lam2 = c*n*pi/L
     
            u2_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
            b2_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0]
    
            u2_2 = b2_n*sp.exp(-lam2**2*t)*sp.sin(n*pi*x/L)
            
            
            
            u2 = u2_1 + u2_2
            sn2.append(u2)
        
        for n in range(1,N3+1):
            lam3 = c*n*pi/L
     
            u3_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
            b3_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] 
    
            u3_2 = b3_n*sp.exp(-lam3**2*t)*sp.sin(n*pi*x/L)
            
            
            
            u3 = u3_1 + u3_2
            sn3.append(u3)
        
    
        sol1 = sumr(sn1)
        sol2 = sumr(sn2)
        sol3 = sumr(sn3)
        
        Sol1.append(sol1)
        Sol2.append(sol2)
        Sol3.append(sol3)