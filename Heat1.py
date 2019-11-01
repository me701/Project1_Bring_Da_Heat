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
    


def Heat_1D_Derivation(c, L, t, T1, T2, fun = lambda x: 100): #maybe don't do this? not sure yet
    """
    This function goes through the derivation of the 1D heat equation
    
    Inputs:
        c: the diffusivity, k/(s*rho)
       
        L: the length of the member
                
        T1: boundary condition: u(0,t) = T1
        
        T2: boundary condition: u(L,t) = T2
        
        fun: boundary condition: u(x,0) = f(x) 
    """
    import sympy as sy
    from sympy.abc import x,t,c
    from sympy import Function, Derivative as D
    sy.init_printing()
    
    u, X, T = map(Function, 'uXT') #creates u(x,t), X(x), T(t)
    
    def initial(L): #this is so we can continuously output the sympy equations within the entire function
       
        print("Here we have the 1D heat equation where u is your temperature as a function of position (x) and time (t): \n")
        eq = sy.Eq(D(u(x,t),t), D(u(x,t), x,2)) 
        return eq
    def seperation(L):
        print('Now through a method know as seperation of variables we say u is a function of x times a funtion of t: \n')
        eq1 = sy.Eq(u(x,t), T(t)*X(x))
        eq1.doit()
        print("Now substituting this expression for the one in our initial equation we get: \n")
        eq2 = sy.Eq(D(T,t)/(c**2*T),D(X,x,2)/X)
        
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
            
            fun: boundary condition: u(x,0) = f(x)
                
        b_n to be called with Heat_1D and ..... 
        """
        return (2/L)*(fun(x) - (x*(T2 - T1)/L + T1))*sp.sin(x*n*pi/L) 
        
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
    
    
    if L <= 10: # this is actually probably too many points since evaluating and trying to run for x and n values
        point = int(10*L)
    else:
        point = 100
        
    if t<= 60:
        times = int(5*t)
    else:
        times = 350
    
    
    T = sp.linspace(0,t,times) # the time increments needed
    X = sp.linspace(0,L,point) # the points on the bar
    N1 = 10
    N2 = 30
    N3 = 80
    Sol1 = []
    Sol2 = []
    Sol3 = []
    
    for ts in T:
        s1 = [T1]
        s2 = [T1]
        s3 = [T1]
        for x in X[1:-1]:
            sn1 = [] #this is list that contains the values of the solution at each n up to N1
            sn2 = [] #this is list that contains the values of the solution at each n up to N2
            sn3 = [] #this is list that contains the values of the solution at each n up to N3
        
            for n in range(1,N1+1):
                lam1 = c*n*pi/L
            
                u1_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
            
                b1_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0]
    
                u1_2 = b1_n*sp.exp(-lam1**2*ts)*sp.sin(n*pi*x/L)
            
            
            
                u1 = u1_1 + u1_2
                sn1.append(u1)
        
            for n in range(1,N2+1):
                lam2 = c*n*pi/L
     
                u2_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
                b2_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0]
    
                u2_2 = b2_n*sp.exp(-lam2**2*ts)*sp.sin(n*pi*x/L)
            
            
            
                u2 = u2_1 + u2_2
                sn2.append(u2)
        
            for n in range(1,N3+1):
                lam3 = c*n*pi/L
     
                u3_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
                b3_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] 
    
                u3_2 = b3_n*sp.exp(-lam3**2*ts)*sp.sin(n*pi*x/L)
            
            
            
                u3 = u3_1 + u3_2
                sn3.append(u3)
        
    
            sol1 = sumr(sn1) #this sums up all the values within the list of the solution evaluated at different n's getting the numerical answer at the position
            sol2 = sumr(sn2)
            sol3 = sumr(sn3)
        
            s1.append(sol1)
            s2.append(sol2)
            s3.append(sol3)
            if x == X[-2]:
                s1.append(T2)
                s2.append(T2)
                s3.append(T2)
        Sol1.append(s1)
        Sol2.append(s2)# fix the appends and then this should be right, you are not taking into account that there are multiple times
        Sol3.append(s3)
        #add 2D plots, need to figure out how we want to do 3D plots, might be best to just have a preset set of data to plot
    #Sols = {'S1': Sol1, 'S2': Sol2, 'S3': Sol3}
    return Sol1,Sol2,Sol3
        
K=Heat_1D(1,pi,1,0,0)
#%% 3D plot example
import numpy as np
def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='plasma')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
ax.view_init(60, 35) #sets view angle

#ax.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#           ) #this code shows how to annotate a graph 