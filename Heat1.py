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
  NOTE: If running in an IDE: go to tools, preferences, iPython console, and then graphics.
  Find graphics backend and switch to automatic inorder to get best results with plots.
"""""""""""
import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits import mplot3d
from scipy import pi
import matplotlib.animation as animation
#%%

def sumr(a):
    """
    a recurrsive way to sum elements with little roundoff or truncation
    """
    if len(a) <= 2:
        return sum(a)
    else:
        return sumr(a[:len(a)//2]) + sumr(a[len(a)//2:])
#%% 



        
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
#%%     
def Heat_1D(c, L, t, T1, T2, fun = lambda x: 100):
    """
    This function gives the solution to the 1D heat equation:
        du(x,t)/dt = c^2 d^2(u(x,t))/dx^2
    Evaluated at multiple times and positions.
    The function only considers members with no internal heat sources. The bar is
    assumed to be homogeneous, made of a single material.
    
    This function's resolution is not very fine, it you want better resolution use
    single_heat as it will be evaluated at a single position and time.
    
    Inputs:
        c: the diffusivity, k/(s*rho)
       
        L: the length of the member
        
        t: the time of interest
        
        T1: boundary condition: u(0,t) = T1
        
        T2: boundary condition: u(L,t) = T2
        
        fun: boundary condition: u(x,0) = f(x) 
        
            fun defaults to 100, this means that the bar is at a uniform temperature of 100 when
            you start observing the system.
    """
    
    
    
    if L <= 10: 
        point = int(10*L)
    else:
        point = 100
        
    if t<= 60:
        times = int(5*t)
    else:
        times = 350
    
    
    T = sp.linspace(0,t,times) # the time increments needed
    X = sp.linspace(0,L,point) # the points on the bar
    N1 = 10 #these are different upperbounds to the series solution to be used as comparison
    N2 = 30
    N3 = 80
    Sol1 = [] # in the end these will each contain a list for every time value of the solution at the different points on the bar
    Sol2 = []
    Sol3 = []
    
    for ts in T: #this loops through the different time increments
        s1 = [T1] #the ends are defined by the boundary conditions so this saves needless calculations
        s2 = [T1]
        s3 = [T1]
        for x in X[1:-1]: # this loops through the defined positions on the bar
            sn1 = [] #this is list that contains the values of the solution at each n up to N1
            sn2 = [] #this is list that contains the values of the solution at each n up to N2
            sn3 = [] #this is list that contains the values of the solution at each n up to N3
        
            for n in range(1,N1+1):# these itterate the n values within the series
                lam1 = c*n*pi/L
            
                u1_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
            
                b1_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] #this finds the fourier coefficient
    
                u1_2 = b1_n*sp.exp(-lam1**2*ts)*sp.sin(n*pi*x/L) #this is the time dependent solution
            
            
            
                u1 = u1_1 + u1_2
                sn1.append(u1) # i think if i had been a little more clever and 
                #had a bit more time i could have made a subfunction so i would
                #not be doing the same calculation 3 times. Use an appropriate bound conditions on for loop as well
        
            for n in range(1,N2+1):# these itterate the n values within the series
                lam2 = c*n*pi/L
     
                u2_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
                b2_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] #this finds the fourier coefficient
    
                u2_2 = b2_n*sp.exp(-lam2**2*ts)*sp.sin(n*pi*x/L) #this is the time dependent solution
            
            
            
                u2 = u2_1 + u2_2
                sn2.append(u2)
        
            for n in range(1,N3+1):# these itterate the n values within the series
                lam3 = c*n*pi/L
     
                u3_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
        
                b3_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] #this finds the fourier coefficient
    
                u3_2 = b3_n*sp.exp(-lam3**2*ts)*sp.sin(n*pi*x/L) #this is the time dependent solution
            
            
            
                u3 = u3_1 + u3_2
                sn3.append(u3)
        
    
            sol1 = sumr(sn1) #this sums up all the values within the list of the solution evaluated at different n's getting the numerical answer at the position
            sol2 = sumr(sn2)
            sol3 = sumr(sn3)
        
            s1.append(sol1)
            s2.append(sol2)
            s3.append(sol3)
            if x == X[-2]:
                s1.append(T2) #the ends are defined by the boundary conditions so this saves needless calculations
                s2.append(T2)
                s3.append(T2)
        Sol1.append(s1)
        Sol2.append(s2)
        Sol3.append(s3)
        
    
    return Sol1,Sol2,Sol3,X,T
        
#%%
def Heat_plot2D(H, x, t, keep=False, show=False):
    """
    This function makes 2D plots of the solution of the heat equation evaluated
    some times t. The idea is that this will be used in with the output of Heat_1D.
    
    Inputs:
        H: the solution to the heat equation at the time specified and the x 
        values inputted, you will want the solutions evaluated at the different
        uppperbounds N of a time t. The H input needs to be in the form:
            [ [N1 solutions at your times, i.e. [t0],[t1],[t2]],[same for N2], [same for N3] ]
        
        x: the positions that the solution was evaluated at. 
        
        t: the times you want to visualize the solution at, you probably do not
        want a lot of time values at once.
        
        keep: set True if you would like to save figures
        
        show: set True if you would like a window to pop up, really only has effect in commandline
        
    Note: 
        The output of Heat_1D contains the x and t values used in the calculations,
        you will definitely want all of the x values. The output also has 3 lists
        of the solutions, 1 for each different upperbound. Usually an easy way
        to define variables is:
            Sol = Heat_1D output
            
            H = Sol[:3]
            
            x = Sol[-2]
            
            t = Sol[-1]
            
        However, for single plots you will need to be a little more vigilant on how you define H.
            ex:
                H = [Heat_1D out[0][2],Heat_1D out[1][2],Heat_1D out[2][2]]
                
            is how you would define H for a given t value.
    """
    
  
        
    if (isinstance(t,sp.ndarray) and len(t)==1) or type(t)==int or type(t)==float or type(t)==sp.float64:
        fig = plt.figure()
        plt.clf()
        title = 'Temperature at Time: ' + str(t)
        plt.plot(x,H[0],'r--',label = 'N = 10')
        plt.plot(x,H[1],'b:',label = 'N = 30')
        plt.plot(x,H[2],'g-',label = 'N = 80')
        plt.xlabel('Position')
        plt.ylabel('Temperature')
        plt.legend(loc=0)
        plt.title(title)
        if keep :
            plt.savefig('Solution_plot', format ='png')
        if show :
            plt.show()
        
    
    elif len(t) > 1: #sorry for very repetitive code Richard/Roberts, did it this way and then didn't want  to make function...
        for i in range(len(t)):
            fig = plt.figure(i)
            plt.clf()
            title = 'Temperature at Time: ' + str(t[i])
            plt.plot(x,H[0][i],'r--',label = 'N = 10')
            plt.plot(x,H[1][i],'b:',label = 'N = 30')
            plt.plot(x,H[2][i],'g-',label = 'N = 80')
            plt.xlabel('Position')
            plt.ylabel('Temperature')
            plt.legend(loc=0)
            plt.title(title)
            if keep:
                fn = 'Solution_plot_' + str(i)
                fn = 'Solution_plot_{}.png'.format(i)
                plt.savefig(fn)
            if show :
                plt.show()

    
        
    else:
        fig = plt.figure()
        plt.clf()
        title = 'Temperature at Time: ' + str(t)
        plt.plot(x,H[0],'r--',label = 'N = 10')
        plt.plot(x,H[1],'b:',label = 'N = 30')
        plt.plot(x,H[2],'g-',label = 'N = 80')
        plt.xlabel('Position')
        plt.ylabel('Temperature')
        plt.legend(loc=0)
        plt.title(title)
        if keep :
            plt.savefig('Solution_plot',format='png')
        if show :
            plt.show()
#%%
def Generate_bar(L,T1,T2,fun=lambda x: 100, keep=False):
    """
    Function to generate informational graphic depicting functions used, 
    starting temperatures at either end, and length of the bar. 
    """
    val = input("Please input the function you have used: ")   #To produce an easy way of identifying the function used.  
    #Configuring table to be visually appealing. 
    fig, ax = plt.subplots(1)
    ax.axis('off')
    rect = patches.Rectangle((5,5), 100, 2.5, linewidth=1, fill = False)
    ax.add_patch(rect)
    ax.margins(.2,5)
    #Adding labels for temperature, length, and function. 
    ax.annotate('T1 = {0:.3f}'.format(T1), xy=(5, 5), xytext=(3, 1.5),
            arrowprops=dict(facecolor='black', shrink=0.05))
    ax.annotate('T2 = {0:.3f}'.format(T2), xy=(105,5), xytext=(107,1.5), 
                arrowprops=dict(facecolor='black', shrink=0.05))
    ax.annotate('Length = {0:.3f}'.format(L), xy=(105, 9), xytext=(-7, 8.5),
                arrowprops=dict(arrowstyle='<->'))
    ax.annotate('u(x,0)={}'.format(str(val)), xy=(50, 7.5), xytext=(50, 10),
                arrowprops=dict(facecolor='black', shrink=0.05)) 
    if keep:
        plt.savefig("Bar",format='png')
    plt.show()
#%%
def Heat_timeplot(h,x, keep=False,show=True):
    """
    Creates an animation of the distrubution of the temperature over the length
    that changes over time    
    
    Inputs: 
        h: one of the solutions outputted from Heat_1D
        
        x: the poistional output of Heat_1D
    """
    
    data = h
    
    fig = plt.figure()

    im = []
    for i in range(len(data)):
        plot, = plt.plot(x, data[i],'b-')
        im.append((plot,))
        plt.xlabel('Position')
        plt.ylabel('Temperature')
        plt.title('Animated Temperature')
    # Each plot is now a frame in the image
    ani = animation.ArtistAnimation(fig, im, interval=200, blit=True)
    # Save the figure 
    if keep:
        ani.save('Temperature.mp4', writer='ffmpeg')
    # Show the outcome
    if show:
        plt.show() #should work for commandline

#%%
def Heat_plot3D(H, x, t, keep=False, show=False, colormap='plasma'):
    """
    This function creates a 3D plot of the 1D heat solution.
    
    Inputs:
        H: the solution, choose one of the three outputted solutions. This is a
        choice between which upperbound you would like to see. WARNING! only 
        give one  of the three solutions
        
        x: the positional output of Heat_1D
        
        t: the time interval you wish to see graphed. WARNING! time cannot exceed
        that of what was inputted to Heat_1D and the step size value must match
        the output as well. i.e. you can enter a t value less than what was inputted
        but it must be a t value that Heat_1D was evaluated at.
        
        ex: t = outputted t [:-1]
    """
    X1, Y1 = np.meshgrid(x,t)
    z = []
    for i in range(len(t)):
        z.append(H[i]) #this is allowed and done because you are giving one of the solutions and not all three
    
    Z1 = np.array(z)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(X1, Y1, Z1, 50, cmap=colormap)
    ax.set_xlabel('Position')
    ax.set_ylabel('Time')
    ax.set_zlabel('Temperature');
    ax.view_init(60, 35) #sets view angle
    
    if keep :
        plt.savefig('3D_Heat_Plot', format ='png') #open this using the folder GUI not gedit
    
    if show :
        plt.show()
#%% 
def single_heat(c, L, x, t, T1, T2, fun = lambda x: 100):
    """
    This function gives the solution to the 1D heat equation:
        du(x,t)/dt = c^2 d^2(u(x,t))/dx^2
    Evaluated at a single time, t, and position, x.
    The function only considers members with no internal heat sources. The bar is
    assumed to be homogeneous, made of a single material.
    
    Inputs:
        c: the diffusivity, k/(s*rho)
       
        L: the length of the member
        
        x: the position of interest
        
        t: the time of interest
        
        T1: boundary condition: u(0,t) = T1
        
        T2: boundary condition: u(L,t) = T2
        
        fun: boundary condition: u(x,0) = f(x) 
        
            fun defaults to 100, this means that the bar is at a uniform temperature of 100 when
            you start observing the system.
    """

    N = 115
    sn = [] #this list will contain the solution at each itteration, n
    if x == 0:
        sol = T1
        
    elif x == L:
        sol = T2
        
    else:
        for n in range(1,N+1):# these itterate the n values within the series
            lam1 = c*n*pi/L
            
            u1_1 = (T2-T1)*x/L + T1 #this is steady state, time independent solution.
            
            b1_n = integrate.quad(b_n,0, L, args=(L,T1,T2,n,fun,) )[0] #this finds the fourier coefficient
    
            u1_2 = b1_n*sp.exp(-lam1**2*t)*sp.sin(n*pi*x/L) #this is the time dependent solution
                   
            u1 = u1_1 + u1_2
            sn.append(u1)
        
    
        sol = sumr(sn) #this sums up all the solutions at the values n getting the solution at x,t
            
    return sol,x,t

#%% An example output
K=Heat_1D(1,pi,1,0,0)
H = K[:3]
x = K[-2]
t = K[-1]

#%% Plot examples
Heat_plot2D(H,x,t, show=True)
h = H[-1]
#%%
Heat_plot3D(h,x,t, show=True)
#%%
Heat_timeplot(h,x, keep=True,show=True)
#%%
Generate_bar(pi,0,0,keep=True)


