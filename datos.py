# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 13:42:37 2023

@author: LENOVO
"""
import math
import random
import numpy as np
from sympy import *
def datos_problema():
    s_0=0.5#Límite izquierdo del parámetro S
    s_1=2*np.pi #Límite derecho del parámetro S
    N=3 #Grado del polinomio de aproximación
    K=100 #Número de divisiones
    d_t=0.000001
    t_f=3000
    d_s=(s_1-s_0)/(K)
    gam=1
    d=1 #Coeficiente de difusión de v
    #d_1=0.026 #Coeficiente de difusión de u
    d_1=0.15
    capt_dat=[500,1000,1500,2000,2500,3000]
    return s_0,s_1,N,K,d_t,t_f,d_s,gam,d,d_1,capt_dat
def func1(u,v):
    #a=0.1
    a=1
    b=-2
    C=1
    H=-2.5
    f=u+a*v-C*u*v-u*v**2
    return f
def func2(u,v):
    #b=0.9
    # a=0.289
    # b=1.49
    a=1
    b=-2
    C=1
    H=-2.5
    g=H*u+b*v+C*u*v+u*v**2
    return g
def phi1_0(x):
    # a=0.1
    # b=0.9
    # a=0.289
    # b=1.49
    #a=0.289
    #b=1.49
    u_0=0
    p=random.randint(-10,10)*0.01+u_0
    #p=random.randint(-10,10)*0.5+1.73447
    return p
def phi2_0(x):
    # a=0.1
    # b=0.9
    # a=0.289
    # b=1.49
    #a=0.289
    #b=1.49
    v_0=0
    q=random.randint(-10,10)*0.01+v_0
    #q=random.randint(-10,10)*0.5+0.4920
    return q
def curve():
    s,t=symbols('s t')
    p=1+0.0006*t
    #p_1=
    #p_2=
    x=p*cos(3*s)/s
    y=p*sin(3*s)/s
    z=0
    #d_3=1 #Coeficiente de difusión de u
    d_3=1 #Coeficiente de difusión de v
    F=1/sqrt(diff(x,s)**2+diff(y,s)**2+diff(z,s)**2)
    F_p=diff(F,s)
    G=F*d_3
    H=F*diff(1/F,t)
    return x,y,z,F,F_p,G,H

def curvevaluation(x,y,z,s_1,t_1):
    s,t=symbols('s t')
    x=lambdify([s,t],x)
    y=lambdify([s,t],y)
    z=lambdify([s,t],z)
    return x(s_1,t_1),y(s_1,t_1),z(s_1,t_1)
def funcevaluation(F,F_p,G,H):
    s,t=symbols('s t')
    F=lambdify([s,t],F)
    F_p=lambdify([s,t],F_p)
    G=lambdify([s,t],G)
    H=lambdify([s,t],H)
    return F,F_p,G,H
    
    
    
    

    
    
    

    
    
        

    
    
    