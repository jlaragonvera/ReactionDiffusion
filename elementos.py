# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 23:04:23 2023

@author: LENOVO
"""

def elemn(x_l,N,K,d_x,L_1):
    x=[]
    for i in range(1,K+1):
        for j in range(N):
            x_j=L_1+(i-1)*d_x+(x_l[j]+1)*d_x/2
            x.append(x_j)
    x.append(L_1+(K)*d_x)
    return x

