# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 13:54:22 2023

@author: LENOVO
"""

def lag(x,nod,M):#Ingresa Nodos(nod) y grado del polinomio de interpolación (M) y el punto de evaluación (x)
 import numpy as np
 w=[]
 l=[]
 D=[]
 for i in range(M+1):
     w_j=1
     for j in range(M+1):
         if nod[j]!=nod[i]:
             w_j=w_j*(nod[i]-nod[j])
     w.append(1/w_j)
 for i in range(M+1):
     l_j=0
     for j in range(M+1):
         if i!=j:
             l_j=l_j+(w[j]/(x-nod[j]))
     l_q=w[i]/((x-nod[i])*l_j)
     l.append(l_q)
 for i in range(M+1):
     for j in range(M+1):
         if i!=j:
             D_1=(w[j]/w[i])*(1/(nod[i]-nod[j]))
             D.append(D_1)
         else:
             D_1=0
             D.append(D_1)
     D[i*(M+1)+i]=-sum(D[i*(M)+i:i*(M+1)+M+1])
 D_v=np.reshape(D,(M+1,M+1))
 return D_v,l
             
             
             
     
             
         
                  
    