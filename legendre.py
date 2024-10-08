# -*- coding: utf-8 -*-
"""lagrange.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/12zfpaSSLslJLFp8HFedxG9FxoyyCDPX2
"""

def evaluacion_polinomios_Legendre(x,M):
  k=2
  L_n_2=1
  L_n_1=x
  dL_n_2=0
  dL_n_1=1
  for i in range(2,M+1):
    L_n=((2*i-1)/(i))*x*L_n_1-((i-1)/(i))*L_n_2
    dL_n=(2*i-1)*L_n_1+dL_n_2
    if i<M:
        L_n_2=L_n_1
        dL_n_2=dL_n_1
        L_n_1=L_n
        dL_n_1=dL_n
  
  k=M+1
  L_n_11=((2*k-1)/(k))*x*L_n-((k-1)/(k))*L_n_1
  dL_n_11=(2*k-1)*L_n+dL_n_1
  q=L_n_11-L_n_1
  dq=dL_n_11-dL_n_1
  return L_n,q,dq