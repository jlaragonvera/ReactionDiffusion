# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 16:50:29 2023

@author: LENOVO
"""

def Mat_G(nod,M,d): #Sólo requiere la lista de nodos y los valores del coeficiente de difusión en cada nodo de c/elemento
    import numpy as np
    from lagrange import lag
    from nodos_pesos import nodos_y_pesos_cuadratura 
    x,w=nodos_y_pesos_cuadratura(M)
    G=np.zeros((M+1)**2).reshape(M+1,M+1)
    for i in range(M+1):
        for j in range(M+1):
            for k in range(M+1):
                D,l=lag(nod[k],nod,M)
                if i<=j:
                    G_1=w[k]*d[k]*D[k][j]*D[k][i]
                    G[i][j]=G[i][j]+G_1
            if i>j:
                G[i][j]=G[j][i]
    return G
     
                
        