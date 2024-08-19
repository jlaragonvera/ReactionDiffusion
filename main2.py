# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 20:00:00 2023

@author: LENOVO
"""
from nodos_pesos import nodos_y_pesos_cuadratura
from matg import Mat_G
import numpy as np
from datos import datos_problema
from datos import func1
from datos import func2
from datos import phi1_0
from datos import phi2_0
from datos import curve
from datos import curvevaluation
from datos import funcevaluation
from elementos import elemn
from matplotlib import pyplot as plt
from lagrange import lag
import time 
s_0,s_1,N,K,d_t,t_f,d_s,gam,d,d_1,capt_dat=datos_problema()
tt=3000#int(t_f/(d_t))
M=np.zeros((K*N+1)*(tt)).reshape(tt,K*N+1)
N_1=np.zeros((K*N+1)*(tt)).reshape(tt,K*N+1)
A=np.zeros((K*N+1)**2).reshape(K*N+1,K*N+1)
B=np.zeros((K*N+1)**2).reshape(K*N+1,K*N+1)
b_1=np.zeros(K*N+1)
b_2=np.zeros(K*N+1)
x_n,w_n=nodos_y_pesos_cuadratura(N)
s_q=np.array(elemn(x_n,N,K,d_s,s_0))
D_v,l=lag(1,x_n,N)
G_1=[]
x,y,z,F,F_p,G,H=curve()
F,F_p,G,H=funcevaluation(F,F_p,G,H)
l_0=0
datos=[]
start=time.time()
for i in range(len(s_q)):
    M[0][i]=phi1_0(s_q[i])
    N_1[0][i]=phi2_0(s_q[i]) 
t_1=[]  
n_1=[1]
n=1
t_2=0
s=np.zeros(len(s_q))
for i in range(len(s)):
    s[i]=s_q[i]

while t_2<=t_f:
    if len(n_1)==0:
        F_0=F(s,n*d_t)*np.ones(K*N+1)
        F_p0=F_p(s,n*d_t)*np.ones(K*N+1)
        G_0=G(s,n*d_t)*np.ones(K*N+1)
        H_0=H(s,n*d_t)*np.ones(K*N+1)
    else:
        F_0=F(s,t_2+d_t)*np.ones(K*N+1)
        F_p0=F_p(s,t_2+d_t)*np.ones(K*N+1)
        G_0=G(s,t_2+d_t)*np.ones(K*N+1)
        H_0=H(s,t_2+d_t)*np.ones(K*N+1)
    
    G_1=[]
    for k in range(1,K+1):
        d_3=F_0[N*(k-1):N*k+1]*G_0[N*(k-1):N*k+1]*np.ones(N+1)
        G_1.append(Mat_G(x_n,N,d_3))
    A=np.zeros((K*N+1)**2).reshape(K*N+1,K*N+1)
    B=np.zeros((K*N+1)**2).reshape(K*N+1,K*N+1)  
    A[0][0]=(d_s/2)*w_n[0]+(d_s/2)*w_n[0]*H_0[0]*d_t
    B[0][0]=(d_s/2)*w_n[0]+(d_s/2)*w_n[0]*H_0[0]*d_t
    A[-1][-1]=(d_s/2)*w_n[-1]+(d_s/2)*w_n[-1]*H_0[-1]*d_t
    B[-1][-1]=(d_s/2)*w_n[-1]+(d_s/2)*w_n[-1]*H_0[-1]*d_t
    A[0][0:N+1]=A[0][0:N+1]+(2/d_s)*d_1*G_1[0][0]*d_t
    A[0][0:N+1]=A[0][0:N+1]+d_1*w_n[0]*F_p0[0]*G_0[0]*D_v[0]*d_t
    B[0][0:N+1]=B[0][0:N+1]+(2/d_s)*d*G_1[0][0]*d_t
    B[0][0:N+1]=B[0][0:N+1]+d*w_n[0]*F_p0[0]*G_0[0]*D_v[0]*d_t
    A[-1][(K-1)*N:K*N+1]=A[-1][(K-1)*N:K*N+1]+(2/d_s)*d_1*G_1[-1][-1]*d_t
    A[-1][(K-1)*N:K*N+1]=A[-1][(K-1)*N:K*N+1]+d_1*w_n[-1]*F_p0[-1]*G_0[-1]*D_v[-1]*d_t
    B[-1][(K-1)*N:K*N+1]=B[-1][(K-1)*N:K*N+1]+(2/d_s)*d*G_1[-1][-1]*d_t
    B[-1][(K-1)*N:K*N+1]=B[-1][(K-1)*N:K*N+1]+d*w_n[-1]*F_p0[-1]*G_0[-1]*D_v[-1]*d_t
    
    for k in range(1,K+1): 
        A[(k-1)*N+1:k*N,(k-1)*N+1:k*N][np.arange(N-1),np.arange(N-1)]=(d_s/2)*w_n[1:N]+(d_s/2)*w_n[1:N]*H_0[(k-1)*N+1:k*N]*d_t
        B[(k-1)*N+1:k*N,(k-1)*N+1:k*N][np.arange(N-1),np.arange(N-1)]=(d_s/2)*w_n[1:N]+(d_s/2)*w_n[1:N]*H_0[(k-1)*N+1:k*N]*d_t
        A[(k-1)*N+1:k*N,(k-1)*N:k*N+1]=A[(k-1)*N+1:k*N,(k-1)*N:k*N+1]+(2/d_s)*d_1*G_1[k-1][1:N]*d_t
        B[(k-1)*N+1:k*N,(k-1)*N:k*N+1]=B[(k-1)*N+1:k*N,(k-1)*N:k*N+1]+(2/d_s)*d*G_1[k-1][1:N]*d_t
        A[(k-1)*N+1:k*N,(k-1)*N:k*N+1]=A[(k-1)*N+1:k*N,(k-1)*N:k*N+1]+np.transpose(d_1*F_p0[(k-1)*N+1:k*N]*G_0[(k-1)*N+1:k*N]*d_t*(w_n[1:N]*np.transpose(D_v[1:N])))
        B[(k-1)*N+1:k*N,(k-1)*N:k*N+1]=B[(k-1)*N+1:k*N,(k-1)*N:k*N+1]+np.transpose(d*F_p0[(k-1)*N+1:k*N]*G_0[(k-1)*N+1:k*N]*d_t*(w_n[1:N]*np.transpose(D_v[1:N])))
        if k>1:
            A[(k-1)*N][(k-1)*N]=(d_s/2)*(w_n[0]+w_n[-1])+(d_s/2)*(w_n[0]+w_n[-1])*H_0[(k-1)*N]*d_t
            B[(k-1)*N][(k-1)*N]=(d_s/2)*(w_n[0]+w_n[-1])+(d_s/2)*(w_n[0]+w_n[-1])*H_0[(k-1)*N]*d_t
            A[(k-1)*N][(k-1)*N:k*N+1]=A[(k-1)*N][(k-1)*N:k*N+1]+(2/d_s)*d_1*G_1[k-1][0]*d_t
            A[(k-1)*N][(k-2)*N:(k-1)*N+1]=A[(k-1)*N][(k-2)*N:(k-1)*N+1]+(2/d_s)*d_1*G_1[k-2][-1]*d_t
            B[(k-1)*N][(k-1)*N:k*N+1]=B[(k-1)*N][(k-1)*N:k*N+1]+(2/d_s)*d*G_1[k-1][0]*d_t
            B[(k-1)*N][(k-2)*N:(k-1)*N+1]=B[(k-1)*N][(k-2)*N:(k-1)*N+1]+(2/d_s)*d*G_1[k-2][-1]*d_t
            A[(k-1)*N][(k-1)*N:k*N+1]=A[(k-1)*N][(k-1)*N:k*N+1]+d_1*w_n[0]*F_p0[(k-1)*N]*G_0[(k-1)*N]*D_v[0]*d_t
            B[(k-1)*N][(k-1)*N:k*N+1]=B[(k-1)*N][(k-1)*N:k*N+1]+d*w_n[0]*F_p0[(k-1)*N]*G_0[(k-1)*N]*D_v[0]*d_t
            A[(k-1)*N][(k-2)*N:(k-1)*N+1]=A[(k-1)*N][(k-2)*N:(k-1)*N+1]+d_1*w_n[-1]*F_p0[(k-1)*N]*G_0[(k-1)*N]*D_v[-1]*d_t
            B[(k-1)*N][(k-2)*N:(k-1)*N+1]=B[(k-1)*N][(k-2)*N:(k-1)*N+1]+d*w_n[-1]*F_p0[(k-1)*N]*G_0[(k-1)*N]*D_v[-1]*d_t
        b_1[(k-1)*N+1:k*N]=(d_s/2)*w_n[1:N]*M[n-1][(k-1)*N+1:k*N]+(d_s/2)*gam*w_n[1:N]*func1(M[n-1][(k-1)*N+1:k*N],N_1[n-1][(k-1)*N+1:k*N])*d_t
        b_2[(k-1)*N+1:k*N]=(d_s/2)*w_n[1:N]*N_1[n-1][(k-1)*N+1:k*N]+(d_s/2)*gam*w_n[1:N]*func2(M[n-1][(k-1)*N+1:k*N],N_1[n-1][(k-1)*N+1:k*N])*d_t
        if k>1:
            b_1[(k-1)*N]=(d_s/2)*(w_n[0]+w_n[-1])*M[n-1][(k-1)*N]+d_t*(d_s/2)*gam*(w_n[0]+w_n[-1])*func1(M[n-1][(k-1)*N],N_1[n-1][(k-1)*N])
            b_2[(k-1)*N]=(d_s/2)*(w_n[0]+w_n[-1])*N_1[n-1][(k-1)*N]+d_t*(d_s/2)*gam*(w_n[0]+w_n[-1])*func2(M[n-1][(k-1)*N],N_1[n-1][(k-1)*N])
    b_1[0]=(d_s/2)*w_n[0]*M[n-1][0]+d_t*(d_s/2)*gam*w_n[0]*func1(M[n-1][0],N_1[n-1][0])
    b_2[0]=(d_s/2)*w_n[0]*N_1[n-1][0]+d_t*(d_s/2)*gam*w_n[0]*func2(M[n-1][0],N_1[n-1][0])
    b_1[-1]=(d_s/2)*w_n[-1]*M[n-1][-1]+d_t*(d_s/2)*gam*w_n[-1]*func1(M[n-1][-1],N_1[n-1][-1])
    b_2[-1]=(d_s/2)*w_n[-1]*N_1[n-1][-1]+d_t*(d_s/2)*gam*w_n[-1]*func2(M[n-1][-1],N_1[n-1][-1])    
    M[n]=np.linalg.solve(A,b_1)
    N_1[n]=np.linalg.solve(B,b_2)
    if max(abs(M[n]-M[n-1]))<1e-3:
        if len(t_1)==0:
            n_1.append(n)
            t_1.append(n*d_t)
        else:
            t_1.append((n-n_1[-1])*d_t)
            n_1.append(n) 
        d_t=2*d_t
    elif max(abs(M[n]-M[n-1]))>1e-1:
        if len(t_1)==0:
            n_1.append(n)
            t_1.append(n*d_t)
        else:
            t_1.append((n-n_1[-1])*d_t)
            n_1.append(n) 
        d_t=d_t/2
    t_2=sum(t_1)+(n-n_1[-1])*d_t
    n=n+1
    if t_2>=capt_dat[l_0]-1 and t_2<=capt_dat[l_0]+1:
        if l_0<len(capt_dat)-1:
            datos.append(n)
            l_0=l_0+1
        else:
            datos.append(n)
end=time.time()           
print(end-start)    
plt.plot(s,M[n-1])
plt.plot(s,N_1[n-1])
(x_k,y_k,z_k)=curvevaluation(x, y, z, s, t_f)
x_k=x_k*np.ones(K*N+1)
y_k=y_k*np.ones(K*N+1)
z_k=z_k*np.ones(K*N+1)
ax = plt.figure().add_subplot(projection='3d')
for i in range(len(s)-1):
    ax.plot(x_k[i:i+2],y_k[i:i+2],z_k[i:i+2],color=[(M[n-1][i]-min(M[n-1]))/(max(M[n-1])-min(M[n-1])),0,0],linewidth=2)     

tim=np.zeros(n)
ml=0
n_1=n_1[1:]

for k in range(1,n):
    if k==1:
        tim[k]=t_1[0]/n_1[0]
        d_t=t_1[0]/(n_1[0])
        ml=ml+1
    elif k-1 in n_1 and (k>1 and ml<len(n_1)):
        d_t=t_1[ml]/(n_1[ml]-n_1[ml-1])
        ml=ml+1
    elif k-1>n_1[len(n_1)-1] and ml>=len(n_1):
        d_t=(t_2-sum(t_1))/(n-1-n_1[-1])
    if k>1:
        tim[k]=tim[k-1]+d_t
    
S,T=np.meshgrid(s,tim)
for k in range(0,n):
    (x_k,y_k,z_k)=curvevaluation(x, y, z, s[-1], tim[k])
    S[k]=S[k]*x_k/max(s)
fig, ax = plt.subplots()
im=ax.pcolormesh(S, T, M[0:n])
fig.colorbar(im)

# for j in range(int((n+1)/10)):
#     plt.clf()
#     plt.plot(s,M[j*10])
#     plt.show()    
    
    
            
            
        
        
