"""part n"""
import math
import numpy as np
import pandas as pd
B=80
K=110
R=300
Nx=1000
r=0.05
sigma=0.2
delta_x = (R-B)/Nx
xk = []
vj = np.zeros(shape=(Nx+1,1))
for i in range(0,1001):
    xk.append(max(B+(i*delta_x),0))

        
vnt = np.asmatrix(vj) 


Nt=100
T=1
delta_t = T/Nt
t=[]
for j in range(0,101):
    t.append(j*delta_t)
    
def CMatrix(j):
    cm1 = np.zeros(shape=(Nx+1,1))
    cm1[Nx] = R - (math.exp(-r*(T-t[j])))*K
    cm1 = np.asmatrix(cm1)
    return cm1


def m_kk(k):
    if k==0 or k==Nx:
        m_kk=0
    else:
        m_kk = 1 - (r*delta_t)-((sigma*sigma*delta_t*(xk[k]*xk[k]))/(delta_x*delta_x))
    return m_kk
    
def mkk_right(k):
    if k==0 or k==Nx:
        mkk1=0
    else:
        mkk1 = delta_t*(((r*xk[k])/(2*delta_x))+((sigma*sigma*xk[k]*xk[k])/(2*delta_x*delta_x)))
    return mkk1

def mkk_left(k):
    if k==0 or k==Nx:
        mkk2=0
    else:
        mkk2 = delta_t*(((-r*xk[k])/(2*delta_x))+((sigma*sigma*xk[k]*xk[k])/(2*delta_x*delta_x)))
    return mkk2
    
    

M = np.zeros(shape=(Nx+1,Nx+1))
for p in range(0,Nx+1):
    if p!=0 and p!=1000:
        M[p][p-1] = mkk_left(p)
        M[p][p] = m_kk(p)
        M[p][p+1] = mkk_right(p)
M = np.asmatrix(M)
        



V_initial = vnt
for j in range(1,40):
    V_initial = CMatrix(Nt-j) + (M*V_initial)
    
df1 = pd.DataFrame(V_initial)
df1.to_csv('output.csv')


        
    



