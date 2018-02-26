import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


Smax = 450
T = 1
r = 0.01
sigma = 0.2 
K = 100
def put_func(M,N):
    dx = Smax/float(M)
    dt = T/float(N)

    
    i=np.arange(1,M,dtype=np.float)
    
    a = -0.5 * dt * sigma ** 2 * ((i * dx) ** 2)/(dx ** 2)
    b = 1 + dt * (r + r * i * dx/dx + (sigma ** 2) * (i * dx) ** 2/(dx ** 2))
    c = dt * (-r * i * dx/dx - 0.5 * sigma ** 2 * ((i * dx) ** 2)/dx ** 2)
    
    A=np.diag(b)+np.diag(a[1:],k=-1)+np.diag(c[0:M-2],k=1)
    A_inverse = np.linalg.inv(A)
    A_inverse = np.matrix(A_inverse)
    put = np.zeros((N+1,M+1))
    put[0,:]=np.maximum(K - np.arange(0,Smax+dx/2.0,dx,dtype=np.float),0)
    put[:,0]=[K * np.exp(-r*(j)*dt) for j in range(N+1)]
    put[:,M]=0
    
    
    tran = np.zeros((0,M-1))
    tran = (np.maximum(K-np.arange(dx,Smax-dx/2.0,dx,dtype=np.float),0))
    d = np.empty_like (put[0,1:M])
    my_list = [put[0,0]]
    for j in range(0,N,1):
        d[:] = put[j,1:M]
        d[0] = d[0] - a[0] * K * np.exp(-r*(j+1)*dt) 
        
        put[j+1, 1:M] = (A_inverse * np.matrix(d).transpose()).transpose()
        put[j+1,1:M] = np.maximum(K-np.arange(dx,Smax-dx/2.0,dx,dtype=np.float),put[j+1,1:M])
        itemindex = min(np.where(np.equal(tran,put[j+1,1:M])==False)[0])
        my_list.append(itemindex*dx)
     
    put = pd.DataFrame(put)
    list = np.arange(0, T+dt, dt)
    index = np.arange(0, Smax+dx, dx)
    put.index = list
    put.columns = index
    return put, my_list

####put values for plotting dx error keeping dt constant

"""
put_3000, mylist_900 = put_func(3000,1000)
put_2000, mylist_900 = put_func(2000,1000)
put_1500, mylist_900 = put_func(1500,1000)
put_1000, mylist_900 = put_func(1000,1000)
put_900, mylist_900 = put_func(900,1000)
put_850, mylist_850 = put_func(850,1000)
put_800, mylist_800 = put_func(800,1000)"""

                              
####put values for plotting dt error keeping dx constant
"""
put_900, my_list = put_func(6000,900)   
put_800, my_list_200 = put_func(6000,800)
put_700, my_list_300 = put_func(6000,700)
put_600, my_list_400 = put_func(6000,600)
put_500, my_list_450 = put_func(6000,500)
put_400, my_list_450 = put_func(6000,400)"""
def free_boundary(put, my_list):
    plt.plot(put.index,my_list)
    plt.show()
    
def american_put_plot(put):
    x = put.index
    y = put.columns
    Y,X = np.meshgrid(y,x)
    Z = put
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)

put_reference, my_list_reference = put_func(6000, 1000)

def error_dt(put):
    dt_error = put.ix[1] - put_reference.ix[1]
    error = np.linalg.norm(np.array(dt_error))
    return error

def error_dx(put):
    p = scipy.interpolate.interp1d(np.array(put.columns),np.array(put.ix[1]),kind='linear')
    list = []
    for i in range(len(put_reference.columns)):
        list.append(np.float(p(put_reference.columns[i])))
    error_dx = put_reference.ix[1] - list
    error = np.linalg.norm(np.array(error_dx))
    return error

free_boundary(put_reference, my_list_reference)
american_put_plot(put_reference)


#### plot for dx error
"""list = [error_dx(put_800),error_dx(put_850),error_dx(put_900),error_dx(put_1000),error_dx(put_1500),error_dx(put_2000),error_dx(put_3000)]
k = [450/800,450/850,450/900,450/1000,450/1500,450/2000,450/3000]
plt.plot(list,k)
plt.show()"""

####plot for dt error
""""
error_dt_1 = error_dt(put_900)
error_dt_2 = error_dt(put_800)
error_dt_3 = error_dt(put_700)
error_dt_4 = error_dt(put_600)
error_dt_5 = error_dt(put_500)
error_dt_6 = error_dt(put_400)

list = [error_dt_1, error_dt_2, error_dt_3, error_dt_4, error_dt_5, error_dt_6]
k = [1/900,1/800,1/700,1/600,1/500,1/400]
plt.plot(list, k)
plt.show()"""
    
    



