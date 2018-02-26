import numpy as np
s = 100
k = 100
up = 2 
down = 1/2
e_tr = 3/2
n = 10
q_up = (e_tr - down)/(up - down)
q_down =  1 - q_up

S = np.zeros((n+1, n+1), dtype = np.float64)
for i in range(n+1):
    m = 0
    while m <= i:
        S[m,i] = s * (up**i) * (down**(m*2))
        m += 1
        
V = np.zeros((n+1,n+1), dtype = np.float64)

for i in range(n+1):
    V[i, n] = max(k-S[i, n],0)

for i in range(n-1,-1,-1):
    m=0
    while m <= i:
        V[m, i] = max(k - S[m, i], (e_tr**-1)*(q_up * V[m, i + 1] + q_down * V[m + 1, i + 1]))
        m += 1

print("The price of the American put is "+str(V[0,0]))