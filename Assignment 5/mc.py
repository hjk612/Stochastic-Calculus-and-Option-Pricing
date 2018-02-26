import numpy as np
from math import sqrt, exp

def mc(nmbrsimul,dt,r,sigma,K,B,T):

    final=[]
    path=[]
    for i in np.arange(nmbrsimul):
        St=[];payoff=[]
        for t in np.linspace(0,T,num=100):
            if t==0:
                St.append(100)
                payoff.append(0)
            elif t!=0:
                St_next=St[-1] + r* St[-1] * dt + sigma * St[-1] * np.random.normal(0,sqrt(dt))
                St.append(St_next)
                if min(St)>B:
                    pf=max(St[-1]-K,0)
                    payoff.append(pf)
                elif min(St)<=B:
                    pf=0
                    payoff.append(pf)
        path.append(payoff)  
    for i in range(nmbrsimul):
        final.append(path[i][-1])
    price=(exp(-r*T)*sum(final))/nmbrsimul
    print(price)
    
mc(10000,1/100,0.05,0.20,110,80,1)

    