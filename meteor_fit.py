import numpy as n
import matplotlib.pyplot as plt

def simple_fit_meteor(t,ea,no,up):
    """
    all in units of km
    """
    tmean=n.mean(t)
    nm=len(t)
    tn=t-tmean
    # r0x, r0y, r0z, v0x, v0y, v0z
    A=n.zeros([3*nm,6])
    m=n.zeros(3*nm)
    # position measurement
    
    # east pos,  r0x, 
    A[0:nm,0]=1
    # west pos r0y
    A[(nm):(2*nm),1]=1
    # up pos r0z
    A[(2*nm):(3*nm),2]=1
      
    # east pos v0x, 
    A[0:nm,3]=tn
    # west pos v0y
    A[(nm):(2*nm),4]=tn
    # up pos v0z
    A[(2*nm):(3*nm),5]=tn

    m[0:nm]=ea
    m[nm:(2*nm)]=no
    m[(2*nm):(3*nm)]=up

    xhat=n.linalg.lstsq(A,m)[0]
    r0=xhat[0:3]
    v0=xhat[3:6]
    model=n.dot(A,xhat)
    r2=n.abs(m-model)**2
    ea_res=n.sqrt(n.mean(r2[0:nm]))
    no_res=n.sqrt(n.mean(r2[nm:(2*nm)]))
    up_res=n.sqrt(n.mean(r2[(2*nm):(3*nm)]))
    
    return(r0,v0,model,ea_res,no_res,up_res)

      

    
    
