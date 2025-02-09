import numpy as n
import matplotlib.pyplot as plt
import pansy_config as pc
import itertools

# index pairs
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
def uv_coverage(N=100):
    u=n.linspace(-1,1,num=N)
    v=n.linspace(-1,1,num=N)
    uu,vv=n.meshgrid(u,v)
    mag=n.sqrt(uu**2.0+vv**2)
    ww=n.sqrt(1-uu**2-vv**2)
    ww[mag>=1]=n.nan
    uu[mag>=1]=n.nan
    vv[mag>=1]=n.nan
    return(uu,vv,ww)

def pair_mat(ch_pairs,antpos):
    dmat=n.zeros([ch_pairs.shape[0],3],dtype=n.float32)
    for i in range(ch_pairs.shape[0]):
        dmat[i,:]=antpos[ch_pairs[i,1],:]-antpos[ch_pairs[i,0],:]
    return(dmat)


def u_phases(dmat,u0):
    k0=2*n.pi/pc.wavelength
    phases0=n.exp(1j*k0*n.sum(dmat[:,:]*u0[None,:],axis=1))
    if False:
        phases=n.zeros(ch_pairs.shape[0],dtype=n.complex64)
        # the above code does the same but faster
        for i in range(ch_pairs.shape[0]):
            d=antpos[ch_pairs[i,1],:]-antpos[ch_pairs[i,0],:]
            phases[i]=n.exp(1j*k0*n.dot(d,u0))  
    return(phases0)

def get_antpos():
    antpos=[]
    for cid in pc.connections:
        if cid != "RFTX":
            p=n.array([pc.module_center[cid][0],pc.module_center[cid][1],pc.module_center[cid][2]])
            antpos.append(p)
    antpos=n.array(antpos)
    return(antpos)

if __name__ =="__main__":
    print(pc.wavelength)
    antpos=get_antpos()
    dmat=pair_mat(ch_pairs,antpos)
    phases=u_phases(dmat,n.array([0,0,-1]))
    print(phases)
    exit(0)

    for i in range(ch_pairs.shape[0]):
        print(ch_pairs[i,:])
        print(n.angle(phases[i]))