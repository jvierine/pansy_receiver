import numpy as n
import matplotlib.pyplot as plt
import pansy_config as pc
import itertools
import h5py

def image_points(phcal,xc,ch_pairs,u,v,w,dmat):
    """
    image all points
    use phcal applied to the phase of each channel
    xc is the cross correlation
    ch_pair is the antenna pair index structure
    u,v,w is the 
    """
    mu=[]
    mv=[]
    mis=[]
    mjs=[]
    mfs=[]
    for i in range(xc.shape[0]):
        z=n.exp(1j*(n.angle(xc[i,:]) + phcal[ch_pairs[:,0]] - phcal[ch_pairs[:,1]]))
        M=n.abs(mf(z,dmat,u,v,w))
        mi,mj=n.unravel_index(n.argmax(M),shape=M.shape)
        if False:
            print(M[mi,mj])
            plt.pcolormesh(u,v,M)
            plt.colorbar()
            plt.plot(u[mi,mj],v[mi,mj],"x")
            plt.show()
        mfs.append(M[mi,mj])
        mu.append(u[mi,mj])
        mv.append(v[mi,mj])
        mis.append(mi)
        mjs.append(mj)
    mu=n.array(mu)
    mv=n.array(mv)
    mfs=n.array(mfs)
    mis=n.array(mis)
    mjs=n.array(mjs)

    return(mu,mv,mfs,mis,mjs)

# index pairs
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
def uv_coverage(N=100,max_zenith_angle=10):
    """
    unit vectors
    """
    max_u=n.arcsin(n.pi*max_zenith_angle/180.0)
    u=n.linspace(-max_u,max_u,num=N)
    v=n.linspace(-max_u,max_u,num=N)
    uu,vv=n.meshgrid(u,v)
    mag=n.sqrt(uu**2.0+vv**2)
    ww=-n.sqrt(1-uu**2-vv**2)
#    ww[mag>=1]=n.nan
 #   uu[mag>=1]=n.nan
  #  vv[mag>=1]=n.nan
    return(uu,vv,ww)

def uv_coverage2(N=100,max_angle=10,az0=0,el0=90.0):
    """
    unit vectors
    az0, el0 is beam pointing direction
    """
    u0_h=n.cos(n.pi*el0/180.0)
    w0=n.sin(n.pi*el0/180.0)
    v0=u0_h*n.cos(-n.pi*az0/180)
    u0=-u0_h*n.sin(-n.pi*az0/180)

    # not entirely correct!
    max_u=n.arcsin(n.pi*max_angle/180.0)
    u=n.linspace(-max_u+u0,max_u+u0,num=N)
    v=n.linspace(-max_u+v0,max_u+v0,num=N)
    uu,vv=n.meshgrid(u,v)
    mag=n.sqrt(uu**2.0+vv**2)
    ww=-n.sqrt(1-uu**2-vv**2)
#    ww[mag>=1]=n.nan
 #   uu[mag>=1]=n.nan
  #  vv[mag>=1]=n.nan
    return(uu,vv,ww)

def zenang(u,v,w):
    return(180*n.arctan(n.sqrt(u**2+v**2)/w)/n.pi)

def pair_mat(ch_pairs,antpos):
    """
    Each row of the matrix is a vector from antenna 0 to antenna 1
    """
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

def mf(meas,dmat,u,v,w):
    """
    given a measurement, calculate the match function (beam forming) response
    for different pointing directions
    """
    M=n.zeros([u.shape[0],u.shape[1]],dtype=n.complex64)
    
    # for each pair
    k0=2*n.pi/pc.wavelength
    for i in range(meas.shape[0]):
        M+=meas[i]*n.exp(-1j*k0*(dmat[i,0]*u+dmat[i,1]*v+dmat[i,2]*w))
        
    return(M)
#    plt.pcolormesh(n.abs(mf))
 #   plt.show()
    
    #dmat[:,0]*u

    
    
def get_antpos():
    antpos=[]
    for cid in pc.connections:
        if cid != "RFTX":
            p=n.array([pc.module_center[cid][0],pc.module_center[cid][1],pc.module_center[cid][2]])
            antpos.append(p)
    antpos=n.array(antpos)
    return(antpos)


def find_meso_zen_cal():
    """
    use zenith direction echoes to figure out the phase delay of each module
    """
    ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
    
    # A bunch of measurements of PMSE in the zenith direction
    h=h5py.File("data/mesocal.h5","r")
    meas_phase=h["cals"][()]
    h.close()
    meas_phase=n.exp(1j*n.angle(meas_phase))
    # the first 7 are self correlations!
    meas_phase=meas_phase[7:len(meas_phase)]
    # these are the antenna module positions
    antpos=get_antpos()
    # matrix of antenna pair distances
    dmat=pair_mat(ch_pairs,antpos)
    
    # expected phase differences for zenith direction (EM plane wave in -z direction)
    # assuming the sky is not falling on us, the up direction should be the same as the zero doppler shift
    # pmse echo direction (on average)
    model_phase=u_phases(dmat,n.array([0,0,-1]))

    # zenith direction phase differences with model zenith direction phases subtracted
    meas=meas_phase*n.conj(model_phase)
    
    nm=len(model_phase)
    model=n.zeros(nm,dtype=n.complex64)
    def ss(x):
        phi=x
        for chi,chp in enumerate(ch_pairs):
            model[chi]=n.exp(1j*(phi[chp[0]]-phi[chp[1]]))
        s=n.sum(n.abs(meas-model)**2.0)
        return(s)
    import scipy.optimize as sio
    xhat=sio.fmin(ss,n.zeros(7))
    xhat=sio.fmin(ss,xhat)
    xhat=sio.fmin(ss,xhat)
    xhat=sio.fmin(ss,xhat)            
    print(xhat)
    
    for chi,chp in enumerate(ch_pairs):
        print("%d %1.2f %1.2f"%(chi,n.angle(meas_phase[chi]*n.conj(model_phase[chi])),n.angle(n.exp(1j*(xhat[chp[0]]-xhat[chp[1]])))))

    h=h5py.File("data/phases.h5","w")
    h["zenith"]=xhat
    h.close()
        

def psf():
    u,v,w=uv_coverage(N=500)
    ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))    
    antpos=get_antpos()
    print(antpos)
    plt.plot(antpos[0:3,0],antpos[0:3,1],".")
    plt.show()
    dmat=pair_mat(ch_pairs[0:3,:],antpos)
    zangs=n.linspace(-30,30,num=100)
    for i in range(len(zangs)):
        zang=zangs[i]
        zenith_point=u_phases(dmat,n.array([0,n.sin(n.pi*zang/180),-n.cos(n.pi*zang/180)]))
        M=mf(zenith_point,dmat,u,v,w)
        za=zenang(u,v,w)
        fig,ax=plt.subplots()
        ax.pcolormesh(u,v,n.abs(M))
        plt.title("zenith angle %1.1f (deg)"%(zangs[i]))
        ax.plot([0],[n.sin(n.pi*zang/180.0)],"x",color="black")
        cs=ax.contour(u,v,za)
        ax.clabel(cs)
        ax.set_xlabel("u")
        ax.set_ylabel("m")        
        plt.gca().set_aspect('equal')        
        plt.show()
#        plt.savefig("psf-%03d.png"%(i))
 #       plt.close()
#        plt.show()

if __name__ =="__main__":
    psf()
    find_meso_zen_cal()#ch_pairs, model_phase, meas_phase)
