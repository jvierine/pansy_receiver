import itertools
import pansy_interferometry as pint
import h5py
import numpy as n
import scipy.optimize as so
import matplotlib.pyplot as plt
import glob

def get_xcs(beam_num=0):
    """
    get n strongest cross-phases from each file (one meteor)
    """
    fl=glob.glob("caldata/meteor-%d-*.h5"%(beam_num))
    xc=[]
    for f in fl:
        h=h5py.File(f,"r")
        mi=n.argsort(n.abs(h["xc"][()][0,:]))[::-1]
        #print(mi)
        for k in range(n.min((len(mi),10))):                
                xc.append([h["xc"][()][:,mi[k]]])
        h.close()
    xc=n.vstack(xc)
    print(xc.shape)
    return(xc)

def image_points(phcal,xc2,ch_pairs,u,v,w,dmat):
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
    for i in range(xc2.shape[0]):
        z=n.exp(1j*(n.angle(xc2[i,:]) + phcal[ch_pairs[:,0]] - phcal[ch_pairs[:,1]]))
        M=n.abs(pint.mf(z,dmat,u,v,w))
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

# create a uniform grid of u,v,w values
u,v,w=pint.uv_coverage(N=120,max_zenith_angle=10.0)
# antenna positions
antpos=pint.get_antpos()
# channel pairs used for interferometry
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
# matrix of antenna direction vector pairs 
# one row is a vector for each interferometry channel pair
dmat=pint.pair_mat(ch_pairs[:,:],antpos)

# find phases that rotate pointing towards phac
def find_phaserot(phac,ch_pairs):
    def ss(x):
        phcal=n.zeros(7)
        phcal[1:7]=x
        zmodel=n.zeros(ch_pairs.shape[0],dtype=n.complex64)
        for i in range(ch_pairs.shape[0]):
            zmodel[i]=n.exp(1j*(phcal[ch_pairs[i][0]]-phcal[ch_pairs[i][1]]))
        s=n.sum(n.abs(zmodel-phac)**2.0)
        #print(s)
        return(s)
    xhat=so.fmin(ss,n.zeros(6))                    
    xhat=so.fmin(ss,xhat)                    
    xhat=so.fmin(ss,xhat)                    
    xhat=so.fmin(ss,xhat)                    
    return(xhat)

xc=get_xcs()
mags=n.abs(xc[:,0])

global xc2
phcal=n.zeros(7,dtype=n.float32)

def ss(x):
    global xc2
    phcal[1:len(phcal)]=x
    mu,mv,mfs,mis,mjs=image_points(phcal,xc2,ch_pairs,u,v,w,dmat)
    mean_off=n.mean(n.sqrt(mu**2.0+mv**2.0))
    mean_mf=n.sum(mfs)/len(mfs)/ch_pairs.shape[0]
    
    wg=mfs/ch_pairs.shape[0]
    wgs=n.sum(wg)
    wmu=n.sum(mu*wg)/wgs
    wmv=n.sum(mv*wg)/wgs
    print("mf %1.2f center %1.2f %1.2f"%(mean_mf,wmu,wmv))
    s=-mean_mf+wmu**2.0+wmv**2.0
    print(s)
    return(s)

# phase calibration guess
xhat=[ 1.14248621,  0.15371235,  2.5528059,  -0.35319514,  3.19109113,  3.58991391]
phcal=n.zeros(7,dtype=n.float32)
phcal[1:len(phcal)]=xhat
#mu,mv,mfs,mis,mjs=image_points(phcal,xc,ch_pairs,u,v,w,dmat)

def direction_shift(dir_orig,
                    dir_new,
                    phcal_orig,
                    dmat,
                    ch_pairs):
    """
     find calibration offset that rotates vector 
     towards (um,vm,-wm) to -> (0,0,-1)
    """
    # we want to shift all interferometry solutions from this direction to zenith
    # using a new calibration xhat+xhat2
    um=0.05
    vm=-0.008
    wm=n.sqrt(1-um**2.0-vm**2.0)
    # interferometer phases pointing towards um, vm, -wm
    phac=pint.u_phases(dmat,n.array([um,vm,-wm]))
    # interferometer phases pointing towards zenith
    phacz=pint.u_phases(dmat,n.array([0,0,-1]))
    dphac=n.conj(phac*n.conj(phacz))
    xhat2=find_phaserot(dphac,ch_pairs)

    # the phase cal now contains a part that
    # focuses everything (xhat), 
    # and the rotation (xhat2) towards wanted direction
    phcal=n.zeros(7,dtype=n.float32)
    phcal[1:len(phcal)]=xhat2
    phcal += phcal_orig
    return(phcal)

phcal=direction_shift(dir_orig=n.array([0.05,-0.008,n.sqrt(1-0.05**2-0.008**2)]),
                      dir_new=n.array([0,0,-1]),
                      phcal_orig=phcal,
                      dmat=dmat,
                      ch_pairs=ch_pairs)

mu,mv,mfs,mis,mjs=image_points(phcal,xc,ch_pairs,u,v,w,dmat)
gidx=n.where(mfs/21.0 > 0.95)[0]
xc=xc[gidx,:]
mu=mu[gidx]
mv=mv[gidx]
mfs=mfs[gidx]
mis=mis[gidx]
mjs=mjs[gidx]
mags=mags[gidx]

# IT WORKS!
#plt.scatter(mu,mv,c=10.0*n.log10(mags[gidx]),s=2)
pwr=n.copy(u)
pwr[:,:]=0.0
npwr=n.copy(u)
npwr[:,:]=0.0
for i in range(len(mis)):
    pwr[mis[i],mjs[i]]+=mags[i]
    npwr[mis[i],mjs[i]]+=1.0#mag[gidx[i]]
dB=10.0*n.log10(pwr/npwr)
dB[npwr==0]=n.nan
nfloor=n.nanmin(dB)

# zenith angles
za=pint.zenang(u,v,w)

fig,ax=plt.subplots()
m=ax.pcolormesh(u,v,dB-nfloor,vmin=0,vmax=n.nanmax(dB-nfloor))
cb=fig.colorbar(m,ax=ax)
cb.set_label("Mean power (dB)")
cs=ax.contour(u,v,za,colors="black")
ax.clabel(cs,fmt=r"%1.0f$^{\circ}$",colors="black")
ax.set_xlabel("u")
ax.set_ylabel("v")
#plt.colorbar()
plt.show()

plt.scatter(mu,mv,c=mfs/21.0,s=2)
cb=plt.colorbar()
cb.set_label("Coherence")
plt.xlabel("u")
plt.ylabel("v")
plt.plot(0,0,"x")
plt.show()
