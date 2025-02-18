import itertools
import pansy_interferometry as pint
import h5py
import numpy as n
import scipy.optimize as so
import matplotlib.pyplot as plt
import glob
import os
import re

def fit_chirp(tx_idx,xc):
    """
    if a linear chirp to phases. can be used to cross-calibrate 
    phases for different transmit pointings
    """
    tv=(tx_idx-tx_idx[0])/1e6
    w=n.abs(xc)
    m=n.exp(1j*n.angle(xc))
    def ss(x):
        om=x[0]
        om2=x[1]
        model=n.exp(1j*(om+om2*tv))
        return(n.sum(w*n.abs(model-m)**2.0))
    guess=n.linspace(-100,100,num=1000)
    best=guess[0]
    best_s=1e99
    for i in range(len(guess)):
        s=ss([n.angle(m[0]),guess[i]])
        if s < best_s:
            best_s=s
            best=guess[i]

    xhat=so.fmin(ss,[n.angle(m[0]),best],disp=False)
    om=xhat[0]
    om2=xhat[1]
    model=n.exp(1j*(om+om2*tv))
    return(xhat,model,tx_idx[0])

def fit_beam_cal(xc0,xcn,weight,ch_pairs):
    """
    xc0 = xc on zenith beam
    xcn = xc on a non-zenith beam
    weight = goodness of xc
    ch_pairs = pairs of antenna channels cross-correlated

    Assume that z_a^1 = z_a^0 e^{i\phi^1_a}.
    """
    def ss(phases):
        phasecal=n.exp(1j*(phases[ch_pairs[:,0]]-phases[ch_pairs[:,1]]))
        model = xc0*phasecal[None,:]
        res=model-xcn
        s=n.sum(res.real**2.0 + res.imag**2.0)
        print(s)
        return(s)
    xhat=so.fmin(ss,n.zeros(7,n.float32))
    print(xhat)

def get_xcs2(max_per_event=30):
    """
    get n strongest cross-phases from each file (one meteor)
    """
    fl=glob.glob("caldata/meteor-0-*.h5")
    fl.sort()
    xc=[]
    ch_pairs=n.array(list(itertools.combinations(n.arange(7,dtype=n.int64),2)))
    #print(ch_pairs)
    #exit(0)
    cal_xc=[]
    for i in range(5):
        cal_xc.append({"xc0":[],"xcn":[]})

    for f in fl:
        h=h5py.File(f,"r")

        pref=re.search("(caldata/meteor-).-(.*.h5)",f).group(1)
        post=re.search("(caldata/meteor-).-(.*.h5)",f).group(2)
        beams=[0]
        beam_xc=[]
        beam_txi=[]
        for i in range(1,5):
            beamf="%s%d-%s"%(pref,i,post)
            print(beamf)
            if os.path.exists(beamf):
                h2=h5py.File(beamf,"r")
                beam_xc.append(h2["xc"][()])
                beam_txi.append(h2["tx_idx"][()])
                beams.append(i)
#                plt.plot(txib,n.angle(xcb[0,:]),".")
                h2.close()
            else:
                print("%s doesn't exist"%(beamf))


        if len(beams)>1:
            print(beams)
            xc=h["xc"][()]
            txi=h["tx_idx"][()]
            # this is where there can be overlap
            mintx=n.min(txi)
            maxtx=n.max(txi)

            pref=re.search("(caldata/meteor-).-(.*.h5)",f).group(1)
            post=re.search("(caldata/meteor-).-(.*.h5)",f).group(2)
            for i in range(1,5):
                beamf="%s%d-%s"%(pref,i,post)
                if os.path.exists(beamf):
                    h2=h5py.File(beamf,"r")
                    xcb=h2["xc"][()]
                    txib=h2["tx_idx"][()]
                    h2.close()

                    overlap_idx=n.where( (txib > mintx) & (txib < maxtx) )[0]

                    n_overlap=len(overlap_idx)
                    xc0=n.zeros([n_overlap,21],dtype=n.complex64)
                    xcn=n.zeros([n_overlap,21],dtype=n.complex64)

                    for ci in range(21):
                        xhat,model,t0=fit_chirp(txi,xc[ci,:])

                        def modelfun(tidx):
                            return(n.exp(1j*xhat[0])*n.exp(1j*xhat[1]*(tidx-t0)/1e6))
                        
                        xc0[:,ci]=modelfun(txib[overlap_idx])
                        xcn[:,ci]=xcb[ci,overlap_idx]# modelfun(txib[overlap_idx])
                    # add to pile of calibrations
                    cal_xc[i]["xc0"].append(xc0)
                    cal_xc[i]["xcn"].append(xcn)
        

                    if False:
                        plt.plot(txi/1e6,n.angle(xc[ci,:]),".")        
                        plt.title("ci %d beam %d"%(ci,i))                
                        plt.plot(txib[overlap_idx]/1e6,n.angle(xc0[:,ci]))
                        plt.plot(txib[overlap_idx]/1e6,n.angle(xcn[:,ci]),".")
                        plt.ylim([-n.pi,n.pi])
                        plt.show()
                    if False:
                        plt.plot(txi/1e6,n.angle(modelfun(txi)))
                        plt.plot(txi/1e6,n.angle(xc[ci,:]),".",label="beam 0",color="C0")

                        plt.plot(txib/1e6,n.angle(xcb[ci,:]),".",label="beam %d"%(i),color="C%d"%(i))

                        plt.legend()
                        plt.ylim([-n.pi,n.pi])
                        plt.title("XC %d"%(ci))
                        plt.xlabel("Time (unix)")
                        plt.ylabel("Cross-phase (rad)")
                        plt.show()

        #if n.random.rand(1)<0.05:    
        ho=h5py.File("caldata/cal_all.h5","w")            
        for i in range(5):
            print(len(cal_xc[i]["xc0"]))
            
            if len(cal_xc[i]["xc0"])>0:
                xc0=n.vstack(cal_xc[i]["xc0"])
                xcn=n.vstack(cal_xc[i]["xcn"])
                ho["beam%d/xc0"%(i)]=xc0
                ho["beam%d/xcn"%(i)]=xcn
        ho.close()

                #for j in range(1):
#                j=0
 #               plt.plot(n.angle(xc0[:,j]*n.conj(xcn[:,j])),".",label="beam, %d (0,1)"%(i))
#                plt.title("beam %d %d-%d"%(i,ch_pairs[j,0],ch_pairs[j,1]))
  #      plt.ylim([-n.pi,n.pi])
   #     plt.legend()
    #    plt.show()

#            cal_xc.append({"xc0":[],"xcn":[]})

        #for i in range(1,5):
         #   if os.path.exists()
        h.close()

#    xc=n.vstack(xc)
 #   print(xc.shape)
  #  return(xc)

def get_xcs(beam_num=0,max_per_event=30):
    """
    get n strongest cross-phases from each file (one meteor)
    """
    fl=glob.glob("caldata/meteor-%d-*.h5"%(beam_num))
    fl.sort()
    xc=[]
    for f in fl[0:200]:
        h=h5py.File(f,"r")
        mi=n.argsort(n.abs(h["xc"][()][0,:]))[::-1]
        #print(mi)
        for k in range(n.min((len(mi),max_per_event))):
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

#
# find phases that rotate pointing towards phac
# phac can be e.g., zenith pointing obtained with:
# phac_zenith=pint.u_phases(dmat,n.array([0,0,-1]))
# 
# tbd: this should go into pint
#
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

def find_initial_phasecal():
    """
    """
    # create a uniform grid of u,v,w values
    u,v,w=pint.uv_coverage(N=100,max_zenith_angle=10.0)
    # antenna positions
    antpos=pint.get_antpos()
    # channel pairs used for interferometry
    ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
    # matrix of antenna direction vector pairs 
    # one row is a vector for each interferometry channel pair
    dmat=pint.pair_mat(ch_pairs[:,:],antpos)
    print("reading xc")
    xc=get_xcs(max_per_event=3)
    mags=n.abs(xc[:,0])
    idx=n.arange(xc.shape[0])
    ridx=n.random.permutation(idx)[0:200]
    phcal=n.zeros(7,dtype=n.float32)


    xc2=get_xcs(max_per_event=30)
    # plot detections and figure out what is the center
    mu,mv,mfs,mis,mjs=image_points(phcal,xc2[0:1000,:],ch_pairs,u,v,w,dmat)
    plt.scatter(mu,mv,c=mfs/21.0,s=2)
    cb=plt.colorbar()
    cb.set_label("Coherence")
    plt.xlabel("u")
    plt.ylabel("v")
    plt.plot(0,0,"x")
    plt.show()


    def ss(x):
        phcal[1:len(phcal)]=x
        mu,mv,mfs,mis,mjs=image_points(phcal,xc[ridx,:],ch_pairs,u,v,w,dmat)
        mean_off=n.mean(n.sqrt(mu**2.0+mv**2.0))
        mean_mf=n.sum(mfs)/len(mfs)/ch_pairs.shape[0]
        wg=mfs/ch_pairs.shape[0]
        wgs=n.sum(wg)
        wmu=n.sum(mu*wg)/wgs
        wmv=n.sum(mv*wg)/wgs
        print("mf %1.2f center %1.2f %1.2f phcal %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f"%(mean_mf,wmu,wmv,x[0],x[1],x[2],x[3],x[4],x[5]))
        s=-mean_mf#+wmu**2.0+wmv**2.0
        print(s)
        return(s)
    print("optimize")
    # find phases
    phcal[1:7]=so.fmin(ss,phcal[1:7])
    # remove outliers
    mu,mv,mfs,mis,mjs=image_points(phcal,xc[ridx,:],ch_pairs,u,v,w,dmat)
    gidx=n.where(mfs>0.9)[0]
    ridx=ridx[gidx]
    phcal[1:7]=so.fmin(ss,phcal[1:7])
    print(len(ridx))

    # pick the center from the plot
    mu,mv,mfs,mis,mjs=image_points(phcal,xc2[:,:],ch_pairs,u,v,w,dmat)
    plt.scatter(mu,mv,c=mfs/21.0,s=2)
    cb=plt.colorbar()
    cb.set_label("Coherence")
    plt.xlabel("u")
    plt.ylabel("v")
    plt.plot(0,0,"x")
    plt.show()

get_xcs2()
exit(0)
#find_initial_phasecal()
#exit(0)



# create a uniform grid of u,v,w values
u,v,w=pint.uv_coverage(N=200,max_zenith_angle=10.0)
# antenna positions
antpos=pint.get_antpos()
# channel pairs used for interferometry
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
# matrix of antenna direction vector pairs 
# one row is a vector for each interferometry channel pair
dmat=pint.pair_mat(ch_pairs[:,:],antpos)



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

# 
# 1. optimize for a phase calibration that provides a good match function
# 2. rotate to correct pointing direction
# 3. optimize for good match function
#

# tbd, add code to spit out the initial phasecal
# phase calibration guess
xhat=[ 1.14248621,  0.15371235,  2.5528059,  -0.35319514,  3.19109113,  3.58991391]
phcal=n.zeros(7,dtype=n.float32)
phcal[1:len(phcal)]=xhat
#mu,mv,mfs,mis,mjs=image_points(phcal,xc,ch_pairs,u,v,w,dmat)



phcal=direction_shift(dir_orig=n.array([0.05,-0.008,n.sqrt(1-0.05**2-0.008**2)]),
                      dir_new=n.array([0,0,-1]),
                      phcal_orig=phcal,
                      dmat=dmat,
                      ch_pairs=ch_pairs)
import time
ho=h5py.File("data/phases.h5","w")
ho["phasecal"]=phcal
ho["beamid"]=0
ho["creation_time"]=time.time()
ho.close()

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
