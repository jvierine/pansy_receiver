import h5py
import matplotlib.pyplot as plt
import numpy as n
import pansy_interferometry as pint
import pansy_config as pc


h=h5py.File("/data1/pansy/metadata/xc-1738210528.h5","r")
print(h.keys())
xc=h["xc"][()]
fv=h["fvec"][()]
rv=h["rvec"][()]
ch_pairs=h["ch_pairs"][()]
h.close()

antpos=pint.get_antpos()
dmat=pint.pair_mat(ch_pairs,antpos)
phasecal=pint.get_phasecal()
#print(dmat)
#exit(0)

#u,v,w=pint.uv_coverage1(N=100,max_angle=30,az0=0,el0=90.0)
#idx=pint.mask_uvw(0,90,10,u,v,w)
#print(len(idx))
#exit(0)

def cal_meas(xc,ch_pairs,phasecal):
    #XC=n.zeros([n_xc,n_beams,n_ipp,ipp],dtype=n.complex64)
    coh=n.copy(xc)

    noise_floors=n.zeros([7,5])
    snr=n.zeros([7,5,xc.shape[2],xc.shape[3]])
    for i in range(7):
        for bi in range(5):
            noise_floors[i,bi]=n.mean(n.abs(xc[i,bi,:,150:]))
            snr[i,bi,:,:]=(n.abs(xc[i,bi,:,:])-noise_floors[i,bi])/noise_floors[i,bi]

            if False:
                plt.pcolormesh(10.0*n.log10(snr[i,bi].T),vmin=0,vmax=40)
                plt.colorbar()
                plt.show()
            
    # for each xc
    for i in range(xc.shape[0]):
        for bi in range(xc.shape[1]):
            ch0=ch_pairs[i,0]
            ch1=ch_pairs[i,1]
            coh[i,bi,:,:]=xc[i,bi]*n.exp(-1j*(phasecal[bi,ch0]-phasecal[bi,ch1]))/n.sqrt(n.abs(xc[ch0,bi,:,:])-noise_floors[ch0,bi])/n.sqrt(n.abs(xc[ch1,bi,:,:])-noise_floors[ch1,bi])
            if False:
                plt.pcolormesh(n.abs(coh[i,bi,:,:].T),vmin=0,vmax=1)
                cb=plt.colorbar()
                cb.set_label("Coherence")
                plt.show()
    return(coh,snr)


coherence,snr=cal_meas(xc,ch_pairs,phasecal)

def forward_model(ch_pairs,dmat,u,v,w,idx):
    
    n_meas=ch_pairs.shape[0]
    n_model_par=len(idx)
    A=n.zeros([n_meas,n_model_par],dtype=n.complex64)

    k0=2*n.pi/pc.wavelength
    
    #for i in range(meas.shape[0]):
    #    M+=meas[i]*n.exp(-1j*k0*(dmat[i,0]*u+dmat[i,1]*v+dmat[i,2]*w))
    
    for i in range(n_meas):
        A[i,:]= n.exp(-1j*k0*(u[idx]*dmat[i,0] + v[idx]*dmat[i,1] + w[idx]*dmat[i,2]))
    print(A.shape)
    U,S,Vh=n.linalg.svd(A)
    return(U,S,Vh)

u,v,w=pint.uv_coverage1(N=100,max_angle=30,az0=0,el0=90.0)





idx=pint.mask_uvw(0,90,9,u,v,w)

u0=n.min(u[idx])
u1=n.max(u[idx])
du=n.min(n.unique(n.diff(n.sort(n.unique(u[idx])))))
Npixu=int(n.ceil((u1-u0)/du))

v0=n.min(v[idx])
v1=n.max(v[idx])
dv=n.min(n.unique(n.diff(n.sort(n.unique(v[idx])))))
Npixv=int(n.ceil((v1-v0)/dv))


def u2pix(uu):
    return(n.array(n.round((uu-u0)/du),dtype=n.int64))

ui=u2pix(u[idx])
vi=u2pix(v[idx])
#print(ui)
#exit(0)



U,S,Vh=forward_model(ch_pairs,dmat,u,v,w,idx)
if False:
    print("U")
    print(U.shape)
    print("Vh")
    print(Vh.shape)
    print("S")
    print(len(S))
    plt.plot(S)
    plt.show()
Sinv=n.zeros([Vh.shape[0],U.shape[0]],dtype=n.complex64)
for i in range(18):
    Sinv[i,i]=1/(S[i]+1)
#    plt.scatter(u[idx],v[idx],c=n.abs(Vh[i,:]))
 #   plt.show()

plt.pcolormesh(10.0*n.log10(snr[0,0,:,:].T),vmin=0,vmax=30)
plt.show()
gidx=n.arange(260,310)#n.where(snr[0,bi,:,rg] > 50)[0]
I=n.zeros([Npixu,Npixv,len(gidx)])
for bi in range(1):
    rg=68
    n_m=len(gidx)
    I[:,:,:]=0.0
    for i in range(n_m):
        xhat=n.dot(n.dot(n.dot(n.conj(Vh.T),Sinv),n.conj(n.transpose(U))),coherence[:,bi,gidx[i],rg])
        if snr[0,bi,gidx[i],rg]>5:
            ki=n.argmax(n.abs(xhat))
            #for ki in range(len(idx)):
            I[ui[ki],vi[ki],i]+=snr[0,bi,gidx[i],rg]*n.abs(xhat[ki])

            if True:
                plt.subplot(121)
                plt.pcolormesh(fv,rv,10.0*n.log10(snr[0,bi,:,:].T),vmin=0,vmax=30)
                plt.axhline(rv[rg],color="red",alpha=1.0)
                plt.axvline(fv[gidx[i]],color="red",alpha=1.0)
                plt.xlabel("Doppler (Hz)")
                plt.ylabel("Range (km)")
                plt.colorbar()
                plt.subplot(122)
                plt.scatter(u[idx],v[idx],c=n.abs(xhat),s=100)
                plt.title("dop %d"%(gidx[i]))
                plt.colorbar()
                plt.show()
    
    import rgb_balance as rb
#              ax0=[0,1],
 #             ax0label="par 1",
  #            ax1=[0,1],
   #           ax1label="par 2",              
    #          cax=[0,1],
     #         peak_fraction=1,  # 1 = all 0 = only maximum
      #        cblabel="par 3",    
    rb.rgb_image(I,peak_fraction=1,ax0=[u0,u1],ax1=[u0,u1],cax=[-1,1],cblabel="Doppler",ax0label="u",ax1label="v")
    plt.show()
    plt.pcolormesh(n.sum(I,axis=2))
    plt.colorbar()
    plt.show()
        # ]phasecal depends on beam pointing direction!
        #        z=n.exp(1j*(n.angle(xc[i0,:]) + phcal[beam_id,ch_pairs[:,0]] - phcal[beam_id,ch_pairs[:,1]]))



