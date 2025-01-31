
import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import scipy.fftpack as fp
import itertools

def fit_obs(tx_idx,rg,dop):
    n_m=len(tx_idx)
    A=n.zeros([2*n_m,3],dtype=n.float32)
    tmean=n.mean(tx_idx)
    t=(tx_idx-tmean)/1e6
    r_std=0.5
    v_std=3
    A[0:n_m,0]=(1.0)/r_std
    A[0:n_m,1]=(t)/r_std
    A[0:n_m,2]=(t**2.0)/r_std
    A[n_m:(2*n_m),1]=1/v_std
    A[n_m:(2*n_m),2]=(2*t)/v_std
    
    m=n.zeros(2*n_m,dtype=n.float32)
    m[0:n_m]=rg/r_std
    m[n_m:(2*n_m)]=dop/1e3/v_std
    xhat=n.linalg.lstsq(A,m)[0]
    print(xhat)
    model=n.dot(A,xhat)
    r_resid=rg-r_std*model[0:n_m]
    dop_resid=dop/1e3-v_std*model[n_m:(2*n_m)]

    v_std=n.std(dop_resid)
    r_std=n.std(r_resid)
    if False:
        plt.subplot(121)
        plt.plot(t,r_std*model[0:n_m])
        plt.title(r_std)
        plt.plot(t,rg,".")
        plt.subplot(122)
        plt.plot(t,v_std*model[n_m:(2*n_m)])
        plt.title(v_std)    
        plt.plot(t,dop/1e3,".")
        plt.show()
    return(r_std,v_std,xhat)

def cluster(tx_idx,
            rg,
            dop,
            snr,
            min_det=5
            ):
    idx=n.argsort(snr)[::-1]
    pairs=[]
    tv=(tx_idx-tx_idx[0])/1e6
    # first pass
    # look for measurement pair that best fits this measurement
    # has to be less than 10 ms apart and fit better than 500 meters
    # also, doppler can't change more than 15 km/s
    while len(idx)>1:
        i = idx[0]
        dt = (tx_idx[idx]-tx_idx[i])/1e6
        dr = rg[i]-(rg[idx]-dop[idx]*dt)
        ddop = dop[i]-dop[idx]

        fit_idx=n.where( (n.abs(dt) < 10e-3) & (n.abs(dr)<500) & (n.abs(ddop)<15e3))[0]
        pair_idx=idx[fit_idx]
        if len(pair_idx) > 1:
            print(pair_idx)
            pairs.append(pair_idx)
            if False:
                plt.subplot(131)
                plt.plot(dt[fit_idx],dr[fit_idx],".")
                plt.subplot(132)
                plt.plot(tx_idx[pair_idx],rg[pair_idx],".")
                plt.subplot(133)
                plt.plot(tx_idx[pair_idx],dop[pair_idx],".")
                plt.tight_layout()
                plt.show()
        else:
            print("no fit. removing")
            print(pair_idx)
            
        idx=n.setdiff1d(idx,pair_idx)
    if len(pairs)>0:
        #    pairs=n.array(pairs,dtype=n.int64)
        
        plt.subplot(121)
        for p in pairs:
            plt.plot(tv[p],rg[p],".")
        plt.subplot(122)
        for p in pairs:
            plt.plot(tv[p],dop[p]/1e3,".")
        plt.show()

        # second pass.
        # find all the 2-combinations of measurement pairs 
        candidates=list(itertools.combinations(pairs,2))
        used_idx=[]
        tuples=[]
        for c in candidates:
            t0=n.mean(tx_idx[c[0]])/1e6
            t1=n.mean(tx_idx[c[1]])/1e6
            print(c[0])
            print(c[1])
            print(n.intersect1d(used_idx,c[0]))
            print(n.intersect1d(used_idx,c[1]))
            
            if (len(n.intersect1d(used_idx,c[0])) == 0) and (len(n.intersect1d(used_idx,c[1])) == 0):
                # at most 5*4*1.6e-3 apart to try merging
                max_dt=5*8*1.6e-3
                if n.abs(t1-t0)<max_dt:
                    try_idx=n.concatenate((c[0],c[1]))
                    print(try_idx)
                    r_resid, v_resid, xhat=fit_obs(tx_idx[try_idx],rg[try_idx],dop[try_idx])
                    if (r_resid < 0.5) and (v_resid < 2):
                        print("merging")
                        used_idx=n.concatenate((used_idx,try_idx))
                        tuples.append(try_idx)
                    else:
                        print("not merging")
                        print(r_resid)
                        print(v_resid)
            else:
                print(used_idx)
                print("already used. skipping")


        plt.subplot(121)
        for p in tuples:
            plt.plot(tv[p],rg[p],".")
        plt.title("Tuples")
        plt.subplot(122)
        for p in tuples:
            plt.plot(tv[p],dop[p]/1e3,".")
        plt.show()
                
        # third pass.
        # find all the 2-combinations of measurement merged tuples
        candidates=list(itertools.combinations(tuples,2))
        used_idx=[]
        tuples2=[]
        for c in candidates:
            t0=n.mean(tx_idx[c[0]])/1e6
            t1=n.mean(tx_idx[c[1]])/1e6

            if (len(n.intersect1d(used_idx,c[0])) == 0) and (len(n.intersect1d(used_idx,c[1]))==0):
                # at most 5*8*1.6e-3 apart to try merging
                max_dt=5*16*1.6e-3
                if n.abs(t1-t0)< max_dt:
                    try_idx=n.concatenate((c[0],c[1]))
                    print(try_idx)
                    r_resid, v_resid, xhat=fit_obs(tx_idx[try_idx],rg[try_idx],dop[try_idx])
                    if (r_resid < 0.5) and  (v_resid < 2):
                        print("merging")
                        used_idx=n.concatenate((used_idx,try_idx))
                        tuples2.append(try_idx)
            else:
                print("already used. skipping")
                
                
        plt.subplot(121)
        plt.title("Tuples 2")        
        for p in tuples2:
            plt.plot(tv[p],rg[p],".")
        plt.subplot(122)
        for p in tuples2:
            plt.plot(tv[p],dop[p]/1e3,".")
        plt.show()


        # Fourth pass.
        # find all the 2-combinations of measurement merged tuples
        candidates=list(itertools.combinations(tuples2,2))
        used_idx=[]
        tuples3=[]
        for c in candidates:
            t0=n.mean(tx_idx[c[0]])/1e6
            t1=n.mean(tx_idx[c[1]])/1e6

            if (len(n.intersect1d(used_idx,c[0])) == 0) and (len(n.intersect1d(used_idx,c[1]))==0):
                # at most 5*8*1.6e-3 apart to try merging
                max_dt=5*32*1.6e-3
                if n.abs(t1-t0)< max_dt:
                    try_idx=n.concatenate((c[0],c[1]))
                    print(try_idx)
                    r_resid, v_resid, xhat=fit_obs(tx_idx[try_idx],rg[try_idx],dop[try_idx])
                    if (r_resid < 0.5) and  (v_resid < 3):
                        print("merging")
                        used_idx=n.concatenate((used_idx,try_idx))
                        tuples3.append(try_idx)
            else:
                print("already used. skipping")
                
                
        plt.subplot(121)
        plt.title("Tuples 3")        
        for p in tuples3:
            plt.plot(tv[p],rg[p],".")
        plt.subplot(122)
        for p in tuples3:
            plt.plot(tv[p],dop[p]/1e3,".")
        plt.show()
        
    
        # now we have everything. look for all measurements that fit
        
        
        
        

def read_mf_output(dm_mf,i0,i1,snr_threshold=7,tx_pwr_threshold=1e9):
    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)
    snra=n.array([],dtype=n.float32)
    beam=n.array([],dtype=n.int32)         
    # read 20 pulse sequences   
    try:
        data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))
    except:
        import traceback
        traceback.print_exc()
    

    if len(data_dict.keys()) == 0:
        print("no data")
    else:
        for k in data_dict.keys():
            data=data_dict[k]
            
            txp=data["tx_pwr"]
#            print(txp)
            # use this threshold for tx power. check that it is okay!!!
            if n.min(txp) > tx_pwr_threshold:
                txpa=n.concatenate((txpa,txp))
                txidxa=n.concatenate((txidxa,data["tx_idxs"]))
                rnga=n.concatenate((rnga,data["max_range"]))
                dopa=n.concatenate((dopa,data["max_dopvel"]))
                snra=n.concatenate((snra,data["max_snr"]))
                beam=n.concatenate((beam,data["beam_pos_idx"]))            
            else:
                print("low txpower. skipping")

        gidx=n.where(snra>snr_threshold)[0]
        txpa=txpa[gidx]
        txidxa=txidxa[gidx]
        rnga=rnga[gidx]
        dopa=dopa[gidx]
        snra=snra[gidx]
        beam=beam[gidx]
    return(txpa,txidxa,rnga,dopa,snra,beam)

mf_metadata_dir = "/media/archive/metadata/mf"
dm_mf = drf.DigitalMetadataReader(mf_metadata_dir)
db_mf = dm_mf.get_bounds()

dt=10000000
#n_min=int(n.floor((db_mf[1]-db_mf[0])/dt))
d=drf.DigitalRFReader("/media/archive/")
# tx channel bounds
b=d.get_bounds("ch007")
#start_idx=b[0]#db_mf[1]-2*60*60*1000000
start_idx=dt*int(n.floor(db_mf[0]/dt))#-2*60*60*1000000
#start_idx=db_mf[0]
n_min=int(n.floor((db_mf[1]-start_idx)/dt))
for i in range(n_min):
    i0=start_idx+i*dt
    i1=start_idx+i*dt+dt 
    txpa,txidxa,rnga,dopa,snra,beam=read_mf_output(dm_mf,i0,i1)

    gidx=n.where(n.abs(dopa)>3e3)[0]
    txpa=txpa[gidx]
    txidxa=txidxa[gidx]
    rnga=rnga[gidx]
    dopa=dopa[gidx]
    snra=snra[gidx]
    beam=beam[gidx]
    if len(txpa)<5:
        print("not enough data")
        continue
    
    cluster_idx=cluster(txidxa,rnga,dopa,snra)
    
