
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
import pansy_config as pc
import glob
import h5py 

def fit_obs(tx_idx,rg,dop,fit_acc=False,return_model=False):
    """
    given range and doppler measurements, find best fit radial trajectory
    """
    # in some cases, there is a lot of PMSE and the doppler sidelobes of these echoes
    # will jam the processing. Limit to maximum 200 of the largest doppler shift detections
    n_m=len(tx_idx)
    
    dur=(n.max(tx_idx)-n.min(tx_idx))/1e6
    tmean=n.mean(tx_idx)
    t=(tx_idx-tmean)/1e6
    r_std=0.5
    v_std=3

#    if dur > 
    if dur < 0.05:
        A=n.zeros([2*n_m,2],dtype=n.float32)
#        AA=n.zeros([2*n_m,2],dtype=n.float32)

        A[0:n_m,0]=(1.0)/r_std
        A[0:n_m,1]=(t)/r_std
        A[n_m:(2*n_m),1]=1/v_std
    if dur < 0.25:
        A=n.zeros([2*n_m,3],dtype=n.float32)
        A[0:n_m,0]=(1.0)/r_std
        A[0:n_m,1]=(t)/r_std
        A[0:n_m,2]=(t**2.0)/r_std
        A[n_m:(2*n_m),1]=1/v_std
        A[n_m:(2*n_m),2]=(2*t)/v_std
    else:
        A=n.zeros([2*n_m,4],dtype=n.float32)
        A[0:n_m,0]=(1.0)/r_std
        A[0:n_m,1]=(t)/r_std
        A[0:n_m,2]=(t**2.0)/r_std
        A[0:n_m,3]=(t**3.0)/r_std
        A[n_m:(2*n_m),1]=1/v_std
        A[n_m:(2*n_m),2]=(2*t)/v_std
        A[n_m:(2*n_m),3]=(3*t**2)/v_std
# r = r0+v*t + at**2 + bt**3
# dr/dt = v + 2*a*t + 3*b*t**2

    m=n.zeros(2*n_m,dtype=n.float32)
    m[0:n_m]=rg/r_std
    m[n_m:(2*n_m)]=dop/1e3/v_std

    xhat=n.linalg.lstsq(A,m)[0]
#    print(xhat)
    model=n.dot(A,xhat)
    r_resid=rg-r_std*model[0:n_m]
    dop_resid=dop/1e3-v_std*model[n_m:(2*n_m)]

    fit_v_std=n.std(dop_resid)
    fit_r_std=n.std(r_resid)
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
    if return_model:
        import scipy.interpolate as sint
        xhatp=n.zeros(4)
        xhatp[0:len(xhat)]=xhat
        def rmodel(tt):
            tp=(tt-tmean)/1e6
            return(xhatp[0]+xhatp[1]*tp+xhatp[2]*tp**2.0+xhatp[3]*tp**3.0)
        def vmodel(tt):
            tp=(tt-tmean)/1e6
            return(xhatp[1]+2*xhatp[2]*tp+3*xhatp[3]*tp**2.0)

#        rmodel=sint.interp1d(t,r_std*model[0:n_m])
 #       vmodel=sint.interp1d(t,v_std*model[n_m:(2*n_m)])
        return(fit_r_std,fit_v_std,xhat,tmean,rmodel,vmodel)
    else:
        return(fit_r_std,fit_v_std,xhat,tmean)

def cluster(tx_idx,
            rg,
            dop,
            snr,
            min_dur=0.06,
            min_det=6,
            plot=True
            ):

    meteor_detections=[]

    idx=n.argsort(snr)[::-1]
    pairs=[]
    pair_time=[]
    tv=(tx_idx-tx_idx[0])/1e6

    # first pass
    # look for measurement pair that fits together with each measurement
    # has to be less than 10 ms apart and fit better than 500 meters
    # also, doppler can't change more than 15 km/s
    # things that don't git together are not considered in the next stage
    while len(idx)>1:
        i = idx[0]
        dt = (tx_idx[idx]-tx_idx[i])/1e6
        dr = rg[i]-(rg[idx]-(dop[i]/1e3)*dt)
        ddop = dop[i]-dop[idx]
        #print(dr)
        # has to be close in time, but also close enough in range and doppler
        fit_idx=n.where( (n.abs(dt) < 10e-3) & (n.abs(dr)<2.0) & (n.abs(ddop)<13e3))[0]
        pair_idx=idx[fit_idx]
        if len(pair_idx) > 1:
#            print(pair_idx)
            pairs.append(pair_idx)
            pair_time.append(n.mean(tv[pair_idx]))
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
            pass
#            print("no fit. removing")
#            print(pair_idx)
            
        idx=n.setdiff1d(idx,pair_idx)
    pair_time=n.array(pair_time)  

    if len(pairs)>0:
        # we could pair some measurements. now see we can merge any of them

        # second pass.
        # find all the 2-combinations of measurement pairs 
        # these are ones that we will try to expand upon later
        candidates=list(itertools.combinations(pairs,2))
        used_idx=[]
        tuples=[]
        for c in candidates:
            t0=n.mean(tx_idx[c[0]])/1e6
            t1=n.mean(tx_idx[c[1]])/1e6
#            print(c[0])
 #           print(c[1])
  #          print(n.intersect1d(used_idx,c[0]))
   #         print(n.intersect1d(used_idx,c[1]))
            
            if (len(n.intersect1d(used_idx,c[0])) == 0) and (len(n.intersect1d(used_idx,c[1])) == 0):
                # at most 5*2*1.6e-3 apart to try merging
                max_dt=3*5*1.6e-3
                if n.abs(t1-t0)<max_dt:
                    try_idx=n.concatenate((c[0],c[1]))
#                    print(try_idx)
                    r_resid, v_resid, xhat, tmean=fit_obs(tx_idx[try_idx],rg[try_idx],dop[try_idx],fit_acc=False)
                    if (r_resid < 0.5) and (v_resid < 3):
#                        print("merging")
                        used_idx=n.concatenate((used_idx,try_idx))
                        tuples.append(try_idx)
                    else:
                        pass
#                        print("not merging")
 #                       print(r_resid)
  #                      print(v_resid)
            else:
                pass
                #print(used_idx)
                #print("already used. skipping")
                
        # third pass.
        # go through all tuple groups. try to add measurements to them. 
        # always add the best measurement, but make sure it is not too far away in time.
        #   
        # Here is one idea
        # go through each tuple list
        # keep adding the best fitting measurement pair into tuple as long as the observations fit
        # once done, call it a cluster
        #
        used_idx=n.array([],dtype=n.int64)
        tuples2=[]
        for tup in tuples:
            # remove already used indices
            idx_this = n.setdiff1d(tup,used_idx)
            if len(idx_this) < 4:
#                print("not enough measurements in tuple. skipping")
                used_idx=n.concatenate((used_idx,idx_this))
                continue

            t0_this=n.mean(tv[idx_this])
            # find the closest one in time to add
            dtidx=n.argsort(n.abs(pair_time-t0_this))

            # try adding each measurement once
            for pi in range(len(pairs)):
                p=pairs[dtidx[pi]]
                # only try adding if this pair is not already in tuple
                if len(n.intersect1d(idx_this,p))==0:
                    # if we are close enough in time
                    if n.min(n.abs(n.mean(tv[p])-tv[idx_this])) < 0.15:
                        try_idx=n.concatenate((idx_this,p))
                        r_resid, v_resid, xhat, tmean=fit_obs(tx_idx[try_idx],rg[try_idx],dop[try_idx],fit_acc=True)
                        if (r_resid < 0.5) and (v_resid < 7):
                            # adding
#                            print("added pair %1.2f %1.2f"%(r_resid,v_resid))
                            idx_this=try_idx
                        else:
                            pass
#                            print("not adding %1.2f %1.2f"%(r_resid,v_resid))
            dur=n.max(tv[idx_this])-n.min(tv[idx_this])
            if dur>min_dur and len(idx_this)>min_det:
                used_idx=n.concatenate((used_idx,idx_this))
                tuples2.append(idx_this)

        # for each cluster, find measurements that fit        
        if plot:
            plt.plot(tv,rg,".",color="gray")
            plt.title("%s %d %d"%(stuffr.unix2datestr(tx_idx[0]/1e6),len(tx_idx),2*len(pairs)))
            for p in pairs:
                plt.plot(tv[p],rg[p],".",color="red")

        for p in tuples2:
            r_resid, v_resid, xhat, tmean, rmodel,vmodel=fit_obs(tx_idx[p],rg[p],dop[p],fit_acc=True,return_model=True)
            if plot:
                plt.plot(tv[p],rg[p],".",color="green")
            #ps=n.sort(p)
            tvlocal=(tx_idx[p]-tmean)/1e6
            rgmodel=rmodel(tvlocal)#xhat[0]+xhat[1]*tvlocal+0.5*xhat[2]*tvlocal**2.0
            if plot:
                plt.plot(tv[p],rgmodel,".",color="blue")
            dur=n.max(tvlocal)-n.min(tvlocal)
            if plot:
                plt.axvline(n.min(tv[p]),color="green")
                plt.axvline(n.max(tv[p]),color="green")
    #                print(xhat)
                plt.text(n.min(tv[p]),xhat[0],"%d %1.2f s\n%1.1f km/s\n%1.1f km/s2\nfr %1.2f"%(len(p),dur,xhat[1],xhat[2],len(p)/dur))

            meteor_data = {"idx":p,"xhat":xhat}
            meteor_detections.append(meteor_data)
            
        if plot:
            plt.ylim([60,140])
            plt.show()
    return(meteor_detections)
    
        
        
        
        

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
        pass
#        print("no data")
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
                pass
#                print("low txpower. skipping")

        gidx=n.where(snra>snr_threshold)[0]
        txpa=txpa[gidx]
        txidxa=txidxa[gidx]
        rnga=rnga[gidx]
        dopa=dopa[gidx]
        snra=snra[gidx]
        beam=beam[gidx]
    return(txpa,txidxa,rnga,dopa,snra,beam)

def find_clusters(txpa,txidxa,rnga,dopa,snra,beam,dmw):
    # at most this many mf outputs
    max_tx_idx=300
    # Any count per range gate exceeding this is counted as PMSE
    pmse_threshold=100
    # doppler shift has to be at least 3 km/s to be analyzed
    dop_thresh=3e3
    
    gidx=n.where(n.abs(dopa)>dop_thresh)[0]
    txpa=txpa[gidx]
    txidxa=txidxa[gidx]
    rnga=rnga[gidx]
    dopa=dopa[gidx]
    snra=snra[gidx]
    beam=beam[gidx]
    if len(txpa)<5:
        #print("not enough data")
        return

    # filtering PMSE based on histogram statistics
    rbins=n.linspace(60,140,num=80)
    h,be=n.histogram(rnga,rbins)
    pmse_rg=n.where(h>100)[0]
    pmse_drg=2.0
    for ri in pmse_rg:
        print("pmse at %1.2f km"%(rbins[ri]))
        pmse_r=rbins[ri]
        gidx=n.where( ((rnga > (pmse_r+pmse_drg))) | ((rnga < (pmse_r-pmse_drg))) )[0]
        txpa=txpa[gidx]
        txidxa=txidxa[gidx]
        rnga=rnga[gidx]
        dopa=dopa[gidx]
        snra=snra[gidx]
        beam=beam[gidx]
        
    n_m=len(txpa)

    if n_m > max_tx_idx:
        #print("too many detections. limiting to %d highest doppler shift detections "%(max_tx_idx))
        dop_idx=n.argsort(n.abs(dopa))[::-1]
        txpa=txpa[dop_idx[0:max_tx_idx]]
        txidxa=txidxa[dop_idx[0:max_tx_idx]]
        rnga=rnga[dop_idx[0:max_tx_idx]]
        dopa=dopa[dop_idx[0:max_tx_idx]]
        snra=snra[dop_idx[0:max_tx_idx]]
        beam=beam[dop_idx[0:max_tx_idx]] 
        n_m=max_tx_idx

    if len(txidxa) > 5:
        meteor_detections=cluster(txidxa,rnga,dopa,snra)
        for det in meteor_detections:
            xhat=det["xhat"]
            print("%s meteor %d detections v=%1.1f km/s r=%1.1f km"%(stuffr.unix2datestr(n.min(txidxa[det["idx"]])/1e6),len(det["idx"]),xhat[1],xhat[0]))

            odata_dict={}
            odata_dict["xhat"]=det["xhat"]
            odata_dict["tx_idx"]=txidxa[det["idx"]]
            odata_dict["range"]=rnga[det["idx"]]
            odata_dict["doppler"]=dopa[det["idx"]]
            odata_dict["snr"]=snra[det["idx"]]
            odata_dict["beam"]=beam[det["idx"]]
            idx0=n.min(txidxa[det["idx"]])
            try:
                dmw.write([idx0],odata_dict)
            except:
                import traceback
                traceback.print_exc()

def analyze_until_now():
    # this is where the existing metadata lives

    det_md_dir=pc.detections_metadata_dir
    #det_md_dir = "/media/archive/metadata/detections"
    b_det=[-1,-1]
    det_md=None
    analysis_end=-1
    if os.path.exists(det_md_dir):
        print("metadata directory exists. searching for last timestamp")
        try:
            det_md = drf.DigitalMetadataReader(det_md_dir)
            b_det = det_md.get_bounds()
            print(b_det)
        except:
            import traceback
            traceback.print_exc()
            print("couldn't read det metadata")

        try:
            # look for mf search output to determine where it has reached. only analyze that far
            fl=glob.glob("/tmp/meteor_mf_*.h5")
            latest_idx=[]
            for f in fl:
                h=h5py.File(f,"r")
                print(f)
                latest_idx.append(h["latest"][()])
                h.close()
            analysis_end=n.min(latest_idx)
        except:
            import traceback
            traceback.print_exc()
            print("couldn't read det metadata")
    else:
        os.system("mkdir -p %s"%(det_md_dir))
    print("%s"%(stuffr.unix2datestr(analysis_end/1e6)))
#    exit(0)
    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 600
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "det"

    dmw = drf.DigitalMetadataWriter(
        det_md_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )

    mf_metadata_dir=pc.mf_metadata_dir
    #mf_metadata_dir = "/media/archive/metadata/mf"
    dm_mf = drf.DigitalMetadataReader(mf_metadata_dir)
    db_mf = dm_mf.get_bounds()

    dt=10000000
    #n_min=int(n.floor((db_mf[1]-db_mf[0])/dt))
    d=drf.DigitalRFReader("/media/archive/")
    # tx channel bounds
    b=d.get_bounds("ch007")
    #start_idx=db_mf[1]-2*60*60*1000000

    # start in the beginning
    start_idx=dt*int(n.floor(db_mf[0]/dt))
    if b_det[1] != -1:
        # start where detections end, if we have already analyzed
        start_idx=dt*int(n.ceil(b_det[1]/dt))

    # only go up to the point where the slowest thread is at (analysis_end)
    if analysis_end == -1:
        # don't know how long to analyze!
        print("don't know how to analyze. quick_search_meteor.py not running?")
        exit(0)
        analysis_end=db_mf[1]

    n_min=int(n.floor((analysis_end-start_idx)/dt))
    for i in range(n_min):
        i0=start_idx+i*dt
        i1=start_idx+i*dt+dt

        txpa,txidxa,rnga,dopa,snra,beam=read_mf_output(dm_mf,i0,i1)
        find_clusters(txpa,txidxa,rnga,dopa,snra,beam,dmw)
        
if __name__ == "__main__":
    while True:
        analyze_until_now()
        time.sleep(3600)