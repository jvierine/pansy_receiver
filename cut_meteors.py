import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import cluster_mf as cmf
import traceback

os.system("mkdir -p %s"%(pc.cut_metadata_dir))
# transmit pulse metadata
dmt = drf.DigitalMetadataReader(pc.tx_metadata_dir)
bt = dmt.get_bounds()

# meteor detections
dm = drf.DigitalMetadataReader(pc.detections_metadata_dir)
bm = dm.get_bounds()

# match function outputs
dmf = drf.DigitalMetadataReader(pc.mf_metadata_dir)
bmf = dmf.get_bounds()

subdirectory_cadence_seconds = 3600
file_cadence_seconds = 60
samples_per_second_numerator = 1000000
samples_per_second_denominator = 1
file_name = "cut"

dmw = drf.DigitalMetadataWriter(
    pc.cut_metadata_dir,
    subdirectory_cadence_seconds,
    file_cadence_seconds,
    samples_per_second_numerator,
    samples_per_second_denominator,
    file_name,
)

# raw voltage
d=drf.DigitalRFReader(pc.raw_voltage_dir)
b=d.get_bounds("ch007")

def cut_raw_voltage(i0,i1,rmodel,n_pad=100000,beams=[0],rx_ch=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],tx_ch="ch007",txlen=132,
                    plot=False,
                    pad=64):

    # one sample in range (km)
    drg=c.c/1e6/2/1e3    
    n_ch=len(rx_ch)

    tx_data_dict = dmt.read(i0, i1, "id")

    # how many 20 ipp segments do we have
    n_seq=len(tx_data_dict.keys())
    if n_seq==0:
        print("no data found!")
        return({})
    
    # this is the beam number
    beam_idx=n.array(n.mod(n.arange(20,dtype=n.int64),5),dtype=n.int8)

    zrx=n.zeros([n_ch,20*1600],dtype=n.complex64)


    #z_tx=n.zeros([20,txlen],dtype=n.complex64)
    txidx=[]
    ztx_pulses_re=[]
    ztx_pulses_im=[]

    zrx_echoes_re=[]
    zrx_echoes_im=[]

    delays=[]
    beam_ids=[]
    kl=n.sort(list(tx_data_dict.keys()))
    print(kl)
    for key in kl:
#        print(key)
        # create new 20 ipps
        zrx_re=n.zeros([n_ch,20*1600],dtype=n.int16)
        zrx_im=n.zeros([n_ch,20*1600],dtype=n.int16)
        ztx_re=n.zeros(20*1600,dtype=n.int16)
        ztx_im=n.zeros(20*1600,dtype=n.int16)

        for chi in range(len(rx_ch)):
            zrx[chi,:]=d.read_vector_c81d(key,1600*20,rx_ch[chi])
        zrx_re[:,:]=n.array(zrx.real,dtype=n.int16)
        zrx_im[:,:]=n.array(zrx.imag,dtype=n.int16)
        z_tx=d.read_vector_c81d(key,1600*20,tx_ch)
        ztx_re[:]=n.array(z_tx.real,dtype=n.int16)
        ztx_im[:]=n.array(z_tx.imag,dtype=n.int16)
#        plt.plot(ztx_re[0:132])
 #       plt.plot(ztx_im[0:132])
  #      plt.show()
        for i in range(20):
            this_beam_idx=i%5
            delay=int(n.round(rmodel(key+i*1600)/drg))
            if this_beam_idx in beams:
                # index of transmit pulse start
                txidx.append(key+i*1600)
                # this is the expected range delay for the meteor
                delays.append(delay-pad)
                         
                ztx_pulses_re.append(ztx_re[(i*1600):(i*1600+txlen)])
                ztx_pulses_im.append(ztx_im[(i*1600):(i*1600+txlen)])

                # echo matrix for all channels
                zrx_echoes_re.append(zrx_re[:,(i*1600+delay-pad):(i*1600+txlen+delay+pad)])
                zrx_echoes_im.append(zrx_im[:,(i*1600+delay-pad):(i*1600+txlen+delay+pad)])

                beam_ids.append(this_beam_idx)
    #plot=True
    if plot:
        n_ipp = len(delays)
        RTI=n.zeros([n_ipp,1600],dtype=n.float32)
        TXI=n.zeros([n_ipp,txlen],dtype=n.complex64)
        for i in range(n_ipp):
            TXI[i,:]=ztx_pulses_re[i]+ztx_pulses_im[i]*1j
            rx_re=zrx_echoes_re[i]
            rx_im=zrx_echoes_im[i]

            for ci in range(len(rx_ch)):
                RTI[i,delays[i]:(delays[i]+2*pad+txlen)]+=n.array(rx_re[ci,:],dtype=n.float32)**2+n.array(rx_im[ci,:],dtype=n.float32)**2.0
        noise=n.median(RTI[RTI!=0.0])
        plt.pcolormesh(10.0*n.log10(RTI.T+1),vmin=10.0*n.log10(noise))
        plt.title("%s beams %s"%(stuffr.unix2datestr(i0/1e6),str(beams)))
        plt.ylim([n.min(delays),n.max(delays)+txlen+2*pad])
        plt.xlabel("IPP")
        plt.ylabel("Range-gate")
        plt.colorbar()
        plt.show()

        #plt.pcolormesh(n.angle(TXI.T))
        #plt.colorbar()
        #plt.show()

    return({"tx_idx":[txidx], # these are the indices of the transmit pulses
            "beam_id":[beam_ids],  # these are the beam directions
            "ztx_pulses_re":[ztx_pulses_re], # 16-bit int real tx
            "ztx_pulses_im":[ztx_pulses_im], # 16-bit int imag tx
            "zrx_echoes_re":[zrx_echoes_re], # 16-bit int real echo
            "zrx_echoes_im":[zrx_echoes_im], # 16-bit int imag eho 
            "delays":[delays], # delay between tx_idx and echo start sample (tx_idx+delays)
            "pad":[pad], # how much did we pad the echo in addition to txlen
            "txlen":[txlen], # what is the tx pulse length
            "channels":[rx_ch], # what channels are stored
            "tx_channel":[tx_ch]}) # what is the tx channel name

dt=60
sr=1000000

n_block=int(n.floor((bm[1]-bm[0])/dt))

# go through all blocks of data
for i in range(n_block):
    i0=bm[0]+i*dt*sr
    i1=bm[0]+i*dt*sr + dt*sr
    
    # find all detections in interval
    data_dict = dm.read(i0, i1, ("tx_idx","xhat","range","doppler","snr","beam"))

    # for each detection in this window
    for k in data_dict.keys():
        data=data_dict[k]

        # radial trajectory model
        xhat=data["xhat"]
        # timestamps of transmit pulses
        tx_idx=data["tx_idx"]
        # detection range
        range_km=data["range"]
        # signal to noise ratio
        snr=data["snr"]
        # transmit beam number
        beam_idx=data["beam"]
        # doppler shift
        doppler_ms=data["doppler"]

        # sort detection related measurements in time
        idxidx=n.argsort(tx_idx)
        tx_idx=tx_idx[idxidx]
        snr=snr[idxidx]
        beam_idx=beam_idx[idxidx]
        range_km=range_km[idxidx]
        doppler_ms=doppler_ms[idxidx]

        # get model for range
        fit_r_std,fit_v_std,xhat,tmean,rmodel,vmodel=cmf.fit_obs(tx_idx,range_km,doppler_ms,return_model=True)

        # which beams are included. we want at least 5 detections in each beam for it to count
        beamnum=n.array(n.mod(beam_idx,5),dtype=int)
        beam_count=n.zeros(5,dtype=int)
        beam_i0=n.zeros(5,dtype=int)
        beam_i1=n.zeros(5,dtype=int)
        beams=[]
        for bi in range(5):
            beam_count[bi]=len(n.where(beamnum==bi)[0])
            if beam_count[bi] >= 5:
                beams.append(bi)
                beam_i0[bi]=n.min(tx_idx[n.where(beamnum==bi)[0]])
                beam_i1[bi]=n.max(tx_idx[n.where(beamnum==bi)[0]])
                print("beam %d i0 %d i1 %d count %d"%(bi,beam_i0[bi],beam_i1[bi],beam_count[bi]))
        
        # store between i0 and i1 plus padding ch000-ch006
        # store +/- 64 samples around the tx pulse of length 128 samples
        # pad by X samples
        cut_res=cut_raw_voltage(n.min(tx_idx)-61*1600,n.max(tx_idx)+61*1600,
                                rmodel,
                                beams=beams,
                                rx_ch=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],
                                tx_ch="ch007",
                                txlen=132,
                                pad=64,
                                plot=False
                                )
        cut_res["c_snr"]=[snr]
        cut_res["c_tx_idx"]=[tx_idx]
        cut_res["c_beam_idx"]=[beam_idx]
        cut_res["c_doppler"]=[doppler_ms]
        cut_res["c_range_km"]=[range_km]
        okey=n.min(tx_idx)
        try:
            dmw.write([okey],cut_res)
        except:
            traceback.print_exc()

#    return({"tx_idx":txidx, # these are the indices of the transmit pulses
#            "beam_id":beam_ids,  # these are the beam directions
#            "ztx_pulses_re":ztx_pulses_re, # 16-bit int real tx
#            "ztx_pulses_im":ztx_pulses_im, # 16-bit int imag tx
#            "zrx_echoes_re":zrx_echoes_re, # 16-bit int real echo
#            "zrx_echoes_im":zrx_echoes_im, # 16-bit int imag eho 
#            "delays":delays, # delay between tx_idx and echo start sample (tx_idx+delays)
#            "pad":pad, # how much did we pad the echo in addition to txlen
#            "txlen":txlen, # what is the tx pulse length
#            "channels":rx_ch, # what channels are stored
#            "tx_channel":tx_ch}) # what is the tx channel name

        if False:

            # number of ipps 
            n_ipp=len(tx_idx)
            RTI=n.zeros([n_ipp,256],dtype=n.float32)
            chs=["ch000","ch002","ch003","ch004","ch005","ch006"]
            drg=c.c/1e6/2/1e3
            w=n.repeat(1/4,4)
            w2=n.repeat(1/16,16)
            for ipp in range(n_ipp):
                rg=int(range_km[ipp]/drg)
                for ch in chs:
                    z=n.convolve(d.read_vector_c81d(tx_idx[ipp]+rg-64,256,ch),w,mode="same")
                    RTI[ipp,:]+=n.convolve(n.abs(z)**2.0,w2,mode="same")

            t0=n.min(tx_idx)/1e6

            if n.max(snr)>10:
                plt.subplot(221)
                plt.plot(tx_idx/1e6-t0,range_km,".")
                plt.ylabel("Range (km)")
                plt.xlabel("Time (s)")        
                plt.subplot(222)
                plt.plot(tx_idx/1e6-t0,doppler_ms/1e3,".")
                plt.ylabel("Doppler (km/s)")
                plt.xlabel("Time (s)")
                plt.subplot(223)
                plt.scatter(tx_idx/1e6-t0,10.0*n.log10(snr),s=1,c=n.mod(beam_idx,5),cmap="berlin",vmin=0,vmax=4)
                plt.ylabel("SNR (dB)")
                plt.xlabel("Time (s)")                    
                cb=plt.colorbar(location="top")
                cb.set_label("Beam index")
                plt.subplot(224)
                plt.pcolormesh(tx_idx/1e6-t0,n.arange(256),10.0*n.log10(RTI.T),cmap="plasma")
                plt.colorbar(location="top")
                plt.xlabel("Time (s)")                
                plt.tight_layout()
                plt.show()
        


