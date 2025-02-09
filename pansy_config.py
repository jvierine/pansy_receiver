import numpy as n

mf_metadata_dir="/media/analysis/metadata/mf"
tx_metadata_dir="/media/analysis/metadata/tx"
detections_metadata_dir="/media/analysis/metadata/detections"
cut_metadata_dir="/media/analysis/metadata/cut"
mesomode_metadata_dir="/media/analysis/metadata/mesomode"

xc_metadata_dir="/media/analysis/metadata/xc2"

raw_voltage_dir="/media/archive/"

# antenna position
lat=-69.010833
lon=39.599722

# parse antenna connections
f=open("cfg/connections.txt","r")
channel_ids=[]
for l in f.readlines():
    channel_ids.append(l.split(",")[0])
f.close()
# parse antenna positions
f=open("cfg/antpos.csv","r")
f.readline()
antenna={}
for l in f.readlines():
    r=l.split(",")
    serial=r[0]
    name=r[1]
    id=r[2]
    ready=r[3]
    x=float(r[4])
    y=float(r[5])
    z=float(r[6])
    antenna[serial]={"serial":serial,"id":id,"name":name,"x":x,"y":y,"z":z,"ready":ready}


if __name__ == "__main__":
    print(channel_ids)
    for cid in channel_ids:
        print(cid)
        if cid != "RFTX":
            print(antenna[cid])