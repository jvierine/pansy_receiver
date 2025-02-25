import numpy as n
import scipy.constants as c

mf_metadata_dir="/media/analysis/metadata/mf"
mf_isr_metadata_dir="/media/analysis/metadata/mf_isr"

tx_metadata_dir="/media/analysis/metadata/tx"
detections_metadata_dir="/media/analysis/metadata/detections"
cut_metadata_dir="/media/analysis/metadata/cut"
mesomode_metadata_dir="/media/analysis/metadata/mesomode"
xc_metadata_dir="/media/analysis/metadata/xc2"
xc_sparse_metadata_dir="/media/analysis/metadata/xc_sparse"
meteor_cal_metadata_dir="/media/analysis/metadata/meteor_cal"
simple_fit_metadata_dir="/media/analysis/metadata/simple_meteor_fit"
gpslock_metadata_dir="/media/analysis/metadata/gpslock"
clock_metadata_dir="/media/analysis/metadata/clock"

raw_voltage_dir="/media/archive/"

# antenna position
lat=-69.010833
lon=39.599722

freq=47.0e6
wavelength=c.c/freq

# parse antenna connections
f=open("cfg/connections.txt","r")
connections=[]
for l in f.readlines():
    connections.append(l.split(",")[0])
f.close()

# parse antenna positions
f=open("cfg/antpos.csv","r")
f.readline()
module_names=[]
for l in f.readlines():
    module_names.append(l.split(",")[1])
f.close()
#print(n.unique(module_names))

f=open("cfg/antpos.csv","r")
f.readline()
antenna={}
modules={}
module_center={}
for mn in module_names:
    modules[mn]=[]

for l in f.readlines():
    r=l.split(",")
    serial=r[0]
    name=r[1]
    id=r[2]
    ready=r[3]
    x=float(r[4])
    y=float(r[5])
    z=float(r[6])
    modules[name].append(n.array([x,y,z]))
    antenna[serial]={"serial":serial,"id":id,"name":name,"x":x,"y":y,"z":z,"ready":ready}
f.close()
for mn in module_names:
    modules[mn]=n.array(modules[mn])
    module_center[mn]=n.mean(modules[mn],axis=0)

def print_antenna():
    import matplotlib.pyplot as plt
    for cid in antenna.keys():
        if cid != "RFTX":
            plt.plot(antenna[cid]["x"],antenna[cid]["y"],".",color="blue")

    for mn in module_names:
        plt.plot(module_center[mn][0],module_center[mn][1],".",color="red")
    for ci,conn in enumerate(connections):
        if conn in module_center.keys():
            t=plt.text(module_center[conn][0],module_center[conn][1],"%d-%s"%(ci,conn),color="black")
            t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='red',linewidth=0))
        
    plt.xlabel("x (m)")
    plt.xlabel("y (m)")
    plt.gca().set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    print_antenna()
