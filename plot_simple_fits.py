import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr



dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
b = dm.get_bounds()
dt=100000000
n_block=int((b[1]-b[0])/dt)
hgts=[]
t0s=[]
vgs=[]
els=[]
azs=[]
for i in range(n_block):
    data=dm.read(b[0]+i*dt,b[0]+i*dt + dt)
    for k in data.keys():
#                    dout["r0"]=r0
 #           dout["v0"]=v0            
  #          dout["std"]=[eres,nres,ures]

        r0=data[k]["r0"]
        v0=data[k]["v0"]

        # radiant direction
        v1=-v0
        v_h=n.sqrt(v1[0]**2.0+v1[1]**2.0)
        el=180*n.arctan(n.abs(v1[2])/v_h)/n.pi
        az=180*n.arccos(v1[1]/v_h)/n.pi
        
        std=data[k]["std"]
        if n.max(std*1e3) < 500.0:
            hgts.append(r0[2])
            vgs.append(n.linalg.norm(v0))
            t0s.append(stuffr.unix2date(k/1e6))
            els.append(el)
            azs.append(az)


plt.scatter(t0s,azs,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("Azimuth (deg)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

plt.scatter(t0s,els,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("El (deg)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

            
plt.scatter(t0s,hgts,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("Height (km)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

plt.plot(vgs,hgts,".")
plt.xlabel("Geocentric velocity (km/s)")
plt.ylabel("Height (km)")
plt.show()


