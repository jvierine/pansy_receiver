import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr
import numpy as n

def plot_latest_fits(save_png=False):
    dm = drf.DigitalMetadataReader("/tmp/simple_fit")
    b = dm.get_bounds()

    dd=dm.read(b[1]-24*3600*1000000,b[1])
    vgs=[]
    slats=[]
    slons=[]

    for k in dd.keys():
        vgs.append(n.linalg.norm(dd[k]["v0"]))
        slats.append(dd[k]["eclat"])
        slons.append(dd[k]["eclon"])

    slons=n.array(slons)
    slats=n.array(slats)
    vgs=n.array(vgs)

    ax = plt.subplot(111, projection="lambert")
    sp=ax.scatter(n.angle(n.exp(1j*n.pi*slons/180.0)*n.exp(1j*n.pi/2)),n.pi*slats/180.0,c=vgs,vmin=10,vmax=72,s=2,cmap="turbo")
    cb=plt.colorbar(sp)
    cb.set_label("Geocentric velocity (km/s)")
    ax.grid(True)
    ax.set_xlabel("Apex-centered ecliptic longitude (deg)")
    ax.set_ylabel("Ecliptic latitude (deg)")
    plt.tight_layout()
    if save_png:
        plt.savefig("/tmp/latest_radiants.png")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    plot_latest_fits()