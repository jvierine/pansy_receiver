"""Inspect stored grouping results without recomputing the analysis."""
import argparse
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument("input")
p.add_argument("output")
a=p.parse_args()
with h5py.File(a.input) as h:
    t=h["sample_idx"][:]
    r=h["range_km"][:]
    v=h["doppler_m_s"][:]/1000
    labels=[h["legacy_track_id"][:],h["track_id"][:]]
t=(t-t.min())/1e6
selected=(labels[0]>=0)|(labels[1]>=0)
fig,axs=plt.subplots(2,2,figsize=(11,6),sharex=True,sharey="row")
for col,title in enumerate(["Existing grouping","Interferometric RANSAC"]):
    axs[0,col].set_title(title)
    for row,y in enumerate([r,v]):
        ax=axs[row,col]
        ax.scatter(t,y,s=3,c="0.75",linewidths=0)
        for label in np.unique(labels[col][labels[col]>=0]):
            use=labels[col]==label
            ax.scatter(t[use],y[use],s=12,linewidths=0,label=str(label))
        if np.any(selected):
            lo,hi=np.min(y[selected]),np.max(y[selected])
            pad=max((hi-lo)*.15,.5)
            ax.set_ylim(lo-pad,hi+pad)
            ax.set_xlim(t[selected].min()-.1,t[selected].max()+.1)
    axs[1,col].set_xlabel("Time (s)")
axs[0,0].set_ylabel("Range (km)")
axs[1,0].set_ylabel("Radial velocity (km/s)")
fig.tight_layout()
fig.savefig(a.output,dpi=160)
