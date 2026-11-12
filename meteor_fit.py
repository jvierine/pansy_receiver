import numpy as n
import matplotlib.pyplot as plt

def robust_polynomial(t, points, *, degree=1, valid=None, floor=.025, threshold=4.5):
    """Deterministic Theil–Sen initialization, followed by joint-echo clipping.

    A single mask rejects complete observations, not unrelated coordinate
    components. Median/MAD scales prevent a handful of large errors from
    inflating their own acceptance threshold. Time is centered/scaled before
    fitting. Observations are retained, with the rejection mask returned.
    """
    t=n.asarray(t,dtype=float);points=n.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[0]!=len(t):
        raise ValueError('Invalid trajectory dimensions')
    finite=n.isfinite(t)&n.all(n.isfinite(points),axis=1)
    if valid is not None: finite &= n.asarray(valid,dtype=bool)
    if n.count_nonzero(finite)<3: raise ValueError('Fewer than three valid echoes')
    epoch=float(n.mean(t[finite]));span=float(n.ptp(t[finite]))
    if span<=0: raise ValueError('No temporal baseline')
    x=(t-epoch)/span
    ordered=n.flatnonzero(finite)[n.argsort(t[finite])]
    selected=ordered[n.linspace(0,len(ordered)-1,min(32,len(ordered))).astype(int)]
    a,b=n.triu_indices(len(selected),1);dt=x[selected[b]]-x[selected[a]]
    useful=abs(dt)>.1
    if not useful.any(): raise ValueError('No distinct echo times')
    slopes=n.median((points[selected[b[useful]]]-points[selected[a[useful]]])/dt[useful,None],axis=0)
    intercept=n.median(points[finite]-x[finite,None]*slopes,axis=0)
    model=intercept+x[:,None]*slopes
    mask=finite.copy();design=n.vander(x,N=degree+1,increasing=True)
    minimum=max(3,degree+2)
    for _ in range(10):
        residual=points-model
        center=n.median(residual[mask],axis=0)
        scale=n.maximum(1.4826*n.median(abs(residual[mask]-center),axis=0),n.asarray(floor))
        score=n.sqrt(n.sum(((residual-center)/scale)**2,axis=1))
        new=finite&(score<=threshold)
        if n.count_nonzero(new)<minimum:
            raise ValueError('Fewer than required consistent echoes')
        coefficient,_,rank,_=n.linalg.lstsq(design[new],points[new],rcond=None)
        if rank<degree+1: raise ValueError('Rank-deficient trajectory')
        model=design@coefficient
        stable=n.array_equal(mask,new);mask=new
        if stable: break
    residual=points-model
    std=n.sqrt(n.mean(residual[mask]**2,axis=0))
    return dict(epoch=epoch,span=span,coefficient=coefficient,model=model,inlier=mask,
                residual=residual,std=std,scale=scale)


def robust_fit_meteor(t,ea,no,up):
    result=robust_polynomial(t,n.column_stack((ea,no,up)))
    result['r0']=result['coefficient'][0]
    result['v0']=result['coefficient'][1]/result['span']
    return result

def simple_fit_meteor(t,ea,no,up):
    """
    all in units of km
    """
    tmean=n.mean(t)
    nm=len(t)
    tn=t-tmean
    # r0x, r0y, r0z, v0x, v0y, v0z
    A=n.zeros([3*nm,6])
    m=n.zeros(3*nm)
    # position measurement
    
    # east pos,  r0x, 
    A[0:nm,0]=1
    # west pos r0y
    A[(nm):(2*nm),1]=1
    # up pos r0z
    A[(2*nm):(3*nm),2]=1
      
    # east pos v0x, 
    A[0:nm,3]=tn
    # west pos v0y
    A[(nm):(2*nm),4]=tn
    # up pos v0z
    A[(2*nm):(3*nm),5]=tn

    m[0:nm]=ea
    m[nm:(2*nm)]=no
    m[(2*nm):(3*nm)]=up

    xhat=n.linalg.lstsq(A,m)[0]
    r0=xhat[0:3]
    v0=xhat[3:6]
    model=n.dot(A,xhat)
    r2=n.abs(m-model)**2
    ea_res=n.sqrt(n.mean(r2[0:nm]))
    no_res=n.sqrt(n.mean(r2[nm:(2*nm)]))
    up_res=n.sqrt(n.mean(r2[(2*nm):(3*nm)]))
    
    return(r0,v0,model,ea_res,no_res,up_res)

      

    
    
