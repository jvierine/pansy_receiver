import numpy as n
import scipy.constants as sc

# PANSY
G=10**3.7
Grx=(10**3.7/54)
Pt=500e3
f=47e6
lam=sc.c/f
A=Grx*lam**2/(4*n.pi)
B=1/128e-6
tsys=5000
#SNR=Pt*G*A*sigma/((4*n.pi)**2*(100e3**4)*sc.k*4000*B)
SNR=4
sigma=SNR*(((4*n.pi)**2)*(100e3**4)*sc.k*tsys*B)/(Pt*G*A) #=*sigma/
#print(sigma)
print(10*n.log10(sigma))


# MAARSY
G=10**3.35
Grx=10**2.5
Pt=800e3
f=53.5e6
lam=sc.c/f
A=Grx*(lam**2)/(4*n.pi)
B=1/32e-6

#SNR=Pt*G*A*sigma/((4*n.pi)**2*(100e3**4)*sc.k*4000*B)
SNR=4
sigma=SNR*(((4*n.pi)**2)*(100e3**4)*sc.k*tsys*B)/(Pt*G*A) #=*sigma/
#print(sigma)
print(10*n.log10(sigma))
