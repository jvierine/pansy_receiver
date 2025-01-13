import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import scipy.constants as sci
import pandas

mpl.rcParams.update({'font.size': 22})
mpl.rcParams['figure.figsize'] = 20, 10

radar_frequency=47E6
wavelength = sci.c / radar_frequency

def uv_space(ant_posx, ant_posy, wavelength):
    """Generate the U,V space from X,Y,Z"""
    n = len(ant_posx)
    u = np.array([])
    v = np.array([])
    for i in range(n):
        for j in range(i + 1, n):
            u = np.append(u, (ant_posx[i] - ant_posx[j]) / wavelength)
            v = np.append(v, (ant_posy[i] - ant_posy[j]) / wavelength)
            u = np.append(u, (ant_posx[j] - ant_posx[i]) / wavelength)
            v = np.append(v, (ant_posy[j] - ant_posy[i]) / wavelength)
    return u, v

#read in the antenna positions here

antenna_positions = pandas.read_csv('antpos.csv')
grp_lookup = pandas.read_csv('ant_grp_assign2.csv')

#SerialNo. GrpName ModuleID  Ready       X(m)       Y(m)  Z(m)

print(antenna_positions)

print(grp_lookup)

grp_name = np.asarray(antenna_positions.get('GrpName'))
x_pos_original = np.asarray(antenna_positions.get('X(m)'),dtype=np.float32)
y_pos_original = np.asarray(antenna_positions.get('Y(m)'),dtype=np.float32)
z_pos_original = np.asarray(antenna_positions.get('Z(m) '),dtype=np.float32)
data_good = np.asarray(antenna_positions.get('Ready'),dtype=np.float32)
cluster_group = np.asarray(grp_lookup.get('ClusterGrp'),dtype=np.int32)

#grp_name = grp_name[np.where(data_good>0)]
#x_pos_original = x_pos_original[np.where(data_good>0)]
#y_pos_original = y_pos_original[np.where(data_good>0)]
#z_pos_original = z_pos_original[np.where(data_good>0)]

print(cluster_group)
cluster_group_nums = np.unique(cluster_group)
cluster_group_nums.sort()

print(cluster_group_nums)

unique_groups = np.unique(grp_name)
unique_groups.sort()

x_avg_hex = []
y_avg_hex = []
z_avg_hex = []

for x in range(len(unique_groups)):
    x_avg_hex.append(np.mean(x_pos_original[np.where(grp_name==unique_groups[x])]))
    y_avg_hex.append(np.mean(y_pos_original[np.where(grp_name==unique_groups[x])]))
    z_avg_hex.append(np.mean(z_pos_original[np.where(grp_name==unique_groups[x])]))

ant_posx = np.asarray(x_avg_hex)
ant_posy = np.asarray(y_avg_hex)

x_avg_hex = np.asarray(x_avg_hex)
y_avg_hex = np.asarray(y_avg_hex)
z_avg_hex = np.asarray(z_avg_hex)

num_rx_antennas = len(ant_posx)

#end reading in antenna positions

#determine grouping of antennas here

#end of determining antenna groupings

u, v = uv_space(ant_posx, ant_posy, wavelength)

# Pretty plot configuration.
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 22
plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labelsa
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


lbls = [0, 1]
fig, axs = plt.subplots(2, 3)#, gridspec_kw={'width_ratios': [1, 1, 1]})
axs[0,0].set_title("Antenna Positions")
axs[0,0].set_xlabel("East-West[m]")
axs[0,0].set_ylabel("North-South [m]")
axs[0,0].grid(which='both')
axs[0,0].axis('equal')
axs[0,0].scatter(ant_posx[:], ant_posy[:], marker='v', color='k')
#axs[0,0].scatter(ant_posx[0], ant_posy[0], marker='o', color='red',s=100)
#axs[0].set_ylim(-225, 10)
#for i , l in enumerate(lbls):
#axs[0,0].annotate('Existing RX', (-23.0, 3.0),color='red')

axs[0,1].set_title("Sampling Space")
axs[0,1].set_xlabel("u")
axs[0,1].set_ylabel("v")
axs[0,1].grid(which='both')
axs[0,1].axis('equal')
axs[0,1].scatter(u, v, color='k',s=1)
axs[0,1].scatter(0, 0, color='k',s=1)

# axs[2].set_title("Dirty Beam")
# axs[2].set_xlabel("Az [deg]")
# axs[2].set_ylabel("El [deg]")
# axs[2].grid(which='both')
# axs[2].axis('equal')
# axs[2].pcolormesh(AZ, EL, C, cmap='inferno')
# axs[2].colorbar()

#plt.savefig('array_map.pdf')

rad_freq = radar_frequency

rad_lambda = 2.997E8/rad_freq

sphere_radius = 1E5 # in meters

elevation_values = np.radians(90-((np.arange(160))*0.25))
azimuth_values = np.radians(np.arange(721)*0.5)

xx,yy = np.meshgrid(azimuth_values,elevation_values,sparse=True)

l_values = np.cos(xx)*np.cos(yy)
m_values = np.sin(xx)*np.cos(yy)
q_values = np.sin(yy)

u_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
v_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
w_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
temp_ind=0
for first_antenna in range(num_rx_antennas):
    for second_antenna in range(first_antenna+1,num_rx_antennas):
        u_values[temp_ind] = (ant_posx[first_antenna]-ant_posx[second_antenna])
        v_values[temp_ind] = (ant_posy[first_antenna]-ant_posy[second_antenna])
        w_values[temp_ind] = 0
        temp_ind+=1

r_values = np.zeros((len(elevation_values),len(azimuth_values),int(num_rx_antennas*(num_rx_antennas-1)/2)),dtype=np.float32)
r_values = np.einsum('b,lm->lmb',u_values,l_values)+np.einsum('b,lm->lmb',v_values,m_values)

antenna_coherence = np.ones(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.complex64)

A_matrix = np.exp(-1.0j*2.0*np.pi*r_values/rad_lambda)

rx_classic_brightness = np.zeros((len(elevation_values),len(azimuth_values)),dtype=np.complex64)

rx_classic_brightness = np.einsum('b,lmb->lm',antenna_coherence,A_matrix)

axs[0,2].contourf(np.arcsin(l_values)*180/np.pi,np.arcsin(m_values)*180/np.pi,10*np.log10(np.abs(rx_classic_brightness)/(int(num_rx_antennas*(num_rx_antennas-1)/2))),levels=np.arange(-10,0,1)+1.0,extend='both')
CS = axs[0,2].contour(np.arcsin(l_values)*180/np.pi,np.arcsin(m_values)*180/np.pi,10*np.log10(np.abs(rx_classic_brightness)/(int(num_rx_antennas*(num_rx_antennas-1)/2))),levels=[-3],colors='w',zorder=3,linestyles='solid')
axs[0,2].set_xlabel('l direction-cosine (degrees)')
axs[0,2].set_ylabel('m direction-cosine (degrees)')
#axs[0,2].set_ylim(-20,20)
#axs[0,2].set_xlim(-45,45)
axs[0,2].clabel(CS, inline=1, fontsize=12,levels=[-3],colors='w',zorder=3,fmt='%d',inline_spacing=0)

#Now implement the analysis with the modified antenna locations

x_avg_grp = []
y_avg_grp = []
z_avg_grp = []

print(np.where(cluster_group==cluster_group_nums[1])[0])

for x in range(1,len(cluster_group_nums)):
    if len(np.where(cluster_group==cluster_group_nums[x])[0])>1:
        x_avg_grp.append(np.mean(x_avg_hex[np.where(cluster_group==cluster_group_nums[x])[0]]))
        y_avg_grp.append(np.mean(y_avg_hex[np.where(cluster_group==cluster_group_nums[x])[0]]))
        z_avg_grp.append(np.mean(z_avg_hex[np.where(cluster_group==cluster_group_nums[x])[0]]))
    else:
        x_avg_grp.append(np.mean(x_avg_hex[np.where(cluster_group==cluster_group_nums[x])]))
        y_avg_grp.append(np.mean(y_avg_hex[np.where(cluster_group==cluster_group_nums[x])]))
        z_avg_grp.append(np.mean(z_avg_hex[np.where(cluster_group==cluster_group_nums[x])]))    

ant_posx = np.asarray(x_avg_grp)
ant_posy = np.asarray(y_avg_grp)

num_rx_antennas = len(ant_posx)

print(ant_posx)

#end reading in antenna positions

#determine grouping of antennas here

#end of determining antenna groupings

u, v = uv_space(ant_posx, ant_posy, wavelength)

axs[1,0].set_title("Antenna Positions")
axs[1,0].set_xlabel("East-West[m]")
axs[1,0].set_ylabel("North-South [m]")
axs[1,0].grid(which='both')
axs[1,0].axis('equal')
axs[1,0].scatter(ant_posx[:], ant_posy[:], marker='v', color='k')
#axs[0,0].scatter(ant_posx[0], ant_posy[0], marker='o', color='red',s=100)
#axs[0].set_ylim(-225, 10)
#for i , l in enumerate(lbls):
#axs[0,0].annotate('Existing RX', (-23.0, 3.0),color='red')

axs[1,1].set_title("Sampling Space")
axs[1,1].set_xlabel("u")
axs[1,1].set_ylabel("v")
axs[1,1].grid(which='both')
axs[1,1].axis('equal')
axs[1,1].scatter(u, v, color='k',s=1)
axs[1,1].scatter(0, 0, color='k',s=1)

# axs[2].set_title("Dirty Beam")
# axs[2].set_xlabel("Az [deg]")
# axs[2].set_ylabel("El [deg]")
# axs[2].grid(which='both')
# axs[2].axis('equal')
# axs[2].pcolormesh(AZ, EL, C, cmap='inferno')
# axs[2].colorbar()

#plt.savefig('array_map.pdf')

rad_freq = radar_frequency

rad_lambda = 2.997E8/rad_freq

sphere_radius = 1E5 # in meters

elevation_values = np.radians(90-((np.arange(160))*0.25))
azimuth_values = np.radians(np.arange(721)*0.5)

xx,yy = np.meshgrid(azimuth_values,elevation_values,sparse=True)

l_values = np.cos(xx)*np.cos(yy)
m_values = np.sin(xx)*np.cos(yy)
q_values = np.sin(yy)

u_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
v_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
w_values = np.zeros(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.float32)
temp_ind=0
for first_antenna in range(num_rx_antennas):
    for second_antenna in range(first_antenna+1,num_rx_antennas):
        u_values[temp_ind] = (ant_posx[first_antenna]-ant_posx[second_antenna])
        v_values[temp_ind] = (ant_posy[first_antenna]-ant_posy[second_antenna])
        w_values[temp_ind] = 0
        temp_ind+=1

r_values = np.zeros((len(elevation_values),len(azimuth_values),int(num_rx_antennas*(num_rx_antennas-1)/2)),dtype=np.float32)
r_values = np.einsum('b,lm->lmb',u_values,l_values)+np.einsum('b,lm->lmb',v_values,m_values)

antenna_coherence = np.ones(int(num_rx_antennas*(num_rx_antennas-1)/2),dtype=np.complex64)

A_matrix = np.exp(-1.0j*2.0*np.pi*r_values/rad_lambda)

rx_classic_brightness = np.zeros((len(elevation_values),len(azimuth_values)),dtype=np.complex64)

rx_classic_brightness = np.einsum('b,lmb->lm',antenna_coherence,A_matrix)

axs[1,2].contourf(np.arcsin(l_values)*180/np.pi,np.arcsin(m_values)*180/np.pi,10*np.log10(np.abs(rx_classic_brightness)/(int(num_rx_antennas*(num_rx_antennas-1)/2))),levels=np.arange(-10,0,1)+1.0,extend='both')
CS = axs[1,2].contour(np.arcsin(l_values)*180/np.pi,np.arcsin(m_values)*180/np.pi,10*np.log10(np.abs(rx_classic_brightness)/(int(num_rx_antennas*(num_rx_antennas-1)/2))),levels=[-3],colors='w',zorder=3,linestyles='solid')
axs[1,2].set_xlabel('l direction-cosine (degrees)')
axs[1,2].set_ylabel('m direction-cosine (degrees)')
#axs[0,2].set_ylim(-20,20)
#axs[0,2].set_xlim(-45,45)
axs[1,2].clabel(CS, inline=1, fontsize=12,levels=[-3],colors='w',zorder=3,fmt='%d',inline_spacing=0)

plt.tight_layout()
plt.savefig('array_map_nov13_2024.png')
plt.close()

