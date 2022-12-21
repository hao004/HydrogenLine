from cProfile import label
from cmath import pi
from matplotlib.colors import LogNorm
import numpy as np
from numpy import maximum, trapz
import matplotlib.pyplot as plt
from matplotlib import colors

#size=2048
#rest_freq=1420.40575e6
#c=299792.458
#freq=np.load("fft_freq.npy")
#velocity_res = c * ((rest_freq / (rest_freq + freq[0] - freq[1])) - 1)
#x_bound = int(73 // velocity_res)
#velocity=c*((rest_freq/freq)-1)

pixels=[]
# min_velocity=-999
# max_velocity=999
# min_arr=[]
# max_arr=[]
# indexes=[]

# for a in range(13):
#     for i in range(36):
#         velocity = np.load('H1Spectra2/ONDEC' + str(60 - 10*a) + '/' + 'velocity' + str(i * 10) + '.npy')
#         min_index = np.where(velocity==[min(velocity[velocity<0])])
#         max_index = np.where(velocity==[max(velocity[velocity>0])])
#         if velocity[min_index[0][0]]>min_velocity:
#             min_velocity = velocity[min_index[0][0]]
#         if velocity[max_index[0][0]]<max_velocity:
#             max_velocity = velocity[max_index[0][0]]
#     min_arr.append(min_velocity)
#     max_arr.append(max_velocity)
# print(np.round(min_arr),np.round(max_arr))

for a in range(13):
    
    source=[]
    
    for i in range(36):
        
        data = np.load('H1Spectra2/ONDEC' + str(60 - 10*a) + '/' + 'corrected' + str(i * 10) + '.npy')*47.5
        velocity = np.load('H1Spectra2/ONDEC' + str(60 - 10*a) + '/' + 'velocity' + str(i * 10) + '.npy')
        pos_indexes = np.where(np.round(velocity) == 150)
        neg_indexes = np.where(np.round(velocity) == -150)
        area = trapz(data[pos_indexes[0][0]:neg_indexes[0][-1]], velocity[pos_indexes[0][0]:neg_indexes[0][-1]])
        #source.append(area)
        #log_density = np.log10(1.82*10**18*abs(area))
        density=abs(1.82*10**18*area)
        #source.append(log_density)
        source.append(density)

    pixels.append(source)

plt.imshow(pixels, cmap = 'jet', interpolation = 'catrom', extent = [-1/3, 24 + 1/3, -65, 65], aspect = 1/15)
plt.colorbar(orientation = "horizontal", label=( r"Hydrogen Column Density $N_{HI}$" " (cm\u207b\u00b2) "))
plt.title('Hydrogen Line Map')
plt.xlabel('RA(hours)')
plt.xlim([0,24])
plt.xticks([0,6,12,18,24])
plt.ylabel('DEC(degrees)')

plt.show()
