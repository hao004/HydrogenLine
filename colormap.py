from cProfile import label
from cmath import pi
from matplotlib.colors import LogNorm
import numpy as np
from numpy import maximum, trapz
import matplotlib.pyplot as plt
from matplotlib import colors

size=2048
rest_freq=1420.40575e6
c=299792.458
freq=np.load("fft_freq.npy")
velocity_res = c * ((rest_freq / (rest_freq + freq[0] - freq[1])) - 1)
x_bound = int(73 // velocity_res)
velocity=c*((rest_freq/freq)-1)
pixels=[]
#print(velocity[x_bound:size-x_bound])
for a in range(13):
    
    source=[]
    
    for i in range(36):
        
        data = np.load('H1Spectra2/ONDEC' + str(60 - 10*a) + '/' + 'corrected' + str(i * 10) + '.npy')*40
        area = trapz(data, np.flipud(velocity[x_bound:size-x_bound]))
        #source.append(area)
        #log_density = np.log10(1.82*10**18*abs(area))
        density=abs(1.82*10**18*area)
        #source.append(log_density)
        source.append(density)

    pixels.append(source)

plt.imshow(pixels, cmap = 'jet', interpolation = 'catrom', extent = [-5,355,-65,65], aspect = 1)
plt.colorbar(orientation = "horizontal", label=( r"$N_{HI}$" " (cm\u207b\u00b2) "))
plt.xlabel('RA(degrees)')
plt.ylabel('DEC(degrees)')
plt.show()
