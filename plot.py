import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import scipy as sp
import scipy.signal

rest_freq=1420.40575e6
c=299792.458
freq=np.load("fft_freq.npy")
velocity=c*((rest_freq/freq)-1)
velocity_res = c * ((rest_freq / (rest_freq + freq[0] - freq[1])) - 1)
std_arr=[]
size=2048
bound=int(93 // velocity_res)
bound2=int(73 // velocity_res)

for i in range(36):
    index=[]
    dec = -10
    data = np.load('H1Spectra2/ONDEC' + str(dec) + '/' + 'corrected' + str(i * 10) + '.npy')*47.5
    velocity = np.load('H1Spectra2/ONDEC' + str(dec) + '/' + 'velocity' + str(i * 10) + '.npy')

    #data = sp.signal.medfilt(data,5)
    peaks, _ = find_peaks(data, prominence=2.5)
    for a in range(len(peaks)):
        if data[peaks[a]]/0.475 < 10:
            index.append(a)
    peaks = np.delete(peaks,index)
    pos_indexes = np.where(np.round(velocity) == 150)
    neg_indexes = np.where(np.round(velocity) == -150)
    plt.subplot(6,6,i+1)
    plt.title(str(i*10))
    #plt.plot(velocity, data)
    plt.plot(velocity[pos_indexes[0][0]:neg_indexes[0][-1]], data[pos_indexes[0][0]:neg_indexes[0][-1]])
    plt.plot(velocity[peaks],data[peaks],"xr")
    
    # yback=data[size-bound-bound2:]
    # yfront=data[:bound-bound2]

    # std=np.std(yfront)
    # std=np.std(yback)
    # y_comb=np.concatenate((yfront,yback))
    # std=np.std(y_comb)
    # std_arr.append(std)

#print(std_arr)
#nS=np.array(std_arr)
#print(nS*np.sqrt((2.4e6/2048)*120))
plt.show()