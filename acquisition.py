import os
import time
import scipy as sp
import scipy.signal
import numpy as np
from rtlsdr import RtlSdr
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

sdr=RtlSdr()
sdr.sample_rate=2.4e6
sdr.freq_correction=1
sdr.gain=40.2
num_points=2048
integration_time=120
block_duration=num_points/sdr.sample_rate
num_spectra=int(integration_time/block_duration)
rest_frequency=1420.40575e6
c=299792.458
freq=np.load("fft_freq.npy")
velocity=c*(1-(freq/rest_frequency))
loc=EarthLocation(lon=101.6,lat=3.01,height=88)
raw_arr=[]
fullraw_arr=[]

dec=input("Enter dec:")
path="H1Spectra/ONDEC"+dec
rawfiles=os.listdir(path)

for i in range (len(rawfiles)):
    filename=os.path.splitext(rawfiles[i])[0]
    raw_arr.append(filename)
raw_arr=list(map(int,raw_arr))

for i in range(36):
    ra=10*i
    fullraw_arr.append(ra)
    
for ra in raw_arr:
    if ra in fullraw_arr:
        fullraw_arr.remove(ra)

while True:
    fulltime=Time.now()
    t=Time(fulltime,scale='utc',location=loc)
    LST=t.sidereal_time('mean')
    hour=LST.value
    degree=np.round(hour*15,2)
    
    if len(os.listdir(path))==36:
        print("Observations at "+str(dec)+" completed.")
        break
    elif degree in fullraw_arr:
        fullraw_arr.remove(degree)
        print("Observing at\nRA(h,m):"+str(LST)+" DEC(degree):"+str(dec))
        sdr.center_freq=1420.40575e6
        zero_arr=np.load("zero_arr.npy")
        t1=time.perf_counter()
        for _ in range(num_spectra):
            raw_samples=sdr.read_samples(num_points)*np.hamming(num_points)
            spectra=(abs(np.fft.fft(raw_samples)))**2
            zero_arr=zero_arr+spectra
        avgsource=np.fft.fftshift((zero_arr/num_spectra))
        medsource=(sp.signal.medfilt(avgsource,15))
        np.save(os.path.join('H1Spectra/ONDEC'+dec,str(int(degree))),medsource)

        sdr.center_freq=1423e6
        zero_arr=np.load("zero_arr.npy")
        for _ in range(num_spectra):
            raw_samples=sdr.read_samples(num_points)*np.hamming(num_points)
            spectra=(abs(np.fft.fft(raw_samples)))**2
            zero_arr=zero_arr+spectra
        t2=time.perf_counter()
        avgcold=np.fft.fftshift((zero_arr/num_spectra))
        medcold=(sp.signal.medfilt(avgcold,15))
        np.save(os.path.join('H1Spectra/OFFDEC'+dec,str(int(degree))),medcold)
        print("Observation completed. Time taken in seconds: "+str(t2-t1))
        source=medsource-medcold
        np.flip(source)
        plt.clf()
        plt.plot(velocity,source)
        plt.savefig("H1Pic/RA"+str(degree)+"DEC"+str(dec))
    else:
        continue

