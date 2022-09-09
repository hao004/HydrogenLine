import os
from socket import IPV6_RTHDR_TYPE_0
import time
import numpy as np
from rtlsdr import RtlSdr
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

sdr = RtlSdr()
sdr.sample_rate = 2.4e6
sdr.center_freq = 1420.40575e6
sdr.freq_correction = 1
sdr.gain = 40.2
num_points = 2048
integration_time = 120
block_duration = num_points / sdr.sample_rate
num_spectra = int(integration_time / block_duration)
rest_frequency = 1420.40575e6
c = 299792.458
freq = np.load("fft_freq.npy")
velocity = c * (1 - (freq / rest_frequency))
loc = EarthLocation(lon = 101.6, lat = 3.01, height = 88)
raw_arr = []
fullraw_arr = []

dec = input("Enter dec:")
path = "H1Spectra/ONDEC"+dec
rawfiles = os.listdir(path)

#Create array of all RA
for i in range(37):
    ra = 10 * i
    fullraw_arr.append(ra)

#Create array of completed RA
for i in range (len(rawfiles)):
    filename = os.path.splitext(rawfiles[i])[0]
    raw_arr.append(filename)
raw_arr = list(map(int, raw_arr))

#Create array of incompleted RA
for ra in raw_arr:
    if ra in fullraw_arr:
        fullraw_arr.remove(ra)

full = np.array(fullraw_arr)

#Rtlsdr warm up for 30 minutes
for _ in range(15): 
    zero_arr = np.load("zero_arr.npy")
    for _ in range(num_spectra):
        raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
        spectra = (abs(np.fft.fft(raw_samples)))**2
        zero_arr = zero_arr + spectra

#Observation begin
while True:
    fulltime = Time.now()
    t = Time(fulltime, scale = 'utc', location = loc)
    LST = t.sidereal_time('mean')
    hour = LST.value
    degree = np.round(hour * 15,3)

    if len(os.listdir(path)) == 36:
        print("Observations at " + str(dec) + " completed.")
        break

    for i in range(len(full)):
        if abs(degree - full[i]) <= 0.375:
            if (full[i] == 0) or (full[i] == 360):
                index = np.where((full == 0) | (full == 360))
                full = np.delete(full, index)
            else:
                full = np.delete(full,[i])
            
            print("Observing at\nRA(h,m):" + str(LST) + " DEC(degree):" + str(dec))
            zero_arr = np.load("zero_arr.npy")
            t1 = time.perf_counter()
            
            for _ in range(num_spectra):
                raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
                spectra = (abs(np.fft.fft(raw_samples)))**2
                zero_arr = zero_arr + spectra
            
            avgsource = np.fft.fftshift((zero_arr / num_spectra))
            t2 = time.perf_counter()
            print("Observation completed. Time taken in seconds: " + str(t2 - t1))
            
            if (full[i] == 0) or (full[i] == 360):
                np.save(os.path.join('H1Spectra/ONDEC' + dec + str(0), avgsource))
            else:
                np.save(os.path.join('H1Spectra/ONDEC' + dec + str(full[i])), avgsource)
    
    zero_arr = np.load("zero_arr.npy")
    for _ in range(num_spectra):
        raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
        spectra = (abs(np.fft.fft(raw_samples)))**2
        zero_arr = zero_arr + spectra
