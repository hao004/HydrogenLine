import os
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
regions = ['S6',"S7","S8","S9"]
degrees = [233,32,87,268]
full = np.array(degrees)

for _ in range(15): 
    zero_arr = np.load("zero_arr.npy")
    for _ in range(num_spectra):
        raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
        spectra = (abs(np.fft.fft(raw_samples)))**2
        zero_arr = zero_arr + spectra

while True:
    fulltime = Time.now()
    t = Time(fulltime, scale = 'utc', location = loc)
    LST = t.sidereal_time('mean')
    hour = LST.value
    degree = np.round(hour * 15,3)

    for i in range(len(full)):
        
        if abs(degree - full[i]) <= 0.375:
            full = np.delete(full,[i])
            print("Observing at "+regions[i])
            zero_arr = np.load("zero_arr.npy")
            t1 = time.perf_counter()
            
            for _ in range(num_spectra):
                raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
                spectra = (abs(np.fft.fft(raw_samples)))**2
                zero_arr = zero_arr + spectra
            
            avgsource = np.fft.fftshift((zero_arr / num_spectra))
            t2 = time.perf_counter()
            print("Observation completed. Time taken in seconds: " + str(t2 - t1))
            np.save(regions[i] + ".npy", avgsource)
            
    
    zero_arr = np.load("zero_arr.npy")
    for _ in range(num_spectra):
        raw_samples = sdr.read_samples(num_points) * np.hamming(num_points)
        spectra = (abs(np.fft.fft(raw_samples)))**2
        zero_arr = zero_arr + spectra

