import time
import scipy as sp
import numpy as np
from rtlsdr import RtlSdr
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, get_sun

sdr = RtlSdr()
sdr.sample_rate = 2.4e6
sdr.center_freq = 1420.40575e6     
sdr.freq_correction = 1  
sdr.gain = 40.2
num_points = 2048
duration = 1800
loc = EarthLocation(lon = 101.6, lat = 3.01, height = 88)
average_arr = []
time_arr = []

while True:
    t=Time(Time.now(), scale = 'utc', location = loc)
    LST = t.sidereal_time('mean')
    point_hour = np.round(LST.value,2)
    sun_hour = get_sun(t).ra.value/15
    sun_dec = get_sun(t).dec.value
    diff = np.round(abs(point_hour - sun_hour), 2)
    if diff <= 1.4:
        break

print(Time.now())
print('Start')
t1 = time.perf_counter()
while True:
    
    raw_samples = sdr.read_samples(num_points)
    power = sum((abs(raw_samples))**2) / num_points
    average_arr.append(power)
    t2 = time.perf_counter()
    diff = t2 - t1
    time_arr.append(diff)
    if diff > duration:
        break
    time.sleep(1)

print(Time.now())
print("Finished")
x = np.array(time_arr)
y = np.array(average_arr)
np.save("bwtime.npy", x)
np.save("bwpower.npy", y)

mean = sum(x * y) / sum(y)
sigma = np.sqrt(sum(y * (x-mean)**2) / sum(y))

def gauss(x, a, u, sig, offset):
    return a * np.exp(- (x - u)**2 / (2 * sig**2)) + offset

para, covariance = curve_fit(gauss, x, y, p0=[1, mean, sigma, 1])
sig = para[2]
fwhm = abs(2 * sig * np.sqrt(2 * np.log(2)))
print("HPBW: ", fwhm * (0.25 / 60) * np.cos((sun_dec * np.pi) / 180))

np.save("gausspower.npy", gauss(x, *para))

plt.plot(x, gauss(x, *para))
plt.plot(x, y)
plt.show()




        
            
    
