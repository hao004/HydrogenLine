import time
import scipy as sp
import scipy.signal
import numpy as np
from rtlsdr import RtlSdr
import matplotlib.pyplot as plt

sdr = RtlSdr()
sdr.sample_rate = 2.4e6
sdr.center_freq = 1420.40575e6     
sdr.freq_correction = 1  
sdr.gain = 40.2
num_points = 2048
duration = 600
average_arr = []
time_arr = []

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
    elif diff > 300:
        print("Hot")
    time.sleep(1)

x = np.array(time_arr)
y = np.array(average_arr)
np.save("powertime.npy", x)
np.save("power.npy", y)
