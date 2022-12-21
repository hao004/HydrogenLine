import os
import datetime
import numpy as np
import scipy as sp
import scipy.signal
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import ICRS, LSR
from astropy.coordinates import SkyCoord, EarthLocation

c = 299792.458
rest_freq = 1420.40575e6
fft_size = 2048
freq = np.load("fft_freq.npy")
velocity = c * ((rest_freq / freq) - 1)
velocity_res = c * ((rest_freq / (rest_freq + freq[0] - freq[1])) - 1)
x = np.linspace(0, fft_size, num = fft_size)
loc = EarthLocation(lon = 101.6, lat = 3.01, height = 88)
bound = int(93 // velocity_res)
bound2 = int(73 // velocity_res)

for a in range(13):
    count_a = a
    for i in range(36):
        count_i =i
        corrected_vel = velocity * u.km/u.s
        rv = c * ((rest_freq / rest_freq) - 1) * u.km/u.s
        file = "H1Spectra2/ONDEC" + str(60 - 10*count_a) + "/" + str(count_i * 10) + ".npy"
        epoch = os.path.getmtime(file)
        fulltime = datetime.datetime.utcfromtimestamp(epoch)
        t = Time(fulltime, scale='utc', location=loc)
        LST = t.sidereal_time('mean')
        t.format = 'mjd'
        LSThms = LST

        ra = (LST.deg) * u.deg
        dec = (3) * u.deg
        sc = SkyCoord(ra, dec)
        vcorr = sc.radial_velocity_correction(kind='barycentric', obstime=t)  
        rv = rv + vcorr

        my_observation = ICRS(ra, dec, \
                pm_ra_cosdec = 0*u.mas/u.yr, pm_dec = 0*u.mas/u.yr, \
                radial_velocity = rv, distance = 1*u.pc)

        new_rv = my_observation.transform_to(LSR()).radial_velocity
        corrected_vel = corrected_vel + new_rv
        
        raw = np.load("H1Spectra2/ONDEC" + str(60 - 10*count_a) + "/" + str(count_i * 10) + ".npy")
        med_raw = sp.signal.medfilt(raw,5)
        cold = np.load("COLD.npy")
        model = np.poly1d(np.polyfit(x, cold, 40))
        curvefit = model(x)

        index = np.where(raw == np.amax(raw))
        if index[0][0]<295:
            product=sum(curvefit[fft_size-bound2:]*med_raw[fft_size-bound2:])
            square=sum(curvefit[fft_size-bound2:]**2)
        else:
            product=sum(curvefit[:bound2]*med_raw[:bound2])
            square=sum(curvefit[:bound2]**2)
        
        z = product / square
        modified = model(x) * z
        source = raw - modified

        yback = source[fft_size - bound: fft_size - bound2]
        yfront = source[bound2:bound]
        y1 = np.concatenate((yback, yfront))
        xback = velocity[fft_size - bound: fft_size - bound2]
        xfront = velocity[bound2:bound]
        x1 = np.concatenate((xback, xfront))

        mean_y = sum(y1) / len(y1)
        mean_x = sum(x1) / len(x1)
        product = (x1 - mean_x) * (y1 - mean_y)
        square = (x1 - mean_x)**2
        b1 = sum(product) / sum(square)
        b0 = mean_y - b1 * mean_x
        line = b1 * (velocity) + b0

        corrected = source - line
        n = 10//2
    
        yback=corrected[fft_size - bound:fft_size - bound2]
        yfront=corrected[bound2:bound]
    
        std = 0.01
        for c in range(fft_size - bound2*2):
            p = c + bound2
            med = np.median(corrected[p - n:p + n + 1])
            if abs(corrected[p] - med) / std > 3:
                corrected[p] = med

        np.save(os.path.join('H1Spectra2/ONDEC' + str(60 - 10*count_a), 'velocity' + str(count_i * 10)), 
        corrected_vel.value)
        np.save(os.path.join('H1Spectra2/ONDEC' + str(60 - 10*count_a), 'corrected' + str(count_i * 10)),
        corrected)
        

