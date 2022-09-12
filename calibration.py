import os
import datetime
import numpy as np
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
velocity_res = c * ((rest_freq / (freq + freq[0] - freq[1] )) - 1)
x = np.linspace(0, fft_size, num = fft_size)
loc = EarthLocation(lon = 101.6, lat = 3.01, height = 88)
bound = int(60 / velocity_res)
bound2 = int(40 / velocity_res)

for a in range(13):
    
    for i in range(36):
        
        corrected_vel = velocity * u.km/u.s
        rv = c * ((rest_freq / rest_freq) - 1) * u.km/u.s
        file = "H1Spectra/ONDEC" + str(-60 + 10*a) + "/" + str(i * 10) + ".npy"
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
        
        raw = np.load("H1Spectra/ONDEC" + str(-60 + 10*a) + "/" + str(i * 10) + ".npy")
        cold = np.load("COLD.npy")
        model = np.poly1d(np.polyfit(x, cold, 30))
        curvefit = model(x)

        product = sum(curvefit[:bound] * raw[:bound])
        square = sum(curvefit[:bound] ** 2)
        a = product / square

        modified = model(x) * a
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

        source = source - line
        n = 10 // 2
        yback = source[fft_size - bound: fft_size - bound2]
        yfront = source[bound2:bound]
        y = np.concatenate((yback, yfront))
        std = np.std(y)
        for i in range(fft_size - bound*2):
            
            p = i + bound
            med = np.median(source[p - n: p + n + 1])
            if abs(source[p] - med) / std > 12:
                source[p] = med
                print(p)
            if p == 401:
                source[p] = med

        np.save(os.path.join('H1Spectra/CALDEC' + str(-60 + 10*a), 'velocity' + str(i * 10)), \
        corrected_vel[bound:fft_size-bound].value)
        np.save(os.path.join('H1Spectra/CALDEC' + str(-60 + 10*a), 'power' + str(i * 10)), \
        source[bound:fft_size-bound])
        

