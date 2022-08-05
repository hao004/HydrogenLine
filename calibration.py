import os
import datetime
import scipy as sp
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import ICRS, LSR
from astropy.coordinates import SkyCoord, EarthLocation

c=299792.458
rest_freq=1420.40575e6
freq=np.load("fft_freq.npy")
velocity=c*((rest_freq/freq)-1)

x=np.linspace(0,2048,num=2048)
loc=EarthLocation(lon=101.6,lat=3.01,height=88)

for a in range(27):
    for i in range(71):
        corrected_vel=velocity*u.km/u.s
        rv=c*((rest_freq/rest_freq)-1)*u.km/u.s
        file="H1Spectra/ONDEC"+str(-65+5*a)+"/"+str(i)+".npy"
        epoch=os.path.getmtime(file)
        fulltime=datetime.datetime.utcfromtimestamp(epoch)
        t=Time(fulltime,scale='utc',location=loc)
        LST=t.sidereal_time('mean')
        t.format='mjd'
        LSThms=LST

        ra=(LST.deg)*u.deg
        dec=(3)*u.deg
        sc=SkyCoord(ra,dec)
        vcorr=sc.radial_velocity_correction(kind='barycentric', obstime=t)  
        rv=rv+vcorr

        my_observation = ICRS(ra,dec, \
                pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr, \
                radial_velocity=rv, distance = 1*u.pc)

        new_rv = my_observation.transform_to(LSR()).radial_velocity
        corrected_vel=corrected_vel+new_rv
        
        raw=np.load("H1Spectra/ONDEC"+str(-65+5*a)+"/"+str(i)+".npy")
        cold=np.load("H1Spectra/OFFDEC"+str(-65+5*a)+"/"+str(i)+".npy")
        model=np.poly1d(np.polyfit(x,cold,30))
        curvefit=model(x)

        product=sum(curvefit[1700:]*raw[1700:])
        square=sum(curvefit[1700:]**2)
        z=product/square
        print(z)

        modified=model(x)*z
        source=raw-modified

        plt.clf()
        plt.xlim(-250,200)
        plt.ylim(-0.03,1.1)
        plt.plot(corrected_vel,source)
        plt.savefig("H1Spectra/DEC0/"+str(LSThms)+".png")
        np.save(os.path.join('H1Spectra/CALDEC'+str(-65+5*a),'velocity'+str(i)),corrected_vel.value)
        np.save(os.path.join('H1Spectra/CALDEC'+str(-65+5*a),'temp'+str(i)),source)
        




        
