import numpy as np
from numpy import trapz
import matplotlib.pyplot as plt

pixels=[]

for a in range(13):
    source=[]
    for i in range(37):
        data=np.load('H1Spectra/CALDEC'+str(-60+10*a)+'/'+str(i*10))
        area=trapz(data,dx=1)
        source.append(area)
    pixels.append(source)

plt.imshow(pixels,vmin=70,vmax=300,cmap='hot')#,interpolation='gaussian')#extent=[4.33,7.97,0.5,5.5],aspect=0.0667)
plt.colorbar(orientation="horizontal")
plt.xlabel('RA(hours)')
plt.ylabel('DEC(degrees)')
plt.show()
