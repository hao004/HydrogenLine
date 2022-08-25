import os
import datetime

for i in range(13):
    path1=r'H1Spectra/ONDEC'+str(-60+(10*i))
    os.mkdir(path1)
for i in range(13):
    path2=r'H1Spectra/OFFDEC'+str(-60+(10*i))
    os.mkdir(path2)
for i in range(13):
    path3=r'H1Spectra/CALDEC'+str(-60+(10*i))
    os.mkdir(path3)




