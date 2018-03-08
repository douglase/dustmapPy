import dustmap
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

inputfile="dustmap/sample_inputfile.dat"
distance = 10.                #pc
fov = [1.*u.arcsecond]                 #mas
pixelsize=0.05*u.arcsecond
numpixels = 200
datatype = 1
composition = 'astrosil'



fstar= np.array([0.],dtype=np.float32)
#scaling= np.array([0.],dtype=np.double)
wavel = np.array([0.55,5],dtype=np.float32)
Tstar = 5778.

hist=np.empty([numpixels,numpixels],dtype=np.double)
print(hist)
dustmap.dustmap_func(inputfile,
          hist,
          distance,
          fov,
          pixelsize,
          numpixels,
          #au,
          #degrees,
          wavel, #was lambda
          Tstar,
          fstar,
          composition,
          scatteredlight_flag=1,
          thermalemission_flag=1,
          verbose_flag=2,
          inclination = 60., #degrees
          longitude = 90., #degrees
          pa = 45. ,#degrees
          opticaldepth_flag = 1,
          rdust = 1.,
          Lstar = 1.,
          kurucz = 1,
          logg=4.5,
          HG_flag = 1,
          HG_g = 0.1,
          iwa=0.05*u.arcsecond,
          Tsublimate=1e10,
          )
plt.imshow(hist)
plt.show()
          #pfunc,
