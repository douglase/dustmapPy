import dustmap
import numpy as np
import astropy.units as u
import astropy.io.fits as fits
import matplotlib.pyplot as plt

inputfile="dustmap/sample_inputfile.dat"
distance = 10.                #pc
fov = [1.*u.arcsecond]                 #mas
pixelsize=0.005*u.arcsecond
numpixels = 200
datatype = 1
composition = 'astrosil'
inclination = 30. #degrees
longitude = 90. #degrees
pa = 0. #;degrees
opticaldepth = 0
rdust = 1.
thermal = 1
Lstar = 1.
Tstar = 5778.
kurucz = 1
logg = 4.5
scattered = 1
thermal = 1
hg = 1
g = 0.1
#scaling= np.array([0.],dtype=np.double)
wavel = np.array([.5],dtype=np.float32)
Tstar = 5778.

hist=np.zeros([numpixels,numpixels],dtype=np.double)
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
           np.array([0.],dtype=np.float32),
          composition,
          scatteredlight_flag=scattered,
          thermalemission_flag=thermal,
          verbose_flag=2,
          inclination = 40.0, #degrees
          longitude = longitude, #degrees
         # pa = pa ,#degrees
          #opticaldepth_flag = 0,
          #densityhisto_flag=1,
          rdust =rdust,
          Lstar = Lstar,
          kurucz = kurucz,
          effrdust=1,
          logg=4.5,
          HG_flag = hg,
          HG_g = g,
          iwa=0.00*u.arcsecond,
          Tsublimate=1e10,
          ncostheta=500,
          scaling=np.array([1.0],dtype=np.double),
          )
print(hist)
hist[0,:]=0 #dustmap is setting first row to 1 ?
fits.writeto("test.fits",hist,overwrite=True)

plt.imshow(np.array(hist),origin="upper left")

plt.show()
