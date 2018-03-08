import dustmap
import numpy as np
import astropy.units as u
inputfile="dustmap/sample_inputfile.dat"
distance = 10.                #pc
fov = [1.*u.arcsecond]                 #mas
pixelsize=1.0*u.arcsecond
numpixels = 200
datatype = 1
inclination = 60. #degrees
longitude = 90. #degrees
pa = 45. #degrees
opticaldepth = 1
rdust = 1.
thermal = 1
wavel = np.array([1.],dtype=np.double)
Lstar = 1.
Tstar = 5778.
kurucz = 1
logg=4.5
composition = 'astrosil'
hg = 1
g = 0.1
iwa=0.05*u.arcsecond
tsublimate=1e10
fstar= np.array([0.],dtype=np.double)
scaling= np.empty(1,dtype=np.int)

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
          Lstar=Lstar,
          scatteredlight_flag=1,
          verbose_flag=5,
          )
print(hist)
          #pfunc,
