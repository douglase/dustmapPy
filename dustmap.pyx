#using Cython:
""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# import both numpy and the Cython declarations for numpy
#import numpy as np
cimport numpy as np
import cython
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example)
np.import_array()

#@cython.boundscheck(False)
#@cython.wraparound(False)

# cdefine the signature of our c function
cdef extern from "dustmap_c.h":
# void cos_doubles (double * in_array, double * out_array, int size)
    void runRT(double * in_array,double * out_array, int in_col, int in_row)

"""
idl code:

;Now call dustmap_c.c and have it do the dirty work      
x = CALL_EXTERNAL(dustmappath+'dustmap_c.so', 'dustmap_c', $
                  histo, histoxsize, histoysize, inputfile4c, numfiles, datatype, $
                  userdistance, fovx, fovy, userpixelsize, userinclination, userlongitude, userpa, $
                  lambda, nlambda, Lstar, Tstar, logg, kurucz, kuruczfile4c, fstar, $
                  rdust, Tsublimate, lnkfile4c, userscaling, $
                  iwa, aitoff_flag, densityhisto_flag, opticaldepth_flag, $
                  scatteredlight_flag, thermalemission_flag, verbose_flag, $
                  oc_generic_flag, generic_oc_exp, generic_oc_albedo, HG_flag, HG_g, $
                  extinction_flag, xshift, effrdust, $
                  markx, marky, markz, markweight, nmarks, markabs_flag, azavg, distmask, $
                  ncostheta, internal_pfunc, costheta, pfunc, qabs, qsca, nrdust)

"""

# create the wrapper code, with numpy type annotations
def runRT_func(np.ndarray[double, ndim=3, mode="c"] histo not None,
                   diskmask=0):
    """
   
    """
    #this is not recommded way of passing variables:https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
    #but like: https://stackoverflow.com/questions/20182147/declaring-numpy-array-with-int-type-and-pass-array-pointer-to-c-code
    #I couldn't get the recommended way to work
    
    #Make sure variables are defined to be compatible with dustmap_c.c
iflen=strlen(inputfile)
ifn=n_elements(inputfile)
ifn=ifn[0]
inputfile4c=bytarr(total(iflen) + ifn)
strloc=0
for i=0,ifn-1 do begin
   inputfile4c[strloc:strloc+iflen[i]-1] = byte(inputfile[i])
   inputfile4c[strloc+iflen[i]] = byte(0)
   strloc += (iflen[i] + 1)
endfor
numfiles = long(n_elements(inputfile))
datatype = long(datatype)
userdistance = float(userdistance)
fovx = float(userfov[0])
fovy = float(userfov[1])
userpixelsize = float(userpixelsize)
userinclination = float(userinclination)
userlongitude = float(userlongitude)
userpa = float(userpa)
#wavel = float(wavel) #so python lambda isn't renamed
#nlambda = long(wav)
lstar = float(lstar)
Tstar = float(Tstar)
logg = float(logg)
kurucz = long(kurucz)
fstar = np.ndarray[float, ndim=1, mode="c"]#fltarr(nlambda)
rdust = float(rdust)
Tsublimate = float(Tsublimate)
lnkfile4c = bytarr(strlen(lnkfile)+1)
lnkfile4c[0:strlen(lnkfile)-1] = byte(lnkfile)
lnkfile4c[strlen(lnkfile)] = byte(0)
kuruczfile4c = bytarr(strlen(kurucz_file)+1)
kuruczfile4c[0:strlen(kurucz_file)-1] = byte(kurucz_file)
kuruczfile4c[strlen(kurucz_file)] = byte(0)
userscaling = double(userscaling)
iwa = float(iwa)
aitoff_flag = long(aitoff_flag)
densityhisto_flag = long(densityhisto_flag)
opticaldepth_flag = long(opticaldepth_flag)
scatteredlight_flag = long(scatteredlight_flag)
thermalemission_flag = long(thermalemission_flag)
verbose_flag = long(verbose_flag)
generic_oc_exp = float(qexp)
generic_oc_albedo = float(albedo)
HG_g = float(g)
oc_generic_flag = long(oc_generic_flag)
HG_flag = long(HG_flag)
extinction_flag = long(ext)
xshift = float(xshift)
effrdust = float(effrdust)
markx = float(markx)
marky = float(marky)
markz = float(markz)
markweight = float(markweight)
nmarks = long(nmarks)
markabs_flag = long(markabs_flag)
azavg = long(azavg)
distmask = float(distmask)
ncostheta = long(ncostheta)
internal_pfunc = np.ndarray(shape=(nlambda,ncostheta), mode="c")#fltarr(nlambda,ncostheta)
costheta =  np.ndarray(shape=(ncostheta))
nrdust=np.unique(rdust, return_inverse=False, return_counts=False, axis=None).size #nrdust = n_elements(uniq(rdust,sort(rdust)))
pfunc = np.ndarray[shape=(nrdust,nlambda,ncostheta), mode="c"]#pfunc = fltarr(nrdust,nlambda,ncostheta)
qabs  = np.ndarray[shape=(nrdust,nlambda), mode="c"]# fltarr(nrdust,nlambda)
qsca  = np.ndarray[shape=(nrdust,nlambda), mode="c"]# fltarr(nrdust,nlambda)

dustmap_c(<double*> histo.data,  histo.shape[1], histo.shape[0], inputfile4c, numfiles, datatype, \
                  userdistance, fovx, fovy, userpixelsize, userinclination, userlongitude, userpa, \
                  wavel.data, wavel.size, Lstar, Tstar, logg, kurucz, kuruczfile4c, fstar, \
                  rdust, Tsublimate, lnkfile4c, userscaling, \
                  iwa, aitoff_flag, densityhisto_flag, opticaldepth_flag, \
                  scatteredlight_flag, thermalemission_flag, verbose_flag, \
                  oc_generic_flag, generic_oc_exp, generic_oc_albedo, HG_flag, HG_g, \
                  extinction_flag, xshift, effrdust, \
                  markx, marky, markz, markweight, nmarks, markabs_flag, azavg, distmask, \
                  ncostheta, internal_pfunc, costheta, pfunc, qabs, qsca, nrdust)
   )


