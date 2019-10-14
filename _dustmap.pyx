
"""
    edouglas@mit.edu
    Modified from scipy-lectures.org:
    http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id10
    CC Attribution license: https://github.com/scipy-lectures/scipy-lecture-notes/blob/master/LICENSE.rst

    also see:
    https://github.com/cython/cython/blob/master/docs/src/tutorial/array.rst#data-fields
    https://stackoverflow.com/questions/3046305/simple-wrapping-of-c-code-with-cython#3071942
    https://stackoverflow.com/questions/47005382/cython-cannot-convert-python-object-error
    https://stackoverflow.com/questions/17855032/passing-and-returning-numpy-arrays-to-c-methods-via-cython/18176741#18176741
    https://stackoverflow.com/questions/24226001/importerror-dynamic-module-does-not-define-init-function-initfizzbuzz
    https://stackoverflow.com/questions/17014379/cython-cant-convert-python-object-to-double
    https://stackoverflow.com/questions/13669961/convert-python-object-to-cython-pointer#13728234
https://stackoverflow.com/questions/23435756/passing-numpy-integer-array-to-c-code?rq=1
"""

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import astropy.units as u

import cython
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example)
np.import_array()

#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef extern from "dustmap/dmdatatypes.h":
    void cos_doubles (double * in_array,
     double * out_array,
     int size)
#    void runRT(double * in_array,double * out_array, int in_col, int in_row)

# cdefine the signature of our c function
# these variables are cast as if they are C
cdef extern from "dustmap/dmdatatypes.h":
  void dustmap(double * hist,
          int histoxsize,
           int histoysize,
           char * inputfile4c,
           int numfiles,
           int datatype,
           float distance,
           float fovx,
           float fovy,
           float pixelsize,
           float inclination,
           float longitude,
           float pa,
           float * wavel,
           int nlambda,
           float Lstar,
           float Tstar,
           float logg,
           int kurucz,
           char* kuruczfile,
           float * fstar,
           float * rdust,
           float Tsublimate,
           char* lnkfile4c,
           double* scaling,
           float iwa,
           int aitoff_flag,
           int densityhisto_flag,
           int opticaldepth_flag,
           int scatteredlight_flag,
           int thermalemission_flag,
           int verbose_flag,
           float oc_generic_flag,
           float generic_oc_exp,
           float generic_oc_albedo,
           int HG_flag,
           float HG_g,
           int extinction_flag,
           float xshift,
           float effrdust,
           float* markx,
           float* marky,
           float* markz,
           float* markweight,
           int nmarks,
           int markabs_flag,
           int azavg,
           float distmask,
           int ncostheta,
           float* internal_pfunc,
           float* costheta,
           float* pfunc,
           float* qabs,
           float* qsca,
           int nrdust)

# create the wrapper function in python, with numpy type annotations
def cos_doubles_func(np.ndarray[double, ndim=1, mode="c"] in_array not None,
                     np.ndarray[double, ndim=1, mode="c"] out_array not None):
    cos_doubles(<double*> np.PyArray_DATA(in_array),
                <double*> np.PyArray_DATA(out_array),
                in_array.shape[0])
def dustmap_func(          inputfile,
          np.ndarray[double, ndim=2, mode="c"] hist not None,
          distance,
          fov,
          pixelsize,
          numpixels,
          np.ndarray[float, ndim=1, mode="c"] wavel, #was lambda
          Tstar,
          np.ndarray[float, ndim=1, mode="c"] fstar not None,
          composition,
          #np.ndarray[float, ndim=1, mode="c"] scaling not None,
          Tsublimate=1e10,
          datatype=1,
          generic_oc_albedo=0,
          generic_oc_exp=0,
          effrdust=0.0,

          markx= np.ascontiguousarray(np.array([0])),
          marky= np.ascontiguousarray(np.array([0])),
          markz= np.ascontiguousarray(np.array([0])),
          markweight= np.ascontiguousarray(np.array([0])),
          scaling= np.ascontiguousarray(np.array([1.0],dtype=np.double)),
          rdust=0,
          nmarks=0,
          azavg=1,
          distmask=0,
          vga=0,
          hd=0,
          allsky=0,
          xshift=0,
          costheta=0,
          ncostheta= 0,
          inclination=0,
          longitude=0,
          pa=0,
          logg=0,
          Lstar=0,
          kurucz=0,
          kurucz_file = 'dustmap/fp00k2odfnew.pck',
          iwa=0.0*u.arcsecond,
          qexp=0,
          aitoff_flag=0,
          densityhisto_flag=0,
          opticaldepth_flag=0,
          scatteredlight_flag=0,
          thermalemission_flag=0,
          verbose_flag=1,
          noprint=0,
          HG_g=0,
          HG_flag=0,
          extinction_flag=0, #ext,
          markabs_flag=0,
          oc_generic_flag = 0,
          ):
          """
          #this is not recommded way of passing variables:https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
          #but like: https://stackoverflow.com/questions/20182147/declaring-numpy-array-with-int-type-and-pass-array-pointer-to-c-code
          #I couldn't get the recommended way to work

          #Make sure variables are defined to be compatible with dustmap_c.c

          """

          if noprint:
            reduceprint=1

          #kurucz_file = dustmappath+'fp00k2odfnew.pck'
          #TEMPORARY - assume running in same path
          kuruczfile4c = kurucz_file.encode()

          #; Convert fov and pixelsize from degrees to mas if necessary

          #eventuallydo unit conversions with astropy units

          fov = [angle.to(u.milliarcsecond).value for angle in fov]
          pixelsize = int(np.round(pixelsize.to(u.milliarcsecond).value))
          iwa = iwa.to(u.milliarcsecond).value

          #; Convert distance from AU to pc if necessary
          #if keyword_set(au) then distance /= 206265.0 ;convert from AU to parsec
          #if keyword_set(au) and keyword_set(distmask) then distmask /= 206265.0 ;convert from AU to parsec



          #dustmap expects a list of files seperated by \0
          #
          #iflen=len(inputfile)
          #ifn=n_elements(inputfile)
          #ifn=ifn[0]
          #inputfile4c=bytarr(total(iflen) + ifn)
          #3strloc=0
          #for i in range(ifn):
          #  inputfile4c[strloc:strloc+iflen[i]-1] = byte(inputfile[i])
          #    inputfile4c[strloc+iflen[i]] = byte(0)
          #      strloc += (iflen[i] + 1)
          '''if (len(inputfile)>1) & (isintance(inputfile,list)):
            print("Warning, not tested with multiple input yet")
            inputfile4c="\0".join(files)
            symlist=["\'","]","["," ",".",","]
            for sym in symlist:
              #print(sym)
              inputfile4c=inputfile4c.replace(sym,"")
              numfiles = np.int(n_elements(inputfile))
          else:
            numfiles=1
            inputfile4c = inputfile'''
          numfiles=1
          inputfile4c = inputfile.encode()
          #distance = np.float_(distance)
          if len(fov) >1:
            fovx = np.float64(fov[0])
            fovy = np.float64(fov[1])
          else:
            fovx=fov[0]
            fovy=fov[0]
          #pixelsize = np.float_(pixelsize)
          #inclination = np.float_(inclination)
          #longitude = np.float_(longitude)
          #pa = np.float_(pa)
          #wavel = np.float_(wavel) #so python lambda isn't renamed
          nlambda = np.size(wavel)
          #lstar = np.float_(lstar)
          #Tstar = np.float_(Tstar)
          #logg = np.float_(logg)
          #kurucz = np.int(kurucz)
          #fstar = np.ndarray(np.float32, ndim=1, mode="c")#fltarr(nlambda)
          rdust = np.float32(rdust)
          #Tsublimate = np.float_(Tsublimate)
          #lnkfile4c = bytarr(strlen(lnkfile)+1)
          lnkfile4c = str('dustmap/lnkfiles/'+composition+'.lnk').encode()
          #lnkfile4c[0:strlen(lnkfile)-1] = byte(lnkfile)
          #lnkfile4c[strlen(lnkfile)] = byte(0)
          #kuruczfile4c = 'fp00k2odfnew.pck'.encode()# bytarr(strlen(kurucz_file)+1)
          #kuruczfile4c[0:strlen(kurucz_file)-1] = byte(kurucz_file)
          #kuruczfile4c[strlen(kurucz_file)] = byte(0)
          #scaling = np.float__(scaling)
          #iwa = np.float_(iwa)
          #aitoff_flag = np.int(aitoff_flag)
          #densityhisto_flag = np.int(densityhisto_flag)
          #opticaldepth_flag = np.int(opticaldepth_flag)
          #scatteredlight_flag = np.int(scatteredlight_flag)
          #thermalemission_flag = np.int(thermalemission_flag)
          #verbose_flag = np.int(verbose_flag)
          #generic_oc_exp = np.float_(qexp)
          #generic_oc_albedo = np.float_(albedo)
          #HG_g = np.float_(g)
          #oc_generic_flag = np.int(oc_generic_flag)
          #HG_flag = np.int(HG_flag)
          #extinction_flag = np.int(ext)
          #xshift = np.float_(xshift)
          #effrdust = np.float_(effrdust)
          #markx = np.float_(markx)
          #marky = np.float_(marky)
          #markz = np.float_(markz)
          #markweight = np.float_(markweight)
          #nmarks = np.int(nmarks)
          #markabs_flag = np.int(markabs_flag)
          #azavg = np.int(azavg)
          #distmask = np.float_(distmask)
          #ncostheta = np.int(ncostheta)
          internal_pfunc = np.ascontiguousarray(np.ndarray(shape=(nlambda,ncostheta)))#, mode="c"))#fltarr(nlambda,ncostheta)
          costheta =  np.ndarray(shape=(ncostheta))
          nrdust=np.unique(rdust, return_inverse=False, return_counts=False, axis=None).size #nrdust = n_elements(uniq(rdust,sort(rdust)))
          pfunc = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda,ncostheta)))#, mode="c")#pfunc = fltarr(nrdust,nlambda,ncostheta)
          qabs  = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda)))#, mode="c")# fltarr(nrdust,nlambda)
          qsca  = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda)))#, mode="c")# fltarr(nrdust,nlambda)

          #args=["1","2","3"]
          hist  = np.ascontiguousarray(hist)
          fstar  = np.ascontiguousarray(fstar)
          rdust  = np.ascontiguousarray(rdust)
          cos_doubles(<double*> np.PyArray_DATA(hist),
                      <double*> np.PyArray_DATA(hist),
                      hist.shape[0])
          print([Lstar,Tstar,wavel])
          print(["distance",distance])
          dustmap(<double*> np.PyArray_DATA(hist),
                  <int> hist.shape[1],
                  <int> hist.shape[0],
                  <char*> inputfile4c,
                  <int> numfiles,
                  <int> datatype,
                  <float> distance,
                  <float> fovx,
                  <float> fovy,
                  <float> pixelsize,
                  <float> inclination,
                  <float> longitude,
                  <float> pa,
                  <float*> np.PyArray_DATA(wavel),
                  <int> nlambda,
                  <float> Lstar,
                  <float> Tstar,
                  <float> logg,
                  <int> kurucz,
                  <char*> kuruczfile4c,
                  <float*> np.PyArray_DATA(fstar),
                  <float*> np.PyArray_DATA(rdust),
                  <float> Tsublimate,
                  <char*> lnkfile4c,
                  <double*> np.PyArray_DATA(scaling),
                  <float> iwa,
                  <int> aitoff_flag,
                  <int> densityhisto_flag,
                  <int> opticaldepth_flag,
                  <int> scatteredlight_flag,
                  <int> thermalemission_flag,
                  <int> verbose_flag,
                  <int> oc_generic_flag,
                  <float> generic_oc_exp,
                  <float> generic_oc_albedo,
                  <int> HG_flag,
                  <float> HG_g,
                  <int> extinction_flag,
                  <float> xshift,
                  <float> effrdust,
                  <float*> np.PyArray_DATA(markx),
                  <float*> np.PyArray_DATA(marky),
                  <float*> np.PyArray_DATA(markz),
                  <float*> np.PyArray_DATA(markweight),
                  <int> nmarks,
                  <int> markabs_flag,
                  <int> azavg,
                  <float> distmask,
                  <int> ncostheta, #npfunc in dustmap_c.c
                  <float*> np.PyArray_DATA(internal_pfunc),
                  <float*> np.PyArray_DATA(costheta),
                  <float*> np.PyArray_DATA(pfunc),
                  <float*> np.PyArray_DATA(qabs),
                  <float*> np.PyArray_DATA(qsca),
                  <int> nrdust)


#functions to skip writing to disk and save processing time:
cdef extern from "dustmap/dmdatatypes.h":
  void dustmap_memory(double * hist,
            int histoxsize,
            int histoysize,

            #Yinzi's edits for passing data directly
            #char * inputfile4c,

            int * pid, # particle ids
            float * x_coords, #particle x coords
            float * y_coords, #particle y coords
            float * z_coords, #particle z coords
            double * particle_intensity, #particle number density
            int nentries,

            #end Yinzi's edits

            int numfiles,
            int datatype,
            float distance,
            float fovx,
            float fovy,
            float pixelsize,
            float inclination,
            float longitude,
            float pa,
            float * wavel,
            int nlambda,
            float Lstar,
            float Tstar,
            float logg,
            int kurucz,
            char* kuruczfile,
            float * fstar,
            float * rdust,
            float Tsublimate,
            char* lnkfile4c,
            double* scaling,
            float iwa,
            int aitoff_flag,
            int densityhisto_flag,
            int opticaldepth_flag,
            int scatteredlight_flag,
            int thermalemission_flag,
            int verbose_flag,
            float oc_generic_flag,
            float generic_oc_exp,
            float generic_oc_albedo,
            int HG_flag,
            float HG_g,
            int extinction_flag,
            float xshift,
            float effrdust,
            float* markx,
            float* marky,
            float* markz,
            float* markweight,
            int nmarks,
            int markabs_flag,
            int azavg,
            float distmask,
            int ncostheta,
            float* internal_pfunc,
            float* costheta,
            float* pfunc,
            float* qabs,
            float* qsca,
            int nrdust)

# create the wrapper function in python, with numpy type annotations
def cos_doubles_func(np.ndarray[double, ndim=1, mode="c"] in_array not None,
                     np.ndarray[double, ndim=1, mode="c"] out_array not None):
    cos_doubles(<double*> np.PyArray_DATA(in_array),
                <double*> np.PyArray_DATA(out_array),
                in_array.shape[0])
def dustmap_memory_func(#Yinzi's edits
          #inputfile,
          pid,
          x_coords,
          y_coords,
          z_coords,
          particle_intensity,
          nentries,
          #end Yinzi's edits


          np.ndarray[double, ndim=2, mode="c"] hist not None,
          distance,
          fov,
          pixelsize,
          numpixels,
          np.ndarray[float, ndim=1, mode="c"] wavel, #was lambda
          Tstar,
          np.ndarray[float, ndim=1, mode="c"] fstar not None,
          composition,
          #np.ndarray[float, ndim=1, mode="c"] scaling not None,
          Tsublimate=1e10,
          datatype=1,
          generic_oc_albedo=0,
          generic_oc_exp=0,
          effrdust=0.0,

          markx= np.ascontiguousarray(np.array([0])),
          marky= np.ascontiguousarray(np.array([0])),
          markz= np.ascontiguousarray(np.array([0])),
          markweight= np.ascontiguousarray(np.array([0])),
          scaling= np.ascontiguousarray(np.array([1.0],dtype=np.double)),
          rdust=0,
          nmarks=0,
          azavg=1,
          distmask=0,
          vga=0,
          hd=0,
          allsky=0,
          xshift=0,
          costheta=0,
          ncostheta= 0,
          inclination=0,
          longitude=0,
          pa=0,
          logg=0,
          Lstar=0,
          kurucz=0,
          kurucz_file = 'dustmap/fp00k2odfnew.pck',
          iwa=0.0*u.arcsecond,
          qexp=0,
          aitoff_flag=0,
          densityhisto_flag=0,
          opticaldepth_flag=0,
          scatteredlight_flag=0,
          thermalemission_flag=0,
          verbose_flag=1,
          noprint=0,
          HG_g=0,
          HG_flag=0,
          extinction_flag=0, #ext,
          markabs_flag=0,
          oc_generic_flag = 0,
          ):
          """
          #this is not recommded way of passing variables:https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
          #but like: https://stackoverflow.com/questions/20182147/declaring-numpy-array-with-int-type-and-pass-array-pointer-to-c-code
          #I couldn't get the recommended way to work

          #Make sure variables are defined to be compatible with dustmap_c.c

          """

          if noprint:
            reduceprint=1

          #kurucz_file = dustmappath+'fp00k2odfnew.pck'
          #TEMPORARY - assume running in same path
          kuruczfile4c = kurucz_file.encode()

          #; Convert fov and pixelsize from degrees to mas if necessary

          #eventuallydo unit conversions with astropy units

          fov = [angle.to(u.milliarcsecond).value for angle in fov]
          pixelsize = int(np.round(pixelsize.to(u.milliarcsecond).value))
          iwa = iwa.to(u.milliarcsecond).value

          #; Convert distance from AU to pc if necessary
          #if keyword_set(au) then distance /= 206265.0 ;convert from AU to parsec
          #if keyword_set(au) and keyword_set(distmask) then distmask /= 206265.0 ;convert from AU to parsec



          #dustmap expects a list of files seperated by \0
          #
          #iflen=len(inputfile)
          #ifn=n_elements(inputfile)
          #ifn=ifn[0]
          #inputfile4c=bytarr(total(iflen) + ifn)
          #3strloc=0
          #for i in range(ifn):
          #  inputfile4c[strloc:strloc+iflen[i]-1] = byte(inputfile[i])
          #    inputfile4c[strloc+iflen[i]] = byte(0)
          #      strloc += (iflen[i] + 1)
          '''if (len(inputfile)>1) & (isintance(inputfile,list)):
            print("Warning, not tested with multiple input yet")
            inputfile4c="\0".join(files)
            symlist=["\'","]","["," ",".",","]
            for sym in symlist:
              #print(sym)
              inputfile4c=inputfile4c.replace(sym,"")
              numfiles = np.int(n_elements(inputfile))
          else:
            numfiles=1
            inputfile4c = inputfile'''
          numfiles=1
          #inputfile4c = inputfile.encode()
          #distance = np.float_(distance)
          if len(fov) >1:
            fovx = np.float64(fov[0])
            fovy = np.float64(fov[1])
          else:
            fovx=fov[0]
            fovy=fov[0]
          #pixelsize = np.float_(pixelsize)
          #inclination = np.float_(inclination)
          #longitude = np.float_(longitude)
          #pa = np.float_(pa)
          #wavel = np.float_(wavel) #so python lambda isn't renamed
          nlambda = np.size(wavel)
          #lstar = np.float_(lstar)
          #Tstar = np.float_(Tstar)
          #logg = np.float_(logg)
          #kurucz = np.int(kurucz)
          #fstar = np.ndarray(np.float32, ndim=1, mode="c")#fltarr(nlambda)
          rdust = np.float32(rdust)
          #Tsublimate = np.float_(Tsublimate)
          #lnkfile4c = bytarr(strlen(lnkfile)+1)
          lnkfile4c = str('dustmap/lnkfiles/'+composition+'.lnk').encode()
          #lnkfile4c[0:strlen(lnkfile)-1] = byte(lnkfile)
          #lnkfile4c[strlen(lnkfile)] = byte(0)
          #kuruczfile4c = 'fp00k2odfnew.pck'.encode()# bytarr(strlen(kurucz_file)+1)
          #kuruczfile4c[0:strlen(kurucz_file)-1] = byte(kurucz_file)
          #kuruczfile4c[strlen(kurucz_file)] = byte(0)
          #scaling = np.float__(scaling)
          #iwa = np.float_(iwa)
          #aitoff_flag = np.int(aitoff_flag)
          #densityhisto_flag = np.int(densityhisto_flag)
          #opticaldepth_flag = np.int(opticaldepth_flag)
          #scatteredlight_flag = np.int(scatteredlight_flag)
          #thermalemission_flag = np.int(thermalemission_flag)
          #verbose_flag = np.int(verbose_flag)
          #generic_oc_exp = np.float_(qexp)
          #generic_oc_albedo = np.float_(albedo)
          #HG_g = np.float_(g)
          #oc_generic_flag = np.int(oc_generic_flag)
          #HG_flag = np.int(HG_flag)
          #extinction_flag = np.int(ext)
          #xshift = np.float_(xshift)
          #effrdust = np.float_(effrdust)
          #markx = np.float_(markx)
          #marky = np.float_(marky)
          #markz = np.float_(markz)
          #markweight = np.float_(markweight)
          #nmarks = np.int(nmarks)
          #markabs_flag = np.int(markabs_flag)
          #azavg = np.int(azavg)
          #distmask = np.float_(distmask)
          #ncostheta = np.int(ncostheta)
          internal_pfunc = np.ascontiguousarray(np.ndarray(shape=(nlambda,ncostheta)))#, mode="c"))#fltarr(nlambda,ncostheta)
          costheta =  np.ndarray(shape=(ncostheta))
          nrdust=np.unique(rdust, return_inverse=False, return_counts=False, axis=None).size #nrdust = n_elements(uniq(rdust,sort(rdust)))
          pfunc = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda,ncostheta)))#, mode="c")#pfunc = fltarr(nrdust,nlambda,ncostheta)
          qabs  = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda)))#, mode="c")# fltarr(nrdust,nlambda)
          qsca  = np.ascontiguousarray(np.ndarray(shape=(nrdust,nlambda)))#, mode="c")# fltarr(nrdust,nlambda)

          #args=["1","2","3"]
          hist  = np.ascontiguousarray(hist)
          fstar  = np.ascontiguousarray(fstar)
          rdust  = np.ascontiguousarray(rdust)
          # cos_doubles(<double*> np.PyArray_DATA(hist),
          #             <double*> np.PyArray_DATA(hist),
          #             hist.shape[0])
          print([Lstar,Tstar,wavel])
          print(["distance",distance])
          dustmap_memory(<double*> np.PyArray_DATA(hist),
                  <int> hist.shape[1],
                  <int> hist.shape[0],
                  #Yinzi edits
                  #<char*> inputfile4c,
                  <int*> np.PyArray_DATA(pid),
                  <float*> np.PyArray_DATA(x_coords),
                  <float*> np.PyArray_DATA(y_coords),
                  <float*> np.PyArray_DATA(z_coords),
                  <double*> np.PyArray_DATA(particle_intensity),
                  <int> nentries,
                  #end Yinzi edits
                  <int> numfiles,
                  <int> datatype,
                  <float> distance,
                  <float> fovx,
                  <float> fovy,
                  <float> pixelsize,
                  <float> inclination,
                  <float> longitude,
                  <float> pa,
                  <float*> np.PyArray_DATA(wavel),
                  <int> nlambda,
                  <float> Lstar,
                  <float> Tstar,
                  <float> logg,
                  <int> kurucz,
                  <char*> kuruczfile4c,
                  <float*> np.PyArray_DATA(fstar),
                  <float*> np.PyArray_DATA(rdust),
                  <float> Tsublimate,
                  <char*> lnkfile4c,
                  <double*> np.PyArray_DATA(scaling),
                  <float> iwa,
                  <int> aitoff_flag,
                  <int> densityhisto_flag,
                  <int> opticaldepth_flag,
                  <int> scatteredlight_flag,
                  <int> thermalemission_flag,
                  <int> verbose_flag,
                  <int> oc_generic_flag,
                  <float> generic_oc_exp,
                  <float> generic_oc_albedo,
                  <int> HG_flag,
                  <float> HG_g,
                  <int> extinction_flag,
                  <float> xshift,
                  <float> effrdust,
                  <float*> np.PyArray_DATA(markx),
                  <float*> np.PyArray_DATA(marky),
                  <float*> np.PyArray_DATA(markz),
                  <float*> np.PyArray_DATA(markweight),
                  <int> nmarks,
                  <int> markabs_flag,
                  <int> azavg,
                  <float> distmask,
                  <int> ncostheta, #npfunc in dustmap_c.c
                  <float*> np.PyArray_DATA(internal_pfunc),
                  <float*> np.PyArray_DATA(costheta),
                  <float*> np.PyArray_DATA(pfunc),
                  <float*> np.PyArray_DATA(qabs),
                  <float*> np.PyArray_DATA(qsca),
                  <int> nrdust)
