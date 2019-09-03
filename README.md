# dustmapPy
python interface to dustmap

Depends on cython.

Test configuration:

## Notes:

Test system configuraton: 

* Cython: 0.26.1
* Python: Python 3.5.4 |Anaconda custom (64-bit)|

##  setup:

clone this repo:

    git clone https://github.com/douglase/dustmapPy

clone the [dustmap](https://github.com/douglase/dustmap) repo and add it as a submodule (folder) in your dustmappy directory:

    cd dustmapPy
    git submodule init
    git submodule update

build the repo:

    python setup.py build_ext --inplace

## test and example:


```
python test_map.py
```

Compares the results of an exampled dustmapPy call (Left panel of figure below) to an analogous IDL call of dustmap (middle panel) sample_dustmap_call.pro) which can be found in the dustmap repo and shows the residual (Right panel).

![DustmapPy call, python call, residual](IDL-python.png)
