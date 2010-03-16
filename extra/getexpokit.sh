#!/bin/bash

wget http://www.maths.uq.edu.au/expokit/expokit.tar.gz
tar -xzvf expokit.tar.gz
cd expokit/fortran


# Based on: http://sf.anu.edu.au/~mhk900/Python_Workshop/short.pdf
# Steps
# 1) Update references to matve
patch -p 0 < ../../expokit.f.patch

# 2) Generate header
f2py -m expokit -h expokit.pyf expokit.f

# 3) Patch header
patch -p 0 < ../../expokit.pyf.patch

# 4) Compile package
f2py -c expokit.pyf expokit.f blas.f lapack.f

# The shiny new "expokit.so" should be ready for you!
cp expokit.so ../../.
