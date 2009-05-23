#/bin/bash
#Makefile for pySpline

#First make the source files for the real version
cd src

make

cd ../

#Now f2py just the two functions we need
f2py  --fcompiler=intele --f90flags=-r8 -c -m pyspline src/tensbs/b2ink.f src/tensbs/b2val.f src/libtensbs.a

#Now make the source files for the complex version

mv pyspline.so ./python

cd src_cs

make

cd ../

f2py  -c --fcompiler=intele  --f90flags=-r8 -I./src_cs -m pyspline_cs src_cs/tensbs/c_b2ink.f src_cs/tensbs/c_b2val.f src_cs/libtensbs_cs.a

mv pyspline_cs.so ./python