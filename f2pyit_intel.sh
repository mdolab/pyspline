#/bin/bash
#Makefile for pySpline

#First make the source files 
cd src

make
if [ ! $? -eq 0 ]; then
    exit
fi

cd ../

#Now f2py just the functions we need in the pyf
f2py  --fcompiler=intel --f90flags=-r8 -c -m pyspline src/pyspline.pyf src/libspline.a -lblas

mv pyspline.so ./python
