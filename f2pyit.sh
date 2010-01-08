#/bin/bash
#Makefile for pySpline

#First make the source files for the real version
cd src

make
if [ ! $? -eq 0 ]; then
    exit
fi

cd ../

#Now f2py just the functions we need in the pyf
f2py  --fcompiler=intel --f90flags=-r8 -c -m pyspline src/pyspline.pyf src/libspline.a
mv pyspline.so ./python

# #Now make the source files for the complex version

# cd src_cs

# make
# if [ ! $? -eq 0 ]; then
#     exit
# fi

# cd ../

# f2py  --fcompiler=intel --f90flags=-r8 -I./src_cs -c -m pyspline_cs src_cs/pyspline_cs.pyf src_cs/libspline_cs.a 
# mv pyspline_cs.so ./python