#!/bin/bash
set -e
mkdir build/
cd build/
if [[ $COMPILERS == "intel" ]]; then
	cmake -DCMAKE_CXX_COMPILER=icpc \
		  -DCMAKE_C_COMPILER=icc \
		  -DCMAKE_Fortran_COMPILER=ifort ../
else
	cmake ../
fi
make check
make install
cd ../
pip install .
