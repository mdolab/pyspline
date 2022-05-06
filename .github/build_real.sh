#!/bin/bash
set -e
mkdir build/
cd build/
cmake ../
make check
make install
cd ../
pip install .
