#!/bin/bash
set -e
mkdir build/
cd build/
cmake ../
make check
pip install .
