#!/bin/bash
set -e

cd tests
testflo -v --coverage --coverpkg pyspline
