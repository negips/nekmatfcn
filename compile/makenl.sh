#!/bin/sh

./makenek clean
rm SIZE
ln -s SIZE_NL SIZE
./makenek ellipse

#mv nek5000 ../big_cossu/re23_nonlinear/
