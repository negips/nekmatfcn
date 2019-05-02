#!/bin/sh

./makenek clean
rm SIZE
ln -s SIZE_L SIZE
./makenek lu

cp nek5000 ../lu_2016/linear/
cp nek5000 ../lu_2016/nonlinear_test/
