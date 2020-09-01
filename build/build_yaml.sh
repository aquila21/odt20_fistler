#!/bin/bash

cd ../source/yaml/
rm -rf include lib
cd yaml-cpp-release-0.5.2/
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=../..
#cmake .. -DCMAKE_INSTALL_PREFIX:PATH=../.. -DCMAKE_CXX_COMPILER=/usr/bin/g++
make
make install
