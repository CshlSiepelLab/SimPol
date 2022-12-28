#!/bin/sh
git submodule init
git submodule update
apt-get update
apt-get -y install libhdf5-dev
apt-get -y install libboost-all-dev
apt-get -y install libomp-dev