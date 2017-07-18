#!/bin/bash

dir=`pwd`
cd $dir/src
aclocal
autoconf
autoheader
automake -a
./configure
make
