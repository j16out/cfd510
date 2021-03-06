#!/bin/bash
 
###########################################################
#
#
# M9 sim setup
#
#
###########################################################
 
# compile m9 sim
echo "Starting compilation"
echo "Reminder: Make sure root is sourced"
#rm mtest2

rm exlaplace
#macro2
g++ mainl.cpp numerical/numerical.cpp vroot/root.cpp -o2 -o exlaplace `root-config --cflags --glibs` -std=c++0x -pthread

#macro1

#g++ muon.cpp gPT/GPT.cpp mROOT/mroot.cpp -o2 -o mugo `root-config --cflags --glibs` -std=c++0x -pthread

echo "done!"
 
#valgrind --leak-check=full -v ./exlaplace
