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

rm f1
#macro2
g++ fluxtest.cpp numerical/numerical.cpp vroot/root.cpp -o2 -o f1 `root-config --cflags --glibs` -std=c++0x -pthread
#macro1

#g++ muon.cpp gPT/GPT.cpp mROOT/mroot.cpp -o2 -o mugo `root-config --cflags --glibs` -std=c++0x -pthread

echo "done!"
 

