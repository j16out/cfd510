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

rm T1
#macro2
g++ BThomasTest.cpp numerical/numerical.cpp vroot/root.cpp -o2 -o T1 `root-config --cflags --glibs` -std=c++0x -pthread
#macro1

#g++ muon.cpp gPT/GPT.cpp mROOT/mroot.cpp -o2 -o mugo `root-config --cflags --glibs` -std=c++0x -pthread

echo "done!"
 

