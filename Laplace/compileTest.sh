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

rm BCtest
#macro2
g++ BC-test.c -o2 -o BCtest `root-config --cflags --glibs` -std=c++0x -pthread

#macro1

#g++ muon.cpp gPT/GPT.cpp mROOT/mroot.cpp -o2 -o mugo `root-config --cflags --glibs` -std=c++0x -pthread

echo "done!"
 
#valgrind --leak-check=full -v ./exlaplace
