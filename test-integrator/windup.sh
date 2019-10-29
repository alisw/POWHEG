#!/bin/bash

../../Scripts/FindReweightFromCounters.py pwgcounters-st4-*.dat | grep -v examining > counters.dat

mergedata 1 pwgpwhgalone-output[0-9][0-9][0-9][0-9].top -o pwgpwhgalone-output.top
