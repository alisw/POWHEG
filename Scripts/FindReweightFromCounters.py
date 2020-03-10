#!/usr/bin/python

# Usage:
# After a POWHEG parallel run, if ubexcess_correct was set to 1 in the powheg.input file:
# $ ./FindReweightFromCounters.py pwgcounters-st4-????.dat
# It will examine the cross sections and print the correction factor that accounts
# for the cross section estimate performed during the generation of the events.



# regular expressions
import math
import re
import sys
import array
# debugger
#import pdb
#pdb.set_trace()

if(len(sys.argv) < 2):
    print(len(sys.argv))
    print('usage: FindReweightFromCounters.py counter files')
    exit(1)

for which in ['btilde','remnant']:
    tot_numpts=0
    for file in sys.argv[1:]:
        counterfile=open(file,"r")
        print "examining ",file
        fullline=counterfile.read()
        match=re.search('^ *'+which+' cross section used:[ =]*([^$ ]*) *$',fullline,re.M)
        if(match == None):
            print file+' does not contain a "'+which+' cross section used:" line'
            print 'skipping '+which
            break
        sig_used=float(match.group(1))
        match=re.search('^ *'+which+' cross section estimate:[ =]*([^$ ]*) *$',fullline,re.M)
        if(match == None):
            print file+' does not contain a "'+which+' cross section estimate:" line'
            print 'skipping '+which
            break
        sig_estim=float(match.group(1))
        match=re.search('^ *'+which+' cross section error estimate:[ =]*([^$ ]*) *$',fullline,re.M)
        if(match == None):
            print file+' does not contain a "'+which+' cross section error estimate:" line'
            print 'skipping '+which
            break
        sig_estimerr=float(match.group(1))
        match=re.search('^ *'+which+' cross section estimate num. points:[ =]*([^$ ]*) *$',fullline,re.M)
        if(match == None):
            print file+' does not contain a "'+which+' cross section estimate num. points:" line'
            print 'skipping '+which
            break
        sig_numpts=float(match.group(1))
        match=re.search('^ *'+which+' bound violation correction factor:[ =]*([^$ ]*) *$',fullline,re.M)
        if(match != None):
            sig_ubcorr=float(match.group(1))
        else:
            sig_ubcorr=1
            sig_numpts=float(match.group(1))
        if  tot_numpts == 0 :
            tot_numpts=sig_numpts
            tot_estim=sig_estim
            tot_estimerr=sig_estimerr
            tot_used=sig_used
            tot_ubcorr=sig_ubcorr
            tot_nfiles=1
        else:
            if tot_used != sig_used:
                print "the used cross section in the current file",
                "is inconsistent with the previous ones, exiting ..."
                exit -1
    
            sumsq=(tot_estimerr**2*tot_numpts+tot_estim**2)*tot_numpts + \
            (sig_estimerr**2*sig_numpts+sig_estim**2)*sig_numpts
            tot_estim=(tot_estim*tot_numpts+sig_estim*sig_numpts)/(tot_numpts+sig_numpts)
            tot_numpts=tot_numpts+sig_numpts
            tot_estimerr=math.sqrt((sumsq/tot_numpts-tot_estim**2)/tot_numpts)
            if tot_ubcorr != 1 or sig_ubcorr != 1:
                tot_ubcorr=(tot_ubcorr*tot_nfiles+sig_ubcorr)/(tot_nfiles+1)
                tot_nfiles=tot_nfiles+1
            
            counterfile.close()
    else:   # when for file loop is over
        print 'Cross section used for '+which+' event generation:',tot_used
        print 'Cross section computed on the fly for '+which+' event generation:',tot_estim,'+-',tot_estimerr
        if tot_ubcorr != 1:
            print 'Weight increment factor due to corrections for upper bound violation in '+which+' events:',tot_ubcorr
            print 'Total correction factor for '+which+' events:',tot_estim/tot_used/tot_ubcorr
