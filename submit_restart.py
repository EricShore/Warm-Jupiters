#!/usr/bin/env python
import os
from setting import *
from submit_nodes import submit
def main():
    rmin=0
    rmax=max_runs*nnodes
    inpath=rundir
    for i in xrange(rmin,rmax):
        runfile=inpath+"rebound%.4d.txt" % (i+1)
        pklfile=inpath+"runrebound%.4d.pkl" % (i+1)
        fout = open(runfile,mode= "r")
	lines = fout.readlines()
	fout.close()
        j = 1
        while(True):#looks for the last line with text (as opposed to an empty one)
            split_line = lines[-j].split()
            if len(split_line)>0:
                break
            j=j+1
            if (j>=10000):
                print "max loops"
                break
        endtime = float(lines[-j].split()[0])/(2*np.pi) #time at the end of the sim
	if endtime<t_max:#restats all sims where the current time is less than the target one
       #if not os.path.exists(pklfile):
            subfile=basename+"qsubrebound_%d" % i
            submit(basename,subfileauto="qsubrebound_%d" % i, start=1+max_runs*i)
            os.system('qsub %s' % subfile)
        continue
    return


if __name__=='__main__':
    main()
