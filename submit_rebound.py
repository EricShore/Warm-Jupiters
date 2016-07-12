import time as timing
import setting as st
from setting import *
import os
import numpy as np
import scipy as sp
import rebound
import sys
import pickle
from util import check_for_bad_dt,set_hill
from init import init_orbit, read_init,read_init2, init_orbit2, init_orbit3
from submit_nodes import submit
from WJgrow import add_mass
def callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=0,hashes=None):
    #initialize a rebound run
    sim = rebound.Simulation()
    sim.t=t
    sim.add(m=1., r=0.005, name="star")
    for i in range(len(mass_pl)):
        sim.add(m = mass_pl[i], r = r_pl[i], a=a_pl[i], e=e_pl[i], inc=i_pl[i], Omega=Omega_pl[i], omega=omega_pl[i], M = M_pl[i])
        if hashes == None:
            sim.particles[-1].hash = i+1
        else: #i.e. when the simulation gets reloaded
            sim.particles[-1].hash = int(hashes[i])
    sim.move_to_com()
    return sim

def orbit2str(particle):
    #write the orbit elements to a string
    orbit=particle.orbit
    string="%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"%(particle.m,particle.r,orbit.a,orbit.e,orbit.inc,orbit.Omega,orbit.omega,orbit.M)
    #print string
    return string

def saveorbit(outfile,sim):
    #save the orbits elements to a file
    fout=open(outfile,mode="a")
    E=sim.calculate_energy()
    for p in sim.particles:
        if p.hash==sim.particles[0].hash:
            continue
        line=orbit2str(p)
        fout.write("%f %15.12f %d %s\n" % (sim.t,E, p.hash,line))
    fout.close()
    return

def saveorbit2 (outfile, sim): #used for checking the smallest distance between all the particles
    #save the orbits elements to a file
    fout=open(outfile,mode="a")
    E=sim.calculate_energy()
                
    for p1 in sim.particles:
        if p1.hash==sim.particles[0].hash:
            continue
        mind2 = 1.e6
        minid = 0
        for p2 in sim.particles:
            if p1.hash == p2.hash:
                continue
            d2 = (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z)
            if d2 < mind2:
                mind2 = d2
                minid = p2.hash
        line=orbit2str(p1)
        fout.write("%f %15.12f %d %s %d %15.12f \n" % (sim.t,E, p1.hash,line,minid,mind2))
    fout.close()
    return

def saveorbit3 (outfile, sim):
    #save the orbits elements to a file
    fout=open(outfile,mode="a")
    E=sim.calculate_energy()
                
    for p1 in sim.particles:
        if p1.hash==sim.particles[0].hash:
            continue
        mind2 = 1.e6
        minid = 0
        for p2 in sim.particles:
            if p1.hash == p2.hash or p2.hash==sim.particles[0].hash:
                continue
            d2 = (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z)
            mind2 = d2
            minid = p2.hash
        line=orbit2str(p1)
        fout.write("%f %15.12f %d %s %d %15.12f \n" % (sim.t,E, p1.hash,line,minid,mind2))
    fout.close()
    return

def integrate(sim,times,E0,outfile,infofile,small_WJ_radius=False,fast_acc=False):
    #the main integration routine 
    finalstatus=np.zeros(N_pl)
    nstep=np.zeros(N_pl)
    end=np.zeros([N_pl,8])
    
    if frombin or fromtxt:#if reloading a sim, tries to get info from the .pkl file
        if os.path.exists(infofile):
            datadump=pickle.load(open(infofile,"r"))
            init,end,nstep,finalstatus,npcount2,necount2,dE2,bad_dt2,time2 = datadump
        else:
            data = np.loadtxt(outfile,unpack=True)
            n=int(max(data[2]))
            finalstatus=np.zeros(n)
            nstep=np.zeros(n)
            end=np.zeros([n,8])
       
    bad_dts=np.zeros(len(times))
    Ncurrent = sim.N
    dEs=np.zeros(len(times))
    for j,time in enumerate(times):
        try:
            if addmass: 
                add_mass(sim,time,small_WJ_radius=small_WJ_radius,fast_acc=fast_acc)

            sim.integrate(time)
            #deal with Escape
        except rebound.Escape as error:
            #print error
            max_d2 = 0.
            peject=None
            #check distance to be >1000, or (e>1 and distance>100)
            for p in sim.particles:
                if p.hash==sim.particles[0].hash:
                    continue
                d2 = p.x*p.x + p.y*p.y + p.z*p.z
                if d2>max_d2:
                    max_d2 = d2
                    mid = p.hash
                    peject=p
       
            if not peject is None:
                if max_d2>1000:
                    print mid
                    end[mid-1,:]=np.array(list(orbit2str(peject).split()),dtype='f8')
                    sim.remove(hash=mid)
                    nstep[mid-1]=int(sim.t/sim.dt)
                    Ncurrent-=1
                    finalstatus[mid-1]=statuscode['eject']
                    #print "final status",mid,"eject"
                elif max_d2>100:
                    orbit=particle.orbit
                    if orbit.e>1:
                        print mid
                        end[mid-1,:]=np.array(list(orbit2str(peject).split()),dtype='f8')
                        sim.remove(hash=mid)
                        nstep[mid-1]=int(sim.t/sim.dt)
                        Ncurrent-=1
                        finalstatus[mid-1]=statuscode['eject']

        #deal with collision
        if Ncurrent>sim.N:
            #print "collision"
            for l in xrange(len(finalstatus)):

                if finalstatus[l]==0:
                    cflag=True
                    for p in sim.particles:
                        #print p,i+1,cflag
                        if p.hash==sim.particles[0].hash:
                            continue
                        if (p.hash)==(l+1):
                            cflag=False
                            break
                    #print i,cflag
                    if cflag:
                        finalstatus[l]=statuscode['collision']
                        nstep[l]=int(sim.t/sim.dt)
                        #print "final status",i+1,'collision'
            Ncurrent=sim.N
            if check_merger_time:
                saveorbit3(outfile,sim)
                break
        #print orbit2str(sim.particles[1].orbit)
        for p in sim.particles:
            print p.hash
            if  p.hash==sim.particles[0].hash:
                continue

            end[p.hash-1,:]=np.array(list(orbit2str(p).split()),dtype='f8')
           
        if not outfile is None:
            if check_distances:
                saveorbit2(outfile,sim)
            elif check_merger_time:
                saveorbit3(outfile,sim)
            else:
                saveorbit(outfile,sim)#end)#sim)
        if checkpoint:
            checkpointfile=os.path.splitext(outfile)[0]+'.bin'
            sim.save(checkpointfile)
        if  (time/(2*np.pi))%5e4< 2.2e2: #update the .pkl file once every 50,000 years
            dE = np.abs((dEs - E0)/E0)
            #time = timing.time() - start_t
            #total number of planets left
            npcount=len(sim.particles)-1
            #total number of earth type planet left
            necount=0
            init=[]
            for p in sim.particles:
                if p.hash==sim.particles[0].hash:
                    continue
                if p.m<0.5e-3:
                    necount+=1
                nstep[p.hash-1]=int(sim.t/sim.dt)
                parr=np.array(list(orbit2str(p).split()),dtype='f8')
                init.append(parr)
        
            datadump=[init,end,nstep,finalstatus,npcount,necount,dE,bad_dts,time]
            pickle.dump(datadump,open(infofile,"w"))
        #bad_dts[j] = check_for_bad_dt(sim)
    	dEs[j] = sim.calculate_energy()
    return [finalstatus,end,nstep,bad_dts,dEs]

def one_run(runnumber,infile="",HSR=None,small_WJ_radius=False,fast_acc=False,dt=None, EPS=None):
    #this function controls the flow for one run

    #set up the output files and directories (need further modify)

    #initialize the run
    t=0
    outfile=rundir+"rebound%.4d.txt" % runnumber
    if frombin or fromtxt:
        fout = open(outfile,"a")
    else:
        fout=open(outfile,"w")
    fout.close()
    infofile=rundir+"runrebound%.4d.pkl" % runnumber
    if not frombin and not fromtxt:
        if infile=="":
            mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=init_orbit(runnumber,small_WJ_radius=small_WJ_radius)
        else:
            t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl=read_init(infile)

        sim = callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=t)
    elif fromtxt:
        if st.txtfile == "":
            txtfile=rundir+"rebound%.4d.txt" % runnumber
        else:
            txtfile =st.txtfile
        t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,hash_pl=read_init2(txtfile)
        sim = callrebound(mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,t=t,hashes=hash_pl)
    else:
        if st.binfile=="":
            binfile=rundir+"rebound%.4d.bin" % runnumber
        else:
            binfile =st.binfile
        sim=rebound.Simulation.from_file(binfile)
        t=sim.t
    #return
    saveorbit(outfile,sim)#save the initial orbits to output file file
    init=[]
    for p in sim.particles:
        if p.hash==sim.particles[0].hash:
            continue
        parr=np.array(list(orbit2str(p).split()),dtype='f8')
        init.append(parr)
    if debug:
        print init
        print HSR
        return


    # set up integrator (TO BE EDITED)

    if integrator=="hybarid":
    	sim.integrator="hybarid"
    	if not HSR is None:
    	    sim.ri_hybarid.switch_radius = HSR  #units of Hill radii
    	else:
    	    sim.ri_hybarid.switch_radius = 8  #units of Hill radii
    	sim.ri_hybarid.CE_radius = 20.  #X*radius

    	#set up time step
    	if not dt is None:
    	    sim.dt = dt #time step in units of yr/2pi
    	else:
    	    sim.dt = 0.001 #time step in units of yr/2pi
    
    elif integrator == "ias15":
        sim.integrator="ias15"
        if not EPS is None:
            sim.ri_ias15.epsilon = EPS
        else:
            sim.ri_ias15.epsilon = 1e-6

    sim.t=t
    #set up collision options
    #by default, the end result of the collision always
    #keep the small id number.
    sim.testparticle_type = 1
    sim.collision="direct"
    sim.collision_resolve = "merge"
    #sim.collisions_track_dE = 1
    track_energy_offset = 1

    #set up escape options
    sim.exit_max_distance = 100.
    #sim.exit_min_distance = 0.01
    #print sim.collisions[0]

    ####create linear time step arr#
    N_acc = (t_max-t)/acc_time
    if fast_acc:
        N_acc*=10.
    times=np.linspace(t+0,t_max*2*np.pi,int(N_acc))
    #################################
    #times = np.logspace(np.log10(t+1000),np.log10(t+t_max),Noutputs)
    E0 = sim.calculate_energy()
    start_t = timing.time()
    #call integration
    finalstatus,end,nstep,bad_dt,dEs=integrate(sim,times,E0,outfile,infofile,small_WJ_radius=small_WJ_radius,fast_acc=fast_acc)

    if checkpoint:
        checkpointfile=rundir+"rebound%.4d.bin" % runnumber
        sim.save(checkpointfile)


    #Final processing
    #bad_dt = sim.ri_hybarid.timestep_too_large_warning
    dE = np.abs((dEs - E0)/E0)
    time = timing.time() - start_t
    #total number of planets left
    npcount=len(sim.particles)-1
    #total number of earth type planet left
    necount=0
    for p in sim.particles:
        if p.hash==sim.particles[0].hash:
            continue
        if p.m<0.5e-3:
            necount+=1
        nstep[p.hash-1]=int(sim.t/sim.dt)


    datadump=[init,end,nstep,finalstatus,npcount,necount,dE,bad_dt,time]

    def write_outcome(infofile,datadump):
        pickle.dump(datadump,open(infofile,"w"))
        return
    write_outcome(infofile, datadump)
    print 'exit run', runnumber
    return

def main(start=1):
    if gridsearch:
        #minpowdt,maxpowdt,numdt = [-4,-1,1] #min/max limits in logspace, i.e. 10**min - 10**max.
        #minpowHSR,maxpowHSR,numHSR = [2,5,4]      #min/max limits in logspace, i.e. 10**min - 10**
        #dt = P_inner/10./365.*2.*np.pi 
        #HSRarr = np.linspace(minpowHSR,maxpowHSR,numHSR)
        repeat=int(nnodes/4)
        ACCarr = [True,True,False,False]
        WJRarr = [True,False,True,False]
        #print nnodes,repeat
        one_run(start,small_WJ_radius = WJRarr[int((start-1)/repeat)],fast_acc = ACCarr[int((start-1)/repeat)])
    else:
        one_run(start)
    return
if __name__=='__main__':
    if len(sys.argv)==1:
        
        main()
    elif sys.argv[1]=='submit':
        abspath=basename
        submit(abspath)
    elif sys.argv[1]=='restart':
        one_run(1,restartinput)
    else:
        start=eval(sys.argv[1])
        main(start)
