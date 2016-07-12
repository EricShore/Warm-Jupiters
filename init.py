#!/usr/bin/env python
import numpy as np
import scipy as sp
from setting import *
from scipy.stats import rayleigh

def init_orbit(runnumber, small_WJ_radius = True):
    #initial the basic param array for a system
    
    mass_pl=np.zeros(N_pl)
    a_pl=np.zeros(N_pl)
    r_pl=np.zeros(N_pl) #changed to radius rather than density

    #a_inner,a_pl[:3]=set_hill(mass_pl[:3])

    #initial semi-major axis and masses of super earths,
    #in solar units
    M_earth=1./300./1000.
    R_earth = 6371./1.49e8
    density = 3*M_earth/(4*np.pi*R_earth**3)
    loop = True

    seeds=(np.random.rand(500)*1000).astype(int) # this serie is likely to be the same for your paralleal runs (at least for those in the same node, but will be different between different trials [today's run and next week's run])
    
    newseed = seeds[runnumber] # this operation ensure the seed between the runs are different. 
    np.random.seed(newseed) #use the new seed

#draw orbital element below without randomstat argument

    while (loop==True):
    #set up other orbital elements
        e_pl=rayleigh.rvs(scale=sigma_e,size=N_pl)
        i_pl=rayleigh.rvs(scale=sigma_i,size=N_pl)
        
        omega_pl=2.*np.pi*np.random.rand(N_pl)
        Omega_pl=2.*np.pi*np.random.rand(N_pl)
        M_pl=2.*np.pi*np.random.rand(N_pl)
        mass_pl = np.random.normal(loc=5*M_earth,scale = 2*M_earth,size=N_pl)
        if any(mass_pl<0.5*M_earth) or any(mass_pl>10*M_earth): #ensures mass is between 0.5 and 10 Earth Masses
            continue
        r_pl =(3*mass_pl/(4*np.pi*density))**(1./3.)
        K = np.random.normal(loc=k_Hill,scale=2)
        if K <=0: #ensures seperation between planets is positive
            continue
        P_min =rayleigh.rvs(scale=P_inner/(365.))
        a_min = P_min**(2./3.)

        in_range=[]
        a_pl[0] = a_min
        
        for i in range(1,N_pl): #set up ladder of planets 
            m_ratio = 2*(3/(mass_pl[i-1]+mass_pl[i]))**(1./3.)
            a_pl[i] =(K/m_ratio + 1)/(1 - K/m_ratio)*a_pl[i-1]
            if a_pl[i]>0.1 and a_pl[i]<0.4:
                in_range.append(i)
        if len(in_range)>0: #ensures there is always at least 1 planet in the WJ range
            accreting_planet = np.random.choice(in_range)
            if accreting_planet > N_pl-2: #ensures each WJ has at least 2 outer neighbors
                continue
            mass_pl[accreting_planet]=12*M_earth
            loop=False
            if small_WJ_radius: #settings for different Warm Jupiter radii
                r_pl[accreting_planet]=a_pl[accreting_planet]*(1-e_pl[accreting_planet])*(mass_pl[accreting_planet]/3)**(1./3.)/10.
            else:
                r_pl[accreting_planet]=a_pl[accreting_planet]*(1-e_pl[accreting_planet])*(mass_pl[accreting_planet]/3)**(1./3.)
            for i in range(accreting_planet,N_pl): #adjust ladder of planets for WJ
                if i>0:
                    m_ratio = 2*(3/(mass_pl[i-1]+mass_pl[i]))**(1./3.)
                    a_pl[i] =(K/m_ratio + 1)/(1 - K/m_ratio)*a_pl[i-1]
    #print mass_pl, a_pl, r_pl, e_pl, i_pl, omega_pl, Omega_pl, M_pl
    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]
def init_orbit2(randomstat=1):# to be used for manually restarting from a text file
   mass_pl = np.array([0.000002546206, 0.000018003167, 0.000040000000, 0.000013145077, 0.000019332091])
   r_pl = np.array([0.000039086455,0.000075020173, 0.002545837352, 0.000067553750, 0.000076822429])
   a_pl = np.array([0.076640105897, 0.090364905459,0.114113732060,0.143132603639, 0.173436994996])
   e_pl = np.array([ 0.012011911466, 0.017225171277, 0.020216409881, 0.001271180551,0.008610238533])
   i_pl = np.array([0.007262291848, 0.015565181510, 0.008963412498,0.004366647019,0.015973512383])
   omega_pl = np.array([-1.861388419611, -2.018638273672,1.083502338206, 2.377066748780,-0.279289559350])
   Omega_pl = np.array([ 1.166041048373, 1.895378506796, 2.122070122255, -0.062594865598, 2.488691715614])
   M_pl  =np.array([-0.707406825816, 1.772171293919, 2.737586757107, 0.481635488316, -2.806408520830])      
   return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]

def init_orbit3(randomstat=1):
    M_earth=1./300./1000.
    R_earth = 6371./1.49e8
    density = 3*M_earth/(4*np.pi*R_earth**3)
    mass_pl = np.ones(N_pl)*M_earth*5
    r_pl = (3*mass_pl/(4*np.pi*density))**(1./3.)
    
    a_pl = np.zeros(N_pl)
    e_pl = np.zeros(N_pl)
    i_pl=rayleigh.rvs(scale=sigma_i,size=N_pl,random_state=randomstat+1000)
    omega_pl=2.*np.pi*np.random.rand(N_pl)
    Omega_pl=2.*np.pi*np.random.rand(N_pl)
    M_pl=2.*np.pi*np.random.rand(N_pl)
    a_pl[0] = 0.15
    e_pl = rayleigh.rvs(scale=0.15,size=N_pl,random_state=randomstat+1000)
    peri = a_pl[0]*(1+e_pl[0])*(1 - np.random.rand()/10.)
    a_pl[1] = peri/(1.-e_pl[1])

    accreting_planet = np.random.choice([0,1])
    mass_pl[accreting_planet]=12*M_earth
    r_pl[accreting_planet]=a_pl[accreting_planet]*(1-e_pl[accreting_planet])*(mass_pl[accreting_planet]/3)**(1./3.)
    
    return [mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]
    


         

def read_init(infile):
    #need to reload orbit elements from end result of a file.
    data=np.loadtxt(infile, skiprows=1)
    t=data[:,0][0]
    mass_pl=data[:,3]
    r_pl=data[:,4]
    a_pl=data[:,5]
    e_pl=data[:,6]
    i_pl=data[:,7]
    Omega_pl=data[:,8]
    omega_pl=data[:,9]
    M_pl=data[:,10]
    return [t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl]

def read_init2(infile):# used for reading from a text file
    data = np.loadtxt(infile, skiprows=1, unpack = True)
    data2 = data[:,:min(np.where(data[0]== max(data[0]))[0])]
    t = data2[0][-1]
    N_pl = np.size(np.unique(data[2][-10:]))
    
    mass_pl = data2[3][-1*N_pl:]
    r_pl = data2[4][-1*N_pl:]    
    a_pl = data2[5][-1*N_pl:]
    e_pl = data2[6][-1*N_pl:]
    i_pl = data2[7][-1*N_pl:]
    Omega_pl = data2[8][-1*N_pl:]
    omega_pl = data2[9][-1*N_pl:]
    M_pl = data2[10][-1*N_pl:]
    hash_pl = data2[2][-1*N_pl:]
    return [t,mass_pl,a_pl,r_pl,e_pl,i_pl,omega_pl,Omega_pl,M_pl,hash_pl]
