#!/usr/bin/env python
import numpy as np
import scipy as sp
from  setting import *

def mass_growth(time, initial_mass,doubling_time):
    new_mass = initial_mass*(2.0)**(time/doubling_time)
    return new_mass

def add_mass(sim,time, small_WJ_radius=True, fast_acc=False):
    scale = 1.0
    if fast_acc:#adjust accretion timing for faster accretion
        scale = 1.0/10.0
    if time<= 4.5*doubling_time*2*np.pi*scale:
    #if time <= 5e4*2*np.pi:


        sp = sim.particles
        mass_list = []
        M_earth = 1./300./1000.
        threshold = 12*M_earth

        for i in range(1,sim.N):
            mass_list.append(sp[i].m)
        accreting_planet = np.where(np.array(mass_list)>=threshold)[0]+1
        for i in accreting_planet:
            planet_density = (3*sp[i].m/(4*np.pi*sp[i].r**3))
            sp[i].m = mass_growth(acc_time,sp[i].m,doubling_time)
            if small_WJ_radius:#settings for different WJ radii
                sp[i].r = sp[i].a*(1-sp[i].e)*(sp[i].m/(3*sp[0].m))**(1./3.)/10.#(3*sp[i].m/(4*np.pi*sp[i].r))**(1./3.)
            else:
                sp[i].r = sp[i].a*(1-sp[i].e)*(sp[i].m/(3*sp[0].m))**(1./3.)/10.
        sim.move_to_com()
    return

