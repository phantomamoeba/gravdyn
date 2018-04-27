__author__ = 'Dustin Davis'
#AST381 GravDyn HW03 N-Body
#May 11, 2018

#reminder to self: convention here is theta is polar, phi is azimuth
#everything in radians

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

pc2cm = 3.08567758e18
M_sun2grams = 1.989e33
a_kpc=35
a_cm=a_kpc*1e3*pc2cm
M_tot_M_sun = 1e12
M_tot_grams = M_tot_M_sun * M_sun2grams
TotalNumParticles = 1e6 #all particles of the same mass


#fixed seed so repeatable (for testing)
np.random.random(1138)

def rng(seed=None): #simple uniform (0,1)
    return np.random.random()

def random_radius():#todo using Hernquist profile
    pass

def random_on_sphere():
    theta = np.arccos(2*rng()-1)
    phi = 2.0*np.pi*rng()
    return theta, phi

def rtp2xyz(theta,phi,radius=1):
    x = radius * np.cos(phi) * np.sin(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    return x,y,z





def main():

    #for each particle
    #get its radius, then angular position
    #then translate to cartesian
    #get its velocity
    #write out to file (x,y,z,v_x,v_y,v_z)

    pass

if __name__ == '__main__':
    main()



