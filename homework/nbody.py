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
TotalNumParticles = 1e6


def rng(seed=None): #simple uniform (0,1)
    if seed is not None:
        np.random.random(seed)
    return np.random.random()

def random_on_sphere(seed=None):
    theta = np.arccos(2*rng(seed)-1)
    phi = 2.0*np.pi*rng(seed)
    return theta, phi

def rtp2xyz(theta,phi,radius=1):
    x = radius * np.cos(phi) * np.sin(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    return x,y,z

def main():

    pass

if __name__ == '__main__':
    main()



