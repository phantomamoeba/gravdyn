__author__ = 'Dustin Davis'
#AST381 GravDyn HW03 N-Body
#May 11, 2018

#reminder to self: convention here is theta is polar, phi is azimuth
#everything in radians and cgs

import numpy as np
from scipy.integrate import quad
from scipy.integrate import nquad
import matplotlib.pyplot as plt

pc2cm = 3.08567758e18
M_sun2grams = 1.989e33
a_kpc=35
a_cm=a_kpc*1e3*pc2cm #scale length aka r_s
M_tot_M_sun = 1e12
M_tot_grams = M_tot_M_sun * M_sun2grams
TotalNumParticles = int(1e6) #all particles of the same mass
ParticleMass = M_tot_grams / TotalNumParticles
G = 6.67428*10**(-8) #cgs
rho_0 = 1.0 #should be something like rho_200? or rhocrit?

#fixed seed so repeatable (for testing)
np.random.seed(1138)


def rng(seed=None): #simple uniform (0,1)
    return np.random.random()

def random_radius():#todo using Hernquist profile
    pass

def random_on_sphere():
    theta = np.arccos(2*rng()-1)
    phi = 2.0*np.pi*rng()
    return theta, phi

def tpr2xyz(theta,phi,radius=1.0):
    x = radius * np.cos(phi) * np.sin(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    return x,y,z



#Hernquist
def potential(radius):
    return -G*mass_enclosed(radius)/(radius+a_cm)

def density_at(radius):
    return M_tot_grams/(2*np.pi)*a_cm/radius *(1/(radius+a_cm)**3)

def mass_at(radius):
    return 4*np.pi*(radius**2)*density_at(radius)

def dM(radius):
    return 4*np.pi*radius**2 * density_at(radius)

def mass_enclosed(radius):
    #return quad(dM,0,radius)[0]
    #since can analytically integrate, just use that instead
    return M_tot_grams*(radius/(radius + a_cm))**2

def random_radius_cm():
    x = np.random.random()
    #solve eqn(3) from Hernquist paper for r in terms of fractional mass
    #where fractional mass = M(r)/M_tot = x which runs (0,1)
    #the return is the radius at which that fractional mass is the interior mass
    #and is then the equivalent of a random radius drawn from the Hernquist distribution
    return -a_cm*(x**2)/(x**2-1)


class Particle():
    def __init__(self, random=True):
        self.x = None
        self.y = None
        self.z = None
        self.r = None
        self.theta = None
        self.phi = None
        self.v_x = -1
        self.v_y = -1
        self.v_z = -1
        self.mass_enclosed = -1

        if random:
            self.initialize()

    def initialize(self):
        self.r = random_radius_cm()
        self.theta,self.phi = random_on_sphere()
        self.x, self.y, self.z = tpr2xyz(self.theta,self.phi,self.r)
        self.mass_enclosed = mass_enclosed(self.r)

    def __str__(self):
        return "%f %f %f %f %f %f" %(self.x, self.y, self.z, self.v_x, self.v_y, self.v_z)

    def spherical(self):
        return "%f %f %f" %(self.r/pc2cm/1e3, 90.0 - self.theta*180.0/np.pi, self.phi*180.0/np.pi)

def main():


    #for each particle
    #get its radius, then angular position
    #then translate to cartesian
    #get its velocity
    #write out to file (x,y,z,v_x,v_y,v_z)

    #get all particle coords
    print("Building particle locations (~15-20seconds) ...")
    particles = []
    for i in range(int(TotalNumParticles)):
        particles.append(Particle())
       # particles.append(Particle())
        #print(particles[i].spherical())
        #print(particles[i])

    # r_grid = np.linspace(1,1e6*pc2cm,10000)
    # plt.plot(r_grid / pc2cm, mass_enclosed(r_grid) / M_tot_grams)
    # plt.show()
    #
    # me = [p.mass_enclosed/M_sun2grams for p in particles]
    # plt.hist(np.array(me),bins=1000,cumulative=True,normed=True)
    # plt.show()


if __name__ == '__main__':
    main()



