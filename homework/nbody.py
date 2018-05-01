__author__ = 'Dustin Davis'
#AST381 GravDyn HW03 N-Body
#May 11, 2018

#reminder to self: convention here is theta is polar, phi is azimuth
#everything in radians and cgs

import numpy as np
from scipy.integrate import quad
from scipy.integrate import nquad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pc2cm = 3.08567758e18
M_sun2grams = 1.989e33
a_kpc=35
a_cm=a_kpc*1e3*pc2cm #scale length aka r_s
M_tot_M_sun = 1e12
M_tot_grams = M_tot_M_sun * M_sun2grams
TotalNumParticles = int(1e6) #all particles of the same mass (set low for testing)
ParticleMass = M_tot_grams / TotalNumParticles
G = 6.67428*10**(-8) #cgs
rho_0 = 1.0 #should be something like rho_200? or rhocrit?

v_g = np.sqrt(G*M_tot_grams/a_cm)


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

def dM(radius): #basically my pdf
    return 4*np.pi*radius**2 * density_at(radius)

def mass_enclosed(radius): #basically my cdf
    #return quad(dM,0,radius)[0]
    #since can analytically integrate, just use that instead
    return M_tot_grams*(radius/(radius + a_cm))**2

#todo: the translation is okay, but uniform random is wrong ...
#todo: you end up with equal numbers of particles at all radii (because of the uniform random)
def random_radius_cm():
    #make the upper limit a bit less than 1 since that would correspond to r=infinity
    #max_x = mass_enclosed(100*a_cm)/M_tot_grams
    max_x = 1.0
    x = np.random.uniform(0,max_x)
    #solve eqn(3) from Hernquist paper for r in terms of fractional mass
    #where fractional mass = M(r)/M_tot = x which runs (0,1)
    #the return is the radius at which that fractional mass is the interior mass
    #and is then the equivalent of a random radius drawn from the Hernquist distribution
    return a_cm*(np.sqrt(x))/(1-np.sqrt(x))


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


def pdf_sample_accept_reject(pdf, samples=1, x_min=0, x_max=1, pdf_probe=int(1e3)):

    #generate uniform random position in x,y
    #if y is "under" the pdf, accept the position, if not, reject it

    #Hernquist is monatonic, but for future use, the pdf might not be
    #get the boundaries of the pdf (between pdf(x_min) and pdf(x_max)
    pdf_min = 0.
    pdf_max = pdf(np.linspace(x_min, x_max, pdf_probe)).max()

    rvs = np.zeros(samples)
    accepts = 0
    rejects = 0 #not really using this right now, but for future might want ratio with accepts?

    while accepts < samples:
        x = np.random.uniform(x_min, x_max)
        y = np.random.uniform(pdf_min, pdf_max)

        #if y is under the pdf curve ...
        if y < pdf(x):
            rvs[accepts] = x
            accepts += 1
            #print("Samples = %d" %accepts)
        else:
            rejects += 1

    return rvs



#def q(E):
#    return np.sqrt(-1/v_g**2 * E)

def q(radius,v):
    return np.sqrt(-1/v_g**2 * E_tot(radius,v))

def E_tot(radius,v): #per unit mass NOT BINDING E
    return potential(radius) + 0.5*v**2

def check_q(radius,v):
    E =  E_tot(radius,v)
    if E > 0:
        return False

    q = np.sqrt(-1 / v_g ** 2 * E)
    if (q < 0) or (q > 1):
        return False
    else:
        return True

def f(radius,v): #f(E) as f(v)
    #still need to calc v .... we want v back
    #can we solve for v?
    q = np.sqrt(-1/v_g**2 * E_tot(radius,v))

    #if (q < 0) or (q > 1):
    #    return None

    return M_tot_grams/(8*np.sqrt(2)*(np.pi*a_cm*v_g)**3)*1/(1-q**2)**(2.5) * \
    (3 * np.arcsin(q) + q*np.sqrt((1-q**2))*(1-2*q**2)*(8*q**4-8*q**2-3))


def random_velocity(radius):
    v = np.random.uniform(0,5e7) #0 to 500 km/s in cm/2
    while not check_q(radius, v):
        v = np.random.uniform(0,5e7) #0 to 500 km/s in cm/2
    return v

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
        p = Particle() #auto-pops a random position
        v = random_velocity(p.r)
        #so v is okay now
        #get a random direction for v
        theta,phi = random_on_sphere()
        x,y,z = tpr2xyz(theta,phi) #unit sphere
        #vector components
        p.v_x = x*v
        p.v_y = y*v
        p.v_z = z*v

        particles.append(p)

       # particles.append(Particle())
        #print(particles[i].spherical())
        #print(particles[i])


    #sanity check velocity distribution is spherical
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #since particles is random already, just plot the first 1000
    x = [p.v_x/1e5 for p in particles[0:1000]]
    y = [p.v_y/1e5for p in particles[0:1000]]
    z = [p.v_z/1e5 for p in particles[0:1000]]

    ax.scatter(x,y,z)
    plt.show()

    # #sanity check the distributions
    # r_grid = np.linspace(1,1e6*pc2cm,10000)
    # plt.plot(r_grid / pc2cm, mass_enclosed(r_grid) / M_tot_grams) #cdf
    # plt.show()
    # plt.plot(r_grid/pc2cm,dM(r_grid)) #pdf
    # plt.show()
    #
    # re = [p.r/a_cm for p in particles]
    # #limit the range to 100a
    # plt.title("Inverse method")
    # plt.hist(np.array(re),bins=1000,cumulative=True,density=True,range=[0,100])
    # plt.show()
    #
    # rvs = pdf_sample_accept_reject(dM, samples=TotalNumParticles, x_min=1, x_max=1000 * a_cm)
    # plt.title("Accept-Reject method") #check this vs inverse method
    # plt.hist(np.array(rvs)/a_cm, bins=1000, cumulative=True, density=True, range=[0, 100])
    # plt.show()




if __name__ == '__main__':
    main()



