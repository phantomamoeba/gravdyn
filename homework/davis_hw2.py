__author__ = 'Dustin Davis'
#AST381 GravDyn HW02
#March 6, 2018


import numpy as np
from scipy.integrate import quad
from scipy.integrate import nquad
import matplotlib.pyplot as plt


NFW_rho_0 = 4.0 #units of rho_2
#GeV_to_grams = 1.783*10**-24
#mass_DM = 1.0 #GeV
#m_p =  mass_DM * GeV_to_grams
h = 0.71
H0 = 2.2683*10**(-18) #s^-1
H_tdyn = 0.1/H0 #seconds ... the rate would be H0/0.1 or 10*H0
M_sun = 1.989*10**33 #grams
Sigma_T_m = 1.0 #cm^2/g or could use 1/12  that is, Sigma_T / m_p
delta_vir = 200.
rho_crit = 1.8788*10**(-29)*h**2 #g/cm^3
G = 6.67428*10**(-8)
#USE c_200, r_200, etc not



def delta_c_nfw(log_mass):
    c = 10**(log_concentration(log_mass))
    return 200./3. * (c**3.)/(np.log(1+c) - (c/(1+c)) ) #log = ln

def log_concentration(log_mass,a=0.905,b=-0.101): #eqn 8
    """
    assume delta = 200, redshift = 0 for the a, b defaults
    :param mass:
    :param a:
    :param b:
    :return:
    """
    #return a + b*np.log10( h*(10**log_mass)/ 10.**12.) #want log10 here
    return a + b*(np.log10(h)+log_mass-12)



def scale_radius(log_mass): #based on r_200
    """

    :param log_mass: in M_sun
    :return: in cm
    """
    c = 10**(log_concentration(log_mass))
    r_200 = ((M_sun * 10 ** log_mass) / (4 / 3 * np.pi * rho_crit *200.)) ** (1. / 3.)  # yes log = ln
    return r_200/c

def nfw_density(log_mass,r):
    """
    :param r: radius
    :param R_s:  scale radius (where log slope of density is -2) (aka r_2)
    :return: denisty in g/cm^3
    """

    R_s = scale_radius(log_mass)
    rho = rho_crit * delta_c_nfw(log_mass) / ((r / R_s) * (1 + r / R_s) ** 2)

    #if log_mass == 10:
    #    print("%d  %d  %0.2g" % (log_mass, r / 3.086e18, rho))

    return rho


#def alpha(rho,r):
#    return -1* np.log10(rho)/np.log10(r)

def sigma_rms(log_mass, r): #in cm/s
    sigma = dispersion(log_mass, r)
   # print ("%d  %d  %0.2g" %(log_mass, r/3.086e18, sigma))
    return sigma


#scattering rate
def gamma(log_mass,r):
    """
    """
    rho = nfw_density(log_mass,r)
    return rho*Sigma_T_m*sigma_rms(log_mass,r)


def interior_mass_integrand(r,log_mass): #without the 4pi and G ... moved out front of outer integrand
    #expects r in cm
    return (nfw_density(log_mass,r))*r**2

def interior_mass(log_mass,r):
    #expects r in cm
    return quad(interior_mass_integrand, 0, r, args=(log_mass))[0]

#def dispersion_integrand(r_exterior,log_mass):
#    return nfw_density(log_mass,r_exterior)/(r_exterior**2)*quad(interior_mass_integrand, 0, r, args=(log_mass))[0]


def double_integrand(r1,r2,log_mass): #r1 inner upper limit, r2 outer lower limit
    x = nfw_density(log_mass,r2)/(r2**2.) * nfw_density(log_mass,r1)*r1**2.
    return x

def limits_r2(r2,log_mass): #outside
    r_200 = scale_radius(log_mass) * 10**(log_concentration(log_mass))
 #   return [r2, np.inf]
    return [r2, 100*r_200] #can't use np.inf (won't actually converge, so pick a big number)

def limits_r1(r1,dummy): #inside
    return [0, r1]

def dbl_integral(r,log_mass):
    options = {'limit': 100}
    #limits .. inner first, then outer
    x = nquad(double_integrand, [limits_r1, limits_r2(r,log_mass)],opts=[options,options],args=(log_mass,))[0]

    #if log_mass == 10:
    #    print()
    return x

def dispersion(log_mass,r):#cm/s

    nfw_density_at_r = nfw_density(log_mass,r)
    sigma = np.sqrt(4.*np.pi*G/nfw_density_at_r*dbl_integral(r,log_mass))

    if log_mass == 10:
        print("%d  %d  %0.2g" % (log_mass, r / 3.086e18, sigma))
    return sigma

# print("%d  %d  %0.2g" % (log_mass, r / 3.086e18, rho))

def main():
    log_mass = np.arange(10.,16.,1.) #M_vir , virial mass
    #r_grid = np.logspace(-2, 2,num=50) #1/10 to 100x the scale_radius
    r_grid = np.logspace(18,25,num=20) #cm

    norm = plt.Normalize()
    color = plt.cm.jet(norm(np.arange(len(log_mass))))

    plt.figure()
    plt.title("DM Scattering Rate vs Radius\n" + "Concordance Cosmology, h = %g" %(h))
    plt.xlabel("r [kpc]")
    plt.ylabel(r"$\Gamma(r)$ [$s^{-1}$]")
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")


    plt.axhline(y=H0,linestyle="--",color='r')
    plt.gca().annotate(r"$H_{0}$", xy=(10**-3,H0),
                       xytext=(10**-3,H0/5))

    plt.axhline(y=1/H_tdyn, linestyle="--")
    plt.gca().annotate(r"$10H_{0}$", xy=(10**-3, 1/H_tdyn),
                       xytext=(10**-3,1/H_tdyn*2))

    #this is a dumb way to do this, but expedient to code
    #at the least should re-work to just extend mass shell by shell rather than fully recompute each time
    for i in range(len(log_mass)):
        g = np.zeros(r_grid.shape)
        #radii = r_grid*scale_radius(log_mass[i])
        radii = r_grid
        for j in range(len(r_grid)):
            g[j] = gamma(log_mass[i],radii[j])

        plt.plot(r_grid/(3.086e21),g,color=color[i],label="$10^{%d}M_{\odot}$"%log_mass[i])

    plt.legend(loc='lower left', bbox_to_anchor=(0.02, 0.02), borderaxespad=0)

    plt.savefig("dd_hw2.png")
    plt.show()


if __name__ == '__main__':
    main()

