__author__ = 'Dustin Davis'
#AST381 GravDyn HW02
#March 6, 2018


import numpy as np
from scipy.integrate import quad
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
    return rho_crit * delta_c_nfw(log_mass)/((r/R_s)*(1+r/R_s)**2)


def alpha(rho,r):
    return -1* np.log(rho)/np.log(r)

#local velocity dispersion, use Jeans
#todo: **** no ... should be the mass inside radius r
def sigma_rms(log_mass, r, rho):
    v2 = G*M_sun*10**log_mass/r
    s2 = v2/alpha(rho,r)
    return np.sqrt(s2)


#scattering rate
def gamma(log_mass,r):
    """
    """
    rho = nfw_density(log_mass,r)
    return rho*Sigma_T_m*sigma_rms(log_mass,r,rho)



def interior_mass_integrand(r,log_mass):
    return 4.*np.pi*(nfw_density(log_mass,r))*r**2

def interior_mass(log_mass,r):
    return quad(interior_mass_integrand, 0, r,args=(log_mass))


def main():
    log_mass = np.arange(10.,16.,1.) #M_vir , virial mass
   # R_s = scale_radius(log_mass)

   # print(R_s/(3.086e18))

    r_grid = np.logspace(-1, 2,num=100)
    m_grid10 = interior_mass(log_mass[0],r_grid[0])[0]

    norm = plt.Normalize()
    color = plt.cm.jet(norm(np.arange(len(log_mass))))

    plt.figure()
    plt.title("DM Scattering Rate vs Radius\n" + "Concordance Cosmology, h = %g" %(h))
    plt.xlabel("$Log_{10}(r/R_s)$")
    plt.ylabel(r"$\Gamma(r)$ [$s^{-1}$]")
    plt.gca().set_yscale("log")


    plt.axhline(y=H0,linestyle="--",color='r')
    plt.gca().annotate(r"$H_{0}$", xy=(0.8*10**10,H0),
                       xytext=(0.6*10**10,H0/100))

    plt.axhline(y=1/H_tdyn, linestyle="--")
    plt.gca().annotate(r"$10H_{0}$", xy=(0.8*10**10, 1/H_tdyn),
                       xytext=(0.6*10**10,1/H_tdyn *10 ))

    for i in range(len(log_mass)):
        plt.plot(r_grid,gamma(log_mass[i],r_grid),color=color[i],label="$10^{%d}M_{\odot}$"%log_mass[i])

    plt.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

    #plt.savefig("dd_hw2.png")
    plt.show()


if __name__ == '__main__':
    main()

