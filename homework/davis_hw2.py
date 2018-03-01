__author__ = 'Dustin Davis'
#AST381 GravDyn HW02
#March 6, 2018


import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


NFW_rho_0 = 4.0 #units of rho_2
Sigma_T = 6.6524*10**(-25) #cm^2
GeV_to_grams = 1.783*10**-24
mass_DM = 1.0 #GeV
m_p =  mass_DM * GeV_to_grams
h = 0.7
M_sun = 1.989*10**33 #grams


def nfw_density(r, R_s, rho_s):
    """

    :param r: radius
    :param R_s:  scale radius (where log slope of density is -2) (aka r_2)
    :param rho_s:  scale density (at the scale radius) (aka rho_2)
    :return: denisty in g/cm^3
    """

    return NFW_rho_0 / ((r / R_s) * (1. + r / R_s) ** 2) * rho_s #units of rho_2

#local velocity dispersion
def sigma_rms(r):
    return 1.


#scattering rate
def gamma(r,R_s,rho_s):
    """

    :param r:
    :param R_s:
    :param rho_s:
    :return: gamma(r)
    """

    return nfw_density(r,R_s,rho_s)*Sigma_T/m_p*sigma_rms(r)

def log_conc_parm(log_mass,a=0.95,b=-0.101):
    """
    assume delta = 200, redshift = 0 for the a, b defaults
    :param mass:
    :param a:
    :param b:
    :return:
    """
    return a + b*np.log((10**log_mass)/((10**12)/h))

#R_s
def scale_radius(log_c, r_vir):
    return r_vir/(10**log_c)





def main():
    log_mass = np.arange(10.,16.,1.)
    r_vir = np.zeros(log_mass.shape) + 1.
    rho_s = np.zeros(log_mass.shape) + 1. #just for now
    R_s = scale_radius(log_conc_parm(log_mass),r_vir)#   np.zeros(log_mass.shape) + 1. #just for now

    r_grid = np.logspace(-4, 10,num=100)
    norm = plt.Normalize()
    color = plt.cm.hsv(norm(np.arange(len(log_mass))))

    plt.figure()
    plt.title("DM Scattering Rate vs Radius\n" + "h = %g, $\m_p$=%g" %(h,mass_DM))
    plt.xlabel("$Log_{10}(r/R_s)$")
    plt.ylabel(r"$\Gamma(r)$ / ($\sigma_{T}$/m)")


    for i in range(len(log_mass)):
        plt.plot(r_grid,gamma(r_grid,R_s[i],rho_s[i]),color=color[i])



    #plt.savefig("dd_hw2.png")
    plt.show()


if __name__ == '__main__':
    main()

