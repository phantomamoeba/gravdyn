__author__ = 'Dustin Davis'
#AST381 GravDyn HW02
#March 6, 2018


import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


NFW_rho_0 = 4.0
Sigma_T_m  = 1.0  #cross section per unit mass (does not really matter ... just sets a scale)



def nfw_density(r, R_s, rho_s):
    """

    :param r: radius
    :param R_s:  scale radius (where log slope of density is -2) (aka r_2)
    :param rho_s:  scale density (at the scale radius)
    :return: dimensionless density (density / scale density)
    """
    return NFW_rho_0/((r/R_s)*(1.+r/R_s)**2)/rho_s

#local velocity dispersion
def sigma_rms(r):
    return 1.


#scattering rate
def gamma(r,R_s,rho_s):
    """

    :param r:
    :param R_s:
    :param rho_s:
    :return: gamma(r) / Sigma_T_m
    """

    return nfw_density(r,R_s,rho_s)*sigma_rms(r) #*Sigma_T_m ... but will just call it gamma/Sigma_T_m







def main():
    log_mass = np.arange(10.,16.,1.)

    rho_s = np.zeros(log_mass.shape) + 1. #just for now
    R_s = np.zeros(log_mass.shape) + 1. #just for now

    r_grid = np.logspace(-4, 10,num=100)
    norm = plt.Normalize()
    color = plt.cm.hsv(norm(np.arange(len(log_mass))))

    plt.figure()
    plt.title("DM Scattering Rate vs Radius")
    plt.xlabel("$Log_{10}(r/R_s)$")
    plt.ylabel(r"$\Gamma(r)$ / ($\sigma_{T}$/m)")


    for i in range(len(log_mass)):
        plt.plot(r_grid,gamma(r_grid,R_s[i],rho_s[i]),color=color[i])



    #plt.savefig("dd_hw2.png")
    plt.show()


if __name__ == '__main__':
    main()

