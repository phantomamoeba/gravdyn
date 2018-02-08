#test

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt



G = 1
M = 1
a = 1



def sigma_integrand(r,a,gamma):
    '''
    :param r: float radius
    :param a: float scaling radius
    :param g: float gamma
    :return:
    '''
    return (r**(1.-2.*gamma))/((r+a)**(7.-2.*gamma))


def rvd(r,a,gamma):
    '''
    Radial Velocity Dispersion aka sigma
    :param r: float radius
    :param a: float scaling radius
    :param g: float gamma
    :return:
    '''

    return np.sqrt(G*M*(r**gamma)*(r+a)**(4.-gamma)*quad(sigma_integrand,r,np.inf,args=(a,gamma))[0])


def plot_sigma(r_i,r_f,num,a,gamma):
    '''
    :param r:
    :param a:
    :param g:
    :return:
    '''

    #grid = np.arange(r_i,r_f,step)
    grid = np.logspace(np.log10(r_i),np.log10(r_f),num=num,base=10.0)
    sigma = []
    #would be more efficient to use matrix operations but this is easier to read and fast enough
    for r in grid:
        sigma.append(rvd(r,a=1,gamma=gamma))

    plt.title("Radial Dispersion (Prob#1)\n$\gamma$=3/2")
    plt.xlabel("$Log_{10}(r/a)$")
    plt.ylabel(r"$\sigma_{r}$ [$\sqrt{GMa}$]")
    #plt.gca().set_xscale("log")
    plt.plot(np.log10(grid),sigma)

    plt.savefig("dd_hw1p1.png")

    plt.show()



def Phi(r,a,gamma):
    if gamma == 2:
        return G*M/a*np.log(r/(r+a))
    else:
        return G*M/a*(-1/(2.-gamma)*(1-(r/(r+a))**(2.-gamma)))

def Psi(r,a,gamma):
    return -1*Phi(r,a,gamma)/(G*M/a)

#def bE(E,a=1): #dimensionless binding energy
    #return -E/(2.*G*M/(3.*a))

def E_c(r,a,gamma): #E_c but is generally true for any orbit? Not just circular (or,at least for this purpose)

    y_c = r/(r+a)
    return 1./(2.-gamma) - (4.-gamma)/(2.*(2.-gamma))*y_c**(2.-gamma) + 0.5*y_c**(3.-gamma)


def y_func(r,a,gamma): #really feel dumb naming this way ...
    if gamma == 2:
        return np.exp(-1*Psi(a))
    else:
        return (1.-(2.-gamma)*Psi(r,a,gamma))**(1./(2.-gamma))


def f_integrand(psi,E,r,gamma,a):
    y = y_func(r,a,gamma)
    return (1.-y)**2. * (gamma + 2.*y + (4.-gamma)*y**2.)/( y**(4.-gamma) * np.sqrt(E-psi))



def df_f(r,a,gamma):

    psi = Psi(r,a,gamma)
    E = E_c(r,a,gamma)
    return (3.-gamma)*M / (2. * (2.*np.pi*G*M*a)**(1.5)) * quad(f_integrand,0,E,args=(E,r,gamma,a))[0]


def plot_f(r_i,r_f,step,a,gamma):
    '''
    :param r:
    :param a:
    :param g:
    :return:
    '''

    grid = np.arange(r_i,r_f,step)
    psd = []
    #would be more efficient to use matrix operations but this is easier to read and fast enough
    for r in grid:
        psd.append(df_f(r,a=1,gamma=gamma))

    plt.title("Phase Space Density (Prob#2)\n$\gamma$=3/2")
    plt.xlabel("Binding E")
    plt.ylabel("log(density of states)")
    #plt.gca().set_yscale("log")
    plt.plot(-1*E_c(grid,a,gamma),np.log10(psd))

    plt.savefig("dd_hw1p2.png")
    plt.show()



def main():

    #will get a warning about convergence at the 0 limit, but is okay
    #plot_sigma(10**(-5),100,1000,1,1.5)

    plot_f(0.00, 100, step=0.01, a=1, gamma=1.5)

if __name__ == '__main__':
    main()



