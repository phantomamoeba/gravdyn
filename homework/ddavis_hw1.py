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


def plot_sigma(r_i,r_f,step,a,gamma):
    '''
    :param r:
    :param a:
    :param g:
    :return:
    '''

    grid = np.arange(r_i,r_f,step)
    sigma = []
    #would be more efficient to use matrix operations but this is easier to read and fast enough
    for r in grid:
        sigma.append(rvd(r,a=1,gamma=gamma))

    plt.title("Radial Dispersion (Prob#1)")
    plt.xlabel("Radius")
    plt.ylabel(r"$\sigma_{r}$ [$\sqrt{\frac{GM}{r}}$]")
    plt.plot(grid,sigma)

    plt.savefig("dd_hw1p1.png")

    plt.show()



def Phi(r,a,gamma):
    if gamma == 2:
        return G*M/a*np.log(r/(r+a))
    else:
        return G*M/a*(-1/(2.-gamma)*(1-(r/(r+a))**(2.-gamma)))

def Psi(a):
    return -1/Phi(G*M/a)

#def bE(E,a=1): #dimensionless binding energy
    #return -E/(2.*G*M/(3.*a))
def E_c(r,a,gamma):

    y_c = r/(r+a)
    return 1./(2.-gamma) - (4.-gamma)/(2.*(2.-gamma))*y_c**(2.-gamma) + 0.5*y_c**(3.-gamma)


def y_func(gamma,a): #really feel dumb naming this way ...
    if gamma == 2:
        return np.exp(-1*Psi(a))
    else:
        return (1.-(2.-gamma)*Psi(a))**(1./(2.-gamma))


def f_integrand(r,y,gamma,E,a):


    y = y_func(gamma,a)

    return (1.-y)**2. * (gamma+2.*y + (4.-gamma)*y**2.)/( y**(4.-gamma) * np.sqrt(E_c(r,a,gamma)-Psi(a)))



def df_f(r,a,gamma):

    return (3.-gamma)*M / (2. * (2.*np.pi*G*M*a)**(1.5)) * quad(f_integrand,r,args=(None))





def main():

    #will get a warning about convergence at the 0 limit, but is okay
    plot_sigma(0.00,100,0.01,1,1.5)

if __name__ == '__main__':
    main()



