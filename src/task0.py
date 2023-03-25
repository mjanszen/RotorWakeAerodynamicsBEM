# Task 0: calculate the Reynolds number over the span

from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
helper = Helper()


def get_reynolds(U,L,nu):
    """
    Reynolds number calculation
    """
    return np.divide( U*L,nu) 

def get_velocity(r,rmax,tsr,u_0):
    """
    Calculate the absolute velocity array
    """
    u_rot = r * (tsr/rmax) * u_0 
    u = np.sqrt(np.multiply(u_rot,u_rot) + u_0**2) # absolute velocity
    return u 

def get_chord(r, r_max):
    """
        function to calculate chord along the span in m
        r: radial position along the blade in [m]
        r_max: maximum radius of the blade in [m]
    """
    return 3*(1-r/r_max)+1

def task0():
    nu = 1.5 * 10**-5
    u_0 = 10
    inner_radius = 10
    outer_radius = 50
    tsr_range = [6,8,10]
    radii = np.linspace(inner_radius, outer_radius,200) 
    #chords = [1] * len(radii) # np.ones((len(radii),1),) #
    chords = [get_chord(r,outer_radius) for r in radii]
    breakpoint()

    u_6 = get_velocity(radii, outer_radius,tsr_range[0],u_0)
    re_6 = get_reynolds(u_6,chords,nu)
    u_8 = get_velocity(radii, outer_radius,tsr_range[1],u_0)
    re_8 = get_reynolds(u_8,chords,nu)
    u_10 = get_velocity(radii, outer_radius,tsr_range[2],u_0)
    re_10 = get_reynolds(u_10,chords,nu)

    #plt.plot(radii, re_6)
    fig, axs = plt.subplots(3,1) 
    axs[0].plot(radii, chords)
    axs[1].plot(radii, u_6)
    axs[1].plot(radii, u_8)
    axs[1].plot(radii, u_10)
    axs[2].plot(radii, u_6*chords /nu)
    axs[2].plot(radii, u_8*chords /nu)
    axs[2].plot(radii, u_10*chords /nu)
    
    axs[2].set_xlabel("Radial position [m]")
    axs[0].set_ylabel("Chord [m]")
    axs[1].set_ylabel("V [m/s]")
    axs[2].set_ylabel("Reynolds number []")
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    #plt.show()
    plt.savefig("../results/reynolds_number.png",bbox_inches="tight")

    

    





if __name__ =="__main__": 
    task0()
