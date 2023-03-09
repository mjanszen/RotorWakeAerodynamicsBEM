# All Inputs from the Assignment 

import numpy as np
import pandas as pd 



radius = 50                         # radius of the rotor
n_blades = 3                        # number of blades
inner_radius = 0.2 * radius         # inner end of the blade section
pitch_deg = -2                      # pitch in degrees
pitch = np.radians(pitch_deg)       # pitch angle in radian
yaw_angles = np.radians([0,15,30])  # yaw angles to be calculated in radians

### Operational data 
v_0 = 10                            # [m] Wind speed
tsr = [6,8,10]                      # Tip speed ratios  to be calculated


def get_twist(r, r_max):
    """
        function to get the twist along the blade in radians
        r: radial position along the blade in [m]
        r_max: maximum radius of the blade in [m]
        out: twist in radians 
    """
    return np.radians(14*(1 - r/r_max)) 

def get_chord(r, r_max):
    """
        function to calculate chord along the span in m
        r: radial position along the blade in [m]
        r_max: maximum radius of the blade in [m]
    """
    return 3*(1-r/r_max) + 1 


airfoil = pd.read_excel("../data/polar.xlsx",skiprows=3)    # read in the airfoil. Columns [alpha, cl, cd cm]


### Test
print("Done!")
