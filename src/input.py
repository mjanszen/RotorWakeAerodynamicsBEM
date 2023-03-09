# All Inputs from the Assignment 

import numpy as np
import pandas as pd 
import scipy

radius = 50                         # radius of the rotor
n_blades = 3                        # number of blades
inner_radius = 0.2 * radius         # inner end of the blade section
pitch_deg = -2                      # pitch in degrees
pitch = np.radians(pitch_deg)       # pitch angle in radian
yaw_angles = np.radians([0,15,30])  # yaw angles to be calculated in radians

### Operational data 
v_0 = 10                            # [m] Wind speed
tsr = [6,8,10]                      # Tip speed ratios  to be calculated





airfoil = pd.read_excel("../data/polar.xlsx",skiprows=3)    # read in the airfoil. Columns [alpha, cl, cd cm]


### Postprocess 


# with a, r, a' , ft, fn calculate power, thrust 

r_list = np.linspace(inner_radius, radius,10)


print("Done!")
