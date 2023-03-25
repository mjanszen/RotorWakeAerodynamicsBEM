#### Optimizing CP 


# We optimise TSR 6, as here the balde obviously stalls. 
# We optimize Chord and twist 

# Kutta Joukowsky 


from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
from task5 import task5
import matplotlib.pyplot as plt
helper = Helper()

# Choose whicht parts of the code to run 
def task7():

    bem = BEM(data_root="../data",
              file_airfoil="polar.xlsx")

    bem.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
    bem.optimize_TUD2(wind_speed=10, tip_speed_ratio=6, pitch=-2)
    

if __name__=="__main__":
    task7()
