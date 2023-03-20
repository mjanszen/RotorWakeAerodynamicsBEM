from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
from task5 import task5
helper = Helper()

# Choose whicht parts of the code to run 
do = {
    "different_tsr": True,
    "plots": False,
    "c": False,
    "task5": False
}

bem = BEM(data_root="../data",
          file_airfoil="polar.xlsx")


if do["different_tsr"]:
    # Parameters
    bem.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
    # Calculation
    for tsr in np.linspace(6,10,41):
        bem.solve_TUD(wind_speed=10, tip_speed_ratio=tsr, pitch=-2)

if do["task5"]:
    task5()

if do["plots"]:
    pass

if do["c"]:
    pass
