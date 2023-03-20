# Tip loss correction stuff for task 5

from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
helper = Helper()

# Choose whicht parts of the code to run 

bem_tsr6 = BEM(data_root="../data",
          file_airfoil="polar.xlsx")

bem_tsr8 = BEM(data_root="../data",
          file_airfoil="polar.xlsx")
bem_tsr10 = BEM(data_root="../data",
          file_airfoil="polar.xlsx")

# Parameters
bem_tsr6.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
bem_tsr8.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
bem_tsr10.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
# Calculation
bem_tsr6.solve_TUD(wind_speed=10, tip_speed_ratio=6, pitch=-2)
bem_tsr8.solve_TUD(wind_speed=10, tip_speed_ratio=8, pitch=-2)
bem_tsr10.solve_TUD(wind_speed=10, tip_speed_ratio=10, pitch=-2)

plt.subplots(1,1)
plt.plot(bem_tsr6.df_results.r_centre,bem_tsr6.df_results.end_correction,label ="Tsr = 6")
plt.plot(bem_tsr8.df_results.r_centre,bem_tsr8.df_results.end_correction,label ="Tsr = 8")
plt.plot(bem_tsr10.df_results.r_centre,bem_tsr10.df_results.end_correction,label ="Tsr = 10")
plt.legend()
plt.grid()

plt.show()

breakpoint()
print("Done")



