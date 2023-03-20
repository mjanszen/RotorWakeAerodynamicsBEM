# Tip loss correction stuff for task 5

from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
helper = Helper()

# Choose whicht parts of the code to run 
def task5():
    print("Task 5: Start")

    ### Create a bem object for every calculation
    bem_tsr6 = BEM(data_root="../data",
              file_airfoil="polar.xlsx")

    bem_tsr8 = BEM(data_root="../data",
              file_airfoil="polar.xlsx")
    bem_tsr10 = BEM(data_root="../data",
              file_airfoil="polar.xlsx")
    
    # Bem object for calc without tip loss 
    bem_tsr8_no_tip_loss = BEM(data_root="../data",
              file_airfoil="polar.xlsx")
    # Parameters
    bem_tsr6.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
    bem_tsr8.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
    bem_tsr10.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)

    bem_tsr8_no_tip_loss.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)

    ##### Calculation at 3 differnt tsrs
    resolution = 200
    bem_tsr6.solve_TUD(wind_speed=10, tip_speed_ratio=6, pitch=-2, resolution=resolution)
    bem_tsr8.solve_TUD(wind_speed=10, tip_speed_ratio=8, pitch=-2 , resolution=resolution)
    bem_tsr10.solve_TUD(wind_speed=10, tip_speed_ratio=10, pitch=-2, resolution =resolution)


    bem_tsr8_no_tip_loss.solve_TUD(wind_speed=10, tip_speed_ratio=8, pitch=-2 , resolution=resolution,tip_loss_correction=False, root_loss_correction = False)
    #### Grab data that we want


    #breakpoint()

    #### Plots


    ## Plot for the tip loss correction factors
    fig, axs = plt.subplots(1,1, figsize=[6,4])
    axs.plot(bem_tsr6.df_results.r_centre.iloc[-resolution+1:],bem_tsr6.df_results.end_correction.iloc[-resolution+1:],label ="Tsr = 6")
    axs.plot(bem_tsr8.df_results.r_centre.iloc[-resolution+1:],bem_tsr8.df_results.end_correction.iloc[-resolution+1:],label ="Tsr = 8")
    axs.plot(bem_tsr10.df_results.r_centre.iloc[-resolution+1:],bem_tsr10.df_results.end_correction.iloc[-resolution+1:],label ="Tsr = 10")
    axs.legend()
    axs.set_xlabel("Radial position [m]")
    axs.set_ylabel("Tip and root loss correction factor []")
    axs.grid()
    #plt.show()
    fig.savefig("../results/tip_root_loss.png",bbox_inches="tight")


    # Plot to show the difference of adding the tip loss
    fig2, axs2 = plt.subplots(1,1, figsize=[6,4])
    
    axs2.plot(bem_tsr8.df_results.r_centre.iloc[-resolution+1:],bem_tsr8.df_results.a.iloc[-resolution+1:],label ="prandtl corrected")
    axs2.plot(bem_tsr8_no_tip_loss.df_results.r_centre.iloc[-resolution+1:],bem_tsr8_no_tip_loss.df_results.a.iloc[-resolution+1:],label ="no correction")
    
    axs2.legend()
    axs2.set_xlabel("Radial position [m]")
    axs2.set_ylabel(r"Axial induciton factor $a$ []")
    axs2.grid()
    plt.show()
    fig2.savefig("../results/tip_root_loss_induction.png",bbox_inches="tight")
    
    # Final stuff
    print("Task 5: Done")

if __name__ == "__main__":
    task5()

