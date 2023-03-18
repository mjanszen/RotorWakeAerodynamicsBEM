from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
helper = Helper()

# Choose whicht parts of the code to run 
do = {
    "a": True,          # The blade from the assignment
    "DTU_test": False,  # The DTU 10 MW RWT 
    "b": False,
    "c": False
}

bem = BEM(root="../data",
          file_airfoil="polar.xlsx",
          save_dir="results")


if do["a"]:
    # Parameters
    bem.set_constants(rotor_radius=50,
                      root_radius=10,
                      n_blades=3,
                      air_density=1.225)
    # Calculation
    results = bem.solve(wind_speed=10, tip_speed_ratio=8, pitch=np.radians(2))
    
    # Plots
    fig, ax = plt.subplots()
    ax.plot(results["positions"], results["a"], label="a")
    ax.plot(results["positions"], results["a_prime"], label="a'")
    helper.handle_axis(ax, legend=True, grid=True, x_label="radius in m", font_size=15, line_width=3)
    plt.show()

    #### All params for the assignment
    #pitch_deg = -2                      # pitch in degrees
    #pitch = np.radians(pitch_deg)       # pitch angle in radian
    #yaw_angles = np.radians([0,15,30])  # yaw angles to be calculated in radians
    #v_0 = 10                            # [m] Wind speed
    #tsr = [6,8,10]                      # Tip speed ratios  to be calculated

    # Calc omega
    # calc the Cl Cd
    # Calc Cn, Ct
    # calc prandtl correction
    # update a

if do["DTU_test"]:
    # Parameters
    bem.set_constants(rotor_radius=31,
                      root_radius=0,
                      n_blades=3,
                      air_density=1.225)

    # Calculation
    results = bem.solve(wind_speed=8, tip_speed_ratio=2.61*31/8, pitch=np.deg2rad(-3), resolution=1,
                        glauert_correction_type="dtu", blade_end_correction_type="dtu", root_loss_correction=False)
    
    # Print
    print(results)

if do["b"]:
    pass

if do["c"]:
    pass
