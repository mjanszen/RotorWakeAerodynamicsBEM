from BEM import BEM
import numpy as np

do = {
    "a": False,
    "b": False,
    "c": False
}

bem = BEM(root="data",
          file_airfoil="polar.xlsx",
          save_dir="results")


if do["a"]:
    # Parameters
    bem.set_constants(rotor_radius=50,
                      n_blades=3,
                      air_density=1.225)
    # the lines below this will go into the BEM class
    inner_radius = 0.2*radius         # inner end of the blade section
    pitch_deg = -2                      # pitch in degrees
    pitch = np.radians(pitch_deg)       # pitch angle in radian
    yaw_angles = np.radians([0,15,30])  # yaw angles to be calculated in radians
    v_0 = 10                            # [m] Wind speed
    tsr = [6,8,10]                      # Tip speed ratios  to be calculated

    # Calc omega
    # calc the Cl Cd
    # Calc Cn, Ct
    # calc prandtl correction
    # update a 

if do["b"]:
    pass

if do["c"]:
    pass
