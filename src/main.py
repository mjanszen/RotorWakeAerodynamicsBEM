from BEM import BEM

do = {
    "a": False,
    "b": False,
    "c": False
}

bem = BEM()
bem.set_constants(rotor_radius=1,
                  n_blades=3,
                  air_density=1.225)

if do["a"]:
    # Parameters

    radius = 50                         # radius of the rotor
    n_blades = 3                        # number of blades
    inner_radius = 0.2 * radius         # inner end of the blade section
    pitch_deg = -2                      # pitch in degrees
    pitch = np.radians(pitch_deg)       # pitch angle in radian
    yaw_angles = np.radians([0,15,30])  # yaw angles to be calculated in radians
    v_0 = 10                            # [m] Wind speed
    tsr = [6,8,10]                      # Tip speed ratios  to be calculated
    pass

if do["b"]:
    pass

if do["c"]:
    pass
