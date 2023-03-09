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
    pass

if do["b"]:
    pass

if do["c"]:
    pass