from BEM import BEM
import numpy as np
from helper_functions import Helper
import matplotlib.pyplot as plt
from task5 import task5
from task7 import task7
from task10 import task10
helper = Helper()

# Choose whicht parts of the code to run 
do = {
    "different_tsr": False,
    "plots": False,
    "c": False,
    "task5": False,
    "task7": False,
    "task10": True,
    "test": False
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

if do["task7"]:
    #bem.optimize_TUD(wind_speed=10, tip_speed_ratio=6, pitch=-2)
    task7()
    pass

if do["task10"]:
    # Plot circulation for TSR = 6,8,10
    task10()
    pass

if do["plots"]:
    pass

if do["c"]:
    pass

if do["test"]:
    # This is a test/ example for using the bem function and using the outputs 
    # the object property "current_results" is a pandas dataframe with the results from the last / current computation
    bem_test = BEM(data_root="../data",
          file_airfoil="polar.xlsx")
    bem_test.set_constants(rotor_radius=50, root_radius=50*0.2, n_blades=3, air_density=1.225)
    # Calculation
    bem_test.solve_TUD(wind_speed=10, tip_speed_ratio=8, pitch=-2)
    breakpoint()
    fig, axs = plt.subplots(4,1)
    axs[0].plot(bem_test.current_results.r_inner, bem_test.current_results.a)
    axs[1].plot(bem_test.current_results.r_inner, bem_test.current_results.a_prime)
    axs[2].plot(bem_test.current_results.r_inner, bem_test.current_results.alpha)
    axs[3].plot(bem_test.current_results.r_inner, bem_test.current_results.end_correction)
    
    ##### Make plot look a bit nicer
    # the following three lines do the same (and slightly more) as the commented lines below them.
    helper.handle_axis(axs, x_label="Radial position of the blade element centres (m)", grid=True,
                       y_label=["Induction (-)", "Induction (-)", r"$\alpha$", "Blade end correction loss (-)"])
    helper.handle_figure(fig, size=(5,7), show=True)
    # axs[0].set_xlabel("Radial position of the center point[m]")
    # axs[1].set_xlabel("Radial position of the center point[m]")
    # axs[2].set_xlabel("Radial position of the center point[m]")
    # axs[3].set_xlabel("Radial position of the center point[m]")

    # axs[0].set_ylabel("Induction []")
    # axs[1].set_ylabel("Induction []")
    # axs[2].set_ylabel(r"$\alpha$")
    # axs[3].set_ylabel("Tip loss factor []")
    # plt.show()
    print("Done testing")



