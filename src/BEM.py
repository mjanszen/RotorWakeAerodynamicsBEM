import numpy as np
from scipy.optimize import brentq, newton, minimize
import pandas as pd
import scipy.interpolate as interpolate
from helper_functions import Helper
import scipy
import matplotlib.pyplot as plt
helper=Helper()


class BEM:
    def __init__(self,
                 data_root: str,
                 file_airfoil: str):
        self.rotor_radius = None
        self.root_radius = None
        self.n_blades = None
        self.air_density = None
        self.root = data_root
        self.a_prime = 0
        try:
            self.df_results = pd.read_csv(data_root+"/BEM_results.dat")
        except FileNotFoundError:
            self.df_results = pd.DataFrame()

        df_tmp = pd.read_excel(data_root+"/"+file_airfoil, skiprows=3)
        self.interp = {"c_l": interpolate.interp1d(df_tmp["alpha"], df_tmp["cl"]),
                       "c_d": interpolate.interp1d(df_tmp["alpha"], df_tmp["cd"])}
        self.implemented_glauert_correction = ["none", "dtu", "tud"]
        self.implemented_tip_correction = ["none", "dtu", "tud"]
        self.implemented_root_correction = ["none", "tud"]
        self.blade_end_correction = 1

    def set_constants(self,
                      rotor_radius: float,
                      root_radius: float,
                      n_blades: int,
                      air_density: float) -> None:
        self._set(**{param: value for param, value in locals().items() if param != "self"})
        return None

    def solve_TUD(self,
                  wind_speed: float,
                  tip_speed_ratio: float,
                  pitch: float or np.ndarray,
                  start_radius: float = None,
                  resolution: int = 200,
                  max_convergence_error: float=1e-5,
                  max_iterations: int=200,
                  tip_loss_correction: bool=True,
                  root_loss_correction: bool=True) -> None:
        """
        Function to run the BEM loop
        Glauert_correction: either 'tud' (TU Delft) or 'dtu' (Denmark's TU). Same for blade_end_correction
        All angles must be in rad.
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch: IN DEGREE
        :param resolution:
        :return:
        """
        start_radius = start_radius if start_radius is not None else self.root_radius
        # Initialise the result containers
        results = {
            "r_centre": list(),     # radius used for the calculations
            "r_inner": list(),      # inner radius of the blade element
            "r_outer": list(),      # outer radius of the blade element
            "a": list(),            # Axial Induction factor
            "a_prime": list(),      # Tangential induction factor
            "f_n": list(),          # Forces normal to the rotor plane in N/m
            "f_t": list(),          # Forces tangential in the rotor plane in N/m
            "bec": list(),          # blade end correction (depending on 'tip' and 'root')
            "C_T": list(),          # thrust coefficient
            "alpha": list(),        # angle of attack
            "circulation": list(),  # magnitude of the circulation using Kutta-Joukowski
            "v0": list(),           # flow velocity normal to rotor plane
            "tsr": list(),          # tip speed ratio
            "pitch": list(),         # pitch in degree
            "end_correction": list(),# Prandtl tip and root loss factor
            "cl": list(),
            "cd": list()
        }
        # delete data with same wind speed, tip speed ratio and pitch angle.
        try:
            self.df_results = self.df_results.loc[~((self.df_results["tsr"]==tip_speed_ratio) &
                                                    (self.df_results["v0"]==wind_speed) &
                                                    (self.df_results["pitch"]==pitch))]
        except KeyError:
            pass
        pitch = np.deg2rad(pitch)
        # Calculate the rotational speed
        #breakpoint()
        omega = tip_speed_ratio*wind_speed/self.rotor_radius
        radii = np.linspace(start_radius, self.rotor_radius, resolution)
        # Loop along the span of the blade (blade sections)
        print(f"Doing BEM for v0={wind_speed}, tsr={tip_speed_ratio}, pitch={pitch}")
        for r_inside, r_outside in zip(radii[:-1], radii[1:]):      # Take the left and right radius of every element
            r_centre = (r_inside+r_outside)/2                       # representative radius in the middle of the section
            elem_length = r_outside-r_inside                        # length of elemen
            # Get/Set values from the local section
            chord = self._get_chord(r_centre, self.rotor_radius)  # Get the chord
            twist = self._get_twist(r_centre, self.rotor_radius)  # Get the twist
            area_annulus = np.pi*(r_outside**2-r_inside**2)
            a, a_new, a_prime, converged = 1/3, 0, 0, False
            for i in range(max_iterations):
                # get inflow angle
                phi = self._phi(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r_centre)
                # get combined lift and drag coefficient projected into the normal and tangential direction
                _, _, _, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch,
                                                             tip_seed_ratio=tip_speed_ratio, university="tud")
                # get the inflow speed for the airfoil
                inflow_speed = self._inflow_velocity(wind_speed, a, a_prime, omega, r_centre)
                # get thrust force (in N) of the whole turbine at the current radius
                thrust = self._aero_force(inflow_speed, chord, c_n)*self.n_blades*elem_length
                # calculate thrust coefficient that results from the blade element
                C_T = thrust/(1/2*self.air_density*wind_speed**2*area_annulus)
                # get Glauert corrected axial induction factor
                a_new = self._a(C_T=C_T)
                # get the combined (tip and root) correction factor
                blade_end_correction = self._blade_end_correction(which="tud", tip=tip_loss_correction,
                                                                  root=root_loss_correction, radius=r_centre,
                                                                  tip_seed_ratio=tip_speed_ratio, a=a_new)
                # correct the Glauert corrected axial induction factor with the blade end losses
                a_new /= blade_end_correction
                # update the axial induction factor for the next iteration
                a = 0.75*a+0.25*a_new
                # get the tangential force (in N/m) of the whole turbine at the current radius
                f_tangential = self._aero_force(inflow_speed, chord, c_t)*self.n_blades
                # get the tangential induction factor that corresponds to this force AND correct it for tip losses
                a_prime = self._a_prime(f_tangential, r_centre, wind_speed, a, tip_speed_ratio)/blade_end_correction
                # check if the axial induction factor has converged. If it has, the tangential induction factor has too
                if np.abs(a-a_new) < max_convergence_error:
                    converged = True
                    break
            # notify user if loop did not converge, but was stopped by the maximum number of iterations
            if not converged:
                print(f"BEM did not converge for the blade element between {r_inside}m and {r_outside}m. Current "
                      f"change after {max_iterations}: {np.abs(a-a_new)}.")

            # Now that we have the converged axial induction factor, we can get the rest of the values
            phi = self._phi(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r_centre)
            alpha, c_l, c_d, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r_centre,
                                                               tip_seed_ratio=tip_speed_ratio, university="tud")
            inflow_speed = self._inflow_velocity(wind_speed, a, a_prime, omega, r_centre)

            # Assemble the result output structure
            results["r_centre"].append(r_centre)
            results["r_inner"].append(r_inside)
            results["r_outer"].append(r_outside)
            results["a"].append(a)
            results["a_prime"].append(a_prime)
            results["f_n"].append(self._aero_force(inflow_velocity=inflow_speed, chord=chord, force_coefficient=c_n))
            results["f_t"].append(self._aero_force(inflow_velocity=inflow_speed, chord=chord, force_coefficient=c_t))
            results["bec"].append(self._blade_end_correction(which="tud", tip=tip_loss_correction,
                                                             root=root_loss_correction, radius=r_centre,
                                                             tip_seed_ratio=tip_speed_ratio, a=a))
            results["C_T"].append(self._C_T(a))
            results["alpha"].append(alpha)
            results["cl"].append(c_l)
            results["cd"].append(c_d)
            results["circulation"].append(1/2*inflow_speed*c_l*chord)
            results["v0"].append(wind_speed)
            results["tsr"].append(tip_speed_ratio)
            results["pitch"].append(np.rad2deg(pitch))
            results["end_correction"].append(blade_end_correction)
        self.current_results = pd.DataFrame(results)
        self.df_results = pd.concat([self.df_results, pd.DataFrame(results)])
        self.df_results.to_csv(self.root+"/BEM_results.dat", index=False)
        return None
    
    def optimize_TUD(self,
                  wind_speed: float,
                  tip_speed_ratio: float,
                  pitch: float or np.ndarray,
                  start_radius: float = None,
                  resolution: int = 200,
                  max_convergence_error: float=1e-5,
                  max_iterations: int=200,
                  tip_loss_correction: bool=True,
                  root_loss_correction: bool=True) -> None:
        """
        Optimze the chord and twist per blade section 
        Glauert_correction: either 'tud' (TU Delft) or 'dtu' (Denmark's TU). Same for blade_end_correction
        All angles must be in rad.
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch: IN DEGREE
        :param resolution:
        :return:
        """
        start_radius = start_radius if start_radius is not None else self.root_radius
        # Initialise the result containers
        results = {
            "r_centre": list(),     # radius used for the calculations
            "r_inner": list(),      # inner radius of the blade element
            "r_outer": list(),      # outer radius of the blade element
            "a": list(),            # Axial Induction factor
            "a_prime": list(),      # Tangential induction factor
            "f_n": list(),          # Forces normal to the rotor plane in N/m
            "f_t": list(),          # Forces tangential in the rotor plane in N/m
            "bec": list(),          # blade end correction (depending on 'tip' and 'root')
            "C_T": list(),          # thrust coefficient
            "alpha": list(),        # angle of attack
            "circulation": list(),  # magnitude of the circulation using Kutta-Joukowski
            "v0": list(),           # flow velocity normal to rotor plane
            "tsr": list(),          # tip speed ratio
            "pitch": list(),         # pitch in degree
            "end_correction": list(),# Prandtl tip and root loss factor
            "cl": list(),
            "cd": list()
        }
        # delete data with same wind speed, tip speed ratio and pitch angle.
        try:
            self.df_results = self.df_results.loc[~((self.df_results["tsr"]==tip_speed_ratio) &
                                                    (self.df_results["v0"]==wind_speed) &
                                                    (self.df_results["pitch"]==pitch))]
        except KeyError:
            pass
        pitch = np.deg2rad(pitch)
        # Calculate the rotational speed
        #breakpoint()

        omega = tip_speed_ratio*wind_speed/self.rotor_radius
        # go from middle from now 
        #radii = np.linspace(0.55*(self.rotor_radius - start_radius), self.rotor_radius, resolution)
        radii = np.linspace(start_radius, self.rotor_radius, resolution)
        radii_left = radii[:int((resolution/2))][::-1]
        radii_right = radii[int(resolution/2):]
      
        # initialize arrays to store the optimum values
        chord_list = np.zeros(len(radii))
        twist_list = np.zeros(len(radii))
        # Loop along the span of the blade (blade sections)
        print(f"Doing BEM for v0={wind_speed}, tsr={tip_speed_ratio}, pitch={pitch}")
        
        # loop through outer part

        #for r_inside, r_outside in zip(radii[:-1], radii[1:]):      # Take the left and right radius of every element
        iter = 0

        chord_estimate = 4.8
        twist_estimate = 0.15
        for r_inside, r_outside in zip(radii_left[:-1], radii_right[1:]):      # Take the left and right radius of every element
            r_centre = (r_inside+r_outside)/2                       # representative radius in the middle of the section
            elem_length = r_outside-r_inside                        # length of elemen
            # Get/Set values from the local section

            #twist = self._get_twist(r_centre, self.rotor_radius)  # Get the twist
            #chord = self._get_chord(r_centre, self.rotor_radius)  # Get the chord

            def get_force_from_bem(chord, twist): #chord_and_twist):
                """
                Function to return the tangential force for a given chord and twist 
                """

                #chord = chord_and_twist[0]
                #twist = chord_and_twist[1]
                try:
                    area_annulus = np.pi*(r_outside**2-r_inside**2)
                    a, a_new, a_prime, converged = 1/3, 0, 0, False
                    for i in range(max_iterations):
                        # get inflow angle
                        phi = self._phi(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r_centre)
                        # get combined lift and drag coefficient projected into the normal and tangential direction
                        _, _, _, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch,
                                                                     tip_seed_ratio=tip_speed_ratio, university="tud")
                        # get the inflow speed for the airfoil
                        inflow_speed = self._inflow_velocity(wind_speed, a, a_prime, omega, r_centre)
                        # get thrust force (in N) of the whole turbine at the current radius
                        thrust = self._aero_force(inflow_speed, chord, c_n)*self.n_blades*elem_length
                        # calculate thrust coefficient that results from the blade element
                        C_T = thrust/(1/2*self.air_density*wind_speed**2*area_annulus)
                        # get Glauert corrected axial induction factor
                        a_new = self._a(C_T=C_T)
                        # get the combined (tip and root) correction factor
                        blade_end_correction = self._blade_end_correction(which="tud", tip=tip_loss_correction,
                                                                          root=root_loss_correction, radius=r_centre,
                                                                          tip_seed_ratio=tip_speed_ratio, a=a_new)
                        # correct the Glauert corrected axial induction factor with the blade end losses
                        a_new /= blade_end_correction
                        # update the axial induction factor for the next iteration
                        a = 0.75*a+0.25*a_new
                        # get the tangential force (in N/m) of the whole turbine at the current radius
                        f_tangential = self._aero_force(inflow_speed, chord, c_t)*self.n_blades
                        # get the tangential induction factor that corresponds to this force AND correct it for tip losses
                        a_prime = self._a_prime(f_tangential, r_centre, wind_speed, a, tip_speed_ratio)/blade_end_correction
                        # check if the axial induction factor has converged. If it has, the tangential induction factor has too
                        if np.abs(a-a_new) < max_convergence_error:
                            converged = True
                            break
                    phi = self._phi(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r_centre)
                    alpha, c_l, c_d, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r_centre,
                                                                       tip_seed_ratio=tip_speed_ratio, university="tud")
                    return f_tangential, alpha , a ,a_prime, converged, blade_end_correction  # output the tangential force
                except:
                    return np.nan, np.nan, 0, 0 , True
            # Brute force through it ! 
            #chord_estimate = 4.8
            #twist_estimate = 0.15
            #chord_range = np.linspace(chord_estimate*0.5, chord_estimate*1.5,40)
            #twist_range = np.linspace(twist_estimate * 0.5, twist_estimate *1.5, 40) #np.deg2rad(np.linspace(7.4,11.4,100))
            c_l_opt = 1.1
            chord_estimate = self._calc_optimum_chord(r_centre,self.rotor_radius,self.n_blades,c_l_opt,tip_speed_ratio)
            chord_range = np.arange(chord_estimate*0.5, chord_estimate*1.5,0.05) # 5cm steps should be alright resolution
            twist_range = np.arange(twist_estimate * 0.5, twist_estimate *1.5, np.deg2rad(1)) #np.deg2rad(np.linspace(7.4,11.4,100)) # 1 degree steps
            #chord_range = [2,3,4,5]
            #twist_range = [0.13,0.14,0.15, 0.16]
            ft_array =np.empty((len(chord_range), len(twist_range))) # initialize chord x twist array
            alpha_array =np.empty((len(chord_range), len(twist_range))) # initialize chord x twist array
            for j,twist in enumerate(twist_range):
                for i,chord in enumerate(chord_range):
                    ft_array[i,j], alpha_array[i,j] , a , a_prime, converged, blade_end_correction= get_force_from_bem(chord,twist)
               

            # Now we have these optimums. Get the values of chord and twist
            max_indices = np.unravel_index(np.argmax(ft_array, axis=None), ft_array.shape)
            chord_list[iter+int(resolution/2)]  = chord_range[max_indices[0]]
            twist_list[iter+int(resolution/2)]  = twist_range[max_indices[1]]
            chord_estimate  = chord_range[max_indices[0]]
            twist_estimate  = twist_range[max_indices[1]]
            print(f"Estimates:\n Chord: {chord_estimate} \n Twist: {twist_estimate}")
            
            #import matplotlib.pyplot as plt
            #plt.contourf(np.rad2deg(twist_range), chord_range,ft_array,20)
            #plt.xlabel("Twist [deg]")
            #plt.ylabel("Chord [m]")
            #plt.show()
            
            #initial_guess = 5
            #bounds = (0,1)
            #optimum = minimize(get_force_from_bem,initial_guess, method="TNC")
            #breakpoint()
            # notify user if loop did not converge, but was stopped by the maximum number of iterations
            if not converged:
                print(f"BEM did not converge for the blade element between {r_inside}m and {r_outside}m. Current "
                      f"change after {max_iterations}: {np.abs(a-a_new)}.")

            # Now that we have the converged axial induction factor, we can get the rest of the values
            phi = self._phi(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r_centre)
            alpha, c_l, c_d, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r_centre,
                                                               tip_seed_ratio=tip_speed_ratio, university="tud")
            inflow_speed = self._inflow_velocity(wind_speed, a, a_prime, omega, r_centre)

            # Assemble the result output structure
            results["r_centre"].append(r_centre)
            results["r_inner"].append(r_inside)
            results["r_outer"].append(r_outside)
            results["a"].append(a)
            results["a_prime"].append(a_prime)
            results["f_n"].append(self._aero_force(inflow_velocity=inflow_speed, chord=chord, force_coefficient=c_n))
            results["f_t"].append(self._aero_force(inflow_velocity=inflow_speed, chord=chord, force_coefficient=c_t))
            results["bec"].append(self._blade_end_correction(which="tud", tip=tip_loss_correction,
                                                             root=root_loss_correction, radius=r_centre,
                                                             tip_seed_ratio=tip_speed_ratio, a=a))
            results["C_T"].append(self._C_T(a))
            results["alpha"].append(alpha)
            results["cl"].append(c_l)
            results["cd"].append(c_d)
            results["circulation"].append(1/2*inflow_speed*c_l*chord)
            results["v0"].append(wind_speed)
            results["tsr"].append(tip_speed_ratio)
            results["pitch"].append(np.rad2deg(pitch))
            results["end_correction"].append(blade_end_correction)
            iter +=1
        breakpoint()
        self.current_results = pd.DataFrame(results)
        self.df_results = pd.concat([self.df_results, pd.DataFrame(results)])
        self.df_results.to_csv(self.root+"/BEM_results.dat", index=False)
        return None

    def solve_DTU(self,
                  wind_speed: float,
                  tip_speed_ratio: float,
                  pitch: float or np.ndarray,
                  start_radius: float=None,
                  resolution: int=100,
                  brent_bracket: tuple=(0,0.9),
                  glauert_correction_type: str="tud",
                  blade_end_correction_type: str="tud",
                  tip_loss_correction: bool=True,
                  root_loss_correction: bool=True):
        """
        Function to run the BEM loop
        Glauert_correction: either 'tud' (TU Delft) or 'dtu' (Denmark's TU). Same for blade_end_correction
        All angles must be in rad.
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch:
        :param resolution:
        :return:
        """
        
        start_radius = start_radius if start_radius is not None else self.root_radius

        # Initialize the result containers
        results = {
            "positions": list(),    # Positions along the blade
            "a": list(),            # Axial Induction
            "a_prime": list(),      # Tangential induction
            "f_n": list(),          # Forces normal to the rotor plane
            "f_t": list(),          # Forces tangential in the rotor plane
            "bec": list(),          # blade end correction
            "v0": wind_speed        # V_infinity 
        }
        
        # Calculate the rotational speed
        omega = tip_speed_ratio*wind_speed/self.rotor_radius

        # Loop along the span of the blade (blade sections)
        for r in np.linspace(start_radius, self.rotor_radius, resolution):
           
            # Get/Set values from the local section
            self.a_prime = 0            # Initialize a start value for the induction
            self.blade_end_correction = 1
            chord = self._get_chord(r, self.rotor_radius)   # Get the chord
            twist = self._get_twist(r, self.rotor_radius)   # Get the twist
            local_solidity = self._local_solidity(chord, r, self.n_blades)  # Calculate the solidity (solidity is something like a parameter to show, how much we block the flow by blades/ blockage)

            def residue(a):
                """
                Function that calculates the difference of the C_T values from a given induction
                
                IN: 
                    a: axial induction factor

                OUT:
                    Difference between the CT values
                """
                
                # Calculate the inflow angle phi
                phi = self._phi(a=a, a_prime=self.a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r)
                
                # Calculate the aerodynamic properties that result from the Blade element theory
                aero_values = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r, a=a,
                                                       tip_seed_ratio=tip_speed_ratio,
                                                       blade_end_correction_type=blade_end_correction_type,
                                                       root=root_loss_correction, tip=tip_loss_correction)
                c_n, c_t, blade_end_correction = aero_values[3], aero_values[4], aero_values[5]
                # Update the tangential induction based on the other params
                self._update_a_prime(local_solidity=local_solidity, c_tangential=c_t, inflow_angle=phi,
                                     blade_end_correction=blade_end_correction)
                return self._equate_blade_element_and_momentum(glauert_correction=glauert_correction_type, a=a,
                                                               blade_end_correction=blade_end_correction, phi=phi,
                                                               local_solidity=local_solidity, c_normal=c_n)
            
            # Optimize the difference of the CT values and get the corresponding induction
            try:
                a = brentq(residue, *brent_bracket)     # Brent method, has some problems sometimes
            except ValueError:
                # When the Brent method does not work, we optimize using the Newton method
                print("Brent could not be used for the convergence, using Newton instead.")
                a =  newton(residue, 1/3)

            # Now that we have optimum a, we can get the rest of the values 
            phi = np.arctan((1-a)*wind_speed/((1+self.a_prime)*omega*r))
            aero_values = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r, a=a,
                                                   tip_seed_ratio=tip_speed_ratio,
                                                   blade_end_correction_type=blade_end_correction_type,
                                                   root=root_loss_correction, tip=tip_loss_correction)
            c_n, c_t, blade_end_correction = aero_values[3], aero_values[4], aero_values[5]
            inflow_velocity = np.sqrt((omega*r*(1+self.a_prime))**2+(wind_speed*(1-a))**2)      # The actual wind speed seen by the blade section
           
            # Correct for non-possible values of the tangential induction
            if self.a_prime == "-inf":
                c_n = 0
                c_t = 0

            # Assemble the result output structure
            results["positions"].append(r)
            results["a"].append(a)
            results["a_prime"].append(self.a_prime)
            results["f_n"].append(self._aero_force(inflow_velocity=inflow_velocity, chord=chord, force_coefficient=c_n))
            results["f_t"].append(self._aero_force(inflow_velocity=inflow_velocity, chord=chord, force_coefficient=c_t))
            results["bec"].append(blade_end_correction)
        return results

    def _calculate_thrust(self, f_n, radial_positions):
        """
            Calculate thrust from the normal forces. Account for f_t = 0 at the tip.
        f_n: normal forces
        radial_positions: radial position along the blade matching the positions of f_n
        n_blades:   number of blades
        radius:     max radius
        """
        thrust = self.n_blades*scipy.integrate.simpson([*f_n, 0], [*radial_positions, self.rotor_radius])
        return thrust

    def _calculate_power(self, f_t, radial_positions, omega):
        """
            Calculate power from the normal forces. Account for f_n = 0 at the tip.
        f_t: tangential forces
        radial_positions: radial position along the blade matching the positions of f_n
        n_blades:   number of blades
        radius:     max radius
        omega:      [rad/s] rotational speed
        """
        power = omega*self.n_blades*scipy.integrate.simpson(np.multiply([*radial_positions, self.rotor_radius], [*f_t, 0]),
                                                           [*radial_positions, self.rotor_radius])
        return power

    def _calc_ct(self, thrust, velocity):
        """
            Calculate the thrust coefficient ct
        """
        return thrust/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**2))

    def _calc_cp(self, power, velocity):
        """
            Calculate the power coefficient ct
        """
        return power/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _calc_ct_distribution(self, f_n, velocity):
        """
        Calculate the distribution of the thrust coefficient along the blade
        f_n: normal forces along the blade
        radius: maximum radius of the Blade
        velocity: fluid velocity od V0
        """
        f_n = np.array(f_n)
        return f_n/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _calc_cp_distribution(self, f_t, velocity):
        """
        Calc the distribution of the power coeff. along the blade
        f_t: tangential forces along the blade
        radius: maximum radius of the Blade
        velocity: fluid velocity od V0
        """
        f_t = np.array(f_t)
        return f_t/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _phi_to_aero_values(self, phi: float, twist:float or np.ndarray, pitch: float,tip_seed_ratio: float,
                            university: str,
                            radius: float=None, a: float=None, blade_end_correction_type: str=None, tip: bool=None,
                            root: bool=None) -> tuple:
        alpha = np.rad2deg(phi-(twist+pitch))
        c_l = self.interp["c_l"](alpha)
        c_d = self.interp["c_d"](alpha)
        c_n = self._c_normal(phi, c_l, c_d)
        c_t = self._c_tangent(phi, c_l, c_d)
        if university == "tud":
            return alpha, c_l, c_d, c_n, c_t
        elif university == "dtu":
            return alpha, c_l, c_d, c_n, c_t, self._blade_end_correction(which=blade_end_correction_type, tip=tip,
                                                                         root=root, phi=phi, radius=radius, a=a,
                                                                         tip_seed_ratio=tip_seed_ratio)

    def _set(self, **kwargs) -> None:
        """
        Sets parameters of the instance. Raises an error if a parameter is trying to be set that doesn't exist.
        :param kwargs:
        :return:
        """
        existing_parameters = [*self.__dict__]
        for parameter, value in kwargs.items():
            if parameter not in existing_parameters:
                raise ValueError(f"Parameter {parameter} cannot be set. Settable parameters are {existing_parameters}.")
            self.__dict__[parameter] = value
        return None

    def _assert_values(self):
        not_set = list()
        for variable, value in vars(self).items():
            if value is None:
                not_set.append(variable)
        if len(not_set) != 0:
            raise ValueError(f"Variable(s) {not_set} not set. Set all variables before use.")

    def _update_a_prime(self, local_solidity: float, c_tangential: float, blade_end_correction: float,
                        inflow_angle: float) -> None:
        self.a_prime = local_solidity*c_tangential*(1+self.a_prime)/(4*blade_end_correction*np.sin(inflow_angle)*
                                                                     np.cos(inflow_angle))

    def _aero_force(self, inflow_velocity: float, chord: float, force_coefficient: float):
        """
        Calculates the tangential force per unit span.
        :param inflow_velocity:
        :param chord:
        :param c_normal:
        :return:
        """
        return 1/2*self.air_density*inflow_velocity**2*chord*force_coefficient

    def _equate_blade_element_and_momentum(self, glauert_correction: str, a: float, blade_end_correction: float,
                                           phi: float, local_solidity: float, c_normal: float):
        """
        Function to calculate the difference of the CT values (Blade element vs momentum theory) from a given axial induction
        """
        if glauert_correction == "none":
            return 1/((4*blade_end_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-a
        elif glauert_correction == "dtu":
            if a < 1/3:
                return 1/((4*blade_end_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-a
            else:
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-4*a*blade_end_correction*(1-a/4*(5-3*a))
        elif glauert_correction == "tud":
            CT_1 = 1.816
            if a < 1-np.sqrt(CT_1)/2:
                return 4*a*(1-a)
            else:
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-(CT_1-4*(np.sqrt(CT_1)-1)*(1-a))

    def _a_prime(self, F_tangential: float, radius: float, wind_speed: float, a: float, tip_speed_ratio: float):
        return F_tangential/(4*self.air_density*np.pi*radius**2/self.rotor_radius*wind_speed**2*(1-a)*tip_speed_ratio)

    def _blade_end_correction(self, which: str, tip: bool, root: bool, radius: float,
                              tip_seed_ratio: float, a: float, phi: float=None) -> float:
        """
        Different Prandtl correction methods.
        :param which:
        :param tip:
        :param root:
        :param phi:
        :param radius:
        :param tip_seed_ratio:
        :param a:
        :return:
        """
        F = 1
        if tip:
            if which == "dtu":
                if np.sin(np.abs(phi)) < 0.01:
                    pass
                else:
                    F = 2/np.pi*np.arccos(np.exp(-(self.n_blades*(self.rotor_radius-radius))/(2*radius*np.sin(np.abs(phi)))))
            elif which == "tud":
                d = 2*np.pi/self.n_blades*(1-a)/(np.sqrt(tip_seed_ratio**2+(1-a)**2))
                F = 2/np.pi*np.arccos(np.exp(-np.pi*((self.rotor_radius-radius)/self.rotor_radius)/d))
            else:
                raise ValueError(f"Parameter 'blade_end_correction' must be one of "
                                 f"{self.implemented_tip_correction}, but was {which}.")
        if root:
            if which == "tud":
                d = 2*np.pi/self.n_blades*(1-a)/(np.sqrt(tip_seed_ratio**2+(1-a)**2))
                F *= 2/np.pi*np.arccos(np.exp(-np.pi*((radius-self.root_radius)/self.rotor_radius)/d))
            else:
                raise ValueError(f"Parameter 'blade_end_correction' must be one of "
                                 f"{self.implemented_root_correction}, but was {which}.")
        return F

    @staticmethod
    def _a(C_T: float) -> float:
        C_T1 = 1.816
        CT_2 = 2*np.sqrt(C_T1)-C_T1
        if C_T < CT_2:
            return 1/2-np.sqrt(1-C_T)/2
        else:
            return 1+(C_T-C_T1)/(4*np.sqrt(C_T1)-4)

    @staticmethod
    def _C_T(a: float):
        C_T1 = 1.816
        if a < 1-np.sqrt(C_T1)/2:
            return 4*a*(1-a)
        else:
            return C_T1-4*(np.sqrt(C_T1)-1)*(1-a)


    @staticmethod
    def _c_normal(phi: float, c_lift: float, c_drag: float) -> float:
        """
        Calculates an aerodynamic "lift" coefficient according to a coordinate transformation with phi
        :param phi: angle between flow and rotational direction in rad
        :param c_lift: lift coefficient old coordinate system
        :param c_drag: lift coefficient old coordinate system
        :return: Normal force in Newton
        """
        return c_lift*np.cos(phi)+c_drag*np.sin(phi)

    @staticmethod
    def _c_tangent(phi: float, c_lift: float, c_drag: float) -> float:
        """
        Calculates an aerodynamic "drag" coefficient according to a coordinate transformation with phi
        :param phi: angle between flow and rotational direction in rad
        :param c_lift: lift coefficient old coordinate system
        :param c_drag: lift coefficient old coordinate system
        :return: Normal force in Newton
        """
        return c_lift*np.sin(phi)-c_drag*np.cos(phi)

    @staticmethod
    def _local_solidity(chord: float, radius: float, n_blades: int) -> float:
        """
        Calculates the local solidity
        :param chord: in m
        :param radius: distance from rotor axis to blade element in m
        :param n_blades: number of blades
        :return: local solidity
        """
        return n_blades*chord/(2*np.pi*radius)

    @staticmethod
    def _inflow_velocity(wind_speed: float, a: float, a_prime: float, rotational_speed: float, radius: float):
        return np.sqrt((wind_speed*(1-a))**2+(rotational_speed*radius*(1+a_prime))**2)

    @staticmethod
    def _phi(a: float, a_prime: float, wind_speed: float, rotational_speed: float, radius: float) -> float:
        """
        Function to calculate the inflow angle based on the two induction factors, the inflow velocity, radius and
        angular_velocity
        :param a:
        :param a_prime:
        :return:
        """
        return np.tan((1-a)*wind_speed/((1+a_prime)*rotational_speed*radius))

    @staticmethod
    def _get_twist(r, r_max):
        """
            function to get the twist along the blade in radians
            r: radial position along the blade in [m]
            r_max: maximum radius of the blade in [m]
            out: twist in radians
        """
        return np.radians(14*(1-r/r_max))

    @staticmethod
    def _get_chord(r, r_max):
        """
            function to calculate chord along the span in m
            r: radial position along the blade in [m]
            r_max: maximum radius of the blade in [m]
        """
        return 3*(1-r/r_max)+1

    @staticmethod
    def _calc_optimum_chord(r_pos,r_max,n_blades,c_l,tsr):
        """
        Based on the Wind turbine handbook
        Calculate the optimum chord length based on eq. 3.72 from the "wind energy handbook"
        IN:
            r_pos : [m] radial position 
            r_max : [m] maximum radius of the blade 
            n_blades: [] number of blades
            c_l:     [] Lift coefficient
            tsr: []  Tip-speed ratio

        
        """
        
        
        mu = r_pos/r_max     # relative spanwise position
        # calculate the rhs first
        rhs = (8/9) / np.sqrt( (1- 1/3)**2 + tsr**2 * mu**2 *( (1+ (2/(9 * tsr**2 * mu**2)))**2))
        #equation solved for c 
        c = (rhs * 2 * np.pi * r_max) / (n_blades *tsr *c_l)
        return c

    @staticmethod
    def _calc_optimum_twist(r_pos,r_max,tsr, alpha):
        """
        Based on the Wind turbine handbook
        Calculate the optimum twist angle beta in degrees based on eq. 3.74 from the "wind energy handbook"
        IN:
            r_pos : [m] radial position 
            r_max : [m] maximum radius of the blade 
            tsr:    [-] Tip-speed ratio
            alpha:  [-] Design angle of attack (minimum lift to drag ratio)
        """
        mu = r_pos / r_max
        phi = np.arctan( (1- 1/3) / (tsr * mu * (1+ (2/(9 * tsr**2 * mu**2)))) ) 
        phi_deg = np.rad2deg(phi)
        beta = phi_deg - alpha
        return beta

    @staticmethod
    def _tangential_induction_factor(phi: float, local_solidity: float, c_tangent: float, tip_loss_correction: float)\
            -> float:
        return 1/((4*tip_loss_correction*np.sin(phi)*np.cos(phi))/(local_solidity*c_tangent)-1)
