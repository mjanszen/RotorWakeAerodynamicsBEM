import numpy as np
from scipy.optimize import brentq, newton
import pandas as pd
import scipy.interpolate as interpolate
from helper_functions import Helper
import scipy
helper=Helper()


class BEM:
    def __init__(self,
                 root: str,
                 file_airfoil: str,
                 save_dir: str):
        self.rotor_radius = None
        self.root_radius = None
        self.n_blades = None
        self.air_density = None
        self.dir_save = save_dir
        self.a_prime = 0
        helper.create_dir(self.dir_save)

        df_tmp = pd.read_excel(root+"/"+file_airfoil, skiprows=3)
        self.interp = {"c_l": interpolate.interp1d(df_tmp["alpha"], df_tmp["cl"]),
                       "c_d": interpolate.interp1d(df_tmp["alpha"], df_tmp["cd"])}
        self.implemented_glauert_correction = ["none", "dtu", "tud"]
        self.implemented_tip_correction = ["none", "dtu", "tud"]
        self.implemented_root_correction = ["none", "tud"]

    def set_constants(self,
                      rotor_radius: float,
                      root_radius: float,
                      n_blades: int,
                      air_density: float) -> None:
        self._set(**{param: value for param, value in locals().items() if param != "self"})
        return None

    def solve(self,
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
        Glauert_correction: either 'tud' (TU Delft) or 'dtu' (Denmark's TU). Same for blade_end_correction
        All angles must be in rad.
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch:
        :param resolution:
        :return:
        """
        start_radius = start_radius if start_radius is not None else self.root_radius
        results = {
            "positions": list(),
            "a": list(),
            "a_prime": list(),
            "f_n": list(),
            "f_t": list(),
            "bec": list(), # blade end correction
            "v0": wind_speed
        }
        omega = tip_speed_ratio*wind_speed/self.rotor_radius
        for r in np.linspace(start_radius, self.rotor_radius, resolution):
            self.a_prime = 0
            chord = self._get_chord(r, self.rotor_radius)
            twist = self._get_twist(r, self.rotor_radius)
            local_solidity = self._local_solidity(chord, r, self.n_blades)

            def residue(a):
                phi = self._phi(a=a, a_prime=self.a_prime, wind_speed=wind_speed, rotational_speed=omega, radius=r)
                aero_values = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r, a=a,
                                                       tip_seed_ratio=tip_speed_ratio,
                                                       blade_end_correction_type=blade_end_correction_type,
                                                       root=root_loss_correction, tip=tip_loss_correction)
                c_n, c_t, blade_end_correction = aero_values[3], aero_values[4], aero_values[5]
                self._update_a_prime(local_solidity=local_solidity, c_tangential=c_t, inflow_angle=phi,
                                     blade_end_correction=blade_end_correction)
                return self._equate_blade_element_and_momentum(glauert_correction=glauert_correction_type, a=a,
                                                               blade_end_correction=blade_end_correction, phi=phi,
                                                               local_solidity=local_solidity, c_normal=c_n)

            try:
                a = brentq(residue, *brent_bracket)
            except ValueError:
                print("Brent could not be used for the convergence, using Newton instead.")
                a =  newton(residue, 1/3)
            phi = np.arctan((1-a)*wind_speed/((1+self.a_prime)*omega*r))
            aero_values = self._phi_to_aero_values(phi=phi, twist=twist, pitch=pitch, radius=r, a=a,
                                                   tip_seed_ratio=tip_speed_ratio,
                                                   blade_end_correction_type=blade_end_correction_type,
                                                   root=root_loss_correction, tip=tip_loss_correction)
            c_n, c_t, blade_end_correction = aero_values[3], aero_values[4], aero_values[5]
            inflow_velocity = np.sqrt((omega*r*(1+self.a_prime))**2+(wind_speed*(1-a))**2)
            if self.a_prime == "-inf":
                c_n = 0
                c_t = 0
            results["positions"].append(r)
            results["a"].append(a)
            results["a_prime"].append(self.a_prime)
            results["f_n"].append(self._f_n(inflow_velocity=inflow_velocity, chord=chord, c_normal=c_n))
            results["f_t"].append(self._f_t(inflow_velocity=inflow_velocity, chord=chord, c_tangential=c_t))
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

    def _phi_to_aero_values(self, phi: float, twist:float or np.ndarray, pitch: float, radius: float,
                            tip_seed_ratio: float, a: float, blade_end_correction_type: str, tip: bool, root: bool) -> \
            tuple:
        alpha = np.rad2deg(phi-(twist+pitch))
        c_l = self.interp["c_l"](alpha)
        c_d = self.interp["c_d"](alpha)
        c_n = self._c_normal(phi, c_l, c_d)
        c_t = self._c_tangent(phi, c_l, c_d)
        return alpha, c_l, c_d, c_n, c_t, self._blade_end_correction(which=blade_end_correction_type, tip=tip, root=root,
                                                                     phi=phi, radius=radius, a=a,
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

    def _f_n(self, inflow_velocity: float, chord: float, c_normal: float):
        """
        Calculates the normal force per unit span.
        :param inflow_velocity:
        :param chord:
        :param c_normal:
        :return:
        """
        return 1/2*self.air_density*inflow_velocity**2*chord*c_normal

    def _f_t(self, inflow_velocity: float, chord: float, c_tangential: float):
        """
        Calculates the tangential force per unit span.
        :param inflow_velocity:
        :param chord:
        :param c_normal:
        :return:
        """
        return 1/2*self.air_density*inflow_velocity**2*chord*c_tangential

    def _equate_blade_element_and_momentum(self, glauert_correction: str, a: float, blade_end_correction: float,
                                           phi: float, local_solidity: float, c_normal: float):
        if a < 1/3 or glauert_correction == "none":
            return 1/((4*blade_end_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-a
        else:
            if glauert_correction == "dtu":
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-4*a*blade_end_correction*(1-a/4*(5-3*a))
            elif glauert_correction == "tud":
                CT_1 = 1.816
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-(CT_1-4*(np.sqrt(CT_1)-1)*(1-a))
            else:
                raise ValueError(f"Parameter 'glauert_correction' must be one of {self.implemented_glauert_correction}, "
                                 f"but was {glauert_correction}.")

    def _blade_end_correction(self, which: str, tip: bool, root: bool, phi: float, radius: float,
                              tip_seed_ratio: float, a: float) -> float:
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
    def _tangential_induction_factor(phi: float, local_solidity: float, c_tangent: float, tip_loss_correction: float)\
            -> float:
        return 1/((4*tip_loss_correction*np.sin(phi)*np.cos(phi))/(local_solidity*c_tangent)-1)

    @staticmethod
    def _raise_loop_break_error(control_type: str, counter: int, step_size: float):
        raise ValueError(f"Could not pitch to rated power with {counter*step_size}°change {control_type}. The number "
                         f"of pitch increments (loop_breaker) might have been too low.")
