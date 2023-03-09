import numpy as np
from scipy.optimize import brentq as root, newton
import pandas as pd
import scipy.interpolate as interpolate
from helper_functions import Helper
helper=Helper()


class BEM:
    def __init__(self,
                 root: str,
                 file_airfoil: str,
                 save_dir: str):
        self.rotor_radius = None
        self.n_blades = None
        self.air_density = None
        self.dir_save = save_dir
        helper.create_dir(self.dir_save)

        df_tmp = pd.read_excel(root+"/"+file_airfoil, skiprows=3)
        self.interp = {"c_l": interpolate.interp1d(df_tmp["alpha"], df_tmp["c_l"]),
                       "c_d": interpolate.interp1d(df_tmp["alpha"], df_tmp["c_l"])}

    def set_constants(self,
                      rotor_radius: float,
                      n_blades: int,
                      air_density: float) -> None:
        self._set(**{param: value for param, value in locals().items() if param != "self"})
        return None

    def solve(self,
              wind_speed: float,
              tip_speed_ratio: float,
              pitch: float,
              resolution: int=30):
        """
        All angles must be in rad
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch:
        :param resolution:
        :return:
        """
        results = {
            "r": np.linspace(0.2*self.rotor_radius, self.rotor_radius, resolution),
            "a": list(),
            "a_prime": list(),
            "f_n": list(),
            "f_t": list()
        }
        for r in np.linspace(0.2*self.rotor_radius, self.rotor_radius, resolution):
            twist = self._get_twist(r, self.rotor_radius)
            chord = self._get_chord(r, self.rotor_radius)



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
    def _tip_loss_correction(r: float, phi: float, rotor_radius: float, n_blades: int) -> float:
        """
        Returns the factor F for the tip loss correction according to Prandtl
        :param r: current radius
        :param phi: angle between flow and rotational direction in rad
        :param rotor_radius: total radius of the rotor
        :param n_blades: number of blades
        :return: Prandtl tip loss correction
        """
        if np.sin(np.abs(phi)) < 0.01:
            return 1
        return 2/np.pi*np.arccos(np.exp(-(n_blades*(rotor_radius-r))/(2*r*np.sin(np.abs(phi)))))



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
    def _root_axial_induction_factor(phi: float,
                                     local_solidity: float,
                                     c_normal: float,
                                     tip_loss_correction: float,
                                     bracket: tuple=(-1,1)) -> float:
        def residue(aif):
            if aif <= 1/3:
                return 1/((4*tip_loss_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-aif
            else:
                return local_solidity*((1-aif)/np.sin(phi))**2*c_normal-4*aif*tip_loss_correction*(1-aif/4*(5-3*aif))
        try:
            return root(residue, *bracket)
        except ValueError:
            print("Brent could not be used, using Newton.")
            return newton(residue, 1/3)

    @staticmethod
    def _tangential_induction_factor(phi: float, local_solidity: float, c_tangent: float, tip_loss_correction: float)\
            -> float:
        return 1/((4*tip_loss_correction*np.sin(phi)*np.cos(phi))/(local_solidity*c_tangent)-1)

    @staticmethod
    def _raise_loop_break_error(control_type: str, counter: int, step_size: float):
        raise ValueError(f"Could not pitch to rated power with {counter*step_size}°change {control_type}. The number "
                         f"of pitch increments (loop_breaker) might have been too low.")
