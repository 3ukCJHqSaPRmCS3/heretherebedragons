"""Define the orbit Class
"""

# ======= third party imports ======
from cmath import pi, sqrt
from math import *
from operator import length_hint 
import numpy as np
# ==================================

# ========= program imports ========
# ==================================
# ======== Sphereflow-Wave  ========
from swave.utils.constants import R_EARTH, MU_EARTH, DEG2RAD

class Orbit():
    """ Orbit object characterized using keplerian elements:
    1) Semi-major axis [km]
    2) eccentricity [/]
    3) inclination  [deg]
    4) Right Ascension of the Ascending node [deg]
    5) Argument of perigee [deg]
    6) Central body gravity constant. [km^3/s^2]
    Default constructor gives [sma = 6871 lm, ecc = 0, inc = 45 deg, RAAN = aop = 0 deg, mu =  ]
    """

    def __init__(
        self,
        sma = R_EARTH + 500,
        ecc = 0,
        inc = 45,
        RAAN = 0,
        aop = 0,
        mu = MU_EARTH
    ):
        self.sma = sma
        self.ecc = ecc
        self.inc = inc
        self.raan = RAAN
        self.aop = aop
        self.mu = mu

    def get_kep_elements(self):
        """ Method to visualize the class attributes
        """
        if self.mu == MU_EARTH:
            central_planet = "Earth"

        print("\nsma [km]:", self.sma, "\n ecc", self.ecc, "\ninc [deg]", self.inc, "\nRAAN [deg]", self.raan, "\naop [deg]", self.aop,
        "\n central planet:", central_planet)

    def set_kep_elements(self, kep_elements):
        """ Method to set the attributes of the object Orbit with an input vector:
        kep_elements = [sma[km], ecc, inc[deg], RAAN[deg], aop[deg]]
        """
        print(np.size(kep_elements))
        print(type(kep_elements))
        self.sma, self.ecc, self.inc, self.raan, self.aop = kep_elements

    def set_mu(self, mu):
        self.mu = mu

    def kep2cart(self):
        """ Method to find the cartesian state from an orbit espressed in keplerian elements
        """
       
        rraan  = self.raan * DEG2RAD
        raop  = self.aop * DEG2RAD
        rinc  = self.inc * DEG2RAD

        # Define the rotation matrices R(phi) = R(-RAAN)R(-i)R(-aop)
        # R_raan = np.array([[cos(rraan), -sin(rraan), 0], [sin(rraan), cos(rraan), 0], [0, 0, 1]])
        # R_inc = np.array([ [1, 0, 0], [0, cos(rinc), -sin(rinc)], [0, sin(rinc), cos(rinc)]])
        # R_aop = np.array([[cos(raop), - sin(raop), 0], [sin(raop), cos(raop), 0],[0, 0, 1]])
        # R = R_raan.dot(R_inc).dot(R_aop)
        R = np.array([[cos(rraan)*cos(raop)-sin(rraan)*sin(raop)*cos(rinc), -cos(rraan)*sin(raop)-sin(rraan)*cos(raop)*cos(rinc), sin(rraan)*sin(rinc)],
        [sin(rraan)*cos(raop)+cos(rraan)*sin(raop)*cos(rinc), -sin(rraan)*sin(raop)+cos(rraan)*cos(raop)*cos(rinc), -cos(rraan)*sin(rinc)],
        [sin(raop)*sin(rinc), cos(raop)*sin(rinc), cos(rinc)]])
        # Set a vector of true anomaly starting from 1 deg
        theta = np.linspace(0, 2*pi, 100)

        # Define the ouptut vrctors 
        r_eci = np.empty((np.size(theta), 3))
        v_eci = np.empty((np.size(theta), 3))

        for i  in range(np.size(theta)):

            # Compute the eccentric anomaly
            E = 2 * atan2(sqrt((1 - self.ecc) / (1 + self.ecc)) * sin(theta[i]/ 2), cos(theta[i] / 2))

            # Compute the radius magnitude
            r = self.sma * (1 - self.ecc * cos(E))

            #Compute the position and velocity vector in the orbital frame
            r_of = r * np.array([cos(theta[i]), sin(theta[i]), 0])
            v_of = sqrt(self.mu * self.sma)/r * np.array([-sin(E), sqrt(1 - self.ecc**2) * cos(E), 0])
            # print("\n r is ", r_of)
            r_xyz = R.dot(r_of)
            v_xyz = R.dot(v_of)

            # fill the output vector
            r_eci[i,:] = r_xyz
            v_eci[i,:] = v_xyz


        return r_eci, v_eci










