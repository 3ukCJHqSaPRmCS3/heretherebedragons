""" Definition of the propagator class its methods
"""

# ======== standard imports ========
import pickle
# ==================================

# ======= third party imports ======
import numpy as np
from scipy.integrate import odeint
# ==================================

# ========= program imports ========
# ==================================
# ======== Sphereflow-Wave  ========
from swave.utils.orbit import Orbit
from swave.utils.eoms import *
from swave.utils.constants import *

class Propagator():
    """ Propagate the equation of motion using the scipy odeint library

    """
    def __init__(self, initial_state, int_time, mu=MU_EARTH, T = None):
        self.initial_state = initial_state
        self.int_time = int_time
        self.mu = mu
        self.T = T
        self.atol = 1e-8
        self.rtol = 1e-8

    def set_tol(self, atol, rtol):
        """Setter for tolerances

            Inputs:
                - atol: absolute tolerance value
                - rtol: relative tolerance value
        """
        self.atol = atol
        self.rtol = rtol



    def propagate_2BP(self):
        """ Function to propagate the 6 states [x,y,z,vx,vy,vz] using the 2-body problem equations of motion

            Output:
                - Matrix with 6 states every t time instants
        """

        # Find the time span
        t = np.linspace(0, self.int_time)

        # Define the initial guess
        y0 = self.initial_state

        # Integrate the system of equation of motions
        sol = odeint(eoms_2BP, y0, t, args = (self.mu,), rtol = 1e-8, atol = 1e-8, hmax = 50.0)

        return sol

        
    def propagate_2BPT(self):
        """ Function to propagate the 6 states [x,y,z,vx,vy,vz] using the 2-body problem equations of motion with Thrust vector

            Output:
                - Matrix with 6 states every t time instants
        """

        # Find the time span
        t = np.linspace(0, self.int_time)

        # Define the initial state
        y0 = self.initial_state

        # Integrate the system of equation of motions
        sol = odeint(eoms_2BPT, y0, t, args = (self.mu, self.T), rtol = 1e-8, atol = 1e-8, hmax = 50.0)

        return sol

        
