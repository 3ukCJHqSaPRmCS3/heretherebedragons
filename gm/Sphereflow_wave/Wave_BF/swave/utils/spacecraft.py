""" Define the spacecraft Class
"""

# ======= third party imports ======
import numpy as np
# ==================================


# ========= program imports ========
# ==================================
# == Sphereflow-WaveConstellation ==
# Add the path of the geenral folder
from swave.utils.propagators import *
from swave.utils.orbit import *


class Spacecraft():

    def __init__(self, m_wet, isp, thrust, state_vec):
        self.m_wet = m_wet
        self.isp = isp
        self. thrust= thrust
        self.state_vec = np.empty((1,6))

    def set_position(self, state_vec):
        self.state_vec[:] = np.array(state_vec)

    def set_isp(self, isp):
        self.isp = isp

    def set_thrust(self, thrust):
        self.thrust = thrust

    def get_position(self):
        return (self.state_vec)

    def step(self, tof, T_vec = None):

        #Set initial position
        x0 = self.state_vec

        # Propagate
        if T_vec is not None:
            prop =  Propagator(x0, tof, T_vec)
            states = prop.propagate_2BPT()
        elif T_vec is None:
            prop =  Propagator(x0, tof)
            states = prop.propagate_2BPT()
        
        self.state_vec[:] = states[-1, :]


    

