""" OpenAI environment for Sphereflow Wave. Definition of observation space and elements.
"""
import numpy as np 
import matplotlib.pyplot as plt
from swave.utils.orbit import *
import gym
import random

from gym import Env, spaces
import CentralBody as cb
from swave.utils.spacecraft import *
from swave.constellations.constellation import *

class WaveTransfer(Env):

    def __init__(self):
        super(WaveTransfer, self).__init__()


        # Define the 3D obervation space
        self.observation_size = (1000, 1000, 1000)
        
        self.observation_space = spaces.Box(low = -self.observation_size[0], high = self.observation_size[2],
                                            shape = (self.observation_size[0], self.observation_size[1]), dtype= np.float32)


        # Define the action space ranging from 0 to 1 ##DISCRETE for NOW
        self.action_space = spaces.Discrete(4)

        # Define max s/c fuel 
        self.max_fuel = 1000

        # Define the elements inside the environment
        self.n_sc = 5
        self.sc = Spacecraft()
        self.sc.set_isp()
        self.sc.set_thrust()
        self.orbits = []
        self.central_body = cb("Earth")

        self.constellation = Constellation()
    
    def reset(self):

        # Reset the fuel consumed
        self.sc.m_wet = self.max_fuel

        # Reset the spacecraft inisital state
        self.sc.set_position()

        # Reset the rewards
        self.distance_obs = 100
        self.dv_consumption = 0
        self.time_of_flight = 0


        # Determine the place to initialise the spacecraft in
        theta = np.random.rand(2*pi)

        # self.sc.set_position = 

 