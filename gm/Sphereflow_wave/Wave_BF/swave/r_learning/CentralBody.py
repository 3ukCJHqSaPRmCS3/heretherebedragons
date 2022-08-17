# ======== standard imports ========
import pickle
# ==================================

# ======= third party imports ======
import numpy as np


# ==================================

# ========= program imports ========
# ==================================
# == Sphereflow-WaveConstellation ==
# Add the path of the geenral folder
from swave.utils.constants import *


class CentralBody(object):
    def __init__(self, name):
        self.name = name
        self.mu
        
        match self.name:
            case "Earth":
                self.mu = MU_EARTH
                self.rad = 6371.0 #[km]
            
            case "Mars":
                self.mu = 428283.7 #[km^3/s^2]
                self.rad = 3389.5 #[km]

            case _:
                return "This central body has not been defined yet or is not present"
                