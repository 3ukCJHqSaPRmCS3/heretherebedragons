""" Compute transfer solving the lambert problem 
"""

# ======== standard imports ========
import pickle
# ==================================

# ======= third party imports ======
import numpy as np
from lamberthub import izzo2015
# ==================================

# ========= program imports ========
# ==================================
# == Sphereflow-WaveConstellation ==
from swave.utils.orbit import Orbit


class Lambert():
    """ Defining Lambert problem class
    """
    def __init__(self, orbit1: Orbit, orbit2: Orbit, tof, n_rev = 0):
        
        self.orbit1 = orbit1
        self.orbit2 = orbit2
        self.tof = tof
        self.n_rev = n_rev

    def set_rev(self, n):
        """ Setter for number of revolution to take in to account for transfer solution

            Input:
                - n: number of revolutions :Int
        """
        self.n_rev = n

    def compute_transfer(self):
        """ Method to compute all the possible transfer between orbit 1 and orbit 2
        """

        # Compute the cartesian position and velocity vectors
        r1, v1 = self.orbit1.kep2cart()
        r2, v2 = self.orbit2.kep2cart()
        dv_min = 10000
        i_min = 0
        j_min = 0
        l_min = 0

        # Iterate betreen every point in orbit 1 and evry point of orbit 2 to find the transfer with lowest dV
        for i in range(np.size(r1, axis = 0)):
            for j in range(np.size(r2, axis = 0)):
                for l in range(np.size(self.tof)):
                    
                    # COmpute the solution to the Lambert problem
                    v1_ij, v2_ij = izzo2015(self.orbit1.mu, r1[i,:], r2[j,:], self.tof[l], M = self.n_rev, prograde=True, low_path=True, maxiter=50, atol=1e-7, rtol=1e-7, full_output=False)
                    #add the computed dV to the list
                    dv = np.linalg.norm(np.array([v1[i,:] - v1_ij, v2_ij - v2[j,:]]))
                    # Find the minimum dv transfer index
                    if dv < dv_min:
                        dv_min = dv
                        i_min = i
                        j_min = j
                        l_min = l

        
        # Propagate the minimum dv trajectory
        v1_min, v2_min = izzo2015(self.orbit1.mu, r1[i_min,:], r2[j_min,:], self.tof[l_min], M = self.n_rev,
                prograde=True, low_path=True, maxiter=35, atol=1e-7, rtol=1e-7, full_output=False)
        
        # Assemble te output as a 1x6 vector in the form sol = [x,y,z,vx,vy,vz]
        x0 = np.concatenate((r1[i_min, :], v1_min), axis=0)
        xf = np.concatenate((r2[j_min, :], v2_min), axis=0)

        return x0, xf, dv_min, l_min
                

                

