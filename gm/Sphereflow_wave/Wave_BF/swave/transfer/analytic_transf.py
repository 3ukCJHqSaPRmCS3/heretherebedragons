"""The function in this file computes the vector and magnitude of the maneuvers required to go from a constellation A to a constellation B analytically.
This method is brute force: taking every orbit in constellation A, it computes the best point along the trajectory to perform a maneuver, then it
computes a trasnfer to every orbit in constellation B. At the end of this process, the matrix of Dv required to go from avery orbit in A to every
orbit in B will be analyzed to pick the best transfer for each orbit in terms of fuel, timing and avoiding overlapping transfers.
"""

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
from swave.transfer.lambert import Lambert
from swave.utils.orbit import Orbit
from swave.utils.plot import Plot
from swave.utils.propagators import Propagator


def analytic_tranfer(orbit1, orbit2):
    """ Function to compute impulsive transfer between orbit 1 and orbit 2

    Input: 
        - orbit 1: Obejct of the class Orbit
        - orbit 2: Object of the class Orbit
    Output:
        - Cartesian states of the propagated transfer
        - dV corresponding to the min dV possible for the trasfer
    """

    # Compute the tof as the hohmann transfer tof
    a_h = (orbit1.sma + orbit2.sma)*0.5
    t_h = np.pi*np.sqrt(a_h**3/orbit1.mu)
    tof = 2*t_h*np.array([1/3, 1/2, 3/4, 1, 5/4, 6/4, 7/4, 8/4])

    # Compute Transfer
    l = Lambert(orbit1, orbit2, tof)
    l.set_rev = 1
    x0, xf, dv_min, t_idx = l.compute_transfer()
        # print("The minimum dv for this time of flight is: ", dv_min, " km/s in tof: ", tof[t_idx], " s")

    # Propagate
    prop =  Propagator(x0, tof[t_idx])
    states = prop.propagate_2BP()

    return states, dv_min
    



def main():

    n_orbits = 3
    sma1, ecc1, inc1, aop1 = [10000.0, 0.3, 45.0, 0.0]
    raan1 = np.linspace(80,300,n_orbits)

    sma2, ecc2, inc2, aop2 = [8000.0, 0.2, 50.0, 30.0]
    raan2 = np.linspace(50,330,n_orbits)

    orbit1 = list()
    orbit2 = list()
    for i in range(np.size(raan1)):
        orbit1.append(Orbit(sma1, ecc1, inc1, raan1[i], aop1))
        orbit2.append(Orbit(sma2, ecc2, inc2, raan2[i], aop2))

    dV_ij = np.empty((n_orbits,n_orbits))
    states_ij = np.empty((n_orbits,n_orbits, 50, 6))
    
    idx_min = list()
    for i in range(n_orbits):
        for j in range(n_orbits):

            states, dv = analytic_tranfer(orbit1[i], orbit2[j])

            dV_ij[i,j] = dv
            states_ij[i,j,:, :] = states
            # states_ij.append(states)
        
        idx_min.append(np.argmin(dV_ij[i]))
        
    dv_min = np.amin(dV_ij, axis=1)

    print("The minimum dv for this time of flight is: ", dv_min, " km/s")

    # Check that every elment in the list in unique
    if len(idx_min) > len(set(idx_min)):
        print("Not unique/independent solution") # do something
    
    # Plot
    p = Plot()
    p.plot_constellation(n_orbits, idx_min, states_ij, orbit1, orbit2)



if __name__ == "__main__":
    main()








    





    


