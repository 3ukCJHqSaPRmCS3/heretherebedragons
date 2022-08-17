""" Definition of the Plot class its methods
"""

# ======== standard imports ========
import pickle
# ==================================

# ======= third party imports ======
import matplotlib.pyplot as plt
import numpy as np
# ==================================

# ========= program imports ========
# ==================================
# ======== Sphereflow-Wave  ========
from swave.utils.orbit import Orbit


class Plot():

    def plot_transfer(self, transfer, orbit1: Orbit, orbit2: Orbit):
        """ Function to plot single transfer

        Inputs:
            - transfer: matrix containing the states from the Propagation
            - orbit1: departure orbit as object of the Orbit class (define it in keplerian elements)
            - orbit2: arrival orbit as object of the Orbit class (define it in keplerian elements)
        Output:
            - 3D plot with departure and arrival orbit and the transfer to get from the first to the second
        """
        # Set figure and axes
        fig = plt.figure()
        ax = plt.axes(projection = '3d')

        # Plot the central body
        ax.scatter(0,0,0)

        #Plot the transfer
        ax.plot3D(transfer[:,0], transfer[:,1], transfer[:,2], 'red', label='transfer tranjectory')
        # COnvert the orbit from keplerian elements to cartesian states
        r_orbit1, v_orbit1 = orbit1.kep2cart()
        r_orbit2, v_orbit1 = orbit2.kep2cart()
        #Plot the departure and departure orbit
        ax.plot3D(r_orbit1[:,0], r_orbit1[:,1], r_orbit1[:,2], 'black', label='departing orbit')
        ax.plot3D(r_orbit2[:,0], r_orbit2[:,1], r_orbit2[:,2], 'black', label = 'arriving orbit')
        # Set labels and title
        ax.set_xlabel("x [km]")
        ax.set_ylabel("y [km]")
        ax.set_zlabel("z [km]")
        ax.set_title("Transfer between orbit 1 and orbit 2")
        ax.legend()

        plt.show()
        

    
    def plot_constellation(self, n_orbits_constA, idx_min_constB, transfers, constA, constB):
        """ Function to plot constellation transfers:

        Inputs: 
            - n_orbits_constA: number of orbits in departure constellation 
            - idx_min_constB: list collecting indexes for the min dV in orbit_a[1, .., n] -> orbit_b[index]
            - transfers: Matrix with the propagated states for the min dV transfers
            - constA: List with orbits belonging to departure constellation
            - constB: List with orbits belonging to arrival constellation
        Output:
            - 3D plot with departure and arrival orbits and the transfers to get from the every orbit in A to every orbit in B
        """
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    
        # Plot the orbits from the constellations
        for i in range(n_orbits_constA):
         
            #Load the transfer trajectory
            t = transfers[i,idx_min_constB[i]]
            #Load orbit from const A and orbit from const B
            o1 = constA[i]
            o2 = constB[idx_min_constB[i]]

            #Convert the orbits in cartesian elements
            r_orbit1, v_orbit1 = o1.kep2cart()
            r_orbit2, v_orbit1 = o2.kep2cart()

            #plot orbit1, orbit2 and the transfer from the first to the second 
            ax.plot(t[:,0], t[:,1], t[:,2], 'red', linewidth = 1)
            ax.plot(r_orbit1[:,0], r_orbit1[:,1], r_orbit1[:,2], 'black', linewidth = 1)
            ax.plot(r_orbit2[:,0], r_orbit2[:,1], r_orbit2[:,2], 'blue', linewidth = 1)
        
        # Plot the central Body
        ax.scatter(0,0,0)

        # Set axis labels and title
        ax.set_xlabel("x [km]")
        ax.set_ylabel("y [km]")
        ax.set_zlabel("z [km]")
        ax.set_title("Transfer between Constellation A to Constellation B")
        ax.legend(['Transfer trajectory', 'Departing orbit', 'Arrival orbit'])

        plt.show()

