""" Definition of the equaiton of motion for the orbit's states propagation
"""

# ======== standard imports ========
from importlib.machinery import DEBUG_BYTECODE_SUFFIXES
import pickle
# ==================================

# ======= third party imports ======
import numpy as np
# ==================================


def eoms_2BP(y, t, mu):
    """ Equation of motion for the two body problem represting keplerian dynamics
    y : state vector [x, y, z, vx, vy, vz]. Position elements [km], velocity elements [km/s]
    t : time variable.
    mu: central body gravity constant. [km^3/s^2]
    """
    y = np.array(y)
    vel = y[3:6]
    acc = -mu/(np.linalg.norm(y[0:3])**3)*y[0:3] 
    dydt = vel.tolist() + acc.tolist()

    return dydt


def eoms_2BPT(y, t, mu, T):
    """ Equation of motion for the two body problem represting keplerian dynamics with Trust vector
    y : state vector (x, y, z, vx, vy, vz, m). position elemwnts [km], velocity elements [km/s], mass [kg]
    t : time variable.
    mu: central body gravity constant. [km^3/s^2]
    T : fixed inertial thrust (ux,uy,uz). [N]
    """
    y = np.array(y)
    T = np.array(T)
    vel = y[3:6]
    acc = -mu/np.linalg.norm(y[0:3])*y[0:3]  + T/y[-1]
    dydt = vel.tolist() + acc.tolist()

    return dydt

    