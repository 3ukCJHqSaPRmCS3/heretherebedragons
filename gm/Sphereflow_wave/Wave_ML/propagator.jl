#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: propagators
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file defines all the propagators (low fidelity: lagrangian) and 
# (high fidelity: taylor) that are used to run the low-thrust trajectories using the 
# Sims-Flanagan method.
#  
# ----------------------------------------------------------------------------------------

using LinearAlgebra
using DifferentialEquations

# propagate_2BPt: High fidelity propagator
# r: start position, x,y,z. [km]
# v: start velocity, vx,vy,vz. [km/s]
# m0: starting mass. [kg]
# T: fixed inertial thrust, ux,uy,uz. [N]
# tof: propagation time. [s]
# mu: central body gravity constant. [km^3/s^2]
# veff: the product (Isp g0) defining the engine efficiency. 
# reltol: the relative tolerance passed to taylor propagator.
# log10rtol: the absolute tolerance passed to taylor propagator.
# abstol a tuple (rf, vf, mf) containing the final position, velocity and mass after the propagation.
function propagate_2BPt(x0, tof, T, params, rel_tol=1e-8, abs_tol=1e-8)

    #Find position from state at t0
    r0 = x0[1:3]
    # Find velocity from state at t0
    v0 = x0[4:6]
    # Find mass from state at t0
    m0 = x0[7]
    # Compute the engine efficiency
    veff = params.Isp * params.g0

    # Define the initial guess and the integration time span
    init_guess = vcat(r0, v0, m0)
    t_span = (0.0, tof)

    # Define the additonal arguments for the eoms function
    p = [norm(T), T, veff, params.mu]
    # Define the ODE Problem
    prob = ODEProblem(eoms_2BP!, init_guess, t_span, p)
    # Solve the ODE Problem using Tsitouras 5/4 Runge-Kutta method and user-defined tolerances
    sol = solve(prob, Tsit5(), reltol=rel_tol, abstol=abs_tol)

    return(sol[:, end])
     
end

# eoms_2BP: equations of motion for the 2 body Problem with continous thrust
# u: start position, x,y,z. [km] and start velocity, vx,vy,vz. [km/s] and mass [kg]
# m0: starting mass. [kg]
# thrust: fixed inertial thrust, ux,uy,uz. [N]
# mu: central body gravity constant. [km^3/s^2]
function eoms_2BP!(du, u, p, t)

    R0 = norm(u[1:3])
     
    F, T, veff, mu = p
    
    # Compute vx, vy, vz
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # Compute m_dot
    du[7] = -F/veff
    # Compute ax, ay, az
    du[4] = -mu/R0^3 * u[1] + T[1]/(u[7] - norm(du[7])*t)
    du[5] = -mu/R0^3 * u[2] + T[2]/(u[7] - norm(du[7])*t)
    du[6] = -mu/R0^3 * u[3] + T[3]/(u[7] - norm(du[7])*t)

end