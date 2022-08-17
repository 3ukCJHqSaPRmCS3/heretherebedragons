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

include("utils.jl")
using LinearAlgebra
using DifferentialEquations

# propagate_lagrangian: low-fidelity propagator
# r: start position, x,y,z [km]
# v: start velocity, vx,vy,vz [km/s]
# tof: propagation time from a point A along the trajectory [s]
# mu: central body gravity constant [km^3/s^2]
# Returns a tuple (rf, vf) containing the final position and velocity after the propagation.
function propagate_lagrangian(r0::Vector{Float64},v0::Vector{Float64}, tof::Float64, mu::Float64 )

     # Find the distance, speed, energy, semi major axis, eccentricity, specific angular momentum, true anomaly and mean motion
     R0 = norm(r0) #[km]
     V0 = norm(v0) #[km/s]
     ϵ = (V0^2 / 2.0 - mu / R0) # [km^2/s^2]
     a = -mu / (2.0 * ϵ) #[km]
     e = ((V0^2 - mu / R0) * r0 - (dot(r0,v0)) * v0) / mu
     h = sqrt(mu * a * (1 - norm(e)^2)) #[km^2/s]
     θ_0 = acos(dot(r0,e)/(R0 * norm(e))) #[rad]
     # check the sign of θ
     if dot(r0,v0) < 0 
          θ_0 -= 2*pi
     end

     # Define the tolerance and maximum iterations
     tol = 1e-16
     max_iter = 50

     # Solve the Kepler's equation for elliptical orbits
     if a > 0 
          # Mean motion 
          n= sqrt(mu / a^3)

          # Compute Eccentric Anomaly at t0
          E0 = 2 * atan(sqrt((1 - norm(e)) / (1 + norm(e))) * sin(θ_0 / 2), cos(θ_0 / 2))
          
          # Find time of flight from the pericenter to the position r0
          t0 = (E0 - norm(e) * sin(E0)) / n

          # Compute Mean Anomaly at t0
          M = n * (t0 + tof)
          E_guess = M

          # Find the eccentric anomaly at the end of the tof
          E = solve_keplers_eq_E(E_guess, M, e, tol, max_iter)

          # Compute the new true Anomaly
          θ = 2 * atan(sqrt((1 + norm(e)) / (1 - norm(e))) * sin(E / 2), cos(E / 2))


     # Solve the Kepler's equation for hyperbolic orbits
     elseif a < 0
          # Mean motion 
          n= sqrt(mu / abs(a)^3)

          # Compute hyperbolic Anomaly at t0
          H0 = 2 * atanh(sqrt((1- norm(e)) / (1 + norm(e))) * tan(θ_0 / 2))

          # Find time of flight from the pericenter to the position r0
          t0 = (norm(e) * sin(H0) - H0) / n

          # Compute Mean Anomaly at t0
          M = n * (t0 + tof)
          H_guess = 1

          # Find the eccentric anomaly at the end of the tof
          H = solve_keplers_eq_H(H_guess, M, e, tol, max_iter)

          # Compute the new true Anomaly
          θ = 2 * atan(sqrt((1 + norm(e)) / (1 - norm(e))) * sinh(H / 2), cosh(H / 2))
#           # Compute the change in eccentric anomaly
#           ΔH = H - H0
#           # Compute the position at the new H
#           R = a + (R0 - a) * cos(H) + σ0 * sqrt(-a) * sin(H)

#           # Compute the Lagrange coefficients in terms of the change in hyperbolic anomaly
#           f = 1 - a/R0*(1 - cosh(ΔH))
#           g = a * σ0 / sqrt(mu) * (1 - cosh(ΔH)) + R0 * sqrt(- a / mu) * sinh(ΔH) 
#           ### CHECK with classic formulation
# #         g = tof - sqrt(-a^3/mu)*(ΔH - sinh(ΔH))
#           f_dot = -sqrt(mu * a) / (R0 * R) * sinh(ΔH)
#           g_dot = 1 - a / R * (1 - cosh(ΔH))

     end

     # Find the radius magnitude using the conic equation
     R = h^2 / (mu * (1 + norm(e) * cos(θ)))
     # Compute the change in true anomaly
     Δθ = θ - θ_0
     # Find the lagrange coefficients in terms of the change of true Anomaly
     f = 1 - mu*R/h^2 * (1 - cos(Δθ))
     g = R*R0/h * sin(Δθ)
     f_dot = mu/h * (1 - cos(Δθ))/sin(Δθ) * (mu/h^2 * (1 - cos(Δθ)) - 1/R0 - 1/R )
     g_dot = 1 - mu*R0/h^2 *(1 - cos(Δθ))

     # Compute the new radius and velocity
     rf = f * r0 + g * v0 #[km]
     vf = f_dot * r0 + g_dot * v0 #[km/s]
    return(rf, vf, θ)
end

# propagate_2BPt: High fidelity propagator
# r: start position, x,y,z. [km]
# v: start velocity, vx,vy,vz. [km/s]
# m0: starting mass. [kg]
# thrust: fixed inertial thrust, ux,uy,uz. [N]
# tof: propagation time. [s]
# mu: central body gravity constant. [km^3/s^2]
# veff: the product (Isp g0) defining the engine efficiency. 
# reltol: the relative tolerance passed to taylor propagator.
# log10rtol: the absolute tolerance passed to taylor propagator.
# abstol a tuple (rf, vf, mf) containing the final position, velocity and mass after the propagation.
function propagate_2BPt(l::Leg, rel_tol, abs_tol)

     #Find position from state at t0
     r0 = l.x0[1:3]
     # Find velocity from state at t0
     v0 = l.x0[4:6]
     # Compute the time of flight in seconds
     tof = (l.end_time - l.start_time)*DAY2SEC
     # Compute the engine efficiency
     veff = l.spacecraft.isp * g0

     # Define the initial guess and the integration time span
     init_guess = vcat(r0, v0, m0)
     t_span = (0.0, tof)

     # Define the additonal arguments for the eoms function
     p = [norm(l.throttles), l.throttles, veff, l.mu]
     # Define the ODE Problem
     prob = ODEProblem(eoms_2BP!, init_guess, t_span, p)
     # Solve the ODE Problem using Tsitouras 5/4 Runge-Kutta method and user-defined tolerances
     sol = solve(prob, Tsit5(), reltol=rel_tol, abstol=abs_tol)

     push!(l.trajectory, sol)
     
end

# eoms_2BP: equations of motion for the 2 body Problem with continous thrust
# vec0: start position, x,y,z. [km] and start velocity, vx,vy,vz. [km/s]
# m0: starting mass. [kg]
# thrust: fixed inertial thrust, ux,uy,uz. [N]
# mu: central body gravity constant. [km^3/s^2]
function eoms_2BP!(du, u, p, t)

     R0 = norm(u[1:3])
     
     F, T, veff, mu = p

     du[1] = u[4]
     du[2] = u[5]
     du[3] = u[6]
     du[4] = -mu/R0^3 * u[1] + T[1]/u[7]
     du[5] = -mu/R0^3 * u[2] + T[2]/u[7]
     du[6] = -mu/R0^3 * u[3] + T[3]/u[7]

     du[7] = -F/veff


end