#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: propagators
#  Module: wave_GA
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file defines all the propagators that are used to run the low-thrust 
#  trajectories or to plot orbits in the 2BP
# 
#  
# ----------------------------------------------------------------------------------------

using LinearAlgebra
using DifferentialEquations
include("constants.jl")
# propagate_2BPt: High fidelity propagator
# r: start position, x,y,z. [km]
# v: start velocity, vx,vy,vz. [km/s]
# m0: starting mass. [kg]
# tof: propagation time. [s]
# abstol a tuple (rf, vf, mf) containing the final position, velocity and mass after the propagation.
# reltol: the relative tolerance passed to taylor propagator.
# flag_thrust: true if want to include thrust perturbation in eom, default false
# flag_J2: true if want to include oblateness perturbation in eom, default false
# flag_STM: true if the computation of the State Transition Matrix is required
# args: additional argument ex: (1) T: fixed inertial thrust, ux,uy,uz. [N]
function propagate_2BP(x0, tof, sc::Spacecraft; rel_tol=1e-8, abs_tol=1e-8, flag_thrust=false, flag_J2=false, flag_STM=false, args = [])

    p = []
    push!(p, flag_thrust, flag_J2, flag_STM)

    if flag_thrust == true
        β = args
        # Compute the engine efficiency
        veff = sc.Isp * g0 * 1000.0
        # Define the additonal arguments for the eoms function
        push!(p, sc.Thrust, veff, β)
    end

    # Define the initial guess and the integration time span
    init_guess = x0
    t_span = (0.0, tof)

    # Define the ODE Problem
    prob = ODEProblem(eoms_2BP!, init_guess, t_span, p)
    # Solve the ODE Problem using Tsitouras 5/4 Runge-Kutta method and user-defined tolerances
    sol = solve(prob, Vern7(),adaptive=false, dt = 100, reltol=rel_tol, abstol=abs_tol)
    sc.mass_prop = sol.u[end][7] - sc.mass_dry  
    return(sol)
     
end

# eoms_2BP: equations of motion for the 2 body Problem 
# With flags choice of adding (1) continous thrust, (2) J2, oblateness effect
# u: start position, x,y,z. [km] and start velocity, vx,vy,vz. [km/s] and mass [kg]
# p: additional arguments for the eoms depending on the flag
# mu: central body gravity constant. [km^3/s^2]
function eoms_2BP!(du, u, p, t)

    R0 = norm(u[1:3])
    mu = MU_EARTH
    flag_thrust, flag_J2, flag_STM = p[1:3]

    # Compute vx, vy, vz [km/s]
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # Compute ax, ay, az [km/s^2]
    du[4] = -mu/R0^3 * u[1] 
    du[5] = -mu/R0^3 * u[2] 
    du[6] = -mu/R0^3 * u[3] 

    if flag_thrust == true
        #Thrust, exhaust velocity and thrust direction/throttle
        T, veff, β = p[4:6]
        #Find the unit vector of the velocity
        normed = u[4:6]/norm(u[4:6])
        # Compute ax, ay, az [km/s^2]
        du[4] += β[1] * normed[1] * T / u[7] / 1000.0
        du[5] += β[2] * normed[2] * T / u[7] / 1000.0
        du[6] += β[3] * normed[3] * T / u[7] / 1000.0

        # Compute m_dot
        du[7] = -T/veff
    end

    if flag_J2 == true 
        Re = R_EARTH
        # Compute ax, ay, az [km/s^2]
        tx = 1 - 5 * u[3]^2 / R0^2 
        ty = 1 - 5 * u[3]^2 / R0^2
        tz = 3 - 5 * u[3]^2 / R0^2
        du[4] -= 1.5 * mu * J2 * Re ^2 * u[1] / R0^5 * tx
        du[5] -= 1.5 * mu * J2 * Re ^2 * u[2] / R0^5 * ty
        du[6] -= 1.5 * mu * J2 * Re ^2 * u[3] / R0^5 * tz
    end

    if flag_STM == true
        # Reconstruct the matrix Phi at t0 from the input vector
        Φ = reshape(u[8:43], (6,6))

        # Initialize the matrix A, where Φ̇ = A * Φ
        A = Matrix{Float64}(undef,6,6)

        # Create the 0x0(3x3) matrix
        O = zeros(3,3)
        I_mtx = I + zeros(3,3)
        # Create the G matrix 
        gxx = mu*(3*u[1]^2/R0^5 - 1/R0^3)
        gxy = mu*(3*u[1]*u[2]/R0^5)
        gxz = mu*(3*u[1]*u[3]/R0^5)
        gyx = gxy
        gyy = mu*(3*u[2]^2/R0^5 - 1/R0^3)
        gyz = mu*(3*u[2]*u[3]/R0^5)
        gzx = gxz
        gzy = gyz
        gzz = mu*(3*u[3]^2/R0^5 - 1/R0^3)


        if flag_thrust == true
            cx = β[1] * T / u[7] / 1000.0
            cy = β[2] * T / u[7] / 1000.0
            cz = β[3] * T / u[7] / 1000.0
            gxx += cx * (1/R0 - u[1]^2/R0^3)
            gxy -= cx * (u[1]*u[2]/R0^3)
            gxz -= cx * (u[1]*u[3]/R0^3)
            gyx -= cy * (u[1]*u[2]/R0^3)
            gyy += cy * (1/R0 - u[2]^2/R0^3)
            gyz -= cy * (u[2]*u[3]/R0^3)
            gzx -= cz * (u[1]*u[3]/R0^3)
            gzy -= cz * (u[2]*u[3]/R0^3)
            gzz += cz * (1/R0 - u[3]^2/R0^3)
        end

        # if flag_J2 == true 
        #     c = 1.5 * mu * J2 * R_EARTH^2
        #     gxx -= c * (1/R0^5 - 5*(u[1]^2 + u[3]^2)/R0^7 + 35*u[1]^2*u[3]^2/R0^9)
        #     gxy -= c * (35*u[1]*u[2]*u[3]^2/R0^9 - 5*u[1]*u[2]/R0^7)
        #     gxz -= c * (35*u[1]*u[3]^3/R0^9 - 15*u[1]*u[3]/R0^7)
        #     gyx -= gxy
        #     gyy -= c * (1/R0^5 - 5*(u[2]^2 + u[3]^2)/R0^7 + 35*u[2]^2*u[3]^2/R0^9)
        #     gyz -= c * (35*u[2]*u[3]^3/R0^9 - 15*u[2]*u[3]/R0^7)
        #     gzx -= gxz
        #     gzy -= gyz
        #     gzz -= c * (3/R0^5 - 30*u[3]^2/R0^7 + 35*u[3]^4/R0^9)
        # end

        G = [gxx gxy gxz; gyx gyy gyz; gzx gzy gzz]


        # Fill in the A matrix as A = [O, I; G, O]
        A[1:3,1:3] = O
        A[4:6, 4:6] = O
        A[4:6,1:3] = G
        A[1:3, 4:6] = I_mtx

        # Compute Φ̇
        Φ̇ = A * Φ

        # Compute the ouptut
        du[end-35:end] = reshape(Φ̇ ,(1,36))

    end 


end

function STM(u, t, p)

    flag_thrust, flag_J2 = p[1:2]
    mu = MU_EARTH
    R0 = norm(u[1:3])
    # Initialize the matrix A, where Φ̇ = A * Φ
    A = Matrix{Float64}(undef,6,6)

    # Create the 0x0(3x3) matrix
    O = zeros(3,3)
    I_mtx = I + zeros(3,3)
    # Create the G matrix 
    gxx = mu*(3*u[1]^2/R0^5 - 1/R0^3)
    gxy = mu*(3*u[1]*u[2]/R0^5)
    gxz = mu*(3*u[1]*u[3]/R0^5)
    gyx = gxy
    gyy = mu*(3*u[2]^2/R0^5 - 1/R0^3)
    gyz = mu*(3*u[2]*u[3]/R0^5)
    gzx = gxz
    gzy = gyz
    gzz = mu*(3*u[3]^2/R0^5 - 1/R0^3)


    if flag_thrust == true
        #Thrust, exhaust velocity and thrust direction/throttle
        T, veff, β = p[4:6]
        cx = β[1] * T / u[7] / 1000.0
        cy = β[2] * T / u[7] / 1000.0
        cz = β[3] * T / u[7] / 1000.0
        gxx += cx * (1/R0 - u[1]^2/R0^3)
        gxy -= cx * (u[1]*u[2]/R0^3)
        gxz -= cx * (u[1]*u[3]/R0^3)
        gyx -= cy * (u[1]*u[2]/R0^3)
        gyy += cy * (1/R0 - u[2]^2/R0^3)
        gyz -= cy * (u[2]*u[3]/R0^3)
        gzx -= cz * (u[1]*u[3]/R0^3)
        gzy -= cz * (u[2]*u[3]/R0^3)
        gzz += cz * (1/R0 - u[3]^2/R0^3)
    end

    if flag_J2 == true 
        c = 1.5 * mu * J2 * R_EARTH^2
        gxx -= c * (1/R0^5 - 5*(u[1]^2 + u[3]^2)/R0^7 + 35*u[1]^2*u[3]^2/R0^9)
        gxy -= c * (35*u[1]*u[2]*u[3]^2/R0^9 - 5*u[1]*u[2]/R0^7)
        gxz -= c * (35*u[1]*u[3]^3/R0^9 - 15*u[1]*u[3]/R0^7)
        gyx -= gxy
        gyy -= c * (1/R0^5 - 5*(u[2]^2 + u[3]^2)/R0^7 + 35*u[2]^2*u[3]^2/R0^9)
        gyz -= c * (35*u[2]*u[3]^3/R0^9 - 15*u[2]*u[3]/R0^7)
        gzx -= gxz
        gzy -= gyz
        gzz -= c * (3/R0^5 - 30*u[3]^2/R0^7 + 35*u[3]^4/R0^9)
    end

    G = [gxx gxy gxz; gyx gyy gyz; gzx gzy gzz]

    # Fill in the A matrix as A = [O, I; G, O]
    A[1:3,1:3] = O
    A[4:6, 4:6] = O
    A[4:6,1:3] = G
    A[1:3, 4:6] = I_mtx
    
    # Compute Φ̇
    Φ = exp(A*t)

    return Φ

end