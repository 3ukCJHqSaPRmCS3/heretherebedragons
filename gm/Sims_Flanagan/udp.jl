#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: udp (user defined problem)
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file defines the optimal control problem to be optimazed using 
#  Sims-Flanagan transcription and Ipopt to solve the NLP. 
# This problem works by manipulating the starting epoch t0, the transfer time T the final 
# mass mf and the controls 
    # The decision vector is::
    # ctrl_vec = [t0, T, mf, Vxi, Vyi, Vzi, Vxf, Vyf, Vzf, controls]
#  
# ----------------------------------------------------------------------------------------

include("sims_flanagan_setup.jl")
include("utils.jl")

# Structure for the udp
# planets: vector of planets considered for the trajectory 
# spacecarft: takes the s/c characteristics (mass, T, isp)
# n_seg: number of segments inside the leg
# v_inf_bound: [lower, upper] bounds for the DV magnitude [km/s]
# t0_tf_bound: [dep, arr] bounds for the sate of arrival and departure in phase_free is false
# hf: High-fidelity. Activates a continuous representation for the thrust.
# mu: gravitational parameter
mutable struct Trajectory_lt_p2p
    planets::Vector{String}
    spacecraft::Spacecraft
    n_seg::Int64
    v_inf_bound::Vector{Float64}
    t0_tf_bound::Vector{String}
    mu::Float64
    hf::Bool
    lb::Vector{Float64}
    up::Vector{Float64}
    start_pos::Planet
    target_pos::Planet
    eq_const::Vector{Float64}
    ineq_const::Vector{Float64}
    Trajectory_lt_p2p(planets, spacecraft, n_max, v_inf_bound,phase_free,multi_obj,t0_tf_bound) =  
    if n_leg != 2 
        error("\nCheck that the number of of impulses N is not less than 2")
    elseif length(planets)!=2 
        error("Check that just two planets have been listed")
    # elseif length(t0_tf_bound) > 1
    #     error("Check that the phase_free option is set to false when time window are given")
    else
        new(planets, spacecraft, n_max, v_inf_bound,phase_free,multi_obj,t0_tf_bound, mu=mu_sun)
    end
end


function define_bound(udp::Trajectory_lt_p2p)
# define: function that takes the udp and computes boundaries and objective function

    #Define start date and arrival date
    udp.start_pos = planets(udp.planets[1]) ### CREATE THE JPL ephemeris here
    udp.target_pos = planets(udp.planets[2])

    # Define boundary condition on time of flight
    e_t0 = Epoch(udp.t0_tf_bound[1],"mjd2000")
    e_tof = Epoch(udp.t0_tf_bound[2],"mjd2000")

    # Define boundaries on velocity
    v_dep = udp.v_inf_bound[1] * 1000 #[m]
    v_arr = udp.v_inf_bound[2] * 1000 #[m]

    # Define boundaries
    ups.lb = [udp.t0[1], udp.tof[1], udp.spacecraft.mass * 0.1, 
    ones(3)*(-v_dep), ones(3)*(-v_arr), (-ones(3*udp.n_seg))]
    ups.ub = [udp.t0[2], udp.tof[2], udp.spacecraft.mass, 
    ones(3)*(v_dep), ones(3)*(v_arr), ones(3*udp.n_seg)]

    return v_dep, v_arr
end

function define_objective_fun(udp::Trajectory_lt_p2p, ctrl_vec)

    p0, pf, v_lb, v_ub = define_bound(udp)

    # Epoch in mjd2000
    t0 = epoch_date_converter(Epoch(ctrl_vec[1],"mjd2000"))
    tf = epoch_date_converter(Epoch(ctrl_vec[1] + ctrl_vec[2],"mjd2000"))

    # Final mass
    mf = ctrl_vec[3]

    # Controls
    u = ctrl_vec[10:end]

    # Cartesian states of the planets
    r0, v0 = eph(t0, p0) 
    rf, vf = eph(tf, pf)

    for (vc, vi) in zip(v0, ctrl_vec[4:6])
        v0 = push!(vc + vi)
    end

    for (vc, vi) in zip(vf, ctrl_vec[7:9])
        vf = push!(vc + vi)
    end

    # Spacecraft states
    x0 = Sc_state(r0, v0, udp.spacecraft.mass)
    xf = Sc_state(rf, vf, mf)

    # Set leg
    p2p_leg = Leg()
    setfield!(p2p_leg, :start_time, t0)
    setfield!(p2p_leg, :x0, x0)
    setfield!(p2p_leg, :throttles, u)
    setfield!(p2p_leg, :end_time, tf)
    setfield!(p2p_leg, :xf, xf)

    ## To do : change this to use the sims flanagan function I wrote 
    # Compute the equality constraints
    udp.eq_const = equality_constraint(p2p_leg)

    # adimensionalization
    # eq_const[1:3] /= AU
    # eq_const[4:6] /= Earth_vel
    # eq_const[7] /= udp.spacecraft.mass

    # Compute the inequality constraints
    udp.ineq_const = inequality_constraint(p2p_leg)

    # Compute inequality constraints on departure and arrival velocities
    v_dep_con = (ctrl_vec[3]^2 + ctrl_vec[4]^2 + ctrl_vec[5]^2 - v_lb^2)
    v_arr_con = (ctrl_vec[6]^ 2 + ctrl_vec[7]^2 + ctrl_vec[8]^ 2 - v_ub^2)

    # nondimensionalize inequality constraints
    v_dep_con /= Earth_vel^2
    v_arr_con /= Earth_vel^2

    return hcat(([-mf], eq_const, ineq_const, [v_dep_con, v_arr_con]))

end

function print_transf(udp::Trajectory_lt_p2p, ctrl_vec)

    print("\nLow-thrust NEP transfer from ",
              udp.start_pos.name, " to ", udp.target_pos.name)
    print("\nLaunch epoch:", ctrl_vec[1] ,"MJD2000, a.k.a. ", epoch_date_converter(Epoch(ctrl_vec[1],"string")))
    print("\nArrival epoch:",ctrl_vec[1] + ctrl_vec[2] ,"MJD2000, a.k.a.",
     Epoch_date_converter(Epoch(ctrl_vec[1] + ctrl_vec[2], "string")))
    print("\nTime of flight (days):",ctrl_vec[2])
    print("\nLaunch DV (km/s)", [ctrl_vec[3]^2, ctrl_vec[4]^2 ,ctrl_vec[5]^2]/1000, " - ",
     [ctrl_vec[3] / 1000, ctrl_vec[4] / 1000, ctrl_vec[5] / 1000])
    print("\nArrival DV (km/s)",[ctrl_vec[6]^2, ctrl_vec[7]^2, ctrl_vec[8]^2]/1000 ,"-",
     [ctrl_vec[6] / 1000, ctrl_vec[7] / 1000, ctrl_vec[8] / 1000])
end
