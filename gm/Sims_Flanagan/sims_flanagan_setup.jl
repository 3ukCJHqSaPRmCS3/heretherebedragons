#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: sims_flanagan_setup
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file defines all the type of structures that are used
#  in the computation of the optimal low-thrust trajectory using the Sims-Flanagan method.
#  A print function is also defined to have easier access to the value of the structures.
# ----------------------------------------------------------------------------------------


# Structure for the spacecraft characteristics
# mass: the spacecraft mass [kg]
# thrust: the maximum thrust of the spacecraft propulsion system [N]
# isp: the specific impulse of the spacecraft propulsion system [s]
struct Spacecraft
    mass::Float64
    thrust::Float64
    isp::Float64
end


# Structure for the spacecraft states
# r: triplet containing the position vector in cartesian coordiantes [km]
# v: triplet containing the velocity vector in cartesian coordiantes [km/s]
# m: mass [kg]
struct Sc_state
    r::Vector{Float64}
    v::Vector{Float64}
    m::Float64
end

# Structure for the throttle profile
# start_time: starting epoch in the format "yyy-mm-ddThh:mm:ss""
# end_time: ending epoch in the format "yyy-mm-ddThh:mm:ss""
# value: triplet containing the cartesian values of the throttle 
mutable struct Throttle
    start_time::String
    end_time::String
    value::Vector{Float64}
    Throttle() = new()
end

# Structure for the leg elements
# start_time: starting epoch in the format "yyy-mm-ddThh:mm:ss""
# x0: starting sc_state
# throttles: tuple containing the 3N cartesian components of the throttle
# end_time: final epoch in the format "yyy-mm-ddThh:mm:ss""
# xe: final sc_state
# spacecraft: spacecraft
# mu: central body gravity parameter [km^3/s^2]
# hf: Boolean to consider maneuver as real low thrust (true) or impulsive (false)
# trajectory: state vector in the format [x,y,z,vx,vy,vz,m]
mutable struct Leg
    start_time::String
    x0::Sc_state
    throttles::Vector{Throttle}
    end_time::String
    xf::Sc_state
    spacecraft::Spacecraft
    mu::Float64
    hf::Bool
    trajectory = []
    freemass = false
    freetime = false
    Leg() = new()
end

# Structure fore the leg states
# n_seg: number of segments
# c: constant in the Sundmann transformation dt = cr^(alpha ds)
# alpha: exponent in the Sundmann transformation dt = cr^(alpha ds)
# tol: log 10 tolerance set in the Taylor integration of the leg
struct Leg_s
    n_seg::Int
    c::Float64
    alpha::Float64
    tol::Float64
end

# print_characteristics: Function that prints out the argument of any structure of the module
function print_characteristics(obj)
    if typeof(obj) == Spacecraft
        print("The spacecraft has mass: ", obj.mass,"\nthrust: ",
        obj.thrust, "\nisp: ", obj.isp)
    elseif typeof(obj) == Sc_state
        print("The spacecraft states are\nposition: ", obj.r,"\nvelocity: ",
        obj.v, "\nmass: ", obj.m)
    elseif typeof(obj) == Throttle
        print("The throttle elements are\nstart time: ", obj.start_time,"\nend_time: ",
        obj.end_time, "\nvalue: ", obj.value)
    elseif typeof(obj) == Leg
        print("The throttle elements are\nstart time: ", obj.start_time,"\nend_time: ",
        obj.end_time, "\nthrottle: ", obj.throttles, "\nx0: ", obj.x0, "\nxf: ", obj.xf, "\nspacecraft: ", obj.spacecraft)
    elseif typeof(obj) == Leg_s
        print("The leg elements are\nn segments: ", obj.n_seg,"\nSundmann const: ",
        obj.c, "\nSudmann exponent: ", obj.alpha, "\nIntegration tolerance:", obj.tol)
    end

end
