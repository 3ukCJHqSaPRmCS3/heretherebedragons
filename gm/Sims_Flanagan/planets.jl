#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: propagators
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file contains the definition of the plantes class and the related
#  functions
#  
# ----------------------------------------------------------------------------------------

include("utils.jl")

# Planet structure
# mu_central_body: Gravity parameter of the central body (this is not used to compute the ephemerides)
# mu_self: Gravity parameter of the target
# radius: Radius of target body
# self_radius: Safe radius of target body
# name: Body name
mutable struct Planet
    name::String
    mu_central_body::Float64
    mu_self::Float64
    radius::Float64
    self_radius::Float64
    epoch
    Planet() = new()
end


function eph(e, p::Planet)
    p.epoch = e
end

function get_planet_velocity(p::Planets)
    #Define function to find velocity given ephemeris
end