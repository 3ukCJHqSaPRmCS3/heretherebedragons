""" definition of the structures in the algorithm which are not defined/do not have their own name_stc.jl file
"""

mutable struct Spacecraft
    mass_dry::Float64 #[kg]
    mass_prop::Float64 #[kg]
    Isp::Float64 #[1/s]
    Thrust::Float64
    state::Vector{Float64} #[km,km,km,km/s,km/s,km/s]
    Spacecraft() = new()
end

# Orbit structure:
# elements are:
# semi-major axis [km]
# eccentricity [deg]
# inclination [deg]
# right ascension of the ascending node [deg]
# argument of perigee [deg]
# true anomaly [deg], default value is 0 deg
mutable struct Orbit
    sma::Float64
    ecc::Float64
    inc::Float64
    raan:: Float64
    aop::Float64
    ta::Float64
    mu::Float64
end

struct Constellation
    n_orbit::Int64
    orbits::Vector{Orbit}
    spacecrafts::Vector{Spacecraft}
end

mutable struct Cost_fnx
    value::Vector{Float64}
    type::String
    weight::Float64
    Cost_fnx() = new()
end

# Population structure:
# list_states: vector cointaing the pairs [xAi, xBi]
# J_st: cost value associated with the time of transfer of single tuple 
# J_sf: cost value associated with the fuel remaining in every staellite after transfer 
# J_mt: cost value associated with the mean time of transfer of the Population
# J_mf: cost value associated with the mean fuel consumption of the Population

mutable struct Population
    list_states::Matrix{Float64}
    J_st::Cost_fnx 
    J_sf::Cost_fnx
    J_mt::Cost_fnx
    J_mf::Cost_fnx
    Population() = new()
end



mutable struct Transfer
    V::Vector{Float64}
    Δm::Float64
    states::Matrix{Float64}
    converged::Bool
    Transfer() = new()
end