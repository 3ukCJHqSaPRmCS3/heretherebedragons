""" Initial guess for the multiple shooting correction"""

include("population.jl")
include("structures.jl")

function initial_guess(x0, xf, sc::Spacecraft, flags)
    """ Propagate from the initial state for at least one revolution until the distance
    between the final stgate and the final point to reach is less than 100 km"""


    # Initialize the error in distance
    distance = 1e5
    flag_thrust, flag_J2 = flags
    @show x0, xf
    # Find orbital elements
    a_0, e_0, i_0, Ω_0, ω_0, θ_0, h_0, Ε_0 = cart2kep(x0)
    a_f, e_f, i_f, Ω_f, ω_f, θ_f, h_f, Ε_f = cart2kep(xf)
    @show a_0, e_0, i_0, Ω_0, ω_0, θ_0, h_0, Ε_0
    @show a_f, e_f, i_f, Ω_f, ω_f, θ_f, h_f, Ε_f
 
    # Find the orbits periods
    T_0 = 2*pi*sqrt(a_0^3/MU_EARTH)
    T_f = 2*pi*sqrt(a_f^3/MU_EARTH)
    Δt = 0.1


    # Check if the points belong to the same orbit
    if e_0 == e_f && h_0 == h_f && Ε_0 == Ε_f
        
        # Check on the TA
        if θ_f > θ_0
            print("same orbit different TA")
            thrust_specs = [1,1,1]
            T_ig = T_0 + 3600*0.5
        elseif θ_f < θ_0
            thrust_specs = [1,1,1]
            T_ig = T_0 - 3600*0.5
            Δt = -Δt
        end
        #Check on the RAAN
        if Ω_f > Ω_0
            thrust_specs = [-1,-1,-1]
            T_ig += 3600*10
        elseif Ω_f < Ω_0
            thrust_specs = [1,1,1]
            T_ig += 3600*10
        end

    elseif a_f != a_0

        if a_f > a_0
            thrust_specs = [1,1,1]
        else
            thrust_specs = [-1,-1,-1]
        end 
        T_ig = 2*T_0
    end

    while distance > 1e3

        # propagation
        tspan = T_ig # propagation in seconds

        m0 = sc.mass_dry + sc.mass_prop 
        state0 = vcat(x0, m0)
        sol = propagate_2BP(state0, tspan, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2 , args = thrust_specs)

        # retrieve states from the integration
        xf_end = sol.u[end][1:3]

        distance = norm(xf_end - xf[1:3])
 
        T_ig += Δt*3600
    end

    return T_ig
end