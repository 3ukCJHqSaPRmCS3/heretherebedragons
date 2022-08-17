#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: sims_flanagan_functions
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file defines all the functions that are used to construct low-thrust 
#  trajectories using the Sims-Flanagan method.
#  *************************NOT IN USE**************************************
# ----------------------------------------------------------------------------------------

include("sims_flanagan_setup.jl")
include("propagators.jl")
include("utils.jl")
using LinearAlgebra
using Dates

function get_leg_states(l::Leg, high_fidelity::Bool)
# get_leg_states: return the spacecrfat states (t,r,v,m) at the leg grid points
# Ex: # times,r,v,m = get_leg_states(l, 'low_fidelity)


    # Compute the number of segments for forward and backward propagation
    n_seg = length(l.throttles)
    fwd_seg = floor(Int64,(n_seg + 1)/ 2)
    back_seg = floor(Int64, n_seg /2)
    print("\nfwd_seg=", fwd_seg, "\nback_seg=", back_seg)
    # start_time_jmd = datetime2julian(DateTime(l.start_time)) - 2400000.5
    # end_time_mjd = datetime2julian(DateTime(l.end_time)) - 2400000.5

    # Extract info from spacecarft
    sc = l.Spacecraft
    isp = sc.isp
    max_thrust = sc.thrust

    # Estract informatio on the Leg
    throttles = l.Throttles
    mu = l.mu

    # time grid (the mismatch is repeated twice)
    t_grid = zeros(n_seg * 2 + 2)

    ### Forward propagation ###
 
    # x,y,z contain the cartesian components of all points (grid+midpoints)
    x =  zeros(fwd_seg * 2 + 1)
    y = zeros(fwd_seg * 2 + 1)
    z = zeros(fwd_seg * 2 + 1)
    vx = zeros(fwd_seg * 2 + 1)
    vy = zeros(fwd_seg * 2 + 1)
    vz = zeros(fwd_seg * 2 + 1)
    mass = zeros(fwd_seg * 2 + 1)

    state = l.x0

    # Initial conditions
    r = state.r
    v = state.v
    m = state.m
    x[1] = r[1]
    y[1] = r[2]
    z[1] = r[3]
    vx[1] = v[1]
    vy[1] = v[2]
    vz[1] = v[3]
    mass[1] = m


    for (idx, val) in enumerate(throttles[1:fwd_seg])
    # Compute propagation points

        # Convert the string date into modified julian dates and mjd2000 (ref: https://kelvins.esa.int/gtoc9-kessler-run/constants/)
        start_time_mjd = datetime2julian(DateTime(val.start_time)) - 2400000.5
        start_time_mjd2000 = start_time_mjd - 51544
        end_time_mjd = datetime2julian(DateTime(val.end_time)) - 2400000.5
        end_time_mjd2000 = end_time_mjd - 51544

        t_grid[idx] = start_time_mjd2000
        t_grid[idx+1] = start_time_mjd2000 + 
            (end_time_mjd2000 - start_time_mjd2000) / 2.
        dt = (end_time_mjd - start_time_mjd) * DAY2SEC
        alpha = min(norm(val.value), 1.0)
        

        # Keplerian propagation and dV application
        if high_fidelity == false
            for dumb in val.value 
                global dV = [max_thrust / m * dt * dumb]
            end

            r, v, θ= propagate_lagrangian(r, v, dt/2, mu)
            x[idx+1] = r[1]
            y[idx+1] = r[2]
            z[idx+1] = r[3]
            vx[idx+1] = v[1]
            vy[idx+1] = v[2]
            vz[idx+1] = v[3]
            mass[idx+1] = m

            # Add dV to current v
            for (a,b) in zip(v,dV)
                global v = [a + b]            
            end
            r, v, θ = propagate_lagrangian(r, v, dt / 2, mu)
            m = m * exp(-norm(dV) / isp / g0)

            x[idx+2] = r[1]
            y[idx+2] = r[2]
            z[idx+2] = r[3]
            vx[idx+2] = v[1]
            vy[idx+2] = v[2]
            vz[idx+2] = v[3]
            mass[idx+2] = m

        # Taylor propagation of constant thrust u
        else
            for dumb in val.value 
                u = [max_thrust * dumb]
                return u
            end

            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp *g0, -12, -12)
            x[2*idx+1] = r[1]
            y[2*idx+1] = r[2]
            z[2*idx+1] = r[3]
            vx[2*idx+1] = v[1]
            vy[2*idx+1] = v[2]
            vz[2*idx+1] = v[3]
            mass[2*idx+1] = m

            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp * g0, -12, -12)
            x[2*idx+2] = r[1]
            y[2*idx+2] = r[2]
            z[2*idx+2] = r[3]
            vx[2*idx+2] = v[1]
            vy[2*idx+2] = v[2]
            vz[2*idx+2] = v[3]
            mass[2*idx+2] = m

        end
        t_grid[idx+2] = end_time_mjd2000
    end
    
    
    ### Backward propagation ###

    # Cartesian components of backward propagation
    x_back = 0.123 * ones(back_seg * 2 + 1)
    y_back = 0.123 * ones(back_seg * 2 + 1)
    z_back = 0.123 * ones(back_seg * 2 + 1)
    vx_back = zeros(back_seg * 2 + 1)
    vy_back = zeros(back_seg * 2 + 1)
    vz_back = zeros(back_seg * 2 + 1)
    mass_back = zeros(back_seg * 2 + 1)
    state = l.xf

    # Final conditions
    r = state.r
    v = state.v
    m = state.m
    x_back[end] =  r[1]
    y_back[end] = r[2]
    z_back[end] = r[3]
    vx_back[end] = v[1]
    vy_back[end] = v[2]
    vz_back[end] = v[3]
    mass_back[end] = m


    for (idx, val) in enumerate(reverse(throttles)[1: back_seg])
    # Compute propagation points
            
        # Convert the string date into modified julian dates and mjd2000 (ref: https://kelvins.esa.int/gtoc9-kessler-run/constants/)
        start_time_mjd = datetime2julian(DateTime(val.start_time)) - 2400000.5
        start_time_mjd2000 = start_time_mjd - 51544
        end_time_mjd = datetime2julian(DateTime(val.end_time)) - 2400000.5
        end_time_mjd2000 = end_time_mjd - 51544
        
        t_grid[end-2*idx+1] = end_time_mjd2000 - 
            (end_time_mjd2000 - start_time_mjd2000) / 2.
        t_grid[end-2*idx+2] = end_time_mjd2000
        dt = (end_time_mjd - start_time_mjd) * DAY2SEC
        alpha = min(norm(val.value), 1.0)

        # Keplerian propagation and dV application
        if high_fidelity == false
            for dumb in val.value
                global dV = [max_thrust / m * dt * dumb]
            end
            r, v, θ = propagate_lagrangian(r, v, -dt / 2, mu)
            x_back[end-2*idx+1] = r[1]
            y_back[end-2*idx+1] = r[2]
            z_back[end-2*idx+1] = r[3]
            vx_back[end-2*idx+1] = v[1]
            vy_back[end-2*idx+1] = v[2]
            vz_back[end-2*idx+1] = v[3]
            mass_back[end-2*idx+1] = m
            # Add dV to current v
            for (a,b) in zip(v, dV)
                global v = [a - b]
            end
            r, v, θ = propagate_lagrangian(r, v, -dt / 2, mu)
            m = m * exp(norm(dV) / isp /g0)

            x_back[end-2*idx] = r[1]
            y_back[end-2*idx] = r[2]
            z_back[end-2*idx] = r[3]
            vx_back[end-2*idx] = v[1]
            vy_back[end-2*idx] = v[2]
            vz_back[end-2*idx] = v[3]
            mass_back[end-2*idx] = m

            # Taylor propagation of constant thrust u
        else
            for dumb in val.value 
                global u = [max_thrust * dumb]
            end

            r, v, m = propagate_taylor(r, v, m, u, -dt / 2, mu, isp * g0, -12, -12)
            x[-2*idx-2] = r[1]
            y[-2*idx-2] = r[2]
            z[-2*idx-2] = r[3]
            vx[-2*idx-2] = v[1]
            vy[-2*idx-2] = v[2]
            vz[-2*idx-2] = v[3]
            mass[-2*idx-2] = m

            r, v, m = propagate_taylor(r, v, m, u, dt / 2, mu, isp * g0, -12, -12)
            x[-2*idx-3] = r[1]
            y[-2*idx-3] = r[2]
            z[-2*idx-3] = r[3]
            vx[-2*idx-3] = v[1]
            vy[-2*idx-3] = v[2]
            vz[-2*idx-3] = v[3]
            mass[-2*idx-3] = m

        end
        t_grid[end-idx] = start_time_mjd2000

    end


    x = vcat(x, x_back)
    y = vcat(y, y_back)
    z = vcat(z, z_back)
    vx = vcat(vx, vx_back)
    vy = vcat(vy, vy_back)
    vz = vcat(vz, vz_back)
    mass = vcat(mass, mass_back)

    return(t_grid, collect(zip(x, y, z)), collect(zip(vx, vy, vz)), mass)
end


## TO DO: Definition of a function to compute the ephemeris for high fidelity computations: leg_eph









    












            


        








