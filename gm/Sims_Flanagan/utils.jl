#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: propagators
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file contains useful variables, functions and structures for the 
#  computation of the Sims-Flanagan algorithm
#  
# ----------------------------------------------------------------------------------------

using LinearAlgebra

# Constants
global DAY2SEC = 24*3600 # factor to transfor day into seconds [s]
global SEC2DAY = 1/24/3600 # factor to transfor day into seconds [s]
global AU = 1495978707*1e11 # Astronomical Unit [km]
global g0 = 9.80665 # gravitational acceleration [m/s^2]
global R_e = 6371 # Earth radius [km] 
global mu_earth = 398600.4415 # earth gravity parameter [km^3/s^2]
global mu_sun = 1.327e11




function newton_raphson(f, df, x0, tol, max_iter )
# newton_raphson: iterative method to find the root of a function
# f: function in the form f(x) = 0 to solve to find x
# df: derivate of f(x) = 0 in function of x
# tol: tolerance on the error on the solution
# max_iter: maximum number of iterations allowed to find x

    error = 1
    iter = 1
    while error > tol && iter < max_iter
        x  = x0 - f/df
        error = abs(E-E0)
        x0 = x
        iter +=1
    end

    return x
end



function solve_keplers_eq_E(E0, M, e, tol, max_iter)
# solve_keplers_eq_E: Find the eccentric anomaly E starting from mean anomaly using Newton's method
# E0: initial guess on the eccentric anomaly E [rad]
# M: Mean anomaly value [rad]
# R: magnitude of the position vector [km]
# a: semimajor axis [km]
# σ0: parameter sigma
# e: eccentricity vector 
# tol: tolerance on the error on the solution
# max_iter: maximum number of iterations allowed to find E

    ΔE = 1
    iter = 1
    E=0
    while ΔE > tol && iter < max_iter
        
        E = E0 - (E0 - norm(e) * sin(E0) - M)/(1 - norm(e) * cos(E0))
        ΔE = abs(E-E0)
        E0 = E
        iter +=1
        
        if iter == max_iter
            print("the algorithm did not converge to a solution for E")
        end
    end

    return E
end


function solve_keplers_eq_H(H0, M, e, tol, max_iter)
# Find the hyperbolic anomaly H starting from mean anomaly using Newton's method

    ΔH = 1
    iter = 1
    H = 0

    while ΔH > tol && iter < max_iter
        H  = H0 - (norm(e) * sinh(H0) - H0 - M)/(norm(e) * cosh(H0) - 1)
        ΔH = abs(H - H0)
        H0 = H
        iter +=1

        if iter == max_iter
            print("the algorithm did not converge to a solution for H")
        end
    end
    
    return H
end

# Compute the epoch in different formats starting from a julian date
# date: a julian date, default is 0
# date_type: julian date type, one of “jd”, “mjd” and “mjd2000”, defaults to “mjd2000”
mutable struct Epoch
    date
    date_type::String
    Epoch() = new(0, "mjd") 
end

function epoch_date_converter(e::Epoch)
# Function to convert date or epoch in Julian date into a Julian date type

    if typeof(e.date) == String && e.date_type == "jd"
        return datetime2julian(DateTime(e.date))
    elseif typeof(e.date) == String && e.date_type == "mjd"
        return datetime2julian(DateTime(e.date)) - 2400000.5
    elseif typeof(e.date) == String && e.date_type == "mjd2000"
        return datetime2julian(DateTime(e.date)) - 2400000.5 - 51544
    elseif typeof(e.date) == Float64 || typeof(e.date) == Int64 && e.date_type == "jd"
        print("Warning: The inserted data should already be in Julian Date, no conversion is applied")
        return e.date
    elseif typeof(e.date) == Float64 || typeof(e.date) == Int64 && e.date_type == "mjd"
        return e.date - 2400000.5
    elseif typeof(e.date) == Float64 || typeof(e.date) == Int64 && e.date_type == "mjd2000"
        return e.date - 2400000.5 - 51544
    elseif typeof(e.date) == Float64 || typeof(e.date) == Int64 && e.date_type == "string"
        return julian2datetime(e.date)
    else
        error("The data format or the data type specified are not compliant")
    end
    
end


# TO DO COMPUTE PERIOD FUNCTION