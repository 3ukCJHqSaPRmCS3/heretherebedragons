"""Definition of orbits and states along the orbits for creation of initial population
"""


using LinearAlgebra
include("constants.jl")
include("structures.jl")
include("constants.jl")

function kep2cart(o::Orbit, angle::Bool)
    """ Method to find the cartesian state from an orbit espressed in keplerian elements
    """
    
    rraan  = o.raan * DEG2RAD
    raop  = o.aop * DEG2RAD
    rta = o.ta * DEG2RAD
    rinc  = o.inc * DEG2RAD

    if o.ecc == 0.0 && o.inc == 0.0
        o.aop = 0.0
        o.raan = 0.0
        # o.ta = ta_true
    elseif o.ecc == 0.0 && o.inc != 0.0
        o.aop = 0.0
        o.ta  = o.aop + o.ta
    elseif o.ecc != 0.0 && o.inc == 0.0
        o.raan = 0.0
        # o.aop + aop_true
    end

    R = [cos(rraan)*cos(raop)-sin(rraan)*sin(raop)*cos(rinc) -cos(rraan)*sin(raop)-sin(rraan)*cos(raop)*cos(rinc) sin(rraan)*sin(rinc);
    sin(rraan)*cos(raop)+cos(rraan)*sin(raop)*cos(rinc) -sin(rraan)*sin(raop)+cos(rraan)*cos(raop)*cos(rinc) -cos(rraan)*sin(rinc);
    sin(raop)*sin(rinc) cos(raop)*sin(rinc) cos(rinc)]

    if angle == true
        theta = [rta]
    elseif angle == false
        # Set a vector of true anomaly starting from 1 deg
        theta = LinRange(0, 2*pi, 100)
    end

    # Define the ouptut vectors 
    r_eci = zeros(size(theta,1), 3)
    v_eci = zeros(size(theta,1), 3)

    for i = 1:size(theta,1)

        # Compute the eccentric anomaly
        E = 2 * atan(sqrt((1 - o.ecc) / (1 + o.ecc)) * sin(theta[i]/ 2), cos(theta[i] / 2))

        # Compute the radius magnitude
        p = o.sma * (1 - o.ecc^2)# * cos(E))

        #Compute the position and velocity vector in the orbital frame
        r_of = p/(1 + o.ecc*cos(theta[i])) * [cos(theta[i]); sin(theta[i]); 0]
        v_of = sqrt(o.mu /p) * [-sin(theta[i]); (o.ecc + cos(theta[i])); 0]
        # print("\n r is ", r_of)
        r_xyz = R * r_of
        v_xyz = R * v_of

        # fill the output vector
        r_eci[i,:] = r_xyz
        v_eci[i,:] = v_xyz

    end 

    return r_eci, v_eci

end


function cart2kep(x::Vector{Float64})
    # Find position and velocity from input vector
    r = x[1:3]
    v = x[4:6]

    mu = MU_EARTH

    #Define X,Y,Z
    X = [1,0,0]
    Y = [0,1,0]
    Z = [0,0,1]

    # Find specific angular momentum
    h = cross(r,v)
    h_norm = norm(h)

    # Find line of nodes
    n = cross(Z, h)

    # Find the eccentricity vector and norm
    e = cross(v,h)/mu - r/norm(r)
    e_norm = norm(e)

    # Find the energy
    Ε = norm(v)^2/2 - mu/norm(r)

    if e_norm != 1.0
        a = -mu/(2*Ε)
        p = a*(1-e_norm^2)
    else
        p = h_norm^2/mu
        a = Inf
    end

    # Find the inclination
    i = acos(h[3]/h_norm)

    #Find the right ascension of the ascending node (RAAN)
    Ω = acos(n[1]/norm(n))
    if n[2] < 0 
        Ω = 2*pi - Ω
    end

    #Find the argument of perigee (aop)
    ω = acos(dot(n, e)/(norm(n)*e_norm))
    if e[3] < 0 
        ω = 2*pi - ω
    end

    #Find true anomaly (θ)
    if i == 0.0 && e_norm > 0.0 #equatorial eccentric
        θ = acos(e[1]/norm(e))
        if e[2] < 0 
            θ = 2*pi - θ
        end
    elseif e_norm == 0.0 && i > 0.0 #circular inclined
        θ = acos(dot(n,r)/(norm(n)*e_norm))
        if r[3] < 0 
            θ = 2*pi - θ
        end
    elseif e_norm == 0.0 && i == 0.0 # circular equatorial
        θ = acos(r[1]/norm(r))
        if r[2] < 0 
            θ = 2*pi - θ
        end
    
    elseif e_norm != 0.0 && i != 0.0   # general trajectory
        @show r, e, e_norm
        θ = acos(dot(e, r)/(norm(r)*e_norm))
        if dot(r,v) < 0 
            θ = 2*pi - θ
        end
    end


    return a, e_norm, i*RAD2DEG, Ω*RAD2DEG, ω*RAD2DEG, θ*RAD2DEG, h_norm, Ε

end