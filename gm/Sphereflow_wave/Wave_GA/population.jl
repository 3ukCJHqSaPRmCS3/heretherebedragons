""" creation of the population for the genetic algorithm optimization
"""

include("Orbits.jl")
include("structures.jl")



function population(population::Population, constA::Constellation, constB::Constellation)

    # Plot the orbits and find the initial and final states to create the population
    r_A = zeros(100,3, constA.n_orbit)
    v_A = zeros(100,3, constA.n_orbit)
    r_B = zeros(100,3, constB.n_orbit)
    v_B = zeros(100,3, constB.n_orbit)
    x0_A = zeros(constA.n_orbit, 6)
    xf_B = zeros(constB.n_orbit, 6)


    for i = 1:constA.n_orbit
        # Transforn the orbit in cartesian states to plot it
        orb_A = constA.orbits[i]
        r_A[:,:,i],v_A[:,:,i] = kep2cart(orb_A, false)

        #find the initial state
        r, v = kep2cart(orb_A, true)
        x0_A[i,:] = hcat(r,v)
        constA.spacecrafts[i].state = x0_A[i,:]
    end

    for i = 1:constB.n_orbit
        # Transforn the orbit in cartesian states to plot it
        orb_B = constB.orbits[i]
        r_B[:,:,i],v_B[:,:,i] = kep2cart(orb_B, false)

        #find the initial state
        r, v = kep2cart(orb_B, true)
        xf_B[i,:] = hcat(r,v)
        constB.spacecrafts[i].state = xf_B[i,:]
    end

    
    # Creation of population by associating the fist states of each orbits from const A to the ones in const B
    population.list_states = hcat(x0_A,xf_B)


    return r_A, v_A, r_B, v_B, x0_A, xf_B

end

