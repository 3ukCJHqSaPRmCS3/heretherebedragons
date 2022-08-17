
include("sims_flanagan_setup.jl")
include("propagators.jl")
include("utils.jl")
include("planets.jl")

# Returns the nondimensional mismatch equality constraints of the arrival boundary conditions.
# This method propagates the nondimensional dynamics of the spacecraft
# from the departure time `t0` to the arrival time `tf`, then evaluates
# the nondimensional mismatch between the desired arrival state `xf`,
# arrival mass costate `lmf == 0` (`if freemass == True`), and arrival
# Hamiltonian `H == 0` (`if freetime == True`) and the supplied desired
# arrival boundary conditions.

function equality_constraint(l::Leg, atol = 1e-5, rtol = 1e-5 )

    #Check that atol and rtol have been set
    for tol in [atol, rtol]
        if typeof(tol) != int || typeof(tol) != float
            ErrorException("Both atol and rtol must be supplied as instances of either float or int")
        end
    end

    for atr in [l.start_time, l.x0, l.end_time, l.xf]
        if hasproperty(l, atr) == false
            ErrorException("Cannot propagate dynamics, as boundary conditions t0, x0, l0, tf, and xf have not been set. Use set(t0, x0, l0, tf, xf) to set boundary conditions.")
        else
            atol = float(atol)
            rtol = float(rtol)
        end
    end

    # Propagate the state forward
    states_fwd = propagate_2BPt(l, rtol, atol)

    # Compute the propagated arrival states
    prf = l.trajectory[end, 1:3]
    pvf = l.trajectory[end, 4:6]
    
    # Position and velocity arrival mismatch
    drf = prf - l.xf[1:3] 
    dvf = pvf - l.xf[4:6]

    # free arrival mass
    if l.freemass == true
        lmf = l.trajectory[end, end]
    else
        dmf = l.trajectory[end, 7] - l.spacecraft.mass
    end


    # free arrival time 
    if l.freetime == true
        # Compute the Hamiltonian - not for sims flanagan but present for pontryagin
    end

    # Create equality constraints
    if (l.freemass && l.freetime)
        ceq = hcat(drf, dvf,[lmf], [hf] )
    elseif (l.freemass==true && l.freetime==false)
        ceq = hcat(drf, dvf, [lmf])
    elseif (l.freetime==true && l.freemass==false)
        ceq = hcat(drf, dvf, [dmf], [Hf])
    elseif (l.freetime==false && l.freemass==false)
        ceq = hcat(dfr, dvf, [dmf])
    else
        ErrorException("Could not determine equality constraint vector")
    end

    return ceq
end




function inequality_constraint(l::Leg)

    for 

end