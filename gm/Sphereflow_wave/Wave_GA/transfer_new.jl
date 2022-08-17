""" compute transfer in low thrust using single shooting
"""

using LinearAlgebra
include("structures.jl")
include("propagator.jl")
include("constants.jl")
include("utils.jl")

# x0: vector of the starting state of the transfer 1x6, pos in km and vel in km/s
# xf: vector of the final state of the transfer 1x6, pos in km and vel in km/s
# sc: the Sacecraft of the constelation incvolved in the transfer
# t: the transfer structure defined in the main function
# n_arcs: the number of arcs used to divide the transfer in smaller segment
# flag_constraints = String that indicates if any constraint has been applied to any
    #  of the states of the arcs that represent the trajectory.  Type:
    # "xf_all" if the final state position and velocity of the final arc has been 
    #  constrained to be equal to the desired final state xd. 
    #  "xf_pos" if only the position of ther last state of the transfer is constrained to match the 
    #  position of xd ( this implieas that an impulsive burn is occurring at the end of the transfer).
    #  "xixf_pos" if both the initial and final position of the transfer have to match desiderd states 
    #  expressed in the xd vector.
# i_max: maximum number of iteration in the multiple shooter
# flags for additional arguments in the STM computation, the order is:
    # 1) flag_thrust
    # 2) flag_J2
    # 3) flag_STM
    # default are true
function transfer_ms(x0::Vector{Float64}, xf::Vector{Float64}, sc::Spacecraft, t::Transfer, n_arcs::Int64, flag_constraints::String; i_max = 50, tol= 1e-8, flags = [true, true, true])
    """ Transfer computed using multiple shooting """

   
    #Computation of the altitude from the Earth surface
    r0 = norm([x0[1], x0[2], x0[3]])
    rf = norm([xf[1], xf[2], xf[3]])
    Δ_alt = rf-r0
    # Setting the throttle ## THIS CAN BE MODIFIED TO BECOME A CONTROL VARIABLE FOR OPT CTRL
    if Δ_alt < 0
        β = -1*ones(3)
    else
        β = ones(3)
    end

    T_ig = 3600*119.67 #104.747

    # Divide the trajectory in arcs of lenght dt
    dt = T_ig/n_arcs 
    # Initialize the free variable vectors X, U and T
    X = []
    U = []
    T = []

    # Create the initial states to put in the free variable vector
    m0_prop = sc.mass_prop
    x_init = push!(x0, sc.mass_dry+sc.mass_prop)

    # Fill the free variable vectors X, U and T in the case constraint is set to "xixf_all" meaning constrain initial and final position and velocity

    #Constrain the initial and final state position and velocity
    i = 1
    while i <= n_arcs 
        sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
        push!(X, x_init[1:6])
        push!(T, dt)
        push!(U, β)
        x_init = sol.u[end]
        i+=1
    end

    #Definition of the desired final state
    xd = vcat(x0[1:6],xf)

    # Replace the propellant mass in the sc attributes
    sc.mass_prop = m0_prop

    #Find the transfer via single shooting
    t.V, t.Δm, t.converged = multiple_shooting(X, U, T, xd, n_arcs, sc, flag_constraints, i_max, tol, T_ig, flags)
    #Find the states of the transfer
    t.converged ? transfer_states(t, sc, n_arcs, β, flags) : print("\n No transfer has been found from correction via multiple shooting")
end

# V0: vector of the free variables, usually defined as V = [xi_0, Δt0, xi_1, Δt1, ..., xi_n, Δtn], where xi_j is the initial state 
    # of each arc of the transfer traj. and Δti is the relative propagation time. Pos in km and vel in km/s and time in s
# xd: initial or final constraint for the transfer traj. Can be the desired final state, the desired initial state or both. pos in km and vel in km/s
# n_arcs: the number of arcs used to divide the transfer in smaller segment
# sc: the Sacecraft of the constelation incvolved in the transfer
# flag_constraints = String that indicates if any constraint has been applied to any
    #  of the states of the arcs that represent the trajectory.  Type:
    # "xf_all" if the final state position and velocity of the final arc has been 
    #  constrained to be equal to the desired final state xd. 
    #  "xf_pos" if only the position of ther last state of the transfer is constrained to match the 
    #  position of xd ( this implieas that an impulsive burn is occurring at the end of the transfer).
    #  "xixf_pos" if both the initial and final position of the transfer have to match desiderd states 
    #  expressed in the xd vector.
# i_max: maximum number of iteration in the multiple shooter
# tol: tolerance for the nom of the constraint vector F to terminate the multiple shooting correction. Default to 1e-7
# args: additional arguments coming from the trasfer_ms function
# flags for additional arguments in the STM computation, the order is:
    # 1) flag_thrust
    # 2) flag_J2
    # 3) flag_STM
    # default are true
function multiple_shooting(X0, U0, T0, xd, n_arcs::Int64, sc::Spacecraft, flag_constraints::String, i_max::Int64, tol, T_ig, flags)
    """ Multiple shooting correction method"""

    #Set iterations value
    i = 1

    pert = 1e-6


    # Set the desired initial free variable vector
    global X = X0
    global U = U0
    global T = T0
    global flag_convergence = true
    Xd_f = [xd[7:12], [0,0,0], T_ig]
    # Set initial constraint vector and initial STM
    # Constrain the final state to be equal in position and velocity to xdf
    error = 1
    # Set initial propellant mass
    m_prop0 = sc.mass_prop

    #Compose V as [X,U,T] for every initial point at every arc
    V = []

    while error>tol && i<i_max

        DF= compute_jacobian_fd(X, U, T, F, Xd_f, n_arcs, n_steps, sc, flags, pert)     
        @show size(DF), rank(DF)

        # Optimize QP subproblem
        (x_update, u_update, tf_update, status, cost) =
        optimizeTraj(X, U, T, F, DF, n_states, n_arcs, X0_times, X0_states, Xf_times, Xf_states, flagEnd_temp, β)


        ############ LINE SEARCH
        alpha = 1.0
        #Call line search:
        if iterCount > 10
            alpha = lineSearch(X, x_update, U, u_update, T, n_states, n_nodes, nsteps, sc)
        end

        # Update equation
        X = X + x_update * alpha
        U = U + u_update * alpha
        T = T + tf_update * alpha

        # New step 
        i +=1
        if i == i_max
            print("multiple shooting failed\n")
            flag_convergence = false
            break
        end

        #check Update
        F = compute_defects(X,U, T, n_states, n_arcs, n_steps, sc, flags)

        #print progress report
        error = norm(defect[:], Inf)

        @show error
    end
    #Create the V vector
    V = hcat(X, U, T)

    @show i
    #Compute Δv
    Δm = m_prop0 - sc.mass_prop
    
    # Replace spacecraft mass with the initial value to allow a correct evaluation during the propagation in the transfer_states function
    sc.mass_prop = m_prop0

    return V, Δm, flag_convergence

end



# t: the transfer structure defined in the main function
# sc: the Sacecraft of the constelation incvolved in the transfer
# n_arcs: the number of arcs used to divide the transfer in smaller segment
# args: additional arguments coming from the trasfer_ms function
# flags for additional arguments in the STM computation, the order is:
    # 1) flag_thrust
    # 2) flag_J2
    # 3) flag_STM
    # default are true
function transfer_states(t::Transfer, sc::Spacecraft, n_arcs::Int64, flags)
""" Function that propagates the states of the transfer from the converged multiple shooting free variable vector V"""

    # Initialize the position elements vectors 
    x = []
    y = []
    z = []

    # Define the additonal arguments for the eoms function
    flag_thrust = flags[1]
    flag_J2 = flags[2]

    #Vector for index propagation          
    i_vec = Vector{Int64}(undef,6)
    fill!(i_vec,1)
    six_vec = [1,2,3,4,5,6]

    # Spcecraft mass
    m0 = sc.mass_dry + sc.mass_prop

    for h = 1:n_arcs

        ## Computation of the state with forward propagation from the departing point
        x0 = t.V[h][:]
        c0 = t.V[n_arcs+h][:]
        t_span = t.V[2*n_arcs+h]

        push!(x0,m0)

        sol = propagate_2BP(x0, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, args = c0)

        for i = 1:size(sol.u,1)
            push!(x, sol.u[i][1])
            push!(y, sol.u[i][2])
            push!(z, sol.u[i][3])
        end
    end
    
    states = hcat(x,y,z)

    t.states = states

end


# V: vector of the free variables, usually defined as V = [xi_0, Δt0, xi_1, Δt1, ..., xi_n, Δtn], where xi_j is the initial state 
    # of each arc of the transfer traj. and Δti is the relative propagation time. Pos in km and vel in km/s and time in s
# F: constraint vector, usually defined as F = [xf_0 - xi_1, xf_1 - xi_2, ...,  xf_n - xd], where xi_j is the initial state 
    # of each arc of the transfer traj. and xf_j is the final state of each arc after propagation. Pos in km and vel in km/s and time in s
# xd: initial or final constraint for the transfer traj. Can be the desired final state, the desired initial state or both. pos in km and vel in km/s
# flag_constraints = String that indicates if any constraint has been applied to any
    #  of the states of the arcs that represent the trajectory.  Type:
    # "xf_all" if the final state position and velocity of the final arc has been 
    #  constrained to be equal to the desired final state xd. 
    #  "xf_pos" if only the position of ther last state of the transfer is constrained to match the 
    #  position of xd ( this implieas that an impulsive burn is occurring at the end of the transfer).
    #  "xixf_pos" if both the initial and final position of the transfer have to match desiderd states 
    #  expressed in the xd vector.
# args: additional arguments coming from the trasfer_ms function
# flags for additional arguments in the STM computation, the order is:
    # 1) flag_thrust
    # 2) flag_J2
    # 3) flag_STM
    # default are true

function compute_jacobian_fd(X, U, T, F, Xd_f, n_arcs::Int64, n_steps::Int64, sc::Spacecraft, flags, pert)
""" Computation of the jacobian of the constant vector wrt the free variable vector V and the control vector u using finite difference"""

    #Initialize the jacobian as [size(F), size(V)]
    n_states = size(X[1],1)
    n_controls = size(U[1],1)
    n_var = 2*(n_states+n_controls+1)
    DF_temp = zeros(n_arcs*n_states,n_var)
    xd = Xd_f[1][:]
    ud = Xd_f[2][:]
    td = Xd_f[3]

    for i = 1:n_arcs-1

        V_nom = vcat(X[i][:], X[i+1][:], U[i][:], U[i+1][:], T[i][:], T[i+1][:])

        for  j = 1:n_var

            #Create the modified arrays
            V_mod = copy(V_nom)
            V_mod[j] = V_mod[j] + pert

            #reshape
            X_temp = reshape(V_mod[1:(n_states*2)], n_states, 2)
            U_temp = reshape(V_mod[(n_states*2)+1 : (n_states*2)+3], 3, 2)
            T_temp = reshape(V_mod[(n_states*2)+4 : end], 1, 2)
            # Compute the single defect between the two affected nodes:
            F_temp = compute_defects(X_temp, U_temp, T_temp, n_states, n_nodes, n_steps, sc, flags)

            # compute the index of the row that corresponds to this defect
            idx = (i-1)*n_states+1 : i*n_states
            DF_temp[idx,j] = (F_temp - F[:,i]) / pert
            
        end
    end

    # Compute the jacobian for the last state
    V_nomf = vcat(X[end][:], xd[7:12], U[end][:], ud[4:6], T[end][:], td)

    for  j = 1:n_var

        #Create the modified arrays
        V_modf = copy(V_nomf)
        V_modf[j] = V_modf[j] + pert

        #reshape
        X_temp = reshape(V_modf[1:(n_states)], n_states, 2)
        U_temp = reshape(V_modf[(n_states*2)+1 : (n_states*2)+3], 3, 2)
        T_temp = reshape(V_modf[(n_states*2)+4 : end], 1, 2)
        # Compute the single defect between the two affected nodes:
        F_temp = compute_defects(X_temp, U_temp, T_temp, n_states, n_nodes, n_steps, sc, flags)

        # compute the index of the row that corresponds to this defect
        idx = (n_arcs-1)*n_states+1 : n_arcs*n_states
        DF_temp[idx,j] = (F_temp - F[:,n_arcs]) / pert
        
    end

    # Recreate the full Jacobian
    # band diagonal:
    DF_full = zeros(n_states*(n_nodes), n_arcs*(n_var))

    for ind in axes(DF_full,1) #loop through rows
        #Which node's defect are we on?
        idx = Int32(ceil(ind/(n_states)))

        #offset (in # of columns) for the state-related ones
        idx2 = (idx-1)*(n_states)

        #offset (in # of columns) for the control-related ones
        idx3 = n_states*n_arcs + (idx-1)*(3)

        # offset (in # of columns) for the time-related ones
        idx4 = n_states*n_arcs + n_controls*n_arcs + (idx-1)

        #set the appropriate columns of DF_full from DF_temp:
        #state-related terms:
        DF_full[ind, (idx2+1):(idx2+2*n_states)] = DF_temp[ind,1:2*n_states]
        #control-related terms:
        DF_full[ind, (idx3+1):(idx3+2*3)] = DF_temp[ind,(2*n_states+1):(2*n_states+3)]
        #time-related terms:
        DF_full[ind, (idx4+1)] = DF_temp[ind,(2*n_states+4)]
    end
    

    #outputs:
    return DF_full
end

function compute_defects(X, U, T, n_states::Int64, n_arcs::Int64, n_steps::Int64, sc::Spacecraft, flags)
    #Calculate defect constraint violations

    #pre-allocate defect:
     F = zeros(n_states, n_arcs-1)

    for i = 1:n_arcs-1
        #propagate node i forward from time i to time between (i+1) and (i)

        ########### propagate forward:
        x_init = copy(X[:,i])
        u = copy(U[:,i])
        t_span = T[i]
        sol =  propagate_2BP(x_init, t_span, sc, flag_thrust = flags[1], flag_J2 = flags[2], args = u)


        ########## defect between forward and backward results:
        x_t1 = copy(X[:,i+1])
        F[:,i] = sol.u[end][1:6] - x_t1
    end

    #outputs:
    return F
end



function trajOpt(X, U, T, F, DF, n_states, n_arcs, X0_times, X0_states, Xf_times, Xf_states)
    """ Function that optimize the free variable vector V, the control vector U and the integration time"""
    m = Model(solver=IpoptSolver(print_level=0))

    #State update variables:
    @variable(m, X_jump[1:length(X)] )
    #constrain initial mass:
    if n_states == 7
        @constraint(m, (X_jump[7] + X[7, 1]) == mass)
    end

    #control update variables:
    #limits on u_jump have to be determined by u_all
    @variable(m, U_jump[1:length(U)])

    #time of flight variable:
    d = 0.0 #freeze endpoints

    @variable(m, -d <= tf_jump <= d)
    @constraint(m, (tf_TU + tf_jump) <= 10*24*3600) #absolute max bound
    @constraint(m, (tf_TU + tf_jump) >= 0.0 ) #absolute min bound

    #impulsive maneuvers at endpoints:
    @variable(m, dV1_jump[1:3])
    @variable(m, dV2_jump[1:3])
    if !allowImpulsive
        @constraint(m, dV1_jump .== 0.) #disallow impulsive maneuvers
        @constraint(m, dV2_jump .== 0.) #
    end

    #Thrust magnitude constraints:
    #
    #Note: Thrust magnitude constraint is disabled because it was found
    #more effective to constrain thrust with the indirect method (used after
    #converging on the direct solution).
    #
    #Since u_jump is just the *update* to control, the value we want to actually
    #limit to be within constraints is (u_jump + u_all).
    #Create a vector of symbolic quadratic constraints:
    # umag = (u_jump[1:3:end] + u_all[1:3:end]).^2 +
    #     (u_jump[2:3:end] + u_all[2:3:end]).^2 +
    #     (u_jump[3:3:end] + u_all[3:3:end]).^2
    # @constraint(m, umag .<= thrustLimit^2 )


    #Need 'dt' so that unequal time steps gives the right result. If we
    #don't use 'dt' like this, then high thrust for a long segment is
    #weighted equally with high thrust for a short segment.
    t_TU_fixed = t0_TU + (tau+1)/2*(tf_TU-t0_TU) #transform tau -> t
    dt = diff(t_TU_fixed)
    dt_temp = vcat(dt/2, dt[end]/2) + vcat(0, dt[1:end-1]/2, 0)
    dt_repeated = ones(3,1) * dt_temp' #same size as control

    #replace this by actual mass vector (repeated like dt)
    if nstate == 7
        mass_repeated = repmat(X[7,:], 3, 1)
    else
        mass_repeated = 1000.0 * ones(size(u_all))
    end

    #Dynamics constraints:
    #constrain bringing the defects to zero (according to linearization)
    @constraint(m, -DF*vcat(X_jump, u_jump, tf_jump) .== F[:])

    #Constraints for endpoints:

    #Finite difference derivatives of endpoints wrt parameters
    pert = 0.05 #needs to be relatively large because the endpoint orbits are only specified at 100 points (so, every 0.01 tau)
    (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)
    (state_0_mod1, state_f_mod1) = interpEndStates(τ1+pert, τ2+pert, X0_times, X0_states, Xf_times, Xf_states, MU)
    (state_0_mod2, state_f_mod2) = interpEndStates(τ1-pert, τ2-pert, X0_times, X0_states, Xf_times, Xf_states, MU)
    dstate0_dp1 = (state_0_mod1 - state_0_mod2) / (2*pert)
    dstatef_dp2 = (state_f_mod1 - state_f_mod2) / (2*pert)
    #2nd derivatives:
    ddstate0_dp1 = (state_0_mod1 - 2*state_0 + state_0_mod2) / (pert^2)
    ddstatef_dp2 = (state_f_mod1 - 2*state_f + state_f_mod2) / (pert^2)


    if flagEnd
        #linear:
        linearEnd1 = (state_0 + dstate0_dp1*(p1_jump))
        linearEnd2 = (state_f + dstatef_dp2*(p2_jump))
        @constraint(m, (X_all[1:6, 1]   + X_jump[1:6]                     + vcat(zeros(3), dV1_jump+dV1) ) - linearEnd1 .== 0 )
        @constraint(m, (X_all[1:6, end] + X_jump[(end-nstate+1):end][1:6] + vcat(zeros(3), dV2_jump+dV2) ) - linearEnd2 .== 0 )

        #quadratic:
        quadEnd1 = (state_0 + dstate0_dp1*(p1_jump) + ddstate0_dp1 / 2 * (p1_jump)^2 )
        quadEnd2 = (state_f + dstatef_dp2*(p2_jump) + ddstatef_dp2 / 2 * (p2_jump)^2 )
        c1_quad = (X_all[1:6,1]   + X_jump[1:6])                     - (state_0 + dstate0_dp1*(p1_jump) + ddstate0_dp1 / 2 * (p1_jump)^2 )
        c2_quad = (X_all[1:6,end] + X_jump[(end-nstate+1):end][1:6]) - (state_f + dstatef_dp2*(p2_jump) + ddstatef_dp2 / 2 * (p2_jump)^2 )

        #Difference from linear:
        quadCost = norm(ddstate0_dp1) * 1/2 * (p1_jump)^2 +
                    norm(ddstatef_dp2) * 1/2 * (p2_jump)^2

    else #hard fix endpoints
        quadCost = 0.
        (state_0, state_f) = interpEndStates(τ1, τ2, X0_times, X0_states, Xf_times, Xf_states, MU)

        @constraint(m, (X_all[1:6, 1]   + X_jump[1:6]                     + vcat(zeros(3), dV1_jump+dV1) ) - state_0 .== 0 )
        @constraint(m, (X_all[1:6, end] + X_jump[(end-nstate+1):end][1:6] + vcat(zeros(3), dV2_jump+dV2) ) - state_f .== 0 )
    end

    costEnd = sum( ( (dV1 + dV1_jump) * DU/TU ).^2) + sum( ( (dV2 + dV2_jump) * DU/TU ).^2)

    #quadratic cost, with quadratic endpoint constraint violation cost added:
    cost = sum( ((u_all[:] + u_jump).^2) .* dt_repeated[:] ) + β * quadCost + costEnd

    @objective(m, Min, cost)

    #Run QP solver:
    status = JuMP.solve(m)

    x_update = reshape(getvalue(X_jump), nstate, n_nodes)
    u_update = reshape(getvalue(u_jump), 3, n_nodes)

    p1_update = getvalue(p1_jump)
    p2_update = getvalue(p2_jump)
    tf_update = getvalue(tf_jump)
    dV1_update = getvalue(dV1_jump)
    dV2_update = getvalue(dV2_jump)

    tf_TU = tf_TU + getvalue(tf_jump)

    cost = getvalue(cost)

    #Outputs:
    (x_update, u_update, p1_update, p2_update, tf_update, dV1_update, dV2_update, status, cost, dstate0_dp1, dstatef_dp2, ddstate0_dp1, ddstatef_dp2)
end