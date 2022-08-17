""" compute transfer in low thrust using single shooting
"""

using LinearAlgebra
include("structures.jl")
include("propagator.jl")
include("constants.jl")
include("utils.jl")
include("initial_guess.jl")

function transfer_ss(x0::Vector{Float64}, xf::Vector{Float64}, sc::Spacecraft, t::Transfer)
    """ Transfer computed using single shooting """
    #Computation of the initial guess for the time of transfer
    #Computation of the altitude from the Eatth surface
    r0 = norm([x0[1], x0[2], x0[3]])
    rf = norm([xf[1], xf[2], xf[3]])
    Δ_alt = rf-r0

    if Δ_alt < 0
        β = -1*ones(3)
    else
        β = ones(3)
    end

    t_base = 20
    if abs(Δ_alt) < 40
        T_ig = 3600*t_base
    elseif abs(Δ_alt) > 40 && abs(Δ_alt) <100
        T_ig = 3600*(t_base) #+ abs(Δ_alt))
    elseif abs(Δ_alt) > 100 && abs(Δ_alt) <200
        T_ig = 3600*(2*t_base)
    elseif abs(Δ_alt) > 200 && abs(Δ_alt) <300
        T_ig = 3600*(3.5*t_base)
    elseif abs(Δ_alt) > 300 && abs(Δ_alt) <400
        T_ig = 3600*(5*t_base)
    elseif abs(Δ_alt) > 400 && abs(Δ_alt) < 500
        T_ig = 3600*(6*t_base)
    elseif abs(Δ_alt) > 500 && abs(Δ_alt) < 600
        T_ig = 3600*(1*t_base )
    else
        print("Delta alt =", Δ_alt)
    end

    # creation of the free variable vector
    V0 = [x0; T_ig]
    # V0 = [x0; T_ig/2; xf;  -T_ig/2]

    #Definition of the desired final state
    xdf = xf

    #Find the transfer via single shooting
    t.V, t.Δm, t.states = single_shooting(V0, xdf, sc, args = β)
    # t.V, t.Δm, t.states = single_shooting(V0, sc, args = β)

end

function single_shooting(V0, xdf, sc::Spacecraft; tol=1e-8, args =[])
    """ Single shooting correction method"""

    # set iterations values
    i = 1
    i_max = 50
    #Set the desired initial free variable vector
    global V = V0
  
    # Set initial constraint vector and initial STM
    F = ones(9, i_max)*1e-3
    Phi_0 = I+zeros(6,6)

    # Set initial propellant mass
    m_prop0 = sc.mass_prop
    
    while norm(F[:,i])>tol && i<i_max
        
        #Computation of the new final state vector with new x1
        x0 = V[1:6] 
        t_span = V[7]
        m0 = sc.mass_dry + m_prop0
        push!(x0,m0)
        x_init = vcat(x0, reshape(Phi_0, (36, 1)))
        flag_thrust = true
        flag_J2 = true
        flag_STM = true
        global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)

        xf_dot = zeros(size(x_init,1))
        xf = sol.u[end]
        # Compute the engine efficiency
        veff = sc.Isp * g0 * 1000.0
        p = []

        # Define the additonal arguments for the eoms function
        push!(p, flag_thrust, flag_J2, flag_STM, sc.Thrust, veff, args)
        eoms_2BP!(xf_dot, xf,  p, t_span)
        
        #  Update equation
        Phi = (reshape(sol.u[end][end-35:end],(6,6)))
        
        # Phi = STM(x0, t_span)
        @show rank(Phi)
        # I = Phi_0[1:3, 1:6]
        # O = zeros(3,1)
        DF = Matrix{Float64}(undef,6,7)
        DF[1:6,1:6] = Phi 
        DF[1:6, 7] = xf_dot[1:6]
        # DF[7:9, 1:6] =  I
        # DF[7:9, 7] = O
        @show Phi
        # New step 
        i = i+1;
        if i ==20
            print("single shooting failed")
            break
        end
        # Evaluation of the new constraint vector
        diff_stepf = sol.u[end][1:6]-xdf
        # diff_stepi = xf.u[1][1:3]-xd0;
        F[1:6,i] = diff_stepf
        @show F[1:6,i]
        # Evaluation of the new free variables vector
        @show rank(DF)
        global V = V - transpose(DF)*((DF*transpose(DF))\F[1:6,i])
        @show V
    end

    #Compute Δv
    Δm = m_prop0 - sc.mass_prop
    
    return V, Δm, sol.u 

end

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

    T_ig = initial_guess(x0, xf, sc, flags)

    # T_ig = 3600*2.035 #104.747
    @show T_ig
    # Divide the trajectory in arcs of lenght dt
    dt = T_ig/n_arcs 
    # Initialize the free variable vector
    V0 = Vector{Float64}()

    # Create the initial states to put in the free variable vector
    m0_prop = sc.mass_prop
    x_init = push!(x0, sc.mass_dry+sc.mass_prop)

    # Fill the free variable vector V
    if flag_constraints == "xf_all" ||  flag_constraints == "xf_pos"
        i = 1
        while i <= n_arcs
            sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
            push_range!(V0, x_init[1:6])
            push!(V0, dt)
            x_init = sol.u[end]
            i+=1
        end
        #Definition of the desired final state
        xd = xf
    elseif flag_constraints == "xixf_pos"
        #Constrain the initial and final state position
        i = 1
        sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
        push_range!(V0, x_init[4:6])
        push!(V0, dt)
        x_init = sol.u[end]
        i+=1
        while i <= n_arcs 
            sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
            push_range!(V0, x_init[1:6])
            push!(V0, dt)
            x_init = sol.u[end]
            i+=1
        end
        #Definition of the desired final state
        xd = vcat(x0[1:6],xf)
    elseif flag_constraints == "xixf_all"

        #Constrain the initial and final state position and velocity
        i = 1
        sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
        push!(V0, dt)
        x_init = sol.u[end]
        i+=1
        while i <= n_arcs 
            sol = propagate_2BP(x_init, dt, sc,flag_thrust = flags[1], flag_J2 = flags[2], args = β)
            push_range!(V0, x_init[1:6])
            push!(V0, dt)
            x_init = sol.u[end]
            i+=1
        end

        #Definition of the desired final state
        xd = vcat(x0[1:6],xf)

    end
    # Replace the propellant mass in the sc attributes
    sc.mass_prop = m0_prop

    #Find the transfer via single shooting
    # t.V, t.Δm, t.converged = multiple_shooting(V0, xd, n_arcs, sc, flag_constraints, i_max, tol, β, flags)
    #Find the states of the transfer
    # t.converged ? transfer_states(t, sc, n_arcs, β, flags) : print("\n No transfer has been found from correction via multiple shooting")
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
function multiple_shooting(V0, xd, n_arcs::Int64, sc::Spacecraft, flag_constraints::String, i_max::Int64, tol, args, flags)
    """ Multiple shooting correction method"""

    #Set iterations value
    i = 1

    # Set the desired initial free variable vector
    global V = V0
    global flag_convergence = true
    # Set initial constraint vector and initial STM
    if flag_constraints == "xf_all" || flag_constraints == "xixf_all"
        # Constrain the final state to be equal in position and velocity to xdf
        F = ones(6*n_arcs)*1e-3
    elseif flag_constraints == "xf_pos" || flag_constraints == "xixf_pos"
        # Constrain the final state to be equal just in position to xdf
        F = ones(6*n_arcs-3)*1e-3
    end 

    # Set initial propellant mass
    m_prop0 = sc.mass_prop


    while norm(F)>tol && i<i_max

        DF, F = compute_jacobian(V, F, sc, n_arcs, xd, m_prop0, flag_constraints, args, flags)     
        @show size(DF), rank(DF)
        @show size(F), norm(F)  
        

        # Update equation
        # Evaluation of the new free variables vector
        global V = V - transpose(DF)*((DF*transpose(DF))\F)
        # New step 
        i +=1
        if i == i_max
            print("multiple shooting failed\n")
            flag_convergence = false
            break
        end
    end

    @show i
    #Compute Δv
    Δm = m_prop0 - sc.mass_prop
    
    # Replace spacecraft mass with the initial value to allow a correct evaluation during the propagation in the transfer_states function
    sc.mass_prop = m_prop0
    if flag_constraints =="xixf_pos" 
        pushfirst_range!(V, xd[1:3])
    elseif flag_constraints =="xixf_pos" 
        pushfirst_range!(V, xd[1:6])
    end
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
function transfer_states(t::Transfer, sc::Spacecraft, n_arcs::Int64, args, flags)
""" Function that propagates the states of the transfer from the converged multiple shooting free variable vector V"""

    # Initialize the position elements vectors 
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []

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


        #Indexes for vector computation
        idx_c = six_vec + i_vec*7*(h-1)

        ## Computation of the state with forward propagation from the departing point
        x0 = t.V[idx_c] 
        t_span = t.V[7*h]

        push!(x0,m0)

        sol = propagate_2BP(x0, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, args = args)

        for i = 1:size(sol.u,1)
            push!(x, sol.u[i][1])
            push!(y, sol.u[i][2])
            push!(z, sol.u[i][3])
            push!(vx, sol.u[i][4])
            push!(vy, sol.u[i][5])
            push!(vz, sol.u[i][6])
        end
    end
    
    states = hcat(x,y,z,vx,vy,vz)

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
function compute_jacobian(V, F, sc, n_arcs::Int64, xd, m_prop0, flag_constraints::String, args, flags)
    """ Compute the matrix of the partial derivatives of 
    the constraint vector F (vector of the differences) wrt
    the control variables in the vector V using STM matrix""" 


    # Set the initial STM matrix to be the identitiy matrix 6x6
    Phi_0 = I+zeros(6,6)
    
    #Vector for index propagation          
    i_vec = Vector{Int64}(undef,6)
    fill!(i_vec,1)

    # Replace the initial propellant mass at every iteration
    sc.mass_prop = m_prop0

    # Compute the engine efficiency
    veff = sc.Isp * g0 * 1000.0

    # Define the additonal arguments for the eoms function
    p = []
    flag_thrust = flags[1]
    flag_J2 = flags[2]
    flag_STM = flags[3]
    push!(p, flag_thrust, flag_J2, flag_STM, sc.Thrust, veff, args)

    x0 = zeros(6)


    if flag_constraints == "xf_all"
        
        six_vec = [1,2,3,4,5,6]
        # Initialize DF Matrix
        global DF = zeros(6*n_arcs, 7*n_arcs) 
        
        for h = 1:n_arcs
            # Use the STM of the state vector to fill the jacobian matrix DF
            #Indexes for vector computation
            idx_r = six_vec + i_vec*6*(h-1)
            idx_c = six_vec + i_vec*7*(h-1)
            idx_c2 = six_vec + i_vec*7*h

            ## Computation of the state with forward propagation from the departing point
            x0 = V[idx_c] 
            t_span = V[7*h]
            
            m0 = sc.mass_dry + sc.mass_prop
            push!(x0,m0)
            x_init = vcat(x0, reshape(Phi_0, (36, 1)))
            global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)
            
            ## Computation of the derivative of the final state of the forward propagation
            xf_dot = zeros(size(x_init,1))
            xf = sol.u[end]
            eoms_2BP!(xf_dot, xf,  p, t_span)
        
            # Computation of the STM for forward and backward propagation
            Phi = (reshape(sol.u[end][end-35:end],(6,6)))
            ## Computation of DF matrix (derivative of F wrt V)
            global DF[idx_r, idx_c] = Phi
            global DF[idx_r, 7*h] = xf_dot[1:6]
            if h != n_arcs
                global DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
            else           
                diff_step = sol.u[end][1:6]-xd[1:6]
            end

            #Compute the constraint vector
            global F[idx_r] = diff_step
        end
            
    elseif flag_constraints == "xf_pos"

        six_vec = [1,2,3,4,5,6]
        # Initialize DF Matrix
        DF = zeros(6*(n_arcs)-3, 7*(n_arcs))

        for h = 1:n_arcs
            # Use the STM of the state vector to fill the jacobian matrix DF
            #Indexes for vector computation
            idx_r = six_vec + i_vec*6*(h-1)
            idx_c = six_vec + i_vec*7*(h-1)
            idx_c2 = six_vec + i_vec*7*h

            ## Computation of the state with forward propagation from the departing point
            x0 = V[idx_c] 
            t_span = V[7*h]
            
            m0 = sc.mass_dry + sc.mass_prop
            push!(x0,m0)
            x_init = vcat(x0, reshape(Phi_0, (36, 1)))
            global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)
            
            ## Computation of the derivative of the final state of the forward propagation
            xf_dot = zeros(size(x_init,1))
            xf = sol.u[end]
            eoms_2BP!(xf_dot, xf,  p, t_span)
        
            # Computation of the STM for forward and backward propagation
            Phi = (reshape(sol.u[end][end-35:end],(6,6)))

            # Constrain just the final position
            if h != n_arcs
                DF[idx_r, idx_c] = Phi
                DF[idx_r, 7*h] = xf_dot[1:6]
                DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
                F[idx_r] = diff_step
            else      
                DF[idx_r[1:3], idx_c] = Phi[1:3,1:6]
                DF[idx_r[1:3], 7*h] = xf_dot[1:3]                     
                diff_step = sol.u[end][1:3]-xd[1:3]
                #Compute the constraint vector
                F[idx_r[1:3]] = diff_step
            end


        end

    elseif flag_constraints == "xixf_pos"
        # Initialize DF Matrix
        DF = zeros(6*(n_arcs)-3, 7*(n_arcs)-3) # contrain the initial and final position to be fixed

        six_vec0 = [1,2,3,4,5,6] 
        six_vec1 = [5,6,7,8,9,10] - i_vec*7

        for h = 1:n_arcs
            #Indexes for vector computation
            idx_r = six_vec0 + i_vec*6*(h-1)
            idx_c = six_vec1 + i_vec*7*(h-1)
            idx_c2 = six_vec0 + i_vec*6*(h-1)

            
            ## Computation of the state with forward propagation from the departing point
            if h == 1 
                x0[1:3] = xd[1:3]
                x0[4:6] = V[1:3]
                t_span = V[4]
            else
                x0 = V[idx_c]
                t_span = V[idx_c[end]+1]
            end
            
            m0 = sc.mass_dry + sc.mass_prop
            push!(x0,m0)

            x_init = vcat(x0, reshape(Phi_0, (36, 1)))
            global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)
            
            ## Computation of the derivative of the final state of the forward propagation
            xf_dot = zeros(size(x_init,1))
            xf = sol.u[end]
            eoms_2BP!(xf_dot, xf,  p, t_span)
        
            # Computation of the STM for forward and backward propagation
            Phi = (reshape(sol.u[end][end-35:end],(6,6)))

            ## Computation of DF matrix (derivative of F wrt V)
            # Constraint just the initial and final position and not the velocity
            if h ==1
                DF[idx_r, idx_c2[1:3]] = Phi[1:6, 4:6]
                DF[idx_r, 4*h] = xf_dot[1:6]
                DF[idx_r, idx_c+i_vec*7] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c+i_vec*7]
                F[idx_r] = diff_step
            elseif h > 1 && h < n_arcs
                DF[idx_r, idx_c] = Phi
                DF[idx_r, idx_c[end]+1] = xf_dot[1:6]
                DF[idx_r, idx_c+i_vec*7] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c+i_vec*7]
                F[idx_r] = diff_step
            elseif h == n_arcs
                DF[idx_r[1:3], idx_c] = Phi[1:3,1:6]
                DF[idx_r[1:3], idx_c[end]+1] = xf_dot[1:3]   
                diff_step = sol.u[end][1:3]-xd[7:9]
                F[idx_r[1:3]] = diff_step
            end
            
        end
    
    elseif flag_constraints == "xixf_all"

        six_vec = [1,2,3,4,5,6]
        six_vec1 = [2,3,4,5,6,7] - i_vec*7
        # Initialize DF Matrix
        DF = zeros(6*n_arcs, 7*n_arcs - 6) 
        
        for h = 1:n_arcs
            # Use the STM of the state vector to fill the jacobian matrix DF
            #Indexes for vector computation
            idx_r = six_vec + i_vec*6*(h-1)
            idx_c = six_vec1 + i_vec*7*(h-1)
            idx_c2 = six_vec1 + i_vec*7*h

            ## Computation of the state with forward propagation from the departing point
            if h == 1 
                x0[1:6] = xd[1:6]
                t_span = V[1]
            else
                x0 = V[idx_c]
                t_span = V[idx_c[end]+1]
            end
            
            m0 = sc.mass_dry + sc.mass_prop
            push!(x0,m0)
            x_init = vcat(x0, reshape(Phi_0, (36, 1)))
            global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)
            
            ## Computation of the derivative of the final state of the forward propagation
            xf_dot = zeros(size(x_init,1))
            xf = sol.u[end]
            eoms_2BP!(xf_dot, xf,  p, t_span)
        
            # Computation of the STM for forward and backward propagation
            Phi = (reshape(sol.u[end][end-35:end],(6,6)))
            
            ## Computation of DF matrix (derivative of F wrt V)
            # Constraint just the initial and final position and velocity
            if h ==1
                DF[idx_r, h] = xf_dot[1:6]
                DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
                F[idx_r] = diff_step
            elseif h > 1 && h < n_arcs
                DF[idx_r, idx_c] = Phi
                DF[idx_r, idx_c[end]+1] = xf_dot[1:6]
                DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
                F[idx_r] = diff_step
            elseif h == n_arcs
                DF[idx_r, idx_c] = Phi
                DF[idx_r, idx_c[end]+1] = xf_dot[1:6]   
                diff_step = sol.u[end][1:6]-xd[7:12]
                F[idx_r] = diff_step
            end

        end

    end

    return DF, F
end

# - flag_fd = true if partial derivatives have to be computed using finite difference instead of using
    #  analytical derivates.
function compute_jacobian_fd(V, u, F, n_arcs, sc)
""" Computation of the jacobian of the constant vector wrt the free variable vector V and the control vector u using finite difference"""

    # Replace the initial propellant mass at every iteration
    sc.mass_prop = m_prop0

    # Compute the engine efficiency
    veff = sc.Isp * g0 * 1000.0

    # Define the additonal arguments for the eoms function
    p = []
    flag_thrust = flags[1]
    flag_J2 = flags[2]
    flag_STM = flags[3]
    push!(p, flag_thrust, flag_J2, flag_STM, sc.Thrust, veff, args)

    # Initialize DF Matrix
    DF = zeros(6*n_arcs, 7*n_arcs - 6)
    
    for i = 1:n_arcs

        Vu_nom = vcat(V[i:i+1,], u[:,i:i+1])
    end 

end

function trajOpt()
    """ Function that optimize the free variable vector V, the control vector u and the integration time"""

end