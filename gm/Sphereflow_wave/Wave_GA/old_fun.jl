function multiple_shooting(V0, xdf, n_arcs::Int64, sc::Spacecraft, flag_constraints::String; tol=1e-8, args =[])
    """ Multiple shooting correction method"""

    #Set iterations values
    i = 1
    i_max = 50

    # Set the desired initial free variable vector
    global V = V0
    global flag_convergence = true
    # Set initial constraint vector and initial STM
    if flag_constraints == "xf_all"
        # Constrain the final state to be equal in position and velocity to xdf
        F = ones(6*n_arcs)*1e-3
    elseif flag_constraints == "xf_pos"
        # Constrain the final state to be equal just in position to xdf
        F = ones(6*(n_arcs-1)+3)*1e-3
    end 
    # Set the initial STM matrix to be the identitiy matrix 6x6
    Phi_0 = I+zeros(6,6)


    # Set initial propellant mass
    m_prop0 = sc.mass_prop
    
    #Vector for index propagation          
    i_vec = Vector{Int64}(undef,6)
    fill!(i_vec,1)
    six_vec = [1,2,3,4,5,6]

    while norm(F)>tol && i<i_max

        # Replace the initial propellant mass at every iteration
        sc.mass_prop = m_prop0

        # Compute the engine efficiency
        veff = sc.Isp * g0 * 1000.0
        p = []
        # Define the additonal arguments for the eoms function
        flag_thrust = true
        flag_J2 = true
        flag_STM = true
        push!(p, flag_thrust, flag_J2, flag_STM, sc.Thrust, veff, args)

        # Initialize DF Matrix
        DF = zeros(6*n_arcs, 7*n_arcs)


        for h = 1:n_arcs
            #Indexes for vector computation
            idx_r = six_vec + i_vec*6*(h-1)
            idx_c = six_vec + i_vec*7*(h-1)
            idx_c2 = six_vec + i_vec*7*h

            ## Computation of the state with forward propagation from the departing point
            x0 = V[idx_c] 
            t_span = V[7*h]
            
            m0 = sc.mass_dry + sc.mass_prop
            push!(x0,m0)
            # Phi = STM(x0,t_span,p)
            # @show Phi
            x_init = vcat(x0, reshape(Phi_0, (36, 1)))
            global sol = propagate_2BP(x_init, t_span, sc, flag_thrust = flag_thrust, flag_J2 = flag_J2, flag_STM = flag_STM, args = args)
            
            ## Computation of the derivative of the final state of the forward propagation
            xf_dot = zeros(size(x_init,1))
            xf = sol.u[end]
            eoms_2BP!(xf_dot, xf,  p, t_span)
        
            # Computation of the STM for forward and backward propagation
            Phi = (reshape(sol.u[end][end-35:end],(6,6)))
  
            ## Computation of DF matrix (derivative of F wrt V)

            # # Regular DF and F 
            DF[idx_r, idx_c] = Phi
            DF[idx_r, 7*h] = xf_dot[1:6]
            if h != n_arcs
                DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
            else           
                diff_step = sol.u[end][1:6]-xdf[1:6]
            end

            #Compute the constraint vector
            F[idx_r] = diff_step
             
            # Constrain juat thw final position
            DF[idx_r, idx_c] = Phi
            DF[idx_r, 7*h] = xf_dot[1:6]
            if h != n_arcs
                DF[idx_r, idx_c2] = -I + zeros(6,6)
                # Evaluation of the new constraint vector
                diff_step = sol.u[end][1:6]-V[idx_c2]
            else           
                diff_step = sol.u[end][1:6]-xdf[1:6]
            end

            #Compute the constraint vector
            F[idx_r] = diff_step
             
        end

        # Update equation
        # Evaluation of the new free variables vector
        global V = V - transpose(DF)*((DF*transpose(DF))\F)
           
        # New step 
        i +=1
        if i == i_max
            print("multiple shooting failed")
            flag_convergence = false
            break
        end
    end

    
    #Compute Δv
    Δm = m_prop0 - sc.mass_prop
    # Replace spacecraft mass with the initial value to allow a correct evaluation during the propagation in the transfer_states function
    sc.mass_prop = m_prop0

    return V, Δm, flag_convergence

end


function multiple_shooting_constraints(V0, xd, n_arcs, sc::Spacecraft; tol=1e-6, args =[])
    """ Multiple shooting correction method"""

    #Set iterations values
    i = 1
    i_max = 5000

    #Set the desired initial free variable vector
    global V = V0
    global flag_convergence = true
    # Set initial constraint vector and initial STM
    if flag_constraint == "xf_all"
        # Constrain the final state to be equal in position and velocity to xdf
        F = ones(6*n_arcs)*1e-3
    elseif flag_constraint == "xf_pos"
        # Constrain the final state to be equal just in position to xdf
        F = ones(6*(n_arcs-1)+3)*1e-3
    end 

    Phi_0 = I+zeros(6,6)


    # Set initial propellant mass
    m_prop0 = sc.mass_prop
    
    #Vector for index propagation          
    i_vec = Vector{Int64}(undef,6)
    fill!(i_vec,1)
    six_vec0 = [1,2,3,4,5,6] 
    six_vec1 = [5,6,7,8,9,10] - i_vec*7

    @show n_arcs
    while norm(F)>tol && i<i_max
        
        # Replace the initial propellant mass at every iteration
        sc.mass_prop = m_prop0

        # Compute the engine efficiency
        veff = sc.Isp * g0 * 1000.0
        p = []
        # Define the additonal arguments for the eoms function
        flag_thrust = true
        flag_J2 = true
        flag_STM = true
        push!(p, flag_thrust, flag_J2, flag_STM, sc.Thrust, veff, args)

        # Initialize DF Matrix
        DF = zeros(6*(n_arcs)-3, 7*(n_arcs)-3) # contrain the initial and final position to be fixed

        x0 = zeros(6)

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
            # Constraint just the final position and not the velocity
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

        @show size(DF), rank(DF)
        @show size(F), norm(F)
        ## Update equation
        # Evaluation of the new free variables vector
        global V = V - transpose(DF)*((DF*transpose(DF))\F)
        # New step 
        i +=1
        if i == i_max
            print("multiple shooting failed")
            flag_convergence = false
            break
        end
    end

    
    #Compute Δv
    Δm = m_prop0 - sc.mass_prop
    # Replace spacecraft mass with the initial value to allow a correct evaluation during the propagation in the transfer_states function
    sc.mass_prop = m_prop0

    @show i

    pushfirst_range!(V, xd[1:3])
    return V, Δm, flag_convergence

end