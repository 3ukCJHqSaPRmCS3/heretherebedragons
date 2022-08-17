""" main file for testing and simulation of the Wave Genetic Algorithm Approach
"""

include("population.jl")
include("cost_fun.jl")
include("structures.jl")
include("transfer.jl")
using PlotlyJS

# using PlotlyBase: scatter, scatter3d

# num. of orbits
n_c = 3

# Selection of the orbits for the departure constellation 
orb11 = Orbit(8000.0, 0.2, 45.0, 100.0, 250.0, 0.0, MU_EARTH)
orb21 = Orbit(7700.0, 0.0, 90.0, 15.0, 0.0, 0.0, MU_EARTH)
orb31 = Orbit(7500.0, 0.05, 150, 45.0, 30.0, 0.0, MU_EARTH)

# Selection of the orbits for the arrival constellation 
# orb12 = Orbit(8000.0, 0.2, 45.0, 80.0, 250.0, 21.2, MU_EARTH)
# orb22 = Orbit(7700.0, 0.0, 90.0, 20.0, 0.0, 90.0, MU_EARTH)
# orb32 = Orbit(7500.0, 0.05, 150, 75.0, 30.0, 90.0, MU_EARTH)
orb12 = Orbit(8000.0, 0.2, 45.0, 100.0, 250.0, 15.0, MU_EARTH)
orb22 = Orbit(7700.0, 0.0, 90.0, 15.0, 0.0, 90.0, MU_EARTH)
orb32 = Orbit(7500.0, 0.05, 150, 45.0, 30.0, 90.0, MU_EARTH)

# Creation of the spacecrafts vector inside the constellations
sc = Spacecraft()
sc.Isp = 4000.0
sc.Thrust = 0.0008 # Newton 
sc.mass_dry = 13.84 #kg
sc.mass_prop = 0.16 
spacecrafts_A  = Vector{Spacecraft}(undef,n_c)
spacecrafts_B = Vector{Spacecraft}(undef,n_c)
fill!(spacecrafts_A, sc)
fill!(spacecrafts_B, sc)

cA = Constellation(n_c, [orb11; orb21; orb31], spacecrafts_A)
cB = Constellation(n_c, [orb12; orb22; orb32], spacecrafts_B)

# Plot the orbits and find the initial and final states to create the population
pop1 = Population()
r_A, v_A, r_B, v_B, x0_A, xf_B = population(pop1, cA, cB);


# # propagation
# tspan = 3600*(1*2.035)#1 hour #119.67
# thrust_specs = [1,1,1]
# m0 = sc.mass_dry + sc.mass_prop 
# state0 = vcat(pop1.list_states[1,1:6], m0)
# xf = propagate_2BP(state0, tspan, sc, flag_thrust = true, flag_J2 = true , args = thrust_specs)


# # retrieve states from the integration
# xf_x = []
# xf_y = []
# xf_z = []
# for i = 1:size(xf.u,1)
#     push!(xf_x, xf.u[i][1])
#     push!(xf_y, xf.u[i][2])
#     push!(xf_z, xf.u[i][3])
# end


# Evaluation of the cost function for the constellation 
# eval(pop1, cA, cB)
transf = Transfer()
transfer_ms(pop1.list_states[1,1:6], pop1.list_states[1,7:12], sc, transf, 20, "xixf_pos", i_max = 50, tol=1e-7)
print("\n The total fuel consumption is ", transf.Δm)

# Plot

N = 70
layout = Layout(
    scene=attr(
        xaxis=attr(
            nticks=4,
            range=[-1e4,1e4]
        ),
        yaxis=attr(
            nticks=4,
            range=[-1e4,1e4]
        ),
        zaxis=attr(
            nticks=4,
            range=[-1e4,1e4]
        ),
    ),
    width=700,
    margin=attr(
        r=20,
        l=20,
        b=20,
        t=20
    ),
)

plot([
    scatter3d(x= r_A[:,1,1],y= r_A[:,2,1],z= r_A[:,3,1], mode = "lines"),
    scatter3d(x= r_B[:,1,1],y= r_B[:,2,1],z= r_B[:,3,1], mode = "lines"),
    # scatter3d(x= xf_x,y= xf_y,z= xf_z, mode = "lines"),
    scatter3d(x= transf.states[:,1],y= transf.states[:,2],z= transf.states[:,3], mode = "lines", color = "black"),
    scatter3d(x=  [transf.states[end,1]],y= [transf.states[end,2]],z= [transf.states[end,3]], mode = "marker",marker=attr(size=5)),
    scatter3d(x=  [transf.states[1,1]],y= [transf.states[1,2]],z= [transf.states[1,3]], mode = "marker",marker=attr(size=5)),
    # scatter3d(x= [xf_x[end]],y= [xf_y[end]],z= [xf_z[end]], mode = "marker",marker=attr(size=5)),
    scatter3d(x= [x0_A[1,1]],y= [x0_A[1,2]],z= [x0_A[1,3]], mode= "marker",marker=attr(size=5)),
    scatter3d(x= [xf_B[1,1]],y= [xf_B[1,2]],z= [xf_B[1,3]], mode = "marker", marker=attr(size=5)),
    scatter3d(x= [0.0],y= [0.0],z = [0.0], mode = "marker", marker=attr(size=100))
    
    ],

    layout,
)

# #using Plot
# plot(xf, vars=(1,2,3), linecolor= "blue")
# plot!(r_A[:,1,1],r_A[:,2,1],r_A[:,3,1], linecolor="red")
# plot!(r_B[:,1,1],r_B[:,2,1],r_B[:,3,1], linecolor="red")
# scatter3d!([x0_A[1,1]],[x0_A[1,2]],[x0_A[1,3]])

# scatter3d!([xf_B[1,1]],[xf_B[1,2]],[xf_B[1,3]])
# scatter3d!(0,0,0)