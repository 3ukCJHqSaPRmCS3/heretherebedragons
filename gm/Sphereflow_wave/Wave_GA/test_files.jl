include("constants.jl")
include("structures.jl")
include("propagator.jl")
include("Orbits.jl")
using Plots
# Creation of the spacecrafts vector inside the constellations
sc = Spacecraft()
# USe the 6U full operational range
sc.Isp = 1890.0
sc.Thrust = 0.001007 # Newton 
sc.mass_dry = 10.0
sc.mass_prop = 2.0
orb = Orbit(R_EARTH+1000.0, 0.1, 10.0,0.0,0.0,0.0,MU_EARTH)
#Plot the orbit
x_orb, = kep2cart(orb, false)
r_max = 1000.0
for i = 1:size(x_orb,1)
    r = norm([x_orb[i,1], x_orb[i,2], x_orb[i,3]])
    alt_orb = r - R_EARTH 
    if alt_orb > r_max
        global r_max = alt_orb
    end
end

#Calculate the initial state for propagation
x0, v0 = kep2cart(orb, true)
m0 = sc.mass_dry + sc.mass_prop 
state0 = hcat(x0,v0,m0)

# propagation
tspan = 3600*20 #1 hour
thrust_specs = [1,1,1]
xf = propagate_2BP(state0, tspan, sc, flag_thrust = true, flag_J2 = true , args = thrust_specs)

alt = []
for i = 1:size(xf.u,1)
    r = norm([xf.u[i][1], xf.u[i][2], xf.u[i][3]])
    push!(alt, r - R_EARTH)
end
@show r_max
@show findmax(alt)



plot(xf, vars=(1,2,3), linecolor= "blue")
plot!(x_orb[:,1],x_orb[:,2],x_orb[:,3], linecolor= "red")
scatter!([0,0,0])

# # Plot the altitude variation
plot(xf.t, alt)
plot!(xf.t, r_max*ones(size(xf.t,1)), linecolor = "green")

# #Plot the mass variation 
# plot(xf, vars=(0,7), linecolor= "blue")