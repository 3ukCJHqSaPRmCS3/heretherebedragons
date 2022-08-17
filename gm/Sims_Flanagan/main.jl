#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Copyright: Morpheus Space iNC.
#  Identifier: main
#  Module: sims_flanagan_module (Sims Flanagan low thrust optimizer)
#  Author: Giuliana E. Miceli
#  Language: julia v(1.7.3)
#  Environment:
#  Description: This file calls the function to compute an optimal low thrust trajectory
#  using the Sims-Flanagan method
#  
# ----------------------------------------------------------------------------------------

include("propagators.jl")
include("sims_flanagan_functions.jl")
include("sims_flanagan_setup.jl")
using Plots

### TEST 1: PROPAGATOR ##
# test the propagator
tof = 5000.0
mu = 398600.4415
r0 =[-2.394141003279680e+03,-6.577848345501359e+03,0.000000000000000e+00]
v0 = [5.716923258301372e+00,-2.080789897605850e+00,6.083822658434887e+00]
e = ((norm(v0)^2 - mu / norm(r0)) * r0 - (dot(r0,v0)) * v0) / mu

rl,vl,θl = propagate_lagrangian(r0,v0,tof,mu)
if dot(rl,vl) < 0 || θl <0
  print("\nnegative")
  print("\n", θl)
  θl = 2*pi - θl
end
print("\n", θl*180/pi, "\n",)

sol = propagate_2BPt(r0,v0, 100.0,[0.0,0.0,0], tof, mu, 2070*9.8, 1e-16, 1e-16)
rf = sol[1:3,end]
vf = sol[4:6,end]
mf = sol[7,end]
print("\n", rf, "\n", vf, "\n", mf)
θf = acos(dot(rf,e)/(norm(rf) * norm(e))) #[rad]
if dot(rf,vf) < 0 
    print("\nnegative")
    print("\n", θf)
    θf = 2*pi - θf
end
print("\n", θf*180/pi, "\n",)

r = [4.372172353629636e+03 ,  1.224224800435459e+04,  -7.859731956781913e+01]
v = [-3.097537219068065e+00 ,  1.067190264498255e+00 , -3.275733434587419e+00]
θ = acos(dot(r,e)/(norm(r) * norm(e))) #[rad]
if dot(r,v) < 0 
    print("negative")
    θ -= 2*pi
end
print("\n", θ*180/pi, "\n",)

### TEST 2: set-up and sims-flanagan functions
# start_time = "2001-05-26T12:00:00.00"
# x0 = Sc_state([5.656854249492387e+03, 3.999999999999994e+03, 3.999999999999996e+03],
#  [-4.991245094538029e+00,3.529343252911940e+00,3.529343252911939e+00], 100.0)
# throttles  = [Throttle(start_time,"2001-05-26T12:00:01.00", [0.0,0.0,0.0] ), 
# Throttle("2001-05-26T12:01:41.059","2001-05-26T12:06:41.059", [0.0,0.0,0.0] ),
# Throttle("2001-05-26T12:18:22.529","2001-05-26T12:18:30.529", [0.,0.0,0.0] )]
# end_time = "2001-05-26T12:20:00.000"
# xf = Sc_state([-2.160090746950178e+03,5.446742511122797e+03 , 5.446742511122798e+0],
#   [-6.796506709963766e+00,-1.347692793058743e+00,-1.347692793058744e+00], 100.0)
# sc_test = Spacecraft(100.0, 1.0, 2000.0)
# mu = 398600.4415
# leg_test = Leg(start_time,x0,throttles,end_time,xf,sc_test,mu)

# # print_characteristics(leg_test)

# t_grid, list_pos, list_vel, mass = get_leg_states(leg_test,false);


