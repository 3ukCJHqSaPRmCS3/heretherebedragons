# File to test the cost functions to implement in the genetic algorthm optimization 
using Plots
# Fuel cost function
# fi: vector of initial fuel percentage in spacecraft tank
# ff: vector of final fuel percentage in spacecraft tank
# n: Number of spacecraft
n = 10
h = 50
fi = []
ff = zeros(n)
ff_tot = zeros(h)
ff_diff_p = zeros(n,h)
# print(rand(Float64, n)) 
J_f_tot = zeros(h)
J_f = zeros(n)
J_f_diff_p = zeros(n,h)
for g = 1:n
    for j = 1:h
        # ff = rand(0:0.01:1,n)
        ff = vcat(rand(0:0.01:1,g),rand(0:0.01:1)*ones(n-g))
        for i = 1:n
            J_f[i] = ff[i]/sum(ff)*log10(ff[i]/sum(ff))
        end
        J_f_tot[j] = sum(J_f)
        ff_tot[j] = sum(ff)
    end
    p = sortperm(ff_tot)

    ff_diff_p[g,:] = ff_tot[p]
    J_f_diff_p[g,:] = J_f_tot[p] 
    
end

plot(ff_diff_p[1,:], J_f_diff_p[1,:], label = "same:9, diff:1")

plot!(ff_diff_p[2,:], J_f_diff_p[2,:], label = "same:8, diff:2")

plot!(ff_diff_p[3,:], J_f_diff_p[3,:], label = "same:7, diff:3")

plot!(ff_diff_p[4,:], J_f_diff_p[4,:], label = "same:6, diff:4")

plot!(ff_diff_p[5,:], J_f_diff_p[5,:], label = "same:5, diff:5")

plot!(ff_diff_p[6,:], J_f_diff_p[6,:], label = "same:4, diff:6")

plot!(ff_diff_p[7,:], J_f_diff_p[7,:], label = "same:3, diff:7")

plot!(ff_diff_p[8,:], J_f_diff_p[8,:], label = "same:2, diff:8")

plot!(ff_diff_p[9,:], J_f_diff_p[9,:], label = "same:1, diff:9")

plot!(ff_diff_p[10,:], J_f_diff_p[10,:], label = "same:0, diff:10")