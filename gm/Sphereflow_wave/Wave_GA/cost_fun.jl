""" create a structure to evaluate the cost of the transfers in the population
    if the population is pop = [x_A1, x_B1; ...;  x_An, x_Bn], the cost J is 
    evaluated for each tuple (x_Ai, x_Bi), with i = {1,n}, so to go from state A to state B
"""

using Statistics

include("constants.jl")
include("structures.jl")
include("population.jl")

# Single spacecraft time cost function
J_st = Cost_fnx()
J_st.type = "singleTime"
J_st.weight = 0.5

# Mean constellation time cost function
J_mt = Cost_fnx()
J_mt.type = "meanTime"
J_mt.weight = 0.25

# Single spacecraft fuel cost function
J_sf = Cost_fnx()
J_sf.type = "singleFuel"
J_sf.weight = 0.5

# Mean distributed fuel in constellation cost function
J_mf = Cost_fnx()
J_mf.type = "meanFuel"
J_mf.weight = 0.25

# Total cost with weights
J_tot = Cost_fnx()
J_tot.type = "totalCost"


function eval(population::Population, constA::Constellation, constB::Constellation)
    """ Evaluate the cost function for the current structure of the population
    """
    
    𝒫, 𝒮_a, 𝒮_b = population, constA.spacecrafts, constB.spacecrafts

    𝒫.J_st,  𝒫.J_sf,  𝒫.J_mt,  𝒫.J_mf = J_st, J_sf, J_mt, J_mf

    # Evaluate single transfers cost functions J_st and J_sf

    for i = 1 : size(𝒫.list_states, 1)

        # Compute transfer with low thrust
        transfer(𝒫.list_states[i,1:6], 𝒫.list_states[i,7:12], 𝒮_b[i], transf)
        t0 = transf.V[7]
        push!(𝒫.J_st.value, t0)
        push!(𝒫.J_st.value, transf.Δm)
    end

    # Evaluate mean transfers cost functions J_mt and J_mf
    𝒫.J_mt.value = mean(J_st)

    J_mf_i= []
    for i = 1 : size(J_sf.value,1)
        J_mf_i[i] = J_sf.value[i]/sum(J_sf.value)*log10(J_sf.value[i]/sum(J_sf.value))
    end 
    𝒫.J_mf.value = sum(J_mf_i)

    # Evalutaiton of the global cost function taking in consideration every single cost function
    𝒫.J_tot.value = J_st.weight*max(J_st.value) + J_mt.weight * J_mt.value + J_sf.weight * sum(J_sf.value) + J_mf.weight * J_mf.value

    print("total cost:",𝒫.J_tot.value, "\n" )
    @show 𝒫
end



