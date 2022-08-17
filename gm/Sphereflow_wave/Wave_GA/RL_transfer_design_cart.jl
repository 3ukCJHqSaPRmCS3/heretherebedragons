using LinearAlgebra
using Test
using ReinforcementLearningBase.RLBase
using Random
using Random: AbstractRNG
using ClosedIntervals
include("propagator.jl")
export SpacecraftEnv

struct SpacecraftEnvParams{T}
    max_radius::T
    max_speed::T
    goal_distance::T
    Isp::T
    g0::T
    Thrust::T
    mu::T
    max_steps::Int
    timestep::Int
    final_state::Vector{T}
end

Base.show(io::IO, params::SpacecraftEnvParams) = print(
    io,
    join(["$p=$(getfield(params, p))" for p in fieldnames(SpacecraftEnvParams)], ","),
)

function SpacecraftEnvParams(;
    T = Float64,
    max_radius = 8371, #km (2000 km above Earth surface)
    max_speed = 7.9, #km/s (speed of a circular Earth orbit with radius = 100 km)
    goal_distance = 0.5, #km (goal distance between sc0 and target final state)
    Isp = 2100, # Nominal operational range for Morpheus thrusters
    g0 = 0.00981, #km/s^2
    Thrust = 0.0002, #N = 200 microN
    mu = 398600.4, #km^3/s^2
    max_steps = 5e6,
    timestep = 50,
    final_state = [-5217.95, -4749.20, 792.49 ,3.4938,-4.6301,-4.7431], #km[1:3] and km/s[4:6] # Final state to be reached: x,y,z,vx,vy,vz
)
    
    SpacecraftEnvParams{T}(
        max_radius,
        max_speed,
        goal_distance,
        Isp,
        g0,
        Thrust,
        mu,
        max_steps,
        timestep,
        final_state,

    )
end

mutable struct SpacecraftEnv{A,T,ACT,R<:AbstractRNG} <: AbstractEnv
    params::SpacecraftEnvParams{T}
    action_space::A
    observation_space::Space{Vector{ClosedInterval{T}}}
    state::Vector{T}
    action::ACT
    reward::Float64
    done::Bool
    t::Int
    rng::R
end

function SpacecraftEnv(;
    T = Float64,
    rng = Random.GLOBAL_RNG,
    kwargs...,
)
    
    params = SpacecraftEnvParams(; T = T, kwargs...)

    action_space = Base.OneTo(200)
    env = SpacecraftEnv(
        params,
        action_space,
        Space([-params.max_radius..params.max_radius,-params.max_radius..params.max_radius,-params.max_radius..params.max_radius, 
        -params.max_speed..params.max_speed,-params.max_speed..params.max_speed,-params.max_speed..params.max_speed,0.0..1e4,0.0..1e4]), 
        zeros(T, 8),
        rand(action_space),
        0.0,
        false,
        0,
        rng,
    )

    reset!(env)
    env
end

Random.seed!(env::SpacecraftEnv, seed) = Random.seed!(env.rng, seed)
RLBase.action_space(env::SpacecraftEnv) = env.action_space
RLBase.state_space(env::SpacecraftEnv) = env.observation_space
RLBase.reward(env::SpacecraftEnv{A,T}) where {A,T} = env.done ? SpacecraftFinalReward(env::SpacecraftEnv{A,T}) : SpacecraftEnvReward(env::SpacecraftEnv{A,T})
RLBase.is_terminated(env::SpacecraftEnv) = env.done
RLBase.state(env::SpacecraftEnv) = env.state

function RLBase.reset!(env::SpacecraftEnv{A,T}) where {A,T}
    # Servicing SC states: x,y,z,vx,vy,vz,m,mdot,Δm
    env.state[1] = 1861.38 #km
    env.state[2] = 6494.39 #km
    env.state[3] = 1257.65 #km
    env.state[4] = -7.2648 #km/s
    env.state[5] = 1.8122  #km/s
    env.state[6] = 1.3938  #km/s
    env.state[7] = 7000    #kg
    env.state[8] = 0       #kg

    env.done = false
    env.t = 0
    nothing
end

# function (env::SpacecraftEnv{<:Vector{Float64}})(a::AbstractFloat)
#     @assert a in env.action_space
#     env.action = a
#     _step!(env, a)
# end

function (env::SpacecraftEnv{<:Base.OneTo{Int}})(a::Int)
    @assert a in env.action_space
    env.action = a
    alpha = Base.LinRange(-pi/4, pi/4, 10)
    beta =  Base.LinRange(-pi/6, pi/6, 10)
    thrust = 0:1
    A = collect(Iterators.product(thrust, alpha, beta))
    v = vec(A)
    _step!(env, v[a])
end

function _step!(env::SpacecraftEnv, u)
    sc = env.state
    timestep = env.params.timestep
    env.t += timestep
    
    Tx = env.params.Thrust*u[1] * cos(u[2]) * sin(u[3])
    Ty = env.params.Thrust*u[1] * sin(u[2]) * sin(u[3])
    Tz = env.params.Thrust*u[1] * cos(u[3])

    rx,ry,rz,vx,vy,vz,m,Δm = sc[1:8]
    x0 = [rx,ry,rz,vx,vy,vz,m]
    T = [Tx, Ty, Tz]
    rx_n,ry_n,rz_n,vx_n,vy_n,vz_n,m_n = propagate_2BPt(x0, timestep, T, env.params)
    Δm = m - m_n
    env.state[1:8] = [rx_n,ry_n,rz_n,vx_n,vy_n,vz_n,m_n,Δm]


    env.done = 
        env.t>=env.params.max_steps ||
        env.state[7] <= 500.0 ||#kg out of fuel 
        env.state[1:6] == env.params.final_state # Reach the final state
    nothing
end

function SpacecraftEnvReward(env::SpacecraftEnv{A,T}) where {A,T}
    state = env.state
    rxf,ryf,rzf = env.params.final_state[1:3]
    rx,ry,rz = state[1:3]

    distance_xf = sqrt((rxf - rx)^2 + (ryf - ry)^2 + (rzf - rz)^2)

    if distance_xf <= env.params.goal_distance
        # Check if the goal was approached, but if not reached
        print("Goal approached")
        return(100*one(T)) 
    else
        # The env is not done and the reward is zero (can be changed with the decrease in mass)
        return(zero(T))
    end
end

function SpacecraftFinalReward(env::SpacecraftEnv{A,T}) where {A,T}
    state = env.state
    rxf,ryf,rzf = env.params.final_state[1:3]
    rx,ry,rz = state[1:3]

    distance_xf = sqrt((rxf - rx)^2 + (ryf - ry)^2 + (rzf - rz)^2)

    if distance_xf <= env.params.goal_distance
        # if the environment is done because the goal was reached
        print("Goal reached")
        return(1000*one(T))
        env.done = true
    else
        # if the environment is done but not because it reached one of the two spacecraft 
        return(-distance_xf*one(T))
    end
end

using ReinforcementLearning

# TEST =================================================================================================================================#
# env = SpacecraftEnv()
# RLBase.test_runnable!(env)
# hook  = TotalRewardPerEpisode()
# run(RandomPolicy(action_space(env)), env, StopAfterEpisode(10), hook)
# run(RandomPolicy(action_space(env)),env,StopAfterStep(10_000),TotalRewardPerEpisode())


#=======================================================================================================================================#

using ReinforcementLearning
using StableRNGs
using Flux
using Flux.Losses

function RL.Experiment(
    ::Val{:JuliaRL},
    ::Val{:BasicDQN},
    ::Val{:SpacecraftEnv},
    ::Nothing;
    seed = 123,
)
    rng = StableRNG(seed)
    env = SpacecraftEnv(; T = Float64, max_steps = 5e4, rng = rng)
    ns, na = length(state(env)), length(action_space(env))
    agent = Agent(
        policy = QBasedPolicy(
            learner = DQNLearner(
                approximator = NeuralNetworkApproximator(
                    model = Chain(
                        Dense(ns, 64, relu; init = glorot_uniform(rng)),
                        Dense(64, 64, relu; init = glorot_uniform(rng)),
                        Dense(64, na; init = glorot_uniform(rng)),
                    ) |> gpu,
                    optimizer = ADAM(),
                ),
                target_approximator = NeuralNetworkApproximator(
                    model = Chain(
                        Dense(ns, 64, relu; init = glorot_uniform(rng)),
                        Dense(64, 64, relu; init = glorot_uniform(rng)),
                        Dense(64, na; init = glorot_uniform(rng)),
                    ) |> gpu,
                    optimizer = ADAM(),
                ),
                loss_func = huber_loss,
                stack_size = nothing,
                batch_size = 50,
                update_horizon = 1,
                min_replay_history = 100,
                update_freq = 1,
                target_update_freq = 100,
                rng = rng,
            ),
            explorer = EpsilonGreedyExplorer(
                kind = :linear,
                ϵ_init = 0.99,
                ϵ_stable = 0.01,
                decay_steps = 500,
                rng = rng,
            ),
        ),
        trajectory = CircularArraySARTTrajectory(
            capacity = 80_000,
            state = Vector{Float64} => (ns,),
        ),
    )

    stop_condition = StopAfterStep(500_000, is_show_progress=!haskey(ENV, "CI"))
    hook = TotalRewardPerEpisode()

    Experiment(agent, env, stop_condition, hook, "")
end
using Plots
ex = E`JuliaRL_BasicDQN_SpacecraftEnv`
run(ex)
plot(ex.hook.rewards)