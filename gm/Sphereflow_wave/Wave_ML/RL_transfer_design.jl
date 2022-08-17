using LinearAlgebra
using Test
using ReinforcementLearningBase.RLBase
using Random
using Random: AbstractRNG
using ClosedIntervals
export SpacecraftEnv

struct SpacecraftEnvParams{T}
    min_radius::T
    max_radius::T
    goal_distance::T
    Isp::T
    g0::T
    Thrust::T
    mu::T
    max_steps::Int
    timestep::Int
end

Base.show(io::IO, params::SpacecraftEnvParams) = print(
    io,
    join(["$p=$(getfield(params, p))" for p in fieldnames(SpacecraftEnvParams)], ","),
)

function SpacecraftEnvParams(;
    T = Float64,
    min_radius = 6571, #km (200 km from Earth surface)
    max_radius = 7571, #km (1000 km above Earth surface)
    goal_distance = 0.5, #km (goal distance between sc0 and target final state)
    Isp = 2100, #Nominal operational range for Morpheus thrusters
    g0 = 0.00981, #km/s^2
    Thrust = 0.0002, #N = 200 microN
    mu = 398600.4, #km^3/s^2
    max_steps = 5e6,
    timestep = 10,
)
    
    SpacecraftEnvParams{T}(
        min_radius,
        max_radius,
        goal_distance,
        Isp,
        g0,
        Thrust,
        mu,
        max_steps,
        timestep,
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

    action_space = Base.OneTo(21)
    env = SpacecraftEnv(
        params,
        action_space,
        Space([6000.0..10000.0,0.0..10,-100.0..100.0,0.0..100.0,0.0..1100.0,-1000.0..1000.0,0.0..1000.0,6000.0..10000.0,0.0..10,-100.0..100.0,0.0..100.0]), 
        zeros(T, 11),
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
RLBase.reward(env::SpacecraftEnv{A,T}) where {A,T} = env.done ? SpacecraftFinalRewardTest(env::SpacecraftEnv{A,T}) : SpacecraftEnvReward(env::SpacecraftEnv{A,T})
RLBase.is_terminated(env::SpacecraftEnv) = env.done
RLBase.state(env::SpacecraftEnv) = env.state

function RLBase.reset!(env::SpacecraftEnv{A,T}) where {A,T}
    # Servicing SC
    env.state[1] = 6571 #km
    env.state[2] = 2.5 #rad
    env.state[3] = 0.0 #km/s
    env.state[4] = 7.7885 #km/s
    env.state[5] = 5000.0 #kg
    env.state[6] = 0 #kg/s
    env.state[7] = 0 #kg

    # Final state to be reached
    env.state[8] = 7000 #km
    env.state[9] = 1.7125 #rad
    env.state[10] = 0.0 #km/s
    env.state[11] = 7.5461 #km/s


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
    v = Vector{Float64}(vcat(Base.LinRange(-pi/4, pi/4, 20), 0.0))
    _step!(env, v[a])
end

function _step!(env::SpacecraftEnv, throttle)
    μ = env.params.mu
    sc = env.state
    # Thr = env.params.Thrust*throttle
    timestep = env.params.timestep
    env.t += timestep
    
    if throttle == 0.0 
        Thr_r = Thr_θ = Thr = 0.0
    else
        Thr_θ = env.params.Thrust * cos(throttle)
        Thr_r = env.params.Thrust * sin(throttle)
        Thr = env.params.Thrust
    end

    r,θ,vr,vθ,m,m_dot,Δm = sc[1:7]
    r_dot, = vr
    θ_dot = vθ/r
    vr_dot = vθ^2/r - μ/(r^2) + Thr_r/(m - norm(m_dot)*timestep)
    vθ_dot = -vr*vθ/r + Thr_θ/(m-norm(m_dot)*timestep) 
    r_new = r+r_dot*timestep
    θ_new = (θ+θ_dot*timestep)%(2*pi)
    vr_new = vr+vr_dot*timestep
    vθ_new = vθ+vθ_dot*timestep
    m_dot = -abs(Thr)/(env.params.Isp*env.params.g0)
    m_new = m+m_dot*timestep
    Δm = m - m_new
    env.state[1:7] = [r_new,θ_new,vr_new,vθ_new,m_new,m_dot, Δm]

    # No update of the final state
    # r,θ,vr,vθ, = sc[8:10]
    # r_dot, = vr
    # θ_dot = vθ/r 
    # vr_dot = vθ^2/r - μ/(r^2)
    # vθ_dot = -vr*vθ/r
    # r_new = r+r_dot*timestep
    # θ_new = (θ+θ_dot*timestep)%(2*pi)
    # vr_new = vr+vr_dot*timestep
    # vθ_new = vθ+vθ_dot*timestep
    # env.state[7:10] = [r_new,θ_new,vr_new,vθ_new]


    env.done = 
        env.t>=env.params.max_steps ||
        env.state[5] <= 0.0 ||#out of fuel 
        env.state[1] > env.state[8] # beyond final state orbit

    nothing
end

function SpacecraftEnvReward(env::SpacecraftEnv{A,T}) where {A,T}
    state = env.state
    r0 = state[1]
    θ0 = state[2]
    Δm = state[7]
    r = state[8]
    θ = state[9]
    distance_xf = sqrt((r*cos(θ) - r0*cos(θ0))^2 + (r*sin(θ) - r0*sin(θ0))^2)
    if distance_xf < env.params.goal_distance
        print("Goal reached")
        return(1000*one(T)) 
        env.done = true
    else
        # print("\n",distance_iss)
        return(-Δm*one(T))
    end
end

function SpacecraftFinalRewardTest(env::SpacecraftEnv{A,T}) where {A,T}
    state = env.state
    r0 = state[1]
    θ0 = state[2]
    r = state[8]
    θ = state[9]
    distance_xf = sqrt((r*cos(θ) - r0*cos(θ0))^2 + (r*sin(θ) - r0*sin(θ0))^2)
    if distance_xf < env.params.goal_distance
        print("Goal reached")
        return(1000*one(T))
        env.done = true
    else
        # if the environment is done but not because it reached one of the two apacecraft 
        return(-distance_xf*one(T))
    end
end

using ReinforcementLearning

# TEST =================================================================================================================================#
# env = SpacecraftEnv()
# RLBase.test_runnable!(env)
# hook  = TotalRewardPerEpisode()
# run(RandomPolicy(action_space(env)), env, StopAfterEpisode(10_000), hook)
# # run(RandomPolicy(action_space(env)),env,StopAfterStep(10_000),TotalRewardPerEpisode())


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
    env = SpacecraftEnv(; T = Float64, max_steps = 5e3, rng = rng)
    ns, na = length(state(env)), length(action_space(env))
    agent = Agent(
        policy = QBasedPolicy(
            learner = DQNLearner(
                approximator = NeuralNetworkApproximator(
                    model = Chain(
                        Dense(ns, 80, relu; init = glorot_uniform(rng)),
                        Dense(80, 80, relu; init = glorot_uniform(rng)),
                        Dense(80, na; init = glorot_uniform(rng)),
                    ) |> gpu,
                    optimizer = ADAM(),
                ),
                target_approximator = NeuralNetworkApproximator(
                    model = Chain(
                        Dense(ns, 80, relu; init = glorot_uniform(rng)),
                        Dense(80, 80, relu; init = glorot_uniform(rng)),
                        Dense(80, na; init = glorot_uniform(rng)),
                    ) |> gpu,
                    optimizer = ADAM(),
                ),
                loss_func = huber_loss,
                stack_size = nothing,
                batch_size = 32,
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

    stop_condition = StopAfterStep(1_000_000, is_show_progress=!haskey(ENV, "CI"))
    hook = TotalRewardPerEpisode()

    Experiment(agent, env, stop_condition, hook, "")
end
using Plots
ex = E`JuliaRL_BasicDQN_SpacecraftEnv`
run(ex)
plot(ex.hook.rewards)