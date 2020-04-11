# This file generates


import Pkg
Pkg.activate(".")
Pkg.instantiate()

using LightGraphs
using Statistics
include("SIR-discrete-time_fh.jl")
include("Metropolis_sampling_cuts.jl")

function ρ_epidemic(g; β, γ, tsteps, n_initial, args...)
    exp(-1. * mean([epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:100])/30.)
end


g = barabasi_albert(100,10)
cut_ratio = 0.2
cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
β = 0.02; γ = 0.02; tsteps=100; n_initial=10
ρ_temp = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
mc_step_SIR!(cut_g, cut_EL, ρ_epidemic, ρ_temp; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
metropolis!(cut_g, cut_EL, ρ_epidemic, Int(1E2); β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=5)

sample_interval = Int(1E2)

# Warm up the sampler
sample_array_wu = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10) for i in 1:100]

# take 100 samples with a spacing of 1E4
sample_array = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10) for i in 1:100]

s_ρ_wu = [s[3] for s in sample_array_wu]
s_ρ = [s[3] for s in sample_array]
println((mean(s_ρ_wu), mean(s_ρ_wu)+std(s_ρ_wu), mean(s_ρ_wu)-std(s_ρ_wu)))
println((mean(s_ρ), mean(s_ρ)+std(s_ρ), mean(s_ρ)-std(s_ρ)))
