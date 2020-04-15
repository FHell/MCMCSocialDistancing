# This file generates


import Pkg
Pkg.activate(".")
# Pkg.instantiate()

using LightGraphs
using Statistics
include("SIR-discrete-time_fh.jl")
include("Metropolis_sampling_cuts.jl")

function ρ_epidemic(g; β, γ, tsteps, n_initial, args...)
    exp(-1. * mean([epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:1000])/30.)
end


g = barabasi_albert(100,10)
cut_ratio = 0.8
cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
β = 0.02; γ = 0.02; tsteps=100; n_initial=10
ρ_temp = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)

Nref = 100
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end

#mc_step_SIR!(cut_g, cut_EL, ρ_epidemic, ρ_temp; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
#metropolis!(cut_g, cut_EL, ρ_epidemic, Int(1E2); β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E4)

sample_interval = Int(1E2)
if "long" in ARGS
    println("Long run")
    sample_interval *= 100
end

# Temperature Schedule:
sample_array_wu1 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E1) for i in 1:50]
sample_array_wu2 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E2) for i in 1:50]
sample_array_wu3 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E3) for i in 1:50]
sample_array_wu4 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E4) for i in 1:500]

using BSON

bson("data.bson", g = g, sa = sample_array_wu4, cut_ratio = cut_ratio, β = β, γ = γ, tsteps=tsteps, n_initial=n_initial)

# I2 = [-log(s[3])*30. for s in sample_array_wu4]

println(mean[-log(s[3])*30. for s in sample_array_wu4])

# b_data = BSON.load("data.bson")[:sa]

# take 100 samples with a spacing of sample_interval
# sample_array = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E4) for i in 1:100]
#
# s_ρ_wu = [s[3] for s in sample_array_wu]
# s_ρ = [s[3] for s in sample_array]
# println((mean(s_ρ_rand), mean(s_ρ_rand)+std(s_ρ_rand), mean(s_ρ_rand)-std(s_ρ_rand)))
# println((mean(s_ρ_wu), mean(s_ρ_wu)+std(s_ρ_wu), mean(s_ρ_wu)-std(s_ρ_wu)))
# println((mean(s_ρ), mean(s_ρ)+std(s_ρ), mean(s_ρ)-std(s_ρ)))

# I = [mean([epidemic(s[1], rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:1000]) for s in b_data]
# mean(I)
# cut_g_r, cut_EL_r = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
# I_rand = [mean([epidemic(initialize_graph!(g, Int(round(cut_ratio * ne(g))))[1], rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:1000]) for i in 1:100]
# mean(I_rand)
