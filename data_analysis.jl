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

sample_interval = Int(1E2)

using BSON

sim = BSON.load("data.bson")

γ = sim[:γ]
n_initial = sim[:n_initial]
tsteps = sim[:tsteps]
cut_ratio = sim[:cut_ratio]
g = sim[:g]
sa = sim[:sa]
β = sim[:β]

Nref = length(sa)
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end


I = [-log(s[3])*30. for s in sa]
I_r = [-log(s)*30. for s in s_ρ_rand]

println(mean(I))
println(mean(I_r))
