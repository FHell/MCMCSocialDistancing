# This file generates


import Pkg
Pkg.activate(".")
# Pkg.instantiate()

@everywhere using LightGraphs
@everywhere using Statistics
@everywhere include("SIR-discrete-time_fh.jl")
@everywhere include("Metropolis_sampling_cuts.jl")

# function ρ_epidemic(g; β, γ, tsteps, n_initial, args...)
#     exp(-1. * mean([epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:1000])/30.)
# end
function ρ_epidemic(g; β, γ, tsteps, n_initial, args...)
    es = @distributed (vcat) for i in 1:1000
        epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1]
    end
    exp(-1. * mean(es)/30.)
end
function ρ_epidemic1(g; β, γ, tsteps, n_initial, args...)
    es = zeros(1000)
    Threads.@threads for i in 1:1000
        es[i] = epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1]
    end
    exp(-1. * mean(es)/30.)
end

g = barabasi_albert(100,10)
cut_ratio = 0.8
cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
β = 0.02; γ = 0.02; tsteps=1000; n_initial=10
@time ρ_temp = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
@btime ρ_temp = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
@time ρ_temp = ρ_epidemic1(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
@btime ρ_temp = ρ_epidemic1(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)

Nref = 10
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = ρ_epidemic(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end

#mc_step_SIR!(cut_g, cut_EL, ρ_epidemic, ρ_temp; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
#metropolis!(cut_g, cut_EL, ρ_epidemic, Int(1E2); β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E4)
I_r = mean([-log(s)*30. for s in s_ρ_rand])

sample_interval = Int(1E2)
if "long" in ARGS
    println("Long run")
    sample_interval *= 100
end

# Temperature Schedule:
# sample_array_wu1 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E2) for i in 1:50]
# sample_array_wu2 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E3) for i in 1:50]
# sample_array_wu3 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1E3) for i in 1:50]
# sample_array_wu4 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10000) for i in 1:300]
# sample_array_wu5 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=50000) for i in 1:200]

using BSON

i = 1
while isfile("data_$i.bson")
    global i
    i += 1
end

bson("data_v4_$i.bson", g = g, sa = sample_array_wu4, cut_ratio = cut_ratio, β = β, γ = γ, tsteps=tsteps, n_initial=n_initial, sample_interval=sample_interval)
bson("data_v5_$i.bson", g = g, sa = sample_array_wu5, cut_ratio = cut_ratio, β = β, γ = γ, tsteps=tsteps, n_initial=n_initial, sample_interval=sample_interval)
println(mean([-log(s[3])*30. for s in sample_array_wu4]))
