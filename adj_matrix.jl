using LightGraphs, BenchmarkTools, Random, LinearAlgebra, Statistics

include("SIR-discrete-time_fh.jl")
include("Metropolis_sampling_cuts.jl")

function ρ_epidemic(g; args...)
    exp(-1. * eigmax(Matrix(adjacency_matrix(g,Float64))))
end



g = barabasi_albert(100,10)
cut_ratio = 0.8
cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
β = 0.02; γ = 0.02; tsteps=100; n_initial=10
sample_interval = Int(1E2)

Nref = 10000
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = epidemic(cut_g, rand(1:nv(cut_g),n_initial), tsteps; β=β, γ=γ)[1]
end


sample_array = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1000) for i in 1:200]

s_ρ = zeros(Nref)
for i in 1:Nref
    s_ρ[i] = epidemic(cut_g, rand(1:nv(cut_g),n_initial), tsteps; β=β, γ=γ)[1]
end

mean(s_ρ_rand)
mean(s_ρ)

using GraphPlot, EmbeddedGraphs
gplot(cut_g)
characteristic_length(cut_g)

using Plots
plot!(degree_array(g, bins=100))
plot(degree_array(cut_g))
