using LightGraphs
using LinearAlgebra
using BenchmarkTools
using Statistics
using EmbeddedGraphs
using LaTeXStrings

include("SIR-discrete-time_fh.jl")
include("Metropolis_sampling_cuts.jl")

g = barabasi_albert(100,10)#barabasi_albert(100,10)
function ρ_epidemic(g; args...)
    eigmax(Matrix(adjacency_matrix(g,Float64)))^ -1.
end
en(g) = eigmax(Matrix(adjacency_matrix(g,Float64)))
function ρ_epidemic(g; args...)
    exp(-1. * eigmax(Matrix(adjacency_matrix(g,Float64))))
end
function ρ_epidemic2(g; β, γ, tsteps, n_initial, args...)
    exp(-1. * mean([epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] for i in 1:1])/30.)
end
function ρ_epidemic1(g; β, γ, tsteps, n_initial, args...)
    epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] ^ -1.
end


g = barabasi_albert(100,10)
cut_ratio = 0.8
cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
β = 0.02; γ = 0.02; tsteps=100; n_initial=10
sample_interval = Int(1E3)

metropolis!(cut_g, cut_EL, ρ_epidemic, Int(1E3); β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=100)

sample_array_wu2 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=500) for i in 1:100]
sample_array_wu1 = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=1000) for i in 1:100]



#Reference distribution
Nref = 1000
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = ρ_epidemic1(cut_g, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end

s_ρ1 = zeros(100)
for i in 1:100
    s_ρ1[i] = mean([ρ_epidemic1(sample_array_wu1[i][1], β=β, γ=γ, tsteps=tsteps, n_initial=n_initial) for j in 1:Nref])
end

s_ρ2 = zeros(100)
for i in 1:100
    s_ρ2[i] = mean([ρ_epidemic1(sample_array_wu2[i][1], β=β, γ=γ, tsteps=tsteps, n_initial=n_initial) for j in 1:Nref])
end
std([ρ_epidemic1(sample_array_wu1[i][1], β=β, γ=γ, tsteps=tsteps, n_initial=n_initial) for j in 1:Nref].^-1)
println((mean(s_ρ_rand.^-1), std(s_ρ_rand.^-1)/sqrt(1000)))
println((mean(s_ρ1.^-1), std(s_ρ1.^-1)))
println((mean(s_ρ2.^-1), std(s_ρ2.^-1)))

using GraphPlot, Plots
gplot(largest_component(cut_g))
using EmbeddedGraphs
C_rand_mean = mean([global_clustering_coefficient(initialize_graph!(g, Int(round(cut_ratio * ne(g))))[1]) for i in 1:10000]) 
C_rand_std = std([global_clustering_coefficient(initialize_graph!(g, Int(round(cut_ratio * ne(g))))[1]) for i in 1:10000])/sqrt(10000)
L_rand_mean = mean([characteristic_length(largest_component(initialize_graph!(g, Int(round(cut_ratio * ne(g))))[1])) for i in 1:1000])
L_rand_std = std([characteristic_length(largest_component(initialize_graph!(g, Int(round(cut_ratio * ne(g))))[1])) for i in 1:1000])

L_mean1 = mean([characteristic_length(largest_component(sample_array_wu1[i][1])) for i in 1:100])
L_std1 = std([characteristic_length(largest_component(sample_array_wu1[i][1])) for i in 1:100])
C_mean1 =  mean([global_clustering_coefficient(sample_array_wu1[i][1]) for i in 1:100])
C_std1 =  std([global_clustering_coefficient(sample_array_wu1[i][1]) for i in 1:100])/sqrt(100)
L_mean2 = mean([characteristic_length(largest_component(sample_array_wu2[i][1])) for i in 1:100])
L_std2 = std([characteristic_length(largest_component(sample_array_wu2[i][1])) for i in 1:100])
C_mean2 =  mean([global_clustering_coefficient(sample_array_wu2[i][1]) for i in 1:100])
C_std2 =  std([global_clustering_coefficient(sample_array_wu2[i][1]) for i in 1:100])
degree_array(cut_g)
graphs1 = [sample_array_wu1[i][1] for i in 1:length(sample_array_wu1)]
degree_arr_mean1 = mean([degree_array(graphs1[i]) for i in 1:100])
degree_arr_std1 = std([degree_array(graphs1[i]) for i in 1:100])
graphs2 = [sample_array_wu2[i][1] for i in 1:length(sample_array_wu2)]
degree_arr_mean2 = mean([degree_array(graphs2[i]) for i in 1:100])
degree_arr_std2 = std([degree_array(graphs2[i]) for i in 1:100])
pgfplotsx()
lastdegree=10
scatter(collect(0:lastdegree),degree_arr_mean1[1:lastdegree+1], yerr=degree_arr_std1[1:lastdegree+1], size=[335.,335. /1.6], xlabel="Degree "*L"k", ylabel="Probability"*L"[\%]", label=L" \nu=1000", title="MH with Adjacency matrix eigenvalue", legend=:topright)
scatter!(collect(0:lastdegree),degree_arr_mean2[1:lastdegree+1], yerr=degree_arr_std2[1:lastdegree+1], label=L"\nu=500")
savefig("Manuscript/adj_mat.tex")
degree_array(sample_array_wu2[2][1])
