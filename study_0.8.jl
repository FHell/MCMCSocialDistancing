
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using LightGraphs
using Statistics
using GraphPlot
include("SIR-discrete-time_fh.jl")
include("Metropolis_sampling_cuts.jl")

function ρ_epidemic(g; β, γ, tsteps, n_initial, args...)
    epidemic(g, rand(1:nv(g),n_initial), tsteps; β=β, γ=γ)[1] ^ -1.
end


g = barabasi_albert(100,10)
sl_x, sl_y = spring_layout(g)
gplot(g, sl_x, sl_y)

cut_ratio = 0.8
β = 0.02; γ = 0.02; tsteps=100; n_initial=10

#Reference distribution
Nref = 10000
s_ρ_rand = zeros(Nref)
for i in 1:Nref
    cut_g_r, cut_EL_r = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
    s_ρ_rand[i] = ρ_epidemic(cut_g_r, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end

s_ρ_rand_095 = zeros(Nref)
for i in 1:Nref
    cut_g_095, cut_EL_095 = initialize_graph!(g, Int(round(0.95 * ne(g))))
    s_ρ_rand_095[i] = ρ_epidemic(cut_g_095, β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
end

cut_g_r, cut_EL_r = initialize_graph!(g, Int(round(0.8 * ne(g))))
cut_g_095, cut_EL_095 = initialize_graph!(g, Int(round(0.95 * ne(g))))

I_rand = [1/ρ for ρ in s_ρ_rand]
I_rand_095 = [1/ρ for ρ in s_ρ_rand_095]
println("Average max infected with 80% of links cut: $(mean(I_rand)) ± $(std(I_rand))")
println("Average max infected with 95% of links cut: $(mean(I_rand_095)) ± $(std(I_rand_095))")

gplot(cut_g_r, sl_x, sl_y)
gplot(cut_g_095, sl_x, sl_y)

# cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))
# mc_step_SIR!(cut_g, cut_EL, ρ_epidemic, ρ_temp; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial)
# metropolis!(cut_g, cut_EL, ρ_epidemic, Int(1E4); β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10)

sample_interval = Int(1E4)
if "long" in ARGS
    println("Long run")
    sample_interval *= 100
end

cut_g, cut_EL = initialize_graph!(g, Int(round(cut_ratio * ne(g))))

# Warm up the sampler
sample_array_wu = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10) for i in 1:100]

I = [mean([1/ρ_epidemic(s[1]; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial) for i in 1:1000]) for s in sample_array_wu]

mean(I)


# take 100 samples with a spacing of 1E4
sample_array = [metropolis!(cut_g, cut_EL, ρ_epidemic, sample_interval; β=β, γ=γ, tsteps=tsteps, n_initial=n_initial, ν=10) for i in 1:100]

s_ρ_wu = [s[3] for s in sample_array_wu]
s_ρ = [s[3] for s in sample_array]

s_I_wu = [1/s[3] for s in sample_array_wu]
s_I = [1/s[3] for s in sample_array]

println("Average max infected with 80% of links cut (MCMC wu): $(mean(s_I_wu)) ± $(std(s_I_wu))")
println("Average max infected with 80% of links cut (MCMC): $(mean(s_I)) ± $(std(s_I))")
