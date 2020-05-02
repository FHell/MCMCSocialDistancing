using LightGraphs

"""
Initializes a given graph g by cutting cut_size edges randomly. Returns
edge list of cut and still existing edges.
"""
function initialize_graph!(g, cut_size)
    EL = collect(edges(g))
    is_in_graph = trues(ne(g))          # true if edge i of EL is in graph
    cut_EL = typeof(EL)(undef, cut_size)     # List of edges cut from the graph
    existing_EL = copy(EL)              # List of edges present in the graph

    index_list = zeros(Integer, cut_size)
    for i in 1:cut_size
        index = rand(1:length(EL))
        while !is_in_graph[index]
            index = rand(1:length(EL))
        end
        is_in_graph[index] = 0
        index_list[i] = index
        cut_EL[i] = EL[index]
        rem_edge!(g, EL[index])
    end
    deleteat!(existing_EL, sort!(index_list))
    g, cut_EL, existing_EL
end


"""
MC step function, swaps one edge from the cut_EL with the existing_EL and also
removes/adds the edges in the graph.
"""
function mc_step_SIR!(g, ρ_func, ρ_old, existing_EL, cut_EL; ν=1, args...)
    cut_size = length(cut_EL)
    # choose index of edges to be swapped
    index_existing = rand(1:(ne(g)-cut_size))
    index_cut = rand(1:cut_size)
    # select edges to be swapped
    cut_edge = cut_EL[index_cut]                    # edge that was cut and is added to the graph
    existing_edge = existing_EL[index_existing]     # edge that was in the graph and is cut from it
    # energy = energy_function(g, 1000)[1] # Maybe no need to calculate energy again
    # ρ_old = ρ(g, args...)

    rem_edge!(g, existing_edge)
    add_edge!(g, cut_edge)
    ρ_new = ρ_func(g; args...)
    if (ρ_new/ρ_old)^ν < rand() # if ... reject proposal. exponent ν to flatten the pdf
        add_edge!(g, existing_edge)
        rem_edge!(g, cut_edge)
        return ρ_old            # return ρ so it is not calculated again
    else
        # swap the two edges in the lists
        cut_EL[index_cut] = existing_edge
        existing_EL[index_existing] = cut_edge

        return ρ_new            # return ρ so it is not calculated again
    end
end


function ρ_epidemic(g; β, γ, tsteps, I_initial, args...)
    epidemic2(g, rand(1:nv(g),I_initial), tsteps; β=β, γ=γ)[1] ^ -1.
end



function metropolis!(graph::AbstractGraph, ρ_func::Function, mc_steps::Integer;
    mc_step!::Function=mc_step_SIR!,  existing_EL, cut_EL, args...)
    ρ = ρ_func(g; args...)
    for step in 1:mc_steps
        ρ = mc_step!(g, ρ_func, ρ, existing_EL, cut_EL;args...)
    end
    graph
end


# WORKING EXAMPLE
g = barabasi_albert(100,10)
g, cut_EL, existing_EL = initialize_graph!(g, 10)
β = 0.02; γ = 0.02; tsteps=100; I_initial=10
ρ_temp = ρ_epidemic(g, β=β, γ=γ, tsteps=tsteps, I_initial=I_initial)
mc_step_SIR!(g, ρ_epidemic, ρ_temp, existing_EL, cut_EL; β=β, γ=γ, tsteps=tsteps, I_initial=I_initial)

metropolis!(g, ρ_epidemic, 100000, existing_EL=existing_EL, cut_EL=cut_EL; β=β, γ=γ, tsteps=tsteps, I_initial=I_initial, ν=1)

#Check if number of max_infected_mean decreases after metropolis

N = 10000
max_infected_mean=0
for i in 1:N
    global max_infected_mean += epidemic2(g, rand(1:nv(g),10), 100, β=0.02, γ=0.02)[1]
end
max_infected_mean /= N



#### OLD STUFF

"""
This functions executes a simulated annealing of a given graph.
"""
function anneal!(graph::AbstractGraph, ρ_func::Function, mc_steps::Integer;
    data::AbstractArray=[], data! =collect_data_energy!, annealing_schedule::Function=ann_sched_mult_exp,
    mc_step!::Function=mc_step_SIR!,
    temperature_changes::Integer=100, β_0::Float64=0., existing_EL, cut_EL, args...)
    energy = energy_function(g,100)[1]
    for step in 1:mc_steps
        β = annealing_schedule(round(Integer, temperature_changes * step/mc_steps); args...)
        ρ = mc_step!(g, ρ_func, existing_EL, cut_EL, ρ, β)
        data!(data; graph=graph, energy=energy, step=step, args...)
    end
    graph
end
data = zeros(1000)

g = anneal!(g, energy_SIR_dt, 1000, existing_EL=existing_EL, cut_EL=cut_EL, data=data)
data
energy_SIR_dt(g,1000)

ann_sched_mult_exp(k::Integer; T_0::Real=1., T_n::Real=0., α::Real=0.95, args...) = 1. / (T_0 * α^k + T_n)

function energy_SIR_dt(g, N_simulations)
    energy= 0
    energy2=0
    for i in 1:N_simulations
        temp = epidemic2(g, [rand(1:nv(g))], 20, β=0.08)[1]
        energy+=temp
        energy2+=temp^2
    end
    energy/N_simulations, sqrt(energy2/N_simulations - (energy/N_simulations)^2)
end

""" this function puts the energy in a given data array every collect_every steps"""
function collect_data_energy!(data::AbstractArray; energy::Real, step::Int, collect_every::Int=1, args...)
    if (step % collect_every) == 0
        index = div(step, collect_every)
        data[index] = energy
    end
end
