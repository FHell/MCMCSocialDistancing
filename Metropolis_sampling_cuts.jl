using Random

"""
Initializes a given graph g by cutting cut_size edges randomly. Returns
edge list of cut and still existing edges.
"""
function initialize_graph!(g, cut_size)
    EL = collect(edges(g))
    cut_EL = shuffle(EL)[1:cut_size] # rand(EL, cut_size) can select Edges multiple times.
    cut_g = copy(g)
    for edge in cut_EL
        rem_edge!(cut_g, edge)
    end
    cut_g, cut_EL
end


"""
MC step function, restores an edge and cuts a new one, then evaluates ρ_func
    on the newly cut graph.
"""
function mc_step_SIR!(cut_g, cut_EL, ρ_func, ρ_old; ν=1, reject_counter = [1], resamp_old = false, args...)
    if resamp_old
        ρ_old = ρ_func(cut_g; args...)
    end

    ind = rand(1:length(cut_EL))
    restore_edge = cut_EL[ind] # Pick a random edge from the list of cut edges
    new_cut_edge = rand(collect(edges(cut_g))) # And a random edge from the graph

    rem_edge!(cut_g, new_cut_edge) # cut the latter
    add_edge!(cut_g, restore_edge) # restore the former
    cut_EL[ind] = new_cut_edge # and update the list of cut edges

    ρ_new = ρ_func(cut_g; args...)

    if (ρ_new/ρ_old)^ν < rand() # if ... reject proposal. exponent ν to flatten the pdf

        #undo the swap:
        add_edge!(cut_g, new_cut_edge)
        rem_edge!(cut_g, restore_edge)
        cut_EL[ind] = restore_edge

        reject_counter[1] += 1

        return ρ_old            # return ρ so it is not calculated again
    else
        return ρ_new            # return ρ so it is not calculated again
    end
end

function metropolis!(cut_g::AbstractGraph, cut_EL, ρ_func::Function, mc_steps::Integer;
    mc_step!::Function=mc_step_SIR!, args...)
    ρ = ρ_func(cut_g; args...)
    reject_counter = [0]
    for step in 1:mc_steps
        ρ = mc_step!(cut_g, cut_EL, ρ_func, ρ; reject_counter=reject_counter, args...)
    end
    println("Rejected $(reject_counter[1]) out of $mc_steps Steps")
    cut_g, cut_EL, ρ
end
