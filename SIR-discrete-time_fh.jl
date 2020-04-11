function step!(g, I, has_not_been_I, β, γ)
    N_inf = length(I)
    for v in I[1:N_inf]
        # Infect neighbours
        for n in neighbors(g, v)
            if has_not_been_I[n] && rand() < β
                push!(I, n)
                has_not_been_I[n] = false
            end
        end
    end

    for i in N_inf:-1:1
        if rand() < γ
            deleteat!(I, i)
        end
    end
end

function epidemic(g, I_init, N; β=0.01, γ=0.2)
    I = copy(I_init)
    hnbi = trues(nv(g))
    for v in I
        hnbi[v] = false
    end
    max_I = length(I)
    for i in 1:N
        step!(g, I, hnbi, β, γ)

        if length(I) > max_I
            max_I = length(I)
        end

    end
    # Returns maximum number of simultaneously infected, number of infected after
    # last step and total number of those who have been infected
    max_I, length(I), nv(g) - sum(hnbi)
end
