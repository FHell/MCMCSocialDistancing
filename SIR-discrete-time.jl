using LightGraphs
using BenchmarkTools


# I :: Infected, R :: Recovered

function step2!(g, I, has_been_I, β, γ)
    N_inf = length(I)
    for v in I[1:N_inf]
        # Infect neighbours
        for n in neighbors(g, v)
            if (! has_been_I[n]) && rand() < β
                push!(I, n)
                has_been_I[n] = true
            end
        end
    end
    # for v,i in enumerate(I[N_inf:-1:1])
    recovering = falses(length(I))
    recovering[1:N_inf] = map(x -> rand() < 0.3, 1:N_inf)
    deleteat!(I, recovering)

end

function epidemic2(g, I, N; β=0.01, γ=0.2)
    hbi = falses(nv(g))
    max_I = length(I)
    for i in 1:N
        step2!(g, I, hbi, β, γ)
        max_I = max(length(I), max_I)
    end
    # Returns maximum number of simultaneously infected, number of infected after
    # last step and total number of those who have been infected
    max_I, length(I), sum(hbi)
end


function step!(I, R, g)
    for v in vertices(g)
        if I[v]
            # Infect neighbours
            for n in neighbors(g, v)
                if (! I[n]) && (! R[n]) && rand() < 0.3
                    #println("I node $n")
                    I[n] = true
                end
            end

            # recover
            if rand() < 0.2
                #println("node $v R")
                I[v] = false
                R[v] = true
            end
        end
    end
end

function epidemic(I_initial, R_initial, g, N)
    I = copy(I_initial)
    R = copy(R_initial)
    max_I = length(I)

    for i in 1:N
        step!(I, R, g)
        max_I = max(sum(I), max_I)
    end
    max_I
    # println("Recovered: $(sum(R)); Infected: $(sum(I))")
end


g = barabasi_albert(1000, 10)

I_initial = falses(nv(g))
I_initial[1] = true
R_initial = falses(nv(g))


epidemic(I_initial, R_initial, g, 1)
@btime epidemic(I_initial, R_initial, g, 20)

I_initial = rand(1:nv(g), 10)
epidemic2(g, copy(I_initial), 100, β=0.07, γ=0.02)
@btime epidemic2(g, copy(I_initial), 100, β=0.07, γ=0.02)
