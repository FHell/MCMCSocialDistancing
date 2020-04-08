using LightGraphs
using BenchmarkTools

g = barabasi_albert(1000, 10)

infected_initial = [false for v in vertices(g)]
infected_initial[1] = true
recovered_initial = [false for v in vertices(g)]

function step2!(infected, has_been_infected, g)
    for v in infected
        # Infect neighbours
        for n in neighbors(g, v)
            if (! has_been_infected[n]) && rand() < 0.3
                push!(infected, n)
                has_been_infected[n] = true
            end
        end
    end

    for (i, v) in enumerate(infected)
        # recover and remove all copies of the infection
        if rand() < 0.2
            # recovered[v] = true
            deleteat!(infected, i)
        end
    end
end

function epidemic2(infected_initial, hbi_initial, g, N)
    infected = copy(infected_initial)
    hbi = copy(hbi_initial)
    for i in 1:N
        step2!(infected, hbi, g)
    end
    # println("Have been infected: $(sum(hbi)); Infected: $(length(infected))")
end


function step!(infected, recovered, g)
    for v in vertices(g)
        if infected[v]
            # Infect neighbours
            for n in neighbors(g, v)
                if (! infected[n]) && (! recovered[n]) && rand() < 0.3
                    #println("infected node $n")
                    infected[n] = true
                end
            end

            # recover
            if rand() < 0.2
                #println("node $v recovered")
                infected[v] = false
                recovered[v] = true
            end
        end
    end
end

function epidemic(infected_initial, recovered_initial, g, N)
    infected = copy(infected_initial)
    recovered = copy(recovered_initial)
    for i in 1:N
        step!(infected, recovered, g)
    end
    # println("Recovered: $(sum(recovered)); Infected: $(sum(infected))")
end

epidemic(infected_initial, recovered_initial, g, 20)
@btime epidemic(infected_initial, recovered_initial, g, 20)

epidemic2([1], recovered_initial, g, 20)
@btime epidemic2([1], recovered_initial, g, 20)
