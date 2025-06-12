using PauliOperators
using LinearAlgebra
using Plots
using BenchmarkTools

function check_majorana_weights(Nmax::Int)
    println("Checking Majorana weight consistency and timing for N = 1 to $Nmax")

    # Lists to store performance and mismatches
    simple_w_times = Float64[]
    P2M_w_times = Float64[]
    mismatch_report = Dict{Int, Vector{Tuple{String, Int, Int}}}()


    for N in 5:Nmax
        println("\n--- N = $N ---")
        paulis, total_paulis = PauliOperators.generate_all_pauli_strings(N)

        println(" Compare times:")
        # Time simple_majorana_weight
        simple_time = @btime begin
            for pauli in $paulis
                a = PauliOperators.simple_majorana_weight(pauli)
            end
        end

        # Time pauli_to_majorana_occupation
        p2m_time = @btime begin
            for pauli in $paulis
                PauliOperators.pauli_to_majorana_occupation(pauli)
            end
        end

        # push!(simple_w_times, simple_time)
        # push!(P2M_w_times, p2m_time)

        mismatches = Tuple{String, Int, Int}[]
        for pauli in paulis
            pauli_str = PauliOperators.string(pauli)
            simple_w = PauliOperators.simple_majorana_weight(pauli)
            P2M_w, monomial = PauliOperators.pauli_to_majorana_occupation(pauli)
            if simple_w != P2M_w
            @show PauliOperators.string(pauli) simple_w P2M_w
            #     push!(mismatches, (pauli_str, simple_w, P2M_w))
            end
        end
        continue

        # if isempty(mismatches)
        #     println("✔ All weights match for N = $N")
        # else
        #     println("✘ Mismatches found for N = $N: $(length(mismatches)) cases from $total_paulis in total")
        #     mismatch_report[N] = mismatches
        # end

        # println("Timing results:")
        # println("  simple_majorana_weight: $(round(simple_time * 1000, digits=3)) ms")
        # println("  pauli_to_majorana_occupation: $(round(p2m_time * 1000, digits=3)) ms")
    end

    return
    # Print and Plot performance results
    Ns = 1:Nmax
    println(" N     simple_w time     P2M_w_time")
    println(Ns,  simple_w_times .* 1000, P2M_w_times .* 1000)
    plot(Ns, simple_w_times .* 1000, label="Simple M weight", xlabel="N (qubits)", ylabel="Time (ms)", lw=2, marker=:circle)
    plot!(Ns, P2M_w_times .* 1000, label="P to M weight", lw=2, marker=:square, title="Majorana Weight Computation Time")
    display(current())

    return mismatch_report
end

# # Example usage (argument of check funciton is the number of qubits):
# mismatches = check_majorana_weights(2)

# # You can inspect mismatches like this:
# for (N, issues) in mismatches
#     println("Mismatches at N = $N:")
#     for (pauli, w1, w2) in issues
#         println("  $pauli -> simple: $w1, P2M: $w2")
#     end
# end

