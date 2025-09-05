#Load the Package
using PauliOperators
using LinearAlgebra
using Plots

#=
 In this example we will explore the weight distribution of Majorana operators
 and compare with the corresponding one for Pauli strings. This is done by 1)
 for a given number N of qubits, generating all 4^N Pauli strings, 2) computing
 the weight distribution of these strings, and 3) comparing it with the Majorana
 weight distribution. Results are plotted as an histogram.
=#

#Define function generate_all_pauli_strings
"""
This function generates all possible Pauli strings for a given number of qubits N.
Uses recursion to build the strings, each string can be seen as a number in base 4
where digits 0-3 map to I,X,Y,Z.
 - Each integer from 0 to 4^N-1 corresponds to a unique sequence of Pauli operatos.
 - We map the remainder % 4 to a character in ['I', 'X', 'Y', 'Z'].
 - Finally we reverse the string because we're filling it from the least significant position first.
"""
function generate_all_pauli_strings(N::Int)
    paulis = ['I', 'X', 'Y', 'Z']
    total = 4^N
    all_strings = Vector{String}(undef, total)

    for i in 0:(total - 1)
        s = IOBuffer()
        num = i
        for _ in 1:N
            digit = num % 4
            num ÷= 4
            print(s, paulis[digit + 1])
        end
        # reverse because we generate least-significant "qubit" first
        all_strings[i + 1] = String(reverse(String(take!(s))))
    end

    return all_strings, total
end

#1) Generate all Pauli strings for a given number of qubits
N = 3  # Example: 4 qubits
all_paulis, total_paulis = PauliOperators.generate_all_pauli_strings(N)
println("Total Pauli strings for $N qubits: ", total_paulis)

#2) Compute the weight distribution of these strings
Pw = Vector{Float32}()
Mw = Vector{Float32}()
for ps in all_paulis
    P = Pauli(ps)
    pw = PauliOperators.pauli_weight(P)
    mw = PauliOperators.simple_majorana_weight(P)
    push!(Pw, pw)
    push!(Mw, mw)
end

#3) Plot the weight distributions and save
# Get max bin limit across both datasets so ticks match
max_w = max(maximum(Pw), maximum(Mw))
xticks_vals = 0:Integer(max_w)

# First histogram
histogram(Pw;
    bins = (-0.5):(max_w + 0.5),                      # Center bars
    xticks = (xticks_vals, string.(xticks_vals)),     # Integer labels
    xlims = (-0.5, max_w + 0.5),                  # keep bounds tight
    title = "N = $N",
    xlabel = "Weight",
    ylabel = "Frequency",
    fillcolor = :linear_blue_95_50_c20_n256,
    label = "Pauli",
    guidefontsize = 14,
    tickfontsize = 10,
    legendfontsize = 10
)

# Overlay second histogram
histogram!(Mw;
    bins = (-0.5):(max_w + 0.5),                      # Same bin edges
    xticks = (xticks_vals, string.(xticks_vals)),     # Same labels
    xlims = (-0.5, max_w + 0.5),                  # keep bounds tight
    fillcolor = :cyclic_grey_15_85_c0_n256,
    label = "Majorana",
    guidefontsize = 14,
    tickfontsize = 10,
    legendfontsize = 10
)

savefig("weight_distribution_N_$N.png")