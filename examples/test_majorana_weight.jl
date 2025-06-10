#using Pkg
#Pkg.activate(".")

#Load the Package
using PauliOperators
using LinearAlgebra


#Example Types definitions
Pb = PauliBasis("ZIXXXZZZXXY")
println("- - - PauliBasis Properties - - - ")

println("Example, pauli basis: ", Pb)
println("ishermitian: ", LinearAlgebra.ishermitian(Pb))
println("Coefficient: ", coeff(Pb))
println("Symplectic Phase: ", symplectic_phase(Pb))
println("Y-Matrix Representation : ", LinearAlgebra.Matrix(PauliBasis("Y")))

w = PauliOperators.simple_majorana_weight(Pb)
println("Simple Majorana weight is: ", w)

w_two = PauliOperators.pauli_to_majorana_occupation(Pb)
println("Majorna weight ver 2 is: ", w_two)

# Test fucntions in helpers
function test_helpers()
    N = 2
    n_lower = 0
    mw_list = Int[]
    pw_list = Int[]
    indx_lst = Int[]


    all_paulis, total_paulis = PauliOperators.generate_all_pauli_strings(N)
    for ps in all_paulis
        println(ps)
    end

    println("Test weights")
    println("weight for IXYZ :", PauliOperators.get_pauli_weight("IXYZ"))  # Output: 3
    println("weight for IIII", PauliOperators.get_pauli_weight("IIII"))  # Output: 0
    
    for (i,op) in enumerate(all_paulis)
        # Extract the Pauli string from our custom example 
        pauli_str = string(op)

        majo_weight, majo_str = PauliOperators.pauli_to_majorana_occupation(op)
        pauli_weight = PauliOperators.get_pauli_weight(pauli_str)

        println("$(i): $pauli_str weight $pauli_weight ==> $majo_str with weight $majo_weight")

        if majo_weight < pauli_weight
            push!(mw_list, majo_weight)
            push!(pw_list, pauli_weight)
            push!(indx_lst, i)
            n_lower +=1
        end
    end

    println("Majorana monomials with weight lower than Paulis: $n_lower from $total_paulis total operators")

    # Plot the results
#    bar_width = 0.35
#    x = collect(1:length(indx_lst))  # ensure x is a Vector{Int}
#    println("typeof(x) ",typeof(x))
#    println("pw_list = ", pw_list, ", typeof(pw_list) = ", typeof(pw_list))
#    println("mw_vals = ", mw_list, ", typeof(mw_list) = ", typeof(mw_list))
#
#
#    x1 = [i - bar_width / 2 for i in x]
#    x2 = [i + bar_width / 2 for i in x]
#
#    bar(
#    x1, pw_list,
#    bar_width = bar_width,
#    label = "Pauli Weight",
#    color = :orange,
#    xlabel = "Operator Index",
#    ylabel = "Weight",
#    title = "Comparison of Pauli and Majorana Weights",
#    legend = :topright
#    )
#
#    bar!(
#    x2, mw_list,
#    bar_width=bar_width,
#    label = "Majorana Weight",
#    color = :gray
#    )
#
#    display(current())

    # # Show example cases
    println("Example case: YXXZZZ")
    pauli_basis = PauliBasis("YXXZZZ")
    println(PauliOperators.pauli_to_majorana_occupation(pauli_basis))
    println(PauliOperators.simple_majorana_weight(pauli_basis))

    println("Example case: ZZ")
    pauli_basis = PauliBasis("ZZ")
    println(PauliOperators.pauli_to_majorana_occupation(pauli_basis))
    println(PauliOperators.simple_majorana_weight(pauli_basis))

    println("Example case: IXXI")
    pauli_basis = PauliBasis("IXXI")
    println(PauliOperators.pauli_to_majorana_occupation(pauli_basis))
    println(PauliOperators.simple_majorana_weight(pauli_basis))

    println("Example case: YZZY")
    pauli_basis = PauliBasis("YZZY")
    println(PauliOperators.pauli_to_majorana_occupation(pauli_basis))
    println(PauliOperators.simple_majorana_weight(pauli_basis))
end


# Run tester if this file is executed directly (scrip mode)
test_helpers()