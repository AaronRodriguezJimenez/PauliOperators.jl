using IterTools

function get_on_bits(x::T) where T<:Integer
    N = count_ones(x)
    inds = Vector{Int}(undef, N)
    if N == 0
        return inds
    end

    count = 1
    for i in 1:length(bitstring(x))
        if x >> (i-1) & 1 == 1
            inds[count] = i
            count += 1
        end
        count <= N || break
    end
    return inds
end


#= 
  #                                       #
  # - - - Majorana weight functions - - - #
  #                                       #
Here we incporporate funcions designed to estimate the weight
of the majorana string associated to a Pauli String.
=#


"""
 Compute the Majorana weight without translating the Pauli string to
 Majorana string. 
 This function iterates through the Pauli string and applies the following:
    - If z[i] == true, x[i] == false, and the control flag is true, it adds 2 to the weight.
    - If x[i] == true, it flips the control flag and adds 1 to the weight.
 Returns the total weight as an integer.
"""
function simple_majorana_weight(Pb::Union{PauliBasis{N}, Pauli{N}}) where N
    w = 0
    control = true
    # tmp = Pb.z & ~Pb.x  # Bitwise AND with bitwise NOT
    Ibits = ~(Pb.z|Pb.x)
    Zbits = Pb.z & ~Pb.x

    for i in reverse(1:N)  # Iterate from N down to 1
        xbit = (Pb.x >> (i - 1)) & 1 != 0
        Zbit = (Zbits >> (i - 1)) & 1 != 0
        Ibit = (Ibits >> (i - 1)) & 1 != 0
        #println("i=$i, xbit=$xbit, Zbit=$Zbit, Ibit=$Ibit, control=$control, w=$w")
        if Zbit && control || Ibit && !control
            w += 2
        elseif xbit
            control = !control
            w += 1
        end
    end
    return w
end

function pauli_weight(Pb::Union{PauliBasis{N}, Pauli{N}}) where N
    w = 0
    for i in 1:N
        xbit = (Pb.x >> (i - 1)) & 1
        zbit = (Pb.z >> (i - 1)) & 1

        if xbit != 0 || zbit != 0
            w += 1
        end
    end
    return w
end

"""
 Generates all possible Pauli strings for a given lenght N
"""
function generate_all_pauli_strings(N::Int)
    paulis = ['I', 'X', 'Y', 'Z']
    all_basis = PauliBasis{N}[]  # Array of PauliBasis{N}
    count = 0
    for term in product(fill(paulis, N)...)
        str = join(term)
        push!(all_basis, PauliBasis(str))
        count += 1
    end
    return all_basis, count
end

# Custom pretty print
#function Base.show(io::IO, ps::Pauli)
#    if isempty(ps.ops)
#        print(io, "$(ps.coeff)*I")
#    else
#        sorted_ops = sort(collect(ps.ops))
#        ops_str = join(["$op$idx" for (idx, op) in sorted_ops], " ")
#        print(io, "$(ps.coeff)*$ops_str")
#    end
#end

# Compute Pauli weight
function get_pauli_weight(pauli_str::AbstractString)
    count(c -> c in ['X','Y','Z'], pauli_str)
end

#                                                         #
# Functions to create Majorana strings from Pauli strings  #
#
function get_majorana_vector_for_Z_padded_XY(pauli_str::AbstractString)
    N = length(pauli_str)
    vec_len = 2 * N
    vec = zeros(Int, vec_len)

    for i in 1:N
        if pauli_str[i] in ['X', 'Y']
            num_Z = i - 1  # Julia is 1-based, but logic is 0-based
            if pauli_str[i] == 'X'
                vec[2 * num_Z + 1] = 1  # Julia index: m_{2k+1}
            elseif pauli_str[i] == 'Y'
                vec[2 * num_Z + 2] = 1  # Julia index: m_{2k+2}
            end
            break
        end
    end

    return vec
end

function get_majorana_vector_for_Z_only(pauli_str::AbstractString)
    N = length(pauli_str)
    vec_len = 2 * N
    vec = zeros(Int, vec_len)
    z_idx = 0

    for i in 1:N
        p = pauli_str[i]
        if p == 'Z'
            vec[2 * z_idx + 1] = 1  # Julia index correction
            vec[2 * z_idx + 2] = 1
            z_idx += 1
        elseif p == 'I'
            z_idx += 1
        else
            error("Non-Z/I found in Z-only string at position $i: '$p'")
        end
    end

    return vec
end

"""
 Decompose full Pauli string into Majorana occupation vectors.

    For each X/Y, extract its Z-padded form AND the trailing Z-only shield.
    Also extract any standalone Z/I substrings.

    Returns:
        A list of occupation vectors to be XOR'ed into the full Majorana form.
"""
function decompose_to_direct_terms(pauli_str::AbstractString)
    N = length(pauli_str)
    vectors = Vector{Vector{Int}}()  # list of occupation vectors
    used_mask = falses(N)            # tracks used Zs

    # Step 1: global Z-only term (excluding shielded Zs)
    zi_str = join([pauli_str[j] == 'Z' && !used_mask[j] ? 'Z' : 'I' for j in 1:N])
    if occursin('Z', zi_str)
        push!(vectors, get_majorana_vector_for_Z_only(zi_str))
    end

    # Step 2 & 3: Handle X/Y terms
    for i in 1:N
        p = pauli_str[i]
        if p in ['X', 'Y']
            # ZX-like pattern
            pad = ['Z' for _ in 1:(i-1)]
            push!(pad, p)
            append!(pad, ['I' for _ in (i+1):N])
            zxi_str = join(pad)
            vec1 = get_majorana_vector_for_Z_padded_XY(zxi_str)
            push!(vectors, vec1)

            # trailing Z shield
            z_shield = [j != i && zxi_str[j] == 'Z' ? 'Z' : 'I' for j in 1:N]
            if 'Z' in z_shield
                for j in (i+1):N
                    if pauli_str[j] == 'Z'
                        used_mask[j] = true
                    end
                end
                z_only_str = join(z_shield)
                vec2 = get_majorana_vector_for_Z_only(z_only_str)
                push!(vectors, vec2)
            end
        end
    end

    return vectors
end

"""
 Bitwise XOR all majorana occupations: Is equivalent to get the vector associated with Majorana
"""
function combine_majorana_vectors(vectors::Vector{Vector{Int64}})
    result = zeros(Int, length(vectors[1]))
    for v in vectors
        result .= xor.(result, v)
    end
    return result
end

""" 
 Given an occupation vector, return the list of Majorana operators
"""
function get_majorana_operators_from_vector(vec::Vector{Int})
    return ["m$(i)" for (i, val) in enumerate(vec) if val == 1]
end


"""
  Main function for Pauli to Majorana conversion
"""
function pauli_to_majorana_occupation(Pb::PauliBasis{N}) where N
    pauli_str = string(Pb)
    if all(c -> c == 'I', pauli_str)
        return 0, String[]
    else
        vectors = decompose_to_direct_terms(pauli_str)
        result_vec = combine_majorana_vectors(vectors)
        #w = count_ones(result_vec)
        w = count(x -> x == 1, result_vec)
        return w, get_majorana_operators_from_vector(result_vec)
    end
end
