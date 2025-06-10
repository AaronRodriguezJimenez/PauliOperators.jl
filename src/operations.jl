# JUNE - 2025
# This has been taken and adapted from "linsolv" version of PauliOperators
# It contains the expectation_value missed function.


"""
    Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if two Paulis are equal (same basis and same scalar coefficient)
"""
function Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.s == p2.s && p1.z == p2.z && p1.x == p2.x
end

"""
    Base.isequal(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if two Paulis have the same Pauli basis (ignoring scalar)
"""
function Base.isequal(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.z == p2.z && p1.x == p2.x
end

"""
    Base.:(>)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Compare based on basis only (ignoring scalar)
"""
function Base.:>(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.z > p2.z || (p1.z == p2.z && p1.x > p2.x)
end

"""
    Base.:(<)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Compare based on basis only (ignoring scalar)
"""
function Base.:<(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.z < p2.z || (p1.z == p2.z && p1.x < p2.x)
end


"""
    rotate_phase(p::Pauli{N}, θ::Integer)

Multiply phase by `i^θ`
"""
function rotate_phase(p::Pauli{N}, θ::Integer) where N
    return Pauli{N}(p.s * 1im^θ, p.z, p.x)
end

"""
    negate(p::Pauli{N})

Multiply Pauli by -1
"""
negate(p::Pauli{N}) where N = rotate_phase(p, 2)

"""
    is_diagonal(p::PauliBasis)

True if Pauli basis has only Z and I (i.e., x == 0)
"""
is_diagonal(p::PauliBasis) = p.x == 0
is_diagonal(p::Pauli) = is_diagonal(PauliBasis{length(p)}(p.z, p.x))  # or just: p.x == 0

"""
    expectation_value(p::PauliBasis{N}, ket::Ket{N})

Computes ±1 if diagonal; 0 if non-diagonal
"""
function expectation_value(p::PauliBasis{N}, ket::Ket{N}) where N
    is_diagonal(p) || return 0.0
    return (-1)^count_ones(p.z & ket.v)
end

"""
    expectation_value(p::Pauli{N}, ket::Ket{N})

Same as above but with scalar
"""
function expectation_value(p::Pauli{N}, ket::Ket{N}) where N
    is_diagonal(p) || return 0.0
    return p.s * (-1)^count_ones(p.z & ket.v)
end

"""
    expectation_value(p::PauliSum{N}, ket::Ket{N})

Sum over terms
"""
function expectation_value(p::PauliSum{N}, ket::Ket{N}) where N
    expval = 0.0
    for (op, coeff) in p.ops
        expval += coeff * expectation_value(op, ket)
    end
    return expval
end

"""
    expectation_value_sign(p, ket)

Alias for expectation_value
"""
expectation_value_sign(p::Pauli{N}, ket::Ket{N}) where N = expectation_value(p, ket)
expectation_value_sign(p::PauliBasis{N}, ket::Ket{N}) where N = expectation_value(p, ket)
