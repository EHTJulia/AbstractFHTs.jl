##############################################################################
# This is a generic FHT Plan using R2R transform

struct R2RPlan{T,P,S} <: DHTPlan{T}
    r2rplan::P
    function R2RPlan(r2rplan::Plan)
        T = eltype(r2rplan)
        P = typeof(r2rplan)
        S = fieldtypes(T)[1]
        return new{T,P,S}(r2rplan)
    end
end
size(p::R2RPlan) = size(p.r2rplan)
fftdims(p::R2RPlan) = fftdims(p.r2rplan)

function LinearAlgebra.mul!(
    y::AbstractArray{S}, p::R2RPlan{T,P,S}, x::AbstractArray{S}
) where {T,P,S}
    return LinearAlgebra.mul!(y, p.r2rplan, x)
end

function *(p::R2RPlan{T,P,S}, x::AbstractArray{S}) where {T,P,S}
    return p.r2rplan * x
end

function plan_inv(p::R2RPlan)
    return ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.r2rplan.region))
end

struct R2RPlanInplace{T,P,S} <: DHTPlan{T}
    r2rplan::P
    function R2RPlanInplace(r2rplan::Plan)
        T = eltype(r2rplan)
        P = typeof(r2rplan)
        S = fieldtypes(T)[1]
        return new{T,P,S}(r2rplan)
    end
end
size(p::R2RPlanInplace) = size(p.r2rplan)
fftdims(p::R2RPlanInplace) = fftdims(p.r2rplan)

function LinearAlgebra.mul!(
    y::AbstractArray{S}, p::R2RPlanInplace{T,P,S}, x::AbstractArray{S}
) where {T,P,S}
    return LinearAlgebra.mul!(y, p.r2rplan, x)
end

function *(p::R2RPlanInplace{T,P,S}, x::AbstractArray{S}) where {T,P,S}
    return LinearAlgebra.mul!(x, p, x)
end

function plan_inv(p::R2RPlanInplace)
    return ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.r2rplan.region))
end