##############################################################################
# This is a generic FHT Plan using BFFT (unnormalized backword FFT)

struct BFFTPlan{T,P,S} <: DHTPlan{T}
    bfftplan::P
    function BFFTPlan(bfftplan::Plan)
        T = eltype(bfftplan)
        P = typeof(bfftplan)
        S = fieldtypes(T)[1]
        return new{T,P,S}(bfftplan)
    end
end
size(p::BFFTPlan) = size(p.bfftplan)
fftdims(p::BFFTPlan) = fftdims(p.bfftplan)

function LinearAlgebra.mul!(
    y::AbstractArray{S}, p::BFFTPlan{T,P,S}, x::AbstractArray{S}
) where {T,P,S}
    fx = p.bfftplan * x
    y .= real(fx) .+ imag(fx)
    return nothing
end

function *(p::BFFTPlan{T,P,S}, x::AbstractArray{S}) where {T,P,S}
    z = similar(x)
    LinearAlgebra.mul!(z, p, x)
    return z
end

function plan_inv(p::BFFTPlan)
    return ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.bfftplan.region))
end

struct BFFTPlanInplace{T,P,S} <: DHTPlan{T}
    bfftplan::P
    function BFFTPlanInplace(bfftplan::Plan)
        T = eltype(bfftplan)
        P = typeof(bfftplan)
        S = fieldtypes(T)[1]
        return new{T,P,S}(bfftplan)
    end
end
size(p::BFFTPlanInplace) = size(p.bfftplan)
fftdims(p::BFFTPlanInplace) = fftdims(p.bfftplan)

function LinearAlgebra.mul!(
    y::AbstractArray{S}, p::BFFTPlanInplace{T,P,S}, x::AbstractArray{S}
) where {T,P,S}
    fx = p.bfftplan * x
    y .= real(fx) .+ imag(fx)
    return nothing
end

function *(p::BFFTPlanInplace{T,P,S}, x::AbstractArray{S}) where {T,P,S}
    LinearAlgebra.mul!(x, p, x)
    return x
end

function plan_inv(p::BFFTPlanInplace)
    return ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.bfftplan.region))
end