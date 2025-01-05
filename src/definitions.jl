# DHT plan where the inputs are an array of eltype T
abstract type DHTPlan{T} end


eltype(::Type{<:DHTPlan{T}}) where {T} = T


"""
    size(p::DHTPlan, [dim])

Return the size of the input of a plan `p`, optionally at a specified dimenion `dim`.
"""
size(p::DHTPlan, dim) = size(p)[dim]
ndims(p::DHTPlan) = length(size(p))
length(p::DHTPlan) = prod(size(p))::Int


"""
    scale_type(::FFT/DHTPlan)

Extract the element type of the scaling factor of an FHT/FFT plan.
"""
scale_type(::Type{<:DHTPlan{T}}) where {T} = fieldtypes(T)[1]
scale_type(P::Type{<:Plan}) = fieldtypes(eltype(P))[1]


# copy to a 1-based array, using circular permutation
function copy1(::Type{T}, x) where T
    y = Array{T}(undef, map(length, axes(x)))
    Base.circcopy!(y, x)
end

to1(x::AbstractArray) = _to1(axes(x), x)
_to1(::Tuple{Base.OneTo,Vararg{Base.OneTo}}, x) = x
_to1(::Tuple, x) = copy1(eltype(x), x)

# implementations only need to provide plan_X(x, region)
# for X in (:fft, :bfft, ...):
for f in (:fht, :ifht, :fht!, :ifht!)
    pf = Symbol("plan_", f)
    @eval begin
        $f(x::AbstractArray) = $f(x, 1:ndims(x))
        $f(x::AbstractArray, region) = (y = to1(x); $pf(y, region) * y)
        $pf(x::AbstractArray; kws...) = (y = to1(x); $pf(y, 1:ndims(y); kws...))
    end
end

"""
    plan_fht(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)

Pre-plan an optimized FHT along given dimensions (`dims`) of arrays matching the shape and type of `A`. 
(The first two arguments have the same meaning as for `fht`(@ref)).)
Returns an object `P` which represents the linear operator computed by the FHT, and 
which contains all of the information needed to compute `fht(A, dims)` quickly.
"""
plan_fht


"""
    plan_fht!(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)

Same as [`plan_fht`](@ref), but operates in-place on `A`.
"""
plan_fht!


"""
    plan_ifht(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)

Same as [`plan_fht`](@ref), but produces a plan that performs inverse transforms
[`ifht`](@ref).
"""
plan_ifht


"""
    plan_ifht!(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)

Same as [`plan_ifht`](@ref), but operates in-place on `A`.
"""
plan_ifht!

"""
    ifht!(A [, dims])

Same as [`ifht`](@ref), but operates in-place on `A`.
"""
ifht!

"""
    ifft(A [, dims])

Multidimensional inverse FHT.

A one-dimensional inverse FHT computes the one-dimensional inverse DHT (IDHT) defined by as
```math
\\operatorname{IDHT}(A)[k] = \\frac{1}{\\operatorname{length}(A)} \\operatorname{DHT}(A)[k]
```

A multidimensional inverse FHT simply performs this operation along each transformed dimension of `A`.
"""
ifht

"""
    fht!(A [, dims])

Same as [`fht`](@ref), but operates in-place on `A`.
"""
fht!

"""
    fht(A [, dims])

Performs a multidimensional Fast Hartley Transform (FHT) of the array `A`.
The optional `dims` argument specifies an iterable subset of dimensions 
(e.g. an integer, range, tuple, or array) to transform along.
Most efficient if the size of `A` along the transformed dimensions is a product of small
primes; see `Base.nextprod`. See also [`plan_fht()`](@ref) for even greater efficiency.

A one-dimensional FHT computes the one-dimensional discrete Hartley transform (DHT) as
defined by

```math
\\operatorname{DFT}(A)[k] =
  \\sum_{n=1}^{\\operatorname{length}(A)}
  \\cas\\left(+i\\frac{2\\pi
  (n-1)(k-1)}{\\operatorname{length}(A)} \\right) A[n],
```

where `\\cas` is the cosine-and-sine function, or alternatively called Hartley kernel, defined by

```
\\cas(x) = \\cos(x) + \\sin(x).
```

A multidimensional FHT simply performs this operation along each transformed dimension of `A`.
"""
function fht end


##############################################################################
# only require implementation to provide *(::DHTPlan{T}, ::Array{T})
*(p::DHTPlan{T}, x::AbstractArray) where {T} = p * copy1(T, x)

##############################################################################
# To support inv, \, and ldiv!(y, p, x), we require Plan subtypes
# to have a pinv::DFTPlan field, which caches the inverse plan, and which
# should be initially undefined.  They should also implement
# plan_inv(p) to construct the inverse of a plan p.

# hack from @simonster (in #6193) to compute the return type of plan_inv
# without actually calling it or even constructing the empty arrays.
_pinv_type(p::DHTPlan) = typeof([plan_inv(x) for x in typeof(p)[]])
pinv_type(p::DHTPlan) = eltype(_pinv_type(p))

function plan_inv end

inv(p::DHTPlan) =
    isdefined(p, :pinv) ? p.pinv::pinv_type(p) : (p.pinv = plan_inv(p))
\(p::DHTPlan, x::AbstractArray) = inv(p) * x
LinearAlgebra.ldiv!(y::AbstractArray, p::DHTPlan, x::AbstractArray) = LinearAlgebra.mul!(y, inv(p), x)

##############################################################################
# implementations only need to provide the forward FHT transform.
# ifht can be computed by scaling the forward transform.

struct ScaledDHTPlan{T,P,N} <: DHTPlan{T}
    p::P
    scale::N
    function ScaledDHTPlan(p::DHTPlan, scale)
        T = eltype(p)
        P = typeof(p)
        N = fieldtypes(T)[1]
        return new{T, P, N}(p, convert(N, scale))
    end
end
ScaledDHTPlan(p::ScaledDHTPlan, α::Number) = ScaledDHTPlan(p.p, p.scale * α)

size(p::ScaledDHTPlan) = size(p.p)
fftdims(p::ScaledDHTPlan) = fftdims(p.p)
fftfreq(p::ScaledDHTPlan) = fftfreq(p.p)

*(p::ScaledDHTPlan, x::AbstractArray) = LinearAlgebra.rmul!(p.p * x, p.scale)
*(α::Number, p::DHTPlan) = ScaledDHTPlan(p, α)
*(p::DHTPlan, α::Number) = ScaledDHTPlan(p, α)

function normalization(::Type{T}, sz, region) where T
    one(T) / mapreduce(r -> Int(sz[r])::Int, *, region; init=1)::Int
end
normalization(X, region) = normalization(real(eltype(X)), size(X), region)

plan_ifht(x::AbstractArray, region; kws...) =
    ScaledDHTPlan(plan_fht(x, region; kws...), normalization(x, region))

plan_ifht!(x::AbstractArray, region; kws...) =
    ScaledDHTPlan(plan_fht!(x, region; kws...), normalization(x, region))

plan_inv(p::ScaledDHTPlan) = ScaledDHTPlan(plan_inv(p.p), inv(p.scale))

# Don't cache inverse of scaled plan (only inverse of inner plan)
inv(p::ScaledDHTPlan) = ScaledDHTPlan(inv(p.p), inv(p.scale))

LinearAlgebra.mul!(y::AbstractArray, p::ScaledDHTPlan, x::AbstractArray) =
    LinearAlgebra.lmul!(p.scale, LinearAlgebra.mul!(y, p.p, x))


##############################################################################
# This is a generic FHT Plan

struct FHTPlan{T,P,S} <: DHTPlan{T}
    bfftplan::P
    function FHTPlan(bfftplan::Plan)
        T = eltype(bfftplan)
        P = typeof(bfftplan)
        S = fieldtypes(T)[1]
        new{T,P,S}(bfftplan)
    end
end
size(p::FHTPlan) = size(p.bfftplan)
fftdims(p::FHTPlan) = fftdims(p.bfftplan)
fftfreq(p::FHTPlan) = fftfreq(p.bfftplan)

function LinearAlgebra.mul!(y::Array{S}, p::FHTPlan{T,P,S}, x::Array{S}) where {T,P,S}
    fx = p.bfftplan * x
    y .= real(fx) .+ imag(fx)
end

function *(p::FHTPlan{T,P,S}, x::Array{S}) where {T,P,S}
    z = similar(x)
    LinearAlgebra.mul!(z, p, x)
    return z
end

plan_inv(p::FHTPlan) = ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.bfftplan.region))

struct FHTPlanInplace{T,P,S} <: DHTPlan{T}
    bfftplan::P
    function FHTPlan(bfftplan::Plan)
        T = eltype(bfftplan)
        P = typeof(bfftplan)
        S = fieldtypes(T)[1]
        new{T,P,S}(bfftplan)
    end
end
size(p::FHTPlanInplace) = size(p.bfftplan)
fftdims(p::FHTPlanInplace) = fftdims(p.bfftplan)
fftfreq(p::FHTPlanInplace) = fftfreq(p.bfftplan)

function LinearAlgebra.mul!(y::Array{S}, p::FHTPlanInplace{T,P,S}, x::Array{S}) where {T,P,S}
    fx = p.bfftplan * x
    y .= real(fx) .+ imag(fx)
end

function *(p::FHTPlanInplace{T,P,S}, x::Array{S}) where {T,P,S}
    LinearAlgebra.mul!(x, p, x)
    return x
end

plan_inv(p::FHTPlanInplace) = ScaledDHTPlan(p, normalization(Array{Bool}(undef, size(p)), p.bfftplan.region))