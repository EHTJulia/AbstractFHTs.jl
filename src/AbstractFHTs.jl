module AbstractFHTs
import AbstractFFTs: Plan, eltype, size, fftdims
import Base: size, ndims, length, eltype,
             *, inv, \, size
import LinearAlgebra

include("definitions.jl")

end # module