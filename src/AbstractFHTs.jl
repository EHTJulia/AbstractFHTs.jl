module AbstractFHTs
import AbstractFFTs: Plan, eltype, size, fftdims, fftfreqs
import Base: size, ndims, length, eltype,
             *, inv, \, size
import LinearAlgebra

include("definitions.jl")

end # module