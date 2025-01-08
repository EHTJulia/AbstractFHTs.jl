using StableRNGs
using Test

import FFTW: plan_bfft, ESTIMATE, plan_r2r, DHT
import AbstractFastHartleyTransforms:
    R2RPlan,
    R2RPlanInplace,
    BFFTPlan,
    BFFTPlanInplace,
    plan_fht,
    plan_fht!,
    plan_ifht,
    plan_ifht!,
    fht,
    fht!,
    ifht,
    ifht!,
    *,
    inv,
    \,
    eltype,
    size,
    length,
    ndims,
    fftdims

function plan_fht(
    A, dims; flags=ESTIMATE, timelimit=Inf, use_r2r::Bool=ndims(A) == 1 ? true : false
)
    if use_r2r
        r2rplan = plan_r2r(A, DHT, dims; flags=flags, timelimit=timelimit)
        return R2RPlan(r2rplan)
    else
        bfftplan = plan_bfft(A, dims; flags=flags, timelimit=timelimit)
        return BFFTPlan(bfftplan)
    end
end

function plan_fht!(
    A, dims; flags=ESTIMATE, timelimit=Inf, use_r2r::Bool=ndims(A) == 1 ? true : false
)
    if use_r2r
        r2rplan = plan_r2r(A, DHT, dims; flags=flags, timelimit=timelimit)
        return R2RPlanInplace(r2rplan)
    else
        bfftplan = plan_bfft(A, dims; flags=flags, timelimit=timelimit)
        return BFFTPlanInplace(bfftplan)
    end
end

function check_inv(xinv, xorg, xinvref)
    @test xorg ≈ xinv
    @test xinvref == xinv
end

@testset "AbstractFastHartleyTransforms.jl" begin
    rng = StableRNG(0)
    x1d = rand(rng, Float64, 128)
    x2d = rand(rng, Float64, 128, 128)
    x3d = rand(rng, Float64, 16, 16, 16)
    for x in (x1d, x2d, x3d)
        # make FastHartleyTransform plans
        p1 = plan_fht(x)
        p2 = plan_fht!(x)

        # make inverse FastHartleyTransform plans
        ip1 = plan_ifht(x)
        ip2 = plan_ifht!(x)
        ip3 = inv(p1)
        ip4 = inv(p2)

        # check if the forward FPT plans are consistent
        y1 = p1 * x
        @test y1 == p2 * deepcopy(x)
        @test y1 == fht(x)
        @test y1 == fht!(deepcopy(x))

        # check inverse of inverse is the same
        x1 = ip1 * y1
        check_inv(x1, x, x1)
        check_inv(ip2 * deepcopy(y1), x, x1)
        check_inv(ip3 * y1, x, x1)
        check_inv(ip4 * deepcopy(y1), x, x1)
        check_inv(ifht(y1), x, x1)
        check_inv(ifht!(deepcopy(y1)), x, x1)
        check_inv(p1 \ y1, x, x1)
        check_inv(p2 \ deepcopy(y1), x, x1)

        # check methods
        for p in (p1, p2)
            @test eltype(p) == eltype(p.bfftplan)
            @test ndims(p) == ndims(p.bfftplan)
            @test length(p) == length(p.bfftplan)
            @test fftdims(p) == fftdims(p.bfftplan)
        end

        for p in (ip1, ip2, ip3, ip4)
            @test eltype(p) == eltype(p1)
            @test ndims(p) == ndims(p1)
            @test length(p) == length(p1)
            @test fftdims(p) == fftdims(p1)
        end
    end
end
