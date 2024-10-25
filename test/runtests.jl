using Test

include("../src/analysis.jl")

@testset "Window length" begin
    tmin = -Millisecond(200)
    tmax = Millisecond(200) 
    ee = ZonedDateTime("2023-11-17T11:43:31.452+00:00")
    ts = datetime2unix(ee.utc_datetime)
    timestamps = range(ee - Second(10), stop=ee+Second(10), step=Millisecond(100))
    n = get_window_length(timestamps, ee, tmin, tmax)
    @test n == 5 
end

@testset "Resolution" begin
    t = [1,2,2,3,3,4]
    dt = get_time_resolution(t)
    @test dt == 0.5
end

@testset "Rescale" begin
    t = rescale_time([1,3,3,4,4,7])
    @test t ≈ [1.0, 3.0, 3.5, 4.0, 4.5, 7.0]
end

@testset "Rewarded Pokes" begin
    poke_idx = [1, 2, 4, 7]
    rewarded_idx = [3,5]
    rewarded_poke_idx = get_rewarded_pokes(rewarded_idx, poke_idx)
    @test rewarded_poke_idx == [2,3]
end