using Test
using CalciumImagingAnalyses
using CalciumImagingAnalyses: CSV
using TimesDates
using LinearAlgebra
using Dates
using TimeZones
using DataFrames
const CaImAn = CalciumImagingAnalyses


@testset "Align matrices" begin
    q_1,_ = qr(randn(13,13))
    q_2,_ = qr(randn(13,13))
    w_1 = permutedims(q_1[:,1:4]) 
    w_2 = permutedims(q_2[:,1:7])

    CC = w_1*w_2'
    u,s,v = svd(CC)
    R = u*v'
    w_1 = R*w_2

    R = CaImAn.align_subspaces(w_1,w_2)
    @test norm(w_1-R*w_2) < sqrt(eps())

end
@testset "Window length" begin
    tmin = -Millisecond(200)
    tmax = Millisecond(200) 
    ee = ZonedDateTime("2023-11-17T11:43:31.452+00:00")
    ts = datetime2unix(ee.utc_datetime)
    timestamps = range(ee - Second(10), stop=ee+Second(10), step=Millisecond(100))
    n = CaImAn.get_window_length(timestamps, ee, tmin, tmax)
    @test n == 5 
end

@testset "Resolution" begin
    t = [1,2,2,3,3,4]
    dt = CaImAn.get_time_resolution(t)
    @test dt == 0.5
end

@testset "Rescale" begin
    t = CaImAn.rescale_time([1,3,3,4,4,7])
    @test t ≈ [1.0, 3.0, 3.5, 4.0, 4.5, 7.0]
end

@testset "Rewarded Pokes" begin
    poke_idx = [1, 2, 4, 7]
    rewarded_idx = [3,5]
    rewarded_poke_idx = CaImAn.get_rewarded_pokes(rewarded_idx, poke_idx)
    @test rewarded_poke_idx == [2,3]
end

@testset "Check time" begin 
    
    bdata = CSV.read("$(@__DIR__)/../../data/E10B_6Days/Behavior/aRR5_D3/2023/11/FED002_112423_00.csv",DataFrame;header=1);
    cell_data, timestamps = CaImAn.get_cell_data("$(@__DIR__)/../../data/E10B_6Days/Calcium_Imaging_Data/E10B_6D.csv");
    t0c = timestamps[1]
    t0b = bdata[!, "Pi_Time"][1]
    @show t0c t0b
end