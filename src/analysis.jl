using Makie
using GLMakie
using NanoDates
using TimesDates
using Dates
using TimeZones
using LinearAlgebra
using MultivariateStats
using CSV
using StatsBase
using DataFrames

plot_theme = theme_minimal()
plot_theme.Axis.xticksvisible[] = true
plot_theme.Axis.yticksvisible[] = true


function get_cell_data(cdata::DataFrame)
    nrows,ncols = size(cdata)
    cell_data = zeros(nrows-1, ncols-1)
    k = 1
    for col in eachcol(cdata)
        if col[1] == " accepted"
            cell_data[:,k] = parse.(Float64, col[2:end])
            k += 1
        end
    end
    cell_data[:,1:k-1]
end

function get_cell_data(fname::String)
    cdata = CSV.read(fname,DataFrame;header=1)
    if all(isdigit, cdata[!, " "][2])
        _timestamps = parse.(UInt64, cdata[!, " "][2:end])
    else
        _timestamps = parse.(Float64, cdata[!, " "][2:end])
    end
    # timestamps are recorded at 100ms, but only with 1s precision, so we need to convert
    # them to their actual precision
    timestamps = ZonedDateTime.(unix2datetime.(rescale_time(_timestamps)), tz"UTC")
    # convert to UTC
    cell_data = get_cell_data(cdata)
    cell_data, timestamps
end

function get_rewarded_pokes(reward_idx::Vector{T}, poke_idx::Vector{T}) where T <: Integer
    d = reward_idx .- permutedims(poke_idx)
    d[d.<0] .= T(10000)
    ridx = dropdims(argmin(d,dims=2),dims=2)
    [_r[2] for _r in ridx]
end

function get_behavioural_data(fname::String;do_correct_time_zone=true)
    bdata = CSV.read(fname, DataFrame;header=1)

    event_time_s = convert.(String, bdata[!, "Pi_Time"])
    event_time = [TimeDate(replace(et, " "=>"T")) for et in event_time_s]
    if do_correct_time_zone
        # this is hackish
        ztd = ZonedDateTime.(event_time .+ Hour(8), tz"Asia/Singapore")
        ztd = astimezone.(ztd, tz"UTC")
    else
        ztd = ZonedDateTime(event_time, tz"UTC")
    end
    # get the pokes
    poke_event_idx = findall(x->ifelse(ismissing(x), false, isfinite(x)), bdata[!, "Poke_Time"])
    poke_time = ztd[poke_event_idx]
    poke_type = bdata[!,"Event"][poke_event_idx]
    if eltype(bdata[!, "Retrieval_Time"]) == Float64
        retrieval_idx = findall(isfinite, bdata[!, "Retrieval_Time"])
    else
        retrieval_idx = findall(x->isnothing(x) ? false : isfinite(x), tryparse.(Float64, bdata[!, "Retrieval_Time"]))
    end
    retrieval_time = ztd[retrieval_idx]
    # find the unrewarded pokes
   
    rewarded_poke_idx = get_rewarded_pokes(retrieval_idx, poke_event_idx)
    unrewarded_poke_idx = setdiff(1:length(poke_event_idx), rewarded_poke_idx)
    (rewarded_poke_time=poke_time[rewarded_poke_idx], unrewarded_poke_time=poke_time[unrewarded_poke_idx],
    retrieval_time=retrieval_time, rewarded_poke_type=poke_type[rewarded_poke_idx],
    unrewarded_poke_type=poke_type[unrewarded_poke_idx])
end

function get_time_resolution(timestamps::Vector{T}) where T <: Integer
    t0 = timestamps[1]
    idx0 = findfirst(timestamps .> t0)
    t0 = timestamps[idx0]
    idx1 = findlast(timestamps .<= t0)
    t1 = timestamps[idx1]
    dt = 1.0/(idx1-idx0+1)
end

rescale_time(timestamps::Vector{T}) where T <: Real = timestamps

function rescale_time(timestamps::Vector{T}) where T <: Integer
    dt = get_time_resolution(timestamps)
    n = length(timestamps)
    t = zeros(n)
    t[1] = timestamps[1]
    # hackish
    i = 2
    while i <= n
        if timestamps[i] == timestamps[i-1]
            t[i] = t[i-1] + dt
        else
            t[i] = timestamps[i]
        end
        i+= 1
    end
    if t[2]-t[1] == 1.0
        t[1] = t[2] - dt
    end
    t
end

function get_window_length(timestamps::AbstractVector{ZonedDateTime},event::ZonedDateTime, tmin::Period, tmax::Period)
    idx0 = searchsortedfirst(timestamps, event+tmin)
    idx1 = searchsortedlast(timestamps, event+tmax) 
    n = idx1-idx0+1
end

function alignto(X::Matrix{T}, timestamps, event::Vector{ZonedDateTime};tmin=-Second(1), tmax=Second(1)) where T <: Real
    dt = timestamps[2] - timestamps[1]
    _,nc = size(X)
    nt = length(event)
    n = get_window_length(timestamps, first(event), tmin,tmax)
    Xm = zeros(T, n+1, nt, nc)
    x =range(Millisecond(tmin),length=n, step=Millisecond(dt)) 
    for i in axes(X,2) 
        Xm[:,:,i],x = alignto(view(X,:,i), timestamps, event;tmin=tmin,tmax=tmax)
    end
    Xm,x
end

function alignto(X::AbstractVector{T}, timestamps::Vector{T2}, event::Vector{T2};tmin=-Second(1), tmax=Second(1)) where T <: Real where T2 <: ZonedDateTime
    dt = timestamps[2]-timestamps[1]
    idx0 = searchsortedfirst(timestamps, event[1]+tmin)
    idx1 = searchsortedlast(timestamps, event[1]+tmax) 
    n = idx1-idx0+1
    nt = length(event)
    Xm = zeros(T, n+1, nt)
    for (i,ts) in enumerate(event)
        idx0 = searchsortedfirst(timestamps, ts+tmin)
        idx1 = searchsortedlast(timestamps, ts+tmax)
        Xm[1:idx1-idx0+1,i] = X[idx0:idx1]
    end
    Xm, range(Millisecond(tmin),length=n, step=Millisecond(dt))
end

function get_per_cell_average(X::Array{T,3}, bins::AbstractVector{T2},window::Tuple{T3,T3}) where T2 <: Period where T3 <: Period where T <: Real
    nbins, nt, nc = size(X)
    idx0 = searchsortedfirst(bins, window[1])
    idx1 = searchsortedlast(bins, window[2])
    Xm = permutedims(dropdims(mean(X[idx0:idx1,:,:],dims=1),dims=1))
    # subtract baseline from each cell
    Xm .-= mean(Xm,dims=2)
end

function get_subspace(X::Array{T,3}, bins::AbstractVector{T2},window::Tuple{T3,T3};pratio=0.5) where T2 <: Period where T3 <: Period where T <: Real
    # find a subspace via svd
    Xm = get_per_cell_average(X, bins, window)
    get_subspace(Xm;pratio=pratio)
end

function get_subspace(Xm::Matrix{T};pratio=0.5) where T <: Real
    u,s,v = svd(Xm)
    pr = cumsum(s)/sum(s)
    ii = findfirst(pr .>= pratio) 
    u[:,1:ii], s[1:ii], sum(s[1:ii])/sum(s)
end

function project(w::Matrix{T}, X::Array{T,3}) where T <: Real
    Y = zeros(T, size(X,1), size(X,2), size(w,2))
    for i in axes(X,1)
        Y[i,:,:] = X[i,:,:]*w
    end
    Y
end

function plot_3d_representation(Zp::Matrix{T}, Zn::Matrix{T}, lp::Vector{T2}, ln::Vector{T2}) where T <: Real where T2 <: Integer 
    fig = Figure()
    ax = Axis3(fig[1,1])
    markers =  [:circle, :star4, :diamond]
    scatter!(ax, Zp[1,:], Zp[2,:], Zp[3,:],marker=markers[lp],label="Left poke")
    scatter!(ax, Zn[1,:], Zn[2,:], Zn[3,:],marker=markers[ln], label="Right poke")
    ll = Legend(fig[1,1], [MarkerElement(marker=:circle, color=Cycled(1)),
                           MarkerElement(marker=:circle, color=Cycled(2)),
                           MarkerElement(marker=:circle, color=:gray),
                           MarkerElement(marker=:star4, color=:gray),
                           MarkerElement(marker=:diamond, color=:gray)],
                           ["Left poke","Right poke","Day 1", "Day 2", "Day 3"],
                           tellwidth=false, tellheight=false, valign=:top, halign=:left)
    fig
end

function train_decoder(X::Matrix{T}, label::Vector{T2}) where T2 <: String where T <: Real
    Xp = X[:,label.=="Left"]
    Xn = X[:, label.=="Right"]
    train_decoder(Xp, Xn)
end

function train_decoder(Zp::Matrix{T}, Zn::Matrix{T}) where T <: Real
    lda = fit(LinearDiscriminant, Zp, Zn;covestimator=SimpleCovariance())
end

function decode(X::Matrix{T}, timestamps::Vector{ZonedDateTime}, poke_time::Vector{ZonedDateTime}, poke_type::AbstractVector{T2}) where T <: Real where T2
    # identify breaks
    breaks = findall(diff(timestamps).>Second(5))
    bidx = [0;breaks;length(timestamps)]
    Xd = zeros(T,size(X,2),length(poke_time))
    Y = fill!(similar(X), zero(T))
    Y .= X
    offset = 0
    for (b1,b2) in zip(bidx[1:end-1], bidx[2:end])
        _X = X[b1+1:b2,:]
        pidx0 = searchsortedfirst(poke_time, timestamps[b1+1])
        pidx1 = searchsortedlast(poke_time, timestamps[b2])
        Xm,x = alignto(_X, timestamps[b1+1:b2], poke_time[pidx0:pidx1])
        _X1 = get_per_cell_average(Xm, x, (-Second(1), Second(0)))
        Xd[:,offset+1:offset+(pidx1-pidx0+1)] = _X1
        offset += pidx1-pidx0+1
        # substract mean
        Y[b1+1:b2,:] .-= mean(X[b1+1:b2,:],dims=1)
    end
    # get a subspace spanning 50% (arbitrary)
    w,s,cc = get_subspace(Xd)
    w = w[:,1:10]
    Xp = Xd[:,poke_type.=="Left"]
    Xn = Xd[:, poke_type.=="Right"] 
    Zp = w'*Xp
    Zn = w'*Xn
    lda = fit(LinearDiscriminant, Zp, Zn)
    cq_p = predict(lda, Zp)
    cq_n = predict(lda, Zn)
    @show mean(cq_p) mean(cq_n.==0)
    Z = Y*w
    q = evaluate(lda, permutedims(Z))
end

function plot_decoder(lda, Zp::Matrix{T}, Zn::Matrix{T}) where T <: Real
    qp = evaluate(lda, Zp)
    qn = evaluate(lda, Zn)
    @show mean(qp.>0.0) mean(qn.<0.0)
    hp = fit(Histogram, qp)
    hn = fit(Histogram, qn)
    with_theme(plot_theme) do
        fig = Figure()
        ax1 = Axis(fig[1,1])
        ax2 = Axis(fig[1,2])
        scatter!(ax1, rand(length(qp)), qp,label="Left poke")
        scatter!(ax1, rand(length(qn)), qn, label="Right poke")

        lines!(ax2, hp.weights, hp.edges[1][1:end-1])
        lines!(ax2, hn.weights, hn.edges[1][1:end-1])
        linkyaxes!(ax1, ax2)
        for ax in [ax1, ax2]
            hlines!(ax, 0.0, linestyle=:dot, color=:black)
        end
        axislegend(ax1,valign=:bottom)
        ax1.xticklabelsvisible = false
        ax1.xticksvisible = false
        ax1.bottomspinevisible = false 
        ax2.yticklabelsvisible = false
        ax2.yticksvisible = false
        ax2.leftspinevisible = false
        ax1.ylabel = "Discriminant value"
        ax2.xlabel = "Counts"
        fig
    end

end

function plot_timeseries(X::Vector{T}, timestamps::Vector{ZonedDateTime}, event::Vector{ZonedDateTime}, event_type::AbstractVector{T2};do_animate=false,animation_speed=1,moviefile::Union{String,Nothing}=nothing) where T <: Real where T2
    _etypes = unique(event_type)
    _etypes = sort(_etypes,by=a->length(a))
    # create colors
    _colors = Makie.wong_colors()
    event_color = fill(_colors[1], length(event_type))
    for (i,et) in enumerate(_etypes[2:end])
        event_color[event_type.==et] .= _colors[i+1]
    end
    # identifiy breaks
    breaks = findall(diff(timestamps).>Second(5))
    # shorten the breaks
    _timestamps = copy(timestamps)
    _event = copy(event)
    for b in breaks
        Δt=Second(5) - (_timestamps[b+1]-_timestamps[b])
        _timestamps[b+1:end] .+= Δt
        # also shift the events
        eidx = _event .>= _timestamps[b+1] 
        _event[eidx] .+= Δt
    end
    # check that nothing changed
    evalues = zeros(T, length(event))
    for (ii,ee) in enumerate(event)
        idx = searchsortedfirst(timestamps, ee)
        evalues[ii] = X[idx]
    end
    evalues_shifted = zeros(T, length(event))
    for (ii,ee) in enumerate(_event)
        idx = searchsortedfirst(_timestamps, ee)
        evalues_shifted[ii] = X[idx]
    end
    @assert evalues ≈ evalues_shifted
    bidx = [0;breaks;length(timestamps)]
    # convert to seconds since Makie doesn't really handle dates very well    
    #tt = [datetime2unix(ts.utc_datetime) for ts in _timestamps]
    tt = [ts.utc_datetime for ts in _timestamps]
    tt0 = [ts.utc_datetime for ts in timestamps]
    tt1 = [datetime2unix(ts.utc_datetime) for ts in _timestamps]
    et = [datetime2unix(_et.utc_datetime) for _et in _event]
    t0 = tt1[1]
    t1 = tt1[end]
    window = (tt1[1], tt1[1]+60*10)
    with_theme(plot_theme) do
        fig = Figure(size=(1200,600))
        ax = Axis(fig[1,1],ypanlock=true,yzoomlock=true)
        for (b1,b2) in zip(bidx[1:end-1], bidx[2:end])
            lines!(ax, tt1[b1+1:b2], X[b1+1:b2],color=:gray)
        end
        vlines!(ax, tt1[breaks],color=:black, linewidth=2.0)
        vlines!(ax, et,color=event_color)
        hlines!(ax, 0.0, color=:black, linestyle=:dot)
        minorticks = [(tt1[b1+1]+60):60:tt1[b2] for (b1,b2) in zip(bidx[1:end-1],bidx[2:end])]
        xtickpositions = Float64[]
        xticklabels = String[]
        for (b1,b2) in zip(bidx[1:end-1],bidx[2:end])
            push!(xtickpositions, tt1[b1+1])
            push!(xticklabels, string(tt0[b1+1]))
            minorticks = tt1[b1+1]+60:60:tt1[b2]
            for (k,mt) in enumerate(minorticks)
                push!(xtickpositions, mt)
                push!(xticklabels, "+$(k)min")
            end
        end
        ll = Legend(fig[1,1], [LineElement(color=c) for c in _colors[1:length(_etypes)]],
                              _etypes, tellwidth=false, tellheight=false, valign=:top, 
                              halign=:left,framevisible=true,margin=(10,10,10,10))

        ax.ylabel = "Decision function"
        #ax.xticks = (tt1[bidx[1:end-1].+1], string.(tt0[bidx[1:end-1].+1]))
        ax.xticks = (xtickpositions, xticklabels)
        #minorticklabels = [["+$(k)min" for k in 1:length(_tt)] for _tt in minorticks] 
        #@show minorticklabels
        #ax.xminorticks = ([minorticks...,],[minorticklabels...,])
        ax2 = Axis(fig[1,2])
        linkyaxes!(ax, ax2)
        #lines!(ax2, hp.weights, hp.edges[1][1:end-1],color=_colors[1])
        #lines!(ax2, hn.weights, hn.edges[1][1:end-1],color=_colors[2])
        y = evalues[event_type.=="Left"] 
        rainclouds!(ax2, fill(1.0, length(y)),y,color=_colors[1],markersize=5px,plot_boxplots=false,gap=0.3)
        y = evalues[event_type.=="Right"] 
        rainclouds!(ax2, fill(2.0, length(y)),y,color=_colors[2], markersize=5px, plot_boxplots=false, gap=0.3)
        ax2.xticklabelsvisible = false
        colsize!(fig.layout, 2, Relative(0.05))
        ax2.yticklabelsvisible = false
        if do_animate
            if moviefile !== nothing
                record(fig, moviefile, 1:length(tt1)-600) do i
                    xlims!(ax, tt1[i],tt1[i]+600)
                end
            else
                display(fig)
                xlims!(ax, window)
                while window[2] <= tt1[end]
                    t0,t1 = window
                    t0 += animation_speed 
                    t1 += animation_speed 
                    window = (t0, t1)
                    xlims!(ax,window)
                    sleep(0.01)
                end
            end
        end
        fig
    end
end

function animate_timeseries(X::Matrix{T},tevents::Vector{Vector{T2}}, t=1:size(X,2);npoints::Int64=size(X,2),do_animate=false) where T <: Real where T2 <: Real
    with_theme(plot_theme) do
        fig = Figure()
        ax = Axis(fig[1,1])
        ylims!(ax, extrema(X)...)
        tt = Observable(t[1])
        points = [Observable([Point2f(t[j],X[i,j]) for j in 1:npoints]) for i in 1:size(X,1)]
        elines = [Observable(tevents[i][t[1] .<= tevents[i] .< t[npoints]]) for i in 1:length(tevents)]
        # hack
        for j in 1:length(elines)
            if isempty(elines[j].val)
                elines[j].val = [NaN]
            end
        end
        on(tt) do _t
            idx0 = searchsortedfirst(t, _t)
            idx1 = idx0+npoints-1
            if idx1 <= size(X,2)

                xlims!(ax, t[idx0],t[idx1])
                for j in 1:size(X,1)
                    points[j][] = [Point2f(t[i], X[j,i]) for i in idx0:idx1]
                    elines[j].val = tevents[j][t[idx0] .<= tevents[j].<t[idx1]]
                    if isempty(elines[j].val)
                        elines[j].val = [NaN]
                    end
                    elines[j][] = elines[j].val
                end
            end
        end

        for j in 1:length(points)
            lines!(ax, points[j])
            vlines!(ax, elines[j],color=Cycled(j))
        end

        on(events(fig).scroll, priority=5) do (dx,dy)
            tt[] = tt[] + sign(dx)
            return Consume()
        end
        if do_animate
            display(fig)
            for _t in t[1:end-npoints+1]
                tt[] = _t
                sleep(0.01)
            end
        end
        fig 
    end
end