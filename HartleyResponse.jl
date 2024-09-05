using NeuroAnalysis,FileIO,JLD2,Images,StatsBase,StatsPlots,ProgressMeter,DataFrames,XLSX,CircStats, GaussianMixtures
using Plots,ImageTransformations,Statistics,Distributions,LsqFit,PDFIO,HistogramThresholding,HypothesisTests

ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")
getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function CustomCas(x)
    return sin(x) + cos(x)
end

function cond2hartley(cond)
    kcond = select(cond,[:Ori,:SpatialFreq,:SpatialPhase]=>ByRow((o,s,p)->map(i->round(i,digits=4),sin2cas(deg2rad(o),s,p)))=>AsTable)
    kx = levels(kcond.kx)
    ky = levels(kcond.ky)
    phase = levels(kcond.phase)
    #print("kx = $kx; \n ky = $ky \n \n ")
    kx == ky || error("kx and ky values are not match, may not be a valid hartley subspace")

    hi = [(findfirst(i->i==r.ky,ky),findfirst(i->i==r.kx,kx),findfirst(i->i==r.phase,phase)) for r in eachrow(kcond)]
    (;ky,kx,phase,hi)
end

function joinhartley(indir;ccode=ccode,btw=-100:0,rtw=10:130,datafile="sta.jld2")

    stas=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            push!(stas,load(joinpath(root,datafile)))
            p = joinpath(root,datafile)
            print("Pushed $p. \n")
        end
    end
    ugood1 = stas[1]["ugood"]
    siteid = stas[1]["siteid"]
    hcond = map(i->cond2hartley(i["xcond"]),stas)
    nkx = length(hcond[1].kx)
    nky = length(hcond[1].ky)
    nphase = length(hcond[1].phase)

    # sizepx = stas[1]["sizepx"]
    # sizedeg = stas[1]["sizedeg"]
    delays = stas[1]["delays"]
    # ppd = (sizepx[1]-1)/sizedeg[1]
    # xi = stas[1]["xi"]
    # notxi = setdiff(1:prod(sizepx),xi)
    # cii=CartesianIndices(sizepx)
    bdi = filter!(i->!isnothing(i),indexin(btw,delays))
    rdi = filter!(i->!isnothing(i),indexin(rtw,delays))


    # Saving into dataset
    dataset = Dict("delays"=>delays,"bdi"=>bdi,"rdi"=>rdi,"siteid"=>siteid,"kx"=>hcond[1].kx,"ky"=>hcond[1].ky,"phase"=>hcond[1].phase)
    colors = [i["exenv"]["color"] for i in stas]
    colorcodes = map(i->getccode(i,ccode),colors)
    corder = sortperm(colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    dataset["eye"] = [stas[i]["exenv"]["eye"] for i in corder]
    dataset["minmaxcolor"] = [(stas[i]["exenv"]["mincolor"],stas[i]["exenv"]["maxcolor"]) for i in corder]
    ugood=Dict();uy = Dict();uresponsive=Dict()

    uns = map(i->length(keys(stas[i]["uy"])),corder)
    @info "Number of units for each test: $uns"
    uids = mapreduce(i->keys(stas[i]["uy"]),intersect,corder) # only include units that spikes in all tests
    for u in uids
        cuy = zeros(nky,nkx,nphase,length(delays),length(corder))
        for j in eachindex(corder)
            cy = stas[corder[j]]["uy"][u]
            hi = hcond[corder[j]].hi
            @views for d in eachindex(delays)
                foreach((i,r)->cuy[i...,d,j]=r,hi,cy[d,:])
            end
        end
        ugood[u] = ugood1[u]
        uy[u] = cuy
        uresponsive[u] = false
    end
    dataset["ugood"] = ugood
    dataset["uy"] = uy
    dataset["uresponsive"] = uresponsive
    dataset["stas"] = stas
    return dataset
end

"""
    ArcBand(M, θ1, θ2, d1, d2, plot=false)

    Returns list of indices in MxM matrix within specified θ & d range
        Note: d is in sin space
"""
function ArcBand(M, θ1, θ2, d1, d2, plot=false)
    c = (M/2, M/2)
    d = [sqrt(((i-1/2) - c[1])^2 + ((j-1/2) - c[2])^2) for i in 1:M, j in 1:M]
    a = [atan((j-1/2) - c[2], (i-1/2) - c[1]) for i in 1:M, j in 1:M]
    a .+= (a .< 0) .* (2pi)
    a .= mod.(a - (ones(size(a))*pi/2), 2pi) # norm. zero angle like in blur

    # get indices
    mask = (d1 .< d) .& (d .< d2) .& 
           ((θ1 .<= a .<= θ2) .| (θ1 .+ 2pi .<= a .<= θ2 .+ 2pi) 
           .| (θ1 .- 2pi .<= a .<= θ2 .- 2pi))

        indices = findall(mask)
        indices = [(I[1], I[2]) for I in indices]

    if (plot)
        #append!(indices, [(34, 34), (34, 35), (35, 35), (35, 34)])

        mask2 = zeros(Bool, M, M)
        comp = [(i,j) for i in 1:M, j in 1:M]
        mask2 .|= [any(comp[i, j] == tup for tup in indices) for i in 1:size(mask2, 1), j in 1:size(mask, 2)]

        # Plotting the heatmap of mask
        θ1 = round(θ1, sigdigits=3)
        θ2 = round(θ2, sigdigits=3)
        d1 = round(d1, sigdigits=5)
        d2 = round(d2, sigdigits=5)
        dθ = round(rad2deg(abs(θ2 - θ1)), sigdigits = 3)
        dd = round(abs(d1 - d2), sigdigits = 3)
        p = heatmap(mask2, color=:blues, xlabel="Index j", ylabel="Index i", 
        title="dθ: $dθ, dd: $dd", ratio=1,
        yflip = true)
        vline!([34], c=:black)
        hline!([34], c=:black)
        display(p)
    end

    return indices
end

""" UNUSED:
    getAvgDict(blurRes::Dict, times::Dict, M::Int64, compList::Vector)

    Returns a dictionary corresponding to the averaged response over the two (four total) phases.
"""
function getAvgDict(parameters::Dict)
    compDict = Dict()
    blurRes, times, M, compList = getindex.(Ref(parameters), ["bR_f", "tD", "M", "cV"])
    for key in compList
        τ, i, j, p1 = times[key] 
        p2 = p1 == 1 ? 2 : 1
        m1 = blurRes[key][:,:,p1,τ,1]
        top_flip_m1 = m1[Int(M/2):-1:1, end:-1:1]
        bott_m1 = m1[Int(M/2)+1:1:end, :]
        m2 = blurRes[key][:,:,p2,τ,1]
        top_flip_m2 = m2[Int(M/2):-1:1, end:-1:1]
        bott_m2 = m2[Int(M/2)+1:1:end, :]
        
        # note compDict[key] has upside-down responses
        compDict[key] = (top_flip_m1 + bott_m1 + top_flip_m2 + bott_m2)./4
    end
    parameters["cDict"] = compDict
    return compDict
end

"""
    getArcDict(M, mF, vmtx, fRatio=1/4)

    Returns dictionary of [2D tup]=>list of indices in ROI
        given size of matrix M, maximum frequency F, variance
        matrix vmtx, and frequency ratio fRatio
"""
function getArcDict(parameters::Dict)
    M, mF, aV, fRatio = getindex.(Ref(parameters), ["M", "mF", "aV", "fRatio"])
    ArcDict = Dict()
    for i=1:M, j=1:M
        dist = sqrt((((i-1/2)-(M/2))^2) + (((j-1/2)-(M/2))^2)) # in sine space
        f = (dist/(M/2))*mF # spacial frequency at i,j in Hartley space

        if dist < M/2   # if distance in range 
            angle = atan((j-1/2) - M/2, (i-1/2) - M/2) # get angle at i,j
            if angle < 0 
                angle += 2pi
            end
            θ = mod(angle - pi/2, 2pi) # normalize to zero angle (w/o, 0º is y ax)
            
            # if dist is further than maximum derivative of variance matrix (aV)
            l = aV[:,34]; p = zeros(length(l))
            p = [l[i+1] - l[i] for i=1:length(l)-1]; gd=argmin(p)-34.5 # expand larger
            maxθ = dist < gd+(M/20) ? 22 : 12 # if in closest 10% to origin, adjust angle
            θ1 = θ - deg2rad(maxθ/2); θ2 = θ + deg2rad(maxθ/2)

            # This is given fRatio and totalθ in the arguments.
            # range of d in sine space: with fRatio proportion of bandwith/f
            tolf = (fRatio*f)   # tolerance in both directions in spacial frequency
            tols = (tolf/mF)*(M/2)  # tolerance in both directions in sin space
            d1 = 0; d2 = tols
            d1 = dist < tols/2 ? d1 : d1 + (dist-tols/2)
            d2 = dist < tols/2 ? d2 : d2 + (dist-tols/2)
            if abs(d1-d2) < 3; d2 = d1 + 3; end # set minimum distance to 3 (3*.2 = .6 spacial frequency band in sin space)

            plotBool =  false
            idx = ArcBand(M, θ1, θ2, d1, d2, plotBool)

            if length(idx) != 0
                ArcDict[(i,j)] = idx
            else
                ArcDict[(i,j)] = [(i,j)]
            end
        end
    end
    parameters["ArcDict"] = ArcDict
    return ArcDict
end

"""
    getGrid(tup, M, sz=2)

    Returns list of indices in kernel of size sz around tup in mtx size M

"""
function getGrid(tup, M; sz=2, sf=false)
    if sf==true; k = [(0,0),(1,0),(-1,0),(0,1),(0,-1)]; idxs = [tup.+ke for ke in k]
    else;        k = [(m, n) for m=-sz:sz, n=-sz:sz]; idxs = [tup.+ke for ke in k] end
    w = findall((M >= i[1] >=1 for i in idxs) .& (M >= i[2] >=1 for i in idxs))
    return [idxs[id] for id in w]
end

"""
    displayMask(idxs::Vector, M::Int64)

    Given a vector of tuples (x, y), displays locations on MxM grid. Can use this to show masks in Dicts
    Used to see indexes returned by getArcDict and getSFIdxs
"""
function displayMask(idxs::Vector, M::Int64)
    mask2 = zeros(Bool, M, M)
    comp = [(i,j) for i in 1:M, j in 1:M]
    mask2 .|= [any(comp[i, j] == tup for tup in idxs) for i in 1:size(mask2, 1), j in 1:size(mask2, 2)]
    for i=1:M, j=1:M
        if sqrt(((j-1/2) - M/2)^2 + ((i-1/2) - M/2)^2) < M/2 && mask2[i,j] < 1
            mask2[i,j] = 0.5 
        end
    end
    p = heatmap(mask2, color=:blues, xlabel="Index j", ylabel="Index i", ratio=1, yflip = false)
    display(p)
end

"""
    getVarMtx(M)

    Returns variance of angle matrix given size M with 0 on horizontal between Q1 & Q4

"""
function getVarMtx(parameters::Dict)
    M = parameters["M"]
    av = Array{Float64, 2}(undef, M, M); dv = Array{Float64, 2}(undef, M, M); hm = Int64(M/2)
    c = (hm, hm)
    d = [sqrt(((i-1/2) - c[1])^2 + ((j-1/2) - c[2])^2) for i in 1:M, j in 1:M]
    d = d./6.6
    a = [atan((j-1/2) - c[2], (i-1/2) - c[1]) for i in 1:M, j in 1:M]
    a .+= (a .< 0) .* (2pi)
    a .= mod.(a - (ones(size(a))*pi/2), 2pi)
    # making axes closer in angle
    bh_a = a[hm+1:M,M:-1:1] .- pi; a[hm+1:M, :] = bh_a
    for i=1:M, j=1:M
        ids = getGrid((i,j), M)
        a_lst = [a[id[1],id[2]] for id in ids]; d_lst = [d[id[1],id[2]] for id in ids]
        av[i,j] = var(a_lst); dv[i,j] = var(d_lst)
    end
    parameters["aV"] = av; parameters["dV"] = dv; 
    return av, dv
end

"""
    getSFIdxs(M::Int64)

    Returns a dict [(x, y)]=> vector of tuples to use in SF tuning curve, a straight line from origin
"""
function getSFIdxs(parameters::Dict)
    M = Int(parameters["M"])
    # make dictionary getting cone for each (i,j)
    SF_dict = Dict(); c = (M/2, M/2)
    for i=1:M/2, j=1:M/2 # get one quadrant
        # get angle & dist to i,j
        a = atan((j-1/2) - M/2, (i-1/2) - M/2)
        a = a < 0 ? a + 2pi : a; θ = mod(a - pi/2, 2pi) # θ is angle
        d = sqrt(((i-1/2) - c[1])^2 + ((j-1/2) - c[2])^2); rd = d/(M/2)

        # Straight Line from origin at c
        idx_lst = []
        
        v=collect(.01:.03:1); relpos = ((i-1/2) - c[1], (j-1/2) - c[2]) 
        reldist = sqrt(((i-1/2) - c[1])^2 + ((j-1/2) - c[2])^2)/(c[1])
        relloc = argmin(broadcast(abs, v.-reldist))
        # locations on line are relpos * proportion of distance from relpos to origin
        scv = [Int.(round.((i/v[relloc]).*(relpos).+c)) for i in v]
        filter!(x->(x[1]!=0&&x[2]!=0), scv)
        append!(scv, [(i,j)])
        idx_lst = deepcopy(scv)
        # filter to be less than M/2 away
        idx_lst = filter(x -> sqrt(((x[1]-1/2) - c[1])^2 + ((x[2]-1/2) - c[2])^2) < M/2, idx_lst)
        SF_dict[(Int(i),Int(j))] = idx_lst
    end

    # get other three quadrants as reflections of the first
    for key in keys(SF_dict)
        x,y = key; xy_idx = SF_dict[key]
        SF_dict[(Int(M-(x-1)), Int(y))] =         [(M-(id[1]-1), id[2]) for id in xy_idx]
        SF_dict[(Int(M-(x-1)), Int(M-(y-1)))] =   [(M-(id[1]-1), M-(id[2]-1)) for id in xy_idx]
        SF_dict[(Int(x), Int(M-(y-1)))] =         [(id[1], M-(id[2]-1)) for id in xy_idx] 
    end

    parameters["sfRef"] = SF_dict
    return SF_dict
end

"""
    blurMatrix(matrix::Matrix, M::Int64, ArcDict::Dict)

    Returns blurred matrix of size MxM with windows specified at ArcDict
"""
function blurMatrix(matrix::Matrix, M::Int64, ArcDict::Dict)
    result = zeros(size(matrix)); halfm = Int(M/2)
    # interpolate along middle axes (M/2+1 in col, M/2 in row) in 3x3 region
    for i=2:M-1
        matrix[i, halfm+1] =
            mean([matrix[i-1, halfm+2], matrix[i-1, halfm+2], 
            matrix[i, halfm+2], matrix[i, halfm+2],
            matrix[i+1, halfm+2], matrix[i+1, halfm+2]]) 
    end
    
    for j=2:M-1
        matrix[halfm, j] =
            mean([matrix[halfm+1, j-1], matrix[halfm-1, j-1],
            matrix[halfm+1, j], matrix[halfm-1, j],
            matrix[halfm+1, j+1], matrix[halfm-1, j+1]]) 
    end

    for i=1:M, j=1:M
        idx = (i,j) in keys(ArcDict) ? ArcDict[(i,j)] : []

        c = (M/2, M/2)
        d = sqrt(((i-1/2) - c[1])^2 + ((j-1/2) - c[2])^2)
        if (false)
            # Plotting the heatmap of mask
            mask2 = zeros(Bool, M, M)
            comp = [(i,j) for i in 1:M, j in 1:M]
            mask2 .|= [any(comp[i, j] == tup for tup in idx) for i in 1:size(mask2, 1), j in 1:size(mask2, 2)]
            p = heatmap(mask2, color=:blues, xlabel="Index j", ylabel="Index i", ratio=1,
            yflip = true)
            display(p)
        end

        if length(idx) != 0; result[i,j] = mean([matrix[id[1], id[2]] for id in idx])
        else; result[i,j] = matrix[i,j]; end
    end
    return result
end

"""
    getBlurredResponses()
    
    Returns dictionary of responses taking "uy" from dataset as input, of all colors
"""
function getBlurredResponses(parameters::Dict)
    blurredResp = Dict()

    responses, ArcDict, M, nPhases, nDelays, ccode = begin
        getindex.(Ref(parameters), ["uy", "ArcDict", "M", "nP", "nD", "ccode"]) end

    numcolors = length(ccode)
    
    @showprogress "Filtering..." for key in keys(responses)
    blurredResp[key] = Array{Float64, 5}(undef, M, M, nPhases, nDelays, numcolors)
    for ph=1:nPhases for t=1:nDelays for cl=1:numcolors
        m = responses[key][:,:,ph,t,cl]
        blurredResp[key][:,:,ph,t,cl] = blurMatrix(m, M, ArcDict) end end end
    end
    parameters["bR"] = blurredResp
    return parameters
end

"""
    getConstrainedMaxIdx(matrix, M::Int64)

    Returns [kx, ky, phase, τ] of maximum value with size(matrix) is MxMx...
"""
function getConstrainedMaxIdx(matrix, M::Int64)
    #print("Matrix Size: ", size(matrix), "\n")
    idxs = argmax(matrix)
    # check dist from origin is less than M/2
    dist = sqrt((((idxs[1]-1/2)-(M/2))^2) + (((idxs[2]-1/2)-(M/2))^2))
    d_lim = (M/2)
    while (dist > d_lim)
        matrix[idxs] = -Inf
        idxs = argmax(matrix)
        dist = sqrt((((idxs[1]-1/2)-(M/2))^2) + (((idxs[2]-1/2)-(M/2))^2))
    end
    return idxs
end

"""
    getMaxResponseTimes(responses::Dict, M::Int64, colorSpec)

    Returns dictionary [key]=> tuple(τ, kx, ky, phase) of max response in τ= 10:135 ms
"""
function getMaxResponseTimes(parameters::Dict)
    responses, M, ccode = getindex.(Ref(parameters), ["bR", "M", "ccode"])
    # for each recording, combines all phases together into single 3-D matrix
    dKeys = collect(keys(responses)); numcolors = length(ccode)
    timeDict = Dict()
    for key in dKeys 
        colorDict = Dict()
        for cl=1:numcolors
            #print("Color $cl, Key $key\n")
            currcolor = ccode[cl]
            # contains a vector of [time, kx, ky, phase]
            idx = getConstrainedMaxIdx(responses[key][:,:,:,11:36,cl], M)
            #print("Idxs: $idx \n")
            # adjust offset for tau since idx is taking argmax() from slice 11:36
            n_tup = (τ=idx[4]+10, kx=idx[1], ky=idx[2], ph=idx[3]) 
            colorDict[currcolor] = n_tup
            #print("Currcolor location: $n_tup\n")
            #print("Response:", responses[key][n_tup.kx, n_tup.ky, n_tup.ph, n_tup.τ, cl], "\n")
        end
        timeDict[key] = colorDict
    end
    parameters["tD"] = timeDict
    return timeDict
end

"""
    getPSIDict(maxResponses::Dict, phases::Vector)

    Returns dictionary [key][colorstring]=>[resp, resp..., PSI/circ_var of all responses]
        for phases vector with elements in ∈ [0, 1] for phase 'location' (plot hist if pl1)
    Also implements the Kruskal-Wallis test for significance between
        spacio-temporal window of max response compared to those responses
        in different phases. Extracts a p-value (plot if pl2)
"""
function getPSIDict(parameters; pl1=false, pl2=false)
    
    responses, times, phases, M, ccode = getindex.(Ref(parameters), ["bR", "tD", "phases", "M", "ccode"])

    nPhases = length(phases); dKeys = keys(responses); numcolors = length(ccode)
    responseDict = Dict(); #KWDict = Dict()
    expPhases = [2π*x for ph in phases for x in (ph, ph + 0.25)]

    # getting all responses and phase selectivity
    for key in dKeys

        responseDict[key] = Dict()
        for cl=1:numcolors

            ## PSI calculation
            τ, kx, ky = times[key][ccode[cl]]
            
            responseDict[key][ccode[cl]] = Vector()
            for ph in 1:nPhases
                # collect both phases in a matrix - for PSI
                push!(responseDict[key][ccode[cl]], responses[key][kx, ky, ph, τ, cl])
                push!(responseDict[key][ccode[cl]], responses[key][M-kx+1, M-ky+1, ph, τ, cl])
            end

            # responseDict contains (nphases*2 responses and the PSI), normalize for zero
            minval = minimum(responseDict[key][ccode[cl]])
            w1 = minval < 0 ? responseDict[key][ccode[cl]] .- minval : responseDict[key][ccode[cl]]
            push!(responseDict[key][ccode[cl]], circ_r(expPhases, w=w1))
        end

        ## Kruskal-Wallis Significance Test:
        # collect responses in ArcDict window τ ± 5ms*temps for KW significance test
        #KWGroups = Vector{Vector{Float64}}(); ts = 2
        #for ph in 1:nPhases
        #    kwg1 = Vector(); kwg2 = Vector()
        #    for tup in ArcDict[(μx, μy)]
        #        append!(kwg1, responses[key][tup[1],          tup[2],          ph, τ-ts:τ+ts])
        #        append!(kwg2, responses[key][abs(M-tup[1])+1, abs(M-tup[2])+1, ph, τ-ts:τ+ts])
        #    end
        #    push!(KWGroups, kwg1); push!(KWGroups, kwg2)
        #end
        #pval = pvalue(KruskalWallisTest(KWGroups...))
        #KWDict[key] = pval
    end

    # getting Otsu threshold for histogram
    #step=0.05; nbins = Int(floor(1/step)); edges = 0:step:1
    #psiV = [responseDict[key][end] for key in keys(responseDict)]
    #edges, count = HistogramThresholding.build_histogram(psiV, nbins)
    #t = find_threshold(count[1:end], edges, Otsu())

    # TODO: We don't use in plotHartMetrics
    # Irrelevant because we find averaged tuning curves in all responses
    # get list of keys for low responses (less than t). 
    #keylist = collect(keys(responseDict))
    #keylist = [(i, responseDict[i][5]) for i in keylist]
    #sort!(keylist, by = x -> x[2]);
    #bottomPSI = filter(x -> x[2]<t, keylist); topPSI = filter(x -> x[2] ≥ t, keylist)
    #bottomPSI = [x[1] for x in bottomPSI];    topPSI = [x[1] for x in topPSI]

    # plotting histogram of PSIs if true
    #n = length(keys(responses)); w = 2 * n^(-1/3) * iqr([responseDict[key][5] for key in keys(responseDict)])
    #h = histogram(psiV, bins = edges, title="PSI Distribution, n=$n", legend=false)
    #h2 = fit(Histogram, psiV, nbins = nbins); weights = h2.weights
    
    # TODO: include or not include gaussian fitting 
    #g = GMM(2,1,kind=:full); mtx=psiV[:,:]; g.μ .= [0.1, 0.9]
    #em!(g, mtx)
    
    #if pl1
    #    vline!([t], label = "Otsu Thresh", color = :black, linewidth = 2)
    #    display(h)
    #end

    #if pl2
    #    #scat = [(responseDict[key][end], KWDict[key]) for key in dKeys]
    #    n2 = length(scat)
    #    #scatter(scat, xlims=(-.05,1), ylims=(-.05,0.2), title="PSI vs KWP, n=$n2", xlabel="PSI", ylabel="p-value")
    #    pvalv = [KWDict[key] for key in dKeys]
    #    h = histogram(pvalv, title="p-val Distribution, n=$n2", legend=false, nbins=40)
    #    display(h)
    #end

    #parameters["cV"] = bottomPSI; parameters["sV"] = topPSI; parameters["KWPs"] = KWDict;

    parameters["PSIs"] = responseDict 
    colorPSIs = Dict()
    for cl=1:numcolors
        colorPSIs[ccode[cl]] = [responseDict[key][ccode[cl]][end] for key in keys(responseDict)] end

    parameters["PSIvals"] = colorPSIs
    
    return responseDict, colorPSIs 
end

"""
    getOS(parameters)

    Returns dict [key][color string] => Vector of [(θ, response)...] for all phases, and [(θ, response)...] avg'd over all phases
        Note: sampling from a band of thickness 1
"""
function getOS(parameters::Dict)

    blurRes, times, M, phases, ccode = getindex.(Ref(parameters), ["bR", "tD", "M", "phases", "ccode"])
    OS = Dict(); c = (M/2, M/2); nPhases = length(phases); numcolors = length(ccode)

    @showprogress "Getting Orientations..." for key in keys(blurRes) 
        OS[key] = Dict()

        for cl = 1:numcolors
            τ, x, y, p = times[key][ccode[cl]] # [τ, i, j, phase], loc of max response in sin space
            xylst = [(x, y), (M-(x-1), M-(y-1))]
            OriList = []

            for ph=1:nPhases, tup in xylst
                # get distance from origin
                xp, yp = tup
                dist = sqrt((((xp-1/2)-c[1])^2) + (((yp-1/2)-c[2])^2));
                addAngle = (M/2)-xp < 0 ? pi : 0 
                idxs = ArcBand(M, 0+addAngle, pi+addAngle, dist-0.5, dist+0.5); l = length(idxs)
        
                if l > 0
                    resMatrix = blurRes[key][:,:, ph, τ, cl]

                    # collect angles from idxs
                   a = [atan((id[2]-1/2) - M/2, (id[1]-1/2) - M/2) for id in idxs]
                    a .+= (a .< 0) .* (2pi)
                   a .= mod.(a - (ones(size(a))*pi/2), 2pi)

                    # collect responses from idxs
                    r = [resMatrix[id[1], id[2]] for id in idxs]

                    # combine thetas and responses
                    final = [(mod(a[i], pi), r[i]) for i=1:l]
                else; final = []; end
                # append OriList with specific phase responses
                push!(OriList, sort(final))
            end
        
            # append OriList with averaged responses over all thetas
            respvals = [last.(OriList[i]) for i=1:lastindex(OriList)] # get responses
            meanresp = [mean([respvals[i][idx] for i=1:lastindex(respvals)]) for idx=1:lastindex(respvals[1])] # average them
            avgs = [(OriList[1][i][1], meanresp[i]) for i=1:lastindex(meanresp)] # append to tuples with orientations
            push!(OriList, sort(avgs))
            OS[key][ccode[cl]] = OriList
        end
    end

    parameters["OSs"] = OS
    return OS
end

"""
    getSF(parameters, complex)

    if complex, use those averaged responses found in compDict
    Returns dict [key] => Vector [(sf, response)...] for each phase, and [(sf, response)...] avg'd over all phases
        Note: sampling from line of thickness 1
"""
function getSF(parameters::Dict)
    responses, times, M, refDict, mF, phases, ccode = begin
        getindex.(Ref(parameters), ["bR", "tD", "M", "sfRef", "mF", "phases", "ccode"])
    end
    SF = Dict(); nPhases = length(phases); numcolors = length(ccode)

    @showprogress "Getting SF..." for key in keys(responses)
        SF[key] = Dict()

        for cl=1:numcolors
            τ, x, y, p = times[key][ccode[cl]]
            xylst = [(x, y), (M-(x-1), M-(y-1))]
            SFList = []

            for ph=1:nPhases, tup in xylst
                idx = refDict[tup]
                l = length(idx)
                if l > 0 # collect sf's from idxs
                    resMatrix = responses[key][:,:, ph, τ, cl]
                    d = [mF*(sqrt(((id[1]-1/2) - M/2)^2 + ((id[2]-1/2) - M/2)^2)/(M/2)) for id in idx]
                    r = [resMatrix[id[1], id[2]] for id in idx]
                    final = [(d[i], r[i]) for i=1:l]
                else; final = []; end
                push!(SFList, final)
            end

            # append SFList with averaged responses over all spaical frequencies
            respvals = [last.(SFList[i]) for i=1:lastindex(SFList)] # get responses
            meanresp = [mean([respvals[i][idx] for i=1:lastindex(respvals)]) for idx=1:lastindex(respvals[1])] # average them
            avgs = [(SFList[1][i][1], meanresp[i]) for i=1:lastindex(meanresp)] # append to tuples with orientations
            push!(SFList, sort(avgs))
            SF[key][ccode[cl]] = SFList
        end
    end

    parameters["SFs"] = SF
    return SF
end

"""
    plotHartMetrics(parameters, colorSpec, path)

    parameters is dataset; colorSpec is ASCII character from ccode
    Plots the responses of highest intensity & orientation/spacial frequency
        from areas of interest in scatter plots. R is original responses
    Save figure to path
    if organize, for non-complex, plot in order by PSI
    if show, display plot

"""
function plotHartMetrics(parameters::Dict, colorSpec, path; organize=true, show=false)

    # load parameters to plot all uids of single color
    PSIDict, OSs, SFs, bR, R, timeDict, phases, ccode, pass_dict = 
    getindex.(Ref(parameters), ["PSIs", "OSs", "SFs", "bR", "uy", "tD", "phases", "ccode", "ColorResp"])

    # get ASCII version of color, num phases
    cl = ccode[colorSpec]; nphases = length(phases)
    numplots = (nphases+1)*2 # original & blurred response for each phase, then Ori & SF plot

    # invoke sort on keys by PSI
    if organize 
        PSIkeys = collect(keys(PSIDict))
        keylist = [(i, PSIDict[i][cl][end]) for i in PSIkeys]
        sort!(keylist, by = x -> x[2]);
        keylist = [x[1] for x in keylist]
    else; keylist = collect(keys(PSIDict)); end

    # get colors for each plot: 4 colors for curves, 1 black for averaged
    ex = OSs[keylist[1]][cl]
    r = collect(palette(:lightrainbow, length(ex)-1)); b = collect(palette(:greys, 1))
    colors = vcat(r, b)

    @showprogress "Plotting ..." for key in keylist

        fig_elems = string.([key, cl])
        fig_path = joinpath(path,join(fig_elems, "_"))

        psi_val = round(PSIDict[key][cl][end], sigdigits=3) # PSI value
        ph = timeDict[key][cl].ph                           # phase, τ of max resp
        τ = timeDict[key][cl].τ
        passed_tests = pass_dict[key][colorSpec]            # if color passed resp. testing
        font_color = passed_tests ? RGB(0.0, 0.5, 0.0) : RGB(1.0, 0.0, 0.0)
        plot_grid = plot(
            layout = (nphases+1,2), 
            plot_title = "Rec: $key, τ: $τ, ph: $ph, PSI: $psi_val",
            plot_titlefontsize = 15,
            plot_titlefontcolor = font_color,
            margin = 1Plots.mm, 
            size = (600, 700),
            left_margin = 0.5Plots.mm,
            right_margin = 0.5Plots.mm,
            bottom_margin = 0.5Plots.mm,
            top_margin = 3Plots.mm
        )

        # display highest intensity phase first
        o_resp_1 = R[key][:,:,ph,τ,colorSpec]; maxoR = maximum(o_resp_1) # original response
        b_resp_1 = bR[key][:,:,ph,τ,colorSpec];   maxbR = maximum(b_resp_1)  # blurred response
        heatmap!(plot_grid, o_resp_1, xlims=(0, 68), subplot=1,
        ylims=(0, 68), aspect_ratio=:equal)
        heatmap!(plot_grid, b_resp_1, subplot=2,  # filtered 
        xlims=(0, 68), ylims=(0, 68), aspect_ratio=:equal)

        ph_idxs = [i for i=1:length(phases)]; deleteat!(ph_idxs, ph)


        for ph_idx=1:length(phases)-1 # display original & filtered response for each phase
            ph = ph_idxs[ph_idx]  # extract the phase (index/integer) of other phase responses
            o_resp = R[key][:,:,ph,τ,colorSpec]
            b_resp = bR[key][:,:,ph,τ,colorSpec]
            plotnum = 2 + ph_idx
            heatmap!(plot_grid, o_resp, xlims=(0, 68), ylims=(0,68), subplot=plotnum,
            aspect_ratio=:equal, clims=(0,maxoR))
            heatmap!(plot_grid, b_resp, xlims=(0, 68), ylims=(0,68), subplot=(plotnum+1),
            aspect_ratio=:equal, clims=(0,maxbR))
        end
 
        # plot Orientation tuning curve
        ori_vect = OSs[key][cl]; 
        for i=1:lastindex(ori_vect)
            lw = i == lastindex(ori_vect) ? 3.2 : 2.5 # thicker line for avg
            phase_oris = sort(ori_vect[i])
            Ox = [rad2deg(point[1]) for point in phase_oris] # OS
            Oy = [point[2] for point in phase_oris]
            plot!(plot_grid, Ox, Oy, xlabel="θ", 
                ylabel="Response", subplot=(numplots-1), color=colors[i], 
                linewidth=lw, title="Ori. Tuning", legend=false)
        end

        # plot SF tuning curve
        sf_vect = SFs[key][cl]
        for i=1:lastindex(sf_vect)
            lw = i == lastindex(ori_vect) ? 3 : 2.5 # thicker line for avg
            phase_sfs = sort(sf_vect[i])
            SFx = [point[1] for point in phase_sfs]
            SFy = [point[2] for point in phase_sfs]
            plot!(plot_grid, SFx, SFy, xlabel="f", 
                ylabel="Response", subplot=numplots, color=colors[i], 
                linewidth=lw, title="SF Tuning", legend=false)
        end

        savefig(plot_grid, fig_path)
        if show
            plotdisplay(plot_grid)
        end
    end
    return 0
end

"""
    responsiveTesting(dataDict::Dict; pl=false, bsl=collect(1:9), slice=11:36, sdf=3.5, minmax=8)

    Runs following test on all colors in responses:
        1. max response > minmax
        2. passes isresponsive: baseline τ = -40:0 ms, responsive τ = 10:135 ms, std. dev. factor sdf
        3. welch's t-test (2-sample t-test): if distribution of responses in ROI around max response with τ ± 5ms
            is sampled from the distribution at the baseline

    Returns dict w/ two new keys, "anyColorResp" (if any color passed tests) and "ColorResp" (bool array, T/F for A, L, M, S)
"""
function responsiveTesting(dataDict::Dict; bsl=collect(1:9), slice=11:36, sdf=3.5, minmax=5)

    filt_resp, numDelays, ArcDict, ccode, M = getindex.(Ref(dataDict), ["bR", "nD", "ArcDict", "ccode", "M"])
    numcolors = length(ccode)

    color_resp_dict = Dict(); any_resp_dict = Dict()

    @showprogress "Testing Responsiveness..." for uid in keys(filt_resp)
        color_passing = []
        for cl = 1:numcolors
            
            m = filt_resp[uid]; l = argmax(m[:,:,:,slice,cl]); l += CartesianIndex(0,0,0,10); # l is loc of max response
            idxs = ArcDict[(l[1], l[2])]

            # TEST 1: remove if max reading < minmax (most pass)
            max_bool = maximum(m[l,cl]) ≥ minmax
                
            # TEST 2: check if passes isresponsive test (most pass)
            mtx = [m[id[1], id[2], l[3], :, cl] for id in idxs]
            mtx = transpose(reshape(vcat(mtx...), numDelays, length(idxs))) .+ 1e-5
            result = isresponsive(mtx; fun=mean, bi=bsl, ri=collect(slice), sdfactor=sdf); 
            isresponsive_bool=result.r

            # check if passes 2-sample t-test (most that fail fail here, TODO: make stricter)
            t_idxs = []; for i in getGrid((l[1], l[2]),M,sz=1); append!(t_idxs, ArcDict[(i[1], i[2])]); end
            t_idxs = unique(t_idxs)
            responsive_samples = Vector{Real}(); base_samples = Vector{Real}()
            for t = l[4]-1:l[4]+1 for id in t_idxs
                append!(responsive_samples, m[id[1], id[2], l[3], t, cl]) end end
            for t in bsl for id in t_idxs
                append!(base_samples, m[id[1], id[2], l[3], t, cl]) end end
            result = UnequalVarianceTTest(responsive_samples, base_samples)
            ttest_bool = pvalue(UnequalVarianceTTest(responsive_samples, base_samples)) < 0.05
            
            passing = max_bool && isresponsive_bool && ttest_bool
            append!(color_passing, passing)
        end
        color_resp_dict[uid] = color_passing
        any_resp_dict[uid] = any(color_passing)
    end

    dataDict["anyColorResp"] = any_resp_dict
    dataDict["ColorResp"]    = color_resp_dict
    return dataDict

end

function initParams(dataDict::Dict; fRatio=1/5)
    kyV = dataDict["ky"]
    dataDict["M"] = length(kyV);
    # phase is in Images Base, renamed to avoid erros
    dataDict["phases"] = deepcopy(dataDict["phase"]); delete!(dataDict, "phase") 
    dataDict["nP"] = length(dataDict["phases"]);
    dataDict["nD"] = length(dataDict["delays"])
    dataDict["mF"] = maximum(kyV); 
    dataDict["fRatio"] = fRatio
    return dataDict
end
function initDataset(dataDict::Dict)
    dataset = initParams(dataDict, fRatio=1/5); print("Initialized parameters. \n")
    getVarMtx(dataset)              # aV, dV
    getArcDict(dataset)             # ArcDict
    getSFIdxs(dataset)              # sfRef
    print("Initialized dictionaries. \n")
    getBlurredResponses(dataset)    # bR
    print("Filtered responses. \n")
    return dataset
end

siteresultdir = "/Users/andrewperun/Salk_2024/Data/AG1"
plotdir = joinpath(siteresultdir,"HartleyPlots/")

hartleyinfo = (siteresultdir;figfmt=[".png"]) -> begin

    dataset = joinhartley(siteresultdir)
    dataset = initDataset(dataset)          # initializing dictionaries, filtering responses
    dataset = responsiveTesting(dataset)    # see responsiveTesting
    tD = getMaxResponseTimes(dataset)       # get max response times

    # r_i is a response in filtered dataset
    # PSI1 collects responses from all phases and circular variance: PSI1[uid][A/L/M/S] = [r_1, r_2, r_3, r_4, circ_var]
    # PSI2 pools all PSIs from each color: PSI2[A/L/M/S] = [psi_1, psi_2, ..., psi_n]
    PSI1, PSI2 = getPSIDict(dataset)   

    # OSs collects points for Ori. Tuning curve
    # OSs[uid][A/L/M/S] = [(θ_1, r_1), (θ_2, r_2), ..., (θ_n, r_n)]
    OSs = getOS(dataset)

    # SFs collects points for SF Tuning curve
    # SFs[uid][A/L/M/S] = [(f_1, r_1), (f_2, r_2), ..., (f_n, r_n)]
    SFs = getSF(dataset)
    
    # For each color, plot statistics and write in "UID_A/L/M/S" in plotdir
    num_colors = length(dataset["ccode"])
    for colorSpec = 1:num_colors
        plotHartMetrics(dataset, colorSpec, plotdir)
    end

    jldsave(joinpath(siteresultdir,"hartleydataset.jld2");dataset)
end




## look if band-pass cells have low pass in cone isolating stimuli
##      in 4cβ, you have cells that are achromatic, no cone opponency, but high SF bandpass
##      other low pass cells would be cone opponent all the time (type 2)
##      => see if cone opponency in columns (more cone opponency, less spacial opponency; vs achromatic and band pass)
        ### good STAS (RFs) only in linear cells
        ### only in Hartley we determine those

# choose the o/sf bandwidth based on the distributions from center to end:
#   see distribution of bandwidth (y-axis) versus spacial frequency (x-axis)

# neuron 564: try different angle of filtering to eliminate false positive

## If PSI is close to 1, then we have simple cell
# if PSI is close to 0, we have more complex like, then
    # use more phases to test responses and get an average over all phases
        # to derive our tuning curves
    # average => find peak => trace from new, averaged peak
# derived from F1/F0 --- complex cell = non-simple


# rec 43: check drifting for 1 or 2 neurons? orthogonal phase responses
# also rec 271, 391
# hypothesis 2: if it IS two cells, possibly complex, then they must be extremely close, going against the dogma
#   that orientation columns have same orientation
# hypothesis 1: one cell responds to orthogonal stimuli


# Stripe Responses in AG1_V1_ODL1
# [126, 421, 463, 444, 527, 345, 427, 313, 448, 414, 426, 396, 447, 377, 411, 372, 428, 433]

# very small, 
# STAS: scoring nonlinearity, eliminate the sign. Use a spike triggered change in contrast
# deviation from pixel value from the mean
# any change in contrast causes spike