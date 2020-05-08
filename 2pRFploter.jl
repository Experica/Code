using NeuroAnalysis,Statistics,StatsBase,FileIO,Images,Plots,Interact,ImageSegmentation,LsqFit,FFTW,ProgressMeter

disk = "H:"
subject = "AE7"

# recordSession = "002"
# testId = ["008","009","010","011"]

# recordSession = "003"
# testId = ["000","001","002","004"]

# recordSession = "004"
# testId = ["001","002","003","004"]

# recordSession = "005"
# testId = ["000","001","003","004"]

# recordSession = "006"
# testId = ["000","001","002","003"]

recordSession = "001"
testId = ["004", "005", "006","003"]

recordPlane = "000"

## Prepare data & result path
siteId=[]
for i =1:size(testId,1)
    siteid = join(["$(recordSession)_", testId[i], "_$(recordPlane)"])
    push!(siteId,siteid)
end

dataFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]))
dataExportFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
resultFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", "DataExport")
resultFolderPlot = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", join(["plane_",recordPlane]),"0. Original maps","STA")
isdir(resultFolder) || mkpath(resultFolder)
isdir(resultFolderPlot) || mkpath(resultFolderPlot)

## Load all stas
# testids = ["$(siteId)_HartleySubspace_$i" for i in 1:4]
dataFile=[]
for i=1:size(testId,1)
    datafile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_sta.jld2"), dir=dataExportFolder[i],adddir=true)[1]
    push!(dataFile, datafile)
end

dataset = joinsta(load.(dataFile))
save(joinpath(resultFolder,join([subject,"_plane",recordPlane,"_sta_dataset.jld2"])),"dataset",dataset)

datasetResp = responsivesta!(dataset)
save(joinpath(resultFolder,join([subject,"_plane",recordPlane,"_sta_datasetResp.jld2"])),"dataset",datasetResp)
datasetFit = rffit!(datasetResp)
save(joinpath(resultFolder,join([subject,"_plane",recordPlane,"_sta_datasetFit.jld2"])),"dataset",datasetFit)

dataset = load(joinpath(resultFolder,join([subject,"_plane",recordPlane,"_sta_dataset.jld2"])),"dataset")

# filter.(x->occursin.("sta.jld2",x), readdir.(dataExportFolder))
#
#
# uids = size(stas[4]["usta"],1)
# uids = mapreduce(c->size(stas[c]["usta"],1),intersect,1:cn)
# uids = mapreduce(c->size(stas[c]["usta"],1),intersect,1:cn)

function exd(stas,opidx;sig=true)
    exi = [argmin(stas),argmax(stas)]
    ex = stas[exi]
    absexi = argmax(abs.(ex))
    (ex=sig ? ex[absexi] : 0, d=opidx[exi[absexi][3]])
end

function joinsta(stas)
    cn = length(stas)
    imagesize = stas[1]["imagesize"]
    delays = stas[1]["delays"]
    opidx=filter(i->!isnothing(i),indexin([0.198, 0.231, 0.264, 0.297, 0.33],delays))
    stisize = stas[1]["stisize"]
    ppd = imagesize[1]/stisize[1]
    xi = stas[1]["xi"]
    nxi = setdiff(1:prod(imagesize),xi)
    cii=CartesianIndices(imagesize)
    bdi = filter(i->!isnothing(i),indexin([-0.066, -0.033, 0.0],delays))   # use delays after 0 (no delay) as blanks

    dataset = Dict("stisize"=>stisize,"imagesize"=>imagesize,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi)
    # dataset["log"] = [i["log"] for i in stas]
    dataset["color"]=[i["color"] for i in stas]
    # dataset["minmaxcolor"] = [(i["mincolor"],i["maxcolor"]) for i in stas]
    # dataset["minmaxcg"] = map(i->minmaxcolorgradient(i...,n=256),dataset["minmaxcolor"])
    usta = Dict();ucex=Dict();uresponsive=Dict()
    # uids = mapreduce(c->size(stas[c]["usta"],1),intersect,1:cn)
    uids = size(stas[1]["usta"],1)
    # uids = mapreduce(c->keys(stas[c]["usta"]),intersect,1:cn)
    for u in 1:uids
        csta = Array{Float64}(undef,imagesize...,length(delays),cn)
        bsta = mapreduce(c->stas[c]["usta"][u,bdi,:],vcat,1:cn)  # sta of blank
        bm = mean(bsta);bsd = std(bsta)
        cex = []
        for c in 1:cn
            zsta = zscore(stas[c]["usta"][u,:,:],bm,bsd)  # z-scored sta
            for d in eachindex(delays)
                csta[cii[nxi],d,c].=mean(zsta[d,:])
                csta[cii[xi],d,c] = zsta[d,:]
            end
            push!(cex,exd(csta[:,:,opidx,c], opidx))   #PL: To check the extrema of STA, I set the time window constrain here
        end
        usta[u] = csta
        ucex[u] = cex
        uresponsive[u] = false
    end

    dataset["usta"] = usta
    dataset["ucex"] = ucex
    dataset["uresponsive"] = uresponsive
    return dataset
end
# "Local RMS Contrast of each image"
function localcontrast(csta;ws=0.5,ppd=46)
    dims = size(csta)
    clc = Array{Float64}(undef,dims)
    w = floor(Int,ws*ppd)
    w = iseven(w) ? w+1 : w
    for c in 1:dims[4], d in 1:dims[3]
        clc[:,:,d,c] = mapwindow(std,csta[:,:,d,c],(w,w))
    end
    return clc
end

function peakroi(clc)
    pi = [Tuple(argmax(clc))...]
    plc = clc[:,:,pi[3:end]...]
    segs = seeded_region_growing(plc,[(CartesianIndex(1,1),1),(CartesianIndex(pi[1:2]...),2)])
    idx = findall(labels_map(segs).==2)
    idxlim = dropdims(extrema(mapreduce(i->[Tuple(i)...],hcat,idx),dims=2),dims=2)
    hw = map(i->i[2]-i[1],idxlim)
    radius = floor(Int,maximum(hw)/2)
    center = round.(Int,mean.(idxlim))
    cpd = [Tuple(argmax(clc[:,:,:,c]))[3] for c in 1:size(clc,4)]
    return idx,center,radius,pi[4],cpd
end

function isresponsive(csta,idx,d,bdi,nonopidx;msdfactor=3.0,csdfactor=3.0)
    # nonopidx=[-0.066, -0.033, 0.0, 0.033, 0.066, 0.099, 0.132, 0.165, 0.363, 0.396]
    d in nonopidx && return false
    blc = [std(csta[idx,j]) for j in bdi]
    bcm = mean(blc);bcsd=std(blc)
    blm = [mean(csta[idx,j]) for j in bdi]
    bmm = mean(blm);bmsd=std(blm)
    (std(csta[idx,d]) > bcm+csdfactor*bcsd) || (mean(csta[idx,d]) > bmm+msdfactor*bmsd)
end

function responsivesta!(dataset;ws=0.5,msdfactor=3.5,csdfactor=3.5,roimargin=0.2,roirlim=1.6)
    imagesize = dataset["imagesize"]
    ppd=dataset["ppd"]
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    ucex = dataset["ucex"]
    delays = dataset["delays"]
    uresponsive=dataset["uresponsive"]
    nonopidx=filter(i->!isnothing(i),indexin([-0.066, -0.033, 0.0, 0.033, 0.066, 0.099, 0.132, 0.165, 0.363, 0.396],delays))
    opidx=filter(i->!isnothing(i),indexin([0.198, 0.231, 0.264, 0.297, 0.33],delays))

    ulsta=Dict();ulroi=Dict();uconeresponsive=Dict();uconeresponse=Dict();
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        # u=671
        clc = localcontrast(usta[u],ws=ws,ppd=ppd)
        idx,center,radius,c,cpd = peakroi(clc)
        # if radius > roirlim*ppd
        #     uresponsive[u] = false
        #     next!(p)
        #     continue
        # end

        # PL: Check responsiveness of each cone type and achromatic
        # PL: based on std and mean of STA image (compare with blank, which are STAs from the onset or before stim)
        rs= [isresponsive(usta[u][:,:,:,i],idx,cpd[i],bdi,nonopidx,msdfactor=msdfactor,csdfactor=csdfactor) for i in 1:size(clc,4)]
        # PL: based on threshold and delay range, to filter out low-response and cell 'response' too early or too late
        # rs2 = [(abs(ucex[u][i][:ex])>5.5) && (ucex[u][i][:d] in opidx) for i in 1:size(clc,4)]
        rs2 = [abs(ucex[u][i][:ex])>3.5 for i in 1:size(clc,4)]
        uconeresponsive[u] = rs2  # save cell's responesiveness to each cone & achromatic stimuli

        # PL: if cell response to one of cone stim and pass the threshold, it is responsive.
        if any(rs) && any(rs2)
            uresponsive[u] = true
            center = [115 115]
            # radius = 37#ceil(Int,radius*(1+roimargin))
            radius = 56
            vr = map(i->intersect(i.+(-radius:radius),1:imagesize[1]),center)
            radius = (minimum ∘ mapreduce)((r,c)->abs.([r[1],r[end]].-c),append!,vr,center)
            ulsta[u] = usta[u][map(i->i.+(-radius:radius),center)...,:,:]  # PL: I using stim size
            # uconeresponse[u] =
            ulroi[u] = (center=center,stisize=(2*radius+1)/ppd,d=cpd[c],c=c)
            ucex[u]=map(i->exd(ulsta[u][:,:,opidx,i],opidx,sig=rs[i]),1:size(ulsta[u],4))
        else
            uresponsive[u] = false
        end
        next!(p)
    end
    dataset["ulsta"]=ulsta
    dataset["ulroi"]=ulroi
    dataset["uconeresponsive"] = uconeresponsive
    return dataset
end

rfgabor(x,y,p...) = gaborf(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[3]*p[5],θ=p[6],f=p[7],phase=p[8])
# rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3]*p[5],θₑ=p[6],aᵢ=p[7],μᵢ₁=p[2]+p[8],σᵢ₁=p[3]*p[9],μᵢ₂=p[4]+p[10],σᵢ₂=p[3]*p[9]*p[5],θᵢ=p[6])
rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[3]*p[6],μᵢ₂=p[4],σᵢ₂=p[3]*p[6],θᵢ=0)


function modelfit(data,ppd;model=:gabor)
    alb,aub = abs.(extrema(data))
    ab = max(alb,aub)
    pr = (size(data)[1]-1)/2
    sr = pr/ppd

    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-(pr+1) pr+1])/ppd
    y = vec(data)
    rlt = mfun = nothing
    try
        if model == :dog
            if aub >= alb
                ae = aub + 3.5alb
                ai = 3.5alb
                ew = 0.2sr;ewl=0.15sr;ewu=0.3sr
                ir=2;irl = 1.1;iru = 3
            else
                ai = alb + 3.5aub
                ae = 3.5aub
                ew = 0.4sr;ewl=0.16sr;ewu=0.6sr
                ir=0.5;irl = 0.3;iru = 0.9
            end
            # lb=[0,          -0.4sr,    0.1sr,   -0.4sr,    0.5,    0,     0,       -0.1sr,     0.1,    -0.1sr]
            # ub=[10,         0.4sr,    0.5sr,    0.4sr,    2,      π,     Inf,      0.1sr,     10,       0.1sr]
            # p0=[0,       0,        0.3sr,    0,        1,      π/4,   aei[2],    0,         0.25,       0]
            lb=[0.5ae,   -0.35sr,    ewl,   -0.35sr,     0.5ai,    irl]
            ub=[1.5ae,    0.35sr,    ewu,    0.35sr,     1.5ai,    iru]
            p0=[ae,       0,         ew,     0,          ai,       ir]
            mfun = (x,p) -> rfdog.(x[:,1],x[:,2],p...)
        elseif model == :gabor
            ori,sf = freqimagestats(powerspectrum(data,ppd)...)
            fub = min(1.5sf,8);flb=max(0.5sf,0.2)
            lb=[0.7ab,  -0.36sr,   0.1sr,   -0.36sr,   0.3,    0,    flb,   0]
            ub=[1.3ab,   0.36sr,   0.36sr,   0.36sr,   2.3,    π,    fub,   1]
            p0=[ab,      0,        0.2sr,    0,        1,      ori,  sf,    0]
            mfun = (x,p) -> rfgabor.(x[:,1],x[:,2],p...)
        end
        if !isnothing(mfun)
            mfit = curve_fit(mfun,x,y,p0,lower=lb,upper=ub,
            maxIter=3000,x_tol=1e-11,g_tol=1e-15,min_step_quality=1e-4,good_step_quality=0.25,lambda_increase=5,lambda_decrease=0.2)
            rlt = (param=mfit.param,converged=mfit.converged,resid=mfit.resid,r=cor(y,mfun(x,mfit.param)))
        end
    catch
    end
    return rlt
end
function rffit!(dataset;model=[:gabor])
    ulsta = dataset["ulsta"]
    ulroi = dataset["ulroi"]
    ppd = dataset["ppd"]

    if !haskey(dataset,"urf")
        dataset["urf"]=Dict()
    end
    urf = dataset["urf"]
    p = ProgressMeter.Progress(length(ulsta),desc="Fit RFs ... ")
    for u in keys(ulsta)
        if !haskey(urf,u)
            urf[u] = Dict()
        end
        # u!=96 && continue
        exs=map(i->i.ex,dataset["ucex"][u])
        ds=map(i->i.d,dataset["ucex"][u])
        for m in model
            urf[u][m] = map(i->exs[i]==0 ? nothing : modelfit(ulsta[u][:,:,ds[i],i],ppd,model=m),1:length(exs))
        end
        next!(p)
    end

    return dataset
end

t=dataset["ulsta"][92][:,:,12,2]

@manipulate for i in 1:3
    heatmap(rand(10,10),xlabel="t")
end
## stas
# @manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
@manipulate for u in sort(collect(keys(dataset["ulsta"])))#,d in eachindex(dataset["delays"])
    usta=dataset["usta"][u]
    # delays = dataset["delays"]
    d=8
    imagesize = dataset["imagesize"]
    stisize = dataset["stisize"]
    uresponsive = dataset["uresponsive"]
    clim = maximum(map(i->abs(i.ex),dataset["ucex"][u]))
    xlim = [0,5]#[0,stisize[2]]
    ylim = [0,5]#stisize[1]]
    x = range(xlim...,length=imagesize[2])
    y = range(ylim...,length=imagesize[1])

    p = Plots.plot(layout=(1,4),legend=false,size=(1600,600),title="Unit_$(u)_STA_$(delays[d])_Responsive_$(uresponsive[u])")
    foreach(c->Plots.heatmap!(p,subplot=c,x,y,usta[:,:,d,c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlim,ylims=ylim,xticks=xlim,yticks=ylim,yflip=true,xlabel=dataset["color"][c]),1:4)
    p
    savefig("./temp/Unit_$u.png")
end
## best roi stas
 # for u in sort(collect(keys(dataset["ulsta"])))
for u in sort(collect(keys(datasetResp["ulsta"])))
    # u=11
    ulsta = datasetResp["ulsta"][u]
    delays = datasetResp["delays"]
    imagesize = size(ulsta)[1]
    stisize = datasetResp["ulroi"][u].stisize
    # stisize = datasetResp["stisize"]
    ucex = map(i->i.ex,datasetResp["ucex"][u])
    ucd = map(i->i.d,datasetResp["ucex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(stisize,digits=1)]
    xy = range(xylim...,length=imagesize)

    p = Plots.plot(layout=(1,5),legend=false,size=(1650,600))
    Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
    foreach(c->Plots.heatmap!(p,subplot=c+1,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:grid,
    color=:coolwarm,clims=(-clim,clim),xlims=xylim,ylims=xylim,xticks=[],yticks=[],
    yflip=true,xlabel=string(datasetResp["color"][c]),title="Cell_$(u)_STA_$(delays[ucd[c]])"),1:4)
    foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]), label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"])),1)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]), label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"])),1)

    p

    savefig(joinpath(resultFolderPlot,join([subject,"_U",recordSession,"_Plane",recordPlane, "_Cell",u,".png"])))
end

sort(collect(keys(dataset["urf"])))
collect(keys(first(values(dataset["urf"]))))
## rf fit
@manipulate for u in sort(collect(keys(dataset["urf"]))),m in collect(keys(first(values(dataset["urf"]))))
    nc = length(dataset["color"])
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    imagesize = size(ulsta)
    ppd = dataset["ppd"]
    stisize = dataset["ulroi"][u].stisize
    ucex = map(i->i.ex,dataset["ucex"][u])
    ucd = map(i->i.d,dataset["ucex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(stisize,digits=1)]
    xy = range(xylim...,length=imagesize[2])

    p = Plots.plot(layout=(4,nc),legend=false,size=(1500,750))
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=xylim,yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ucd[c]])"),1:nc)

    fit=dataset["urf"][u]
    pr = (imagesize[1]-1)/2
    x=y=(-pr:pr)./ppd
    sr = round(x[end],digits=1)
    xylim=[-sr,sr]
    if m == :dog
        rfs = map(f->isnothing(f) ? nothing : [rfdog(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    else
        rfs = map(f->isnothing(f) ? nothing : [rfgabor(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    end
    foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+nc,frame=:none) :
    Plots.heatmap!(p,subplot=c+nc,x,y,rfs[c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=xylim,yticks=xylim,yformatter=i->-i,yflip=true,xlabel=string(dataset["color"][c]),
    title="Fit_$(m)"),1:nc)

    # foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+2nc,frame=:none) :
    # Plots.histogram!(p,subplot=c+2nc,fit[m][c].resid,frame=:semi,color=:coolwarm,linecolor=:match,bar_width=1,xlims=[-abs(maximum(fit[m][c].resid)),abs(maximum(fit[m][c].resid))],
    # xlabel="Residual",title="",grid=false),1:nc)
    #
    # foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+3nc,frame=:none) :
    # Plots.scatter!(p,subplot=c+3nc,vec(ulsta[:,:,ucd[c],c]),vec(rfs[c]),frame=:semi,color=:coolwarm,grid=false,
    # xlabel="y",ylabel="predict y",title="r = $(round(fit[m][c].r,digits=3))",markerstrokewidth=0,markersize=1),1:nc)
    p
end





function getbestrf(dataset;rt=0.65)
    urf = dataset["urf"]
    nc = length(dataset["color"])
    ubrf=Dict()
    for u in keys(urf)
        # u!=96 && continue
        mfs=[]
        for m in keys(urf[u])
            push!(mfs, map(f->begin
             if isnothing(f)
                 nothing
             else
                 (param=f.param,r=f.r,m=m)
             end
         end,urf[u][m]))
         end
         # display(mfs)
         crf = []
         for i in 1:nc
             push!(crf,mapreduce(mf->mf[i],(m1,m2)->begin
             if isnothing(m1)
                 if isnothing(m2)
                     nothing
                 else
                     m2.r > rt ? m2 : nothing
                 end
             else
                 if isnothing(m2)
                     m1.r > rt ? m1 : nothing
                 else
                     t = m1.r > m2.r ? m1 : m2
                     t.r > rt ? t : nothing
                 end
             end
         end,mfs))
         end
         ubrf[u]=crf
    end
    return ubrf
end


usrf = getbestrf(dataset)


urf = map(i->i[:gabor][2],values(dataset["urf"]))
urft = map(i->i.converged,urf)
vurf = urf[urft]

vups = mapreduce(i->i.param,hcat,vurf)

using Clustering

r = kmeans(vups[[3,5],:],4,display=:final)

r = kmeans([abs.(vups[3:3,:]./vups[5:5,:].-1);maximum(vups[[3,5],:],dims=1).*vups[7:7,:]],3,display=:final)

r = kmeans(maximum(vups[[3,5],:],dims=1).*vups[7:7,:],2,display=:final)


ci = assignments(r)


cu = map(i->us[findall(i.==ci)],1:3)


sort(cu[1])

Plots.scatter(vups[3,:],vups[7,:],group=ci)



# spike = load(joinpath(resultFolder,testids[1],"spike.jld2"),"spike")
# unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

us = collect(keys(dataset["ulsta"]))
# up = unitposition[indexin(us,unitid),:]
# ul = unitlayer[indexin(us,unitid)]



uc = [argmax(map(j->abs(j.ex),dataset["ucex"][u])) for u in keys(dataset["ulsta"])]
ud = map((u,c)->dataset["ucex"][u][c].d,keys(dataset["ulsta"]),uc)

ucsta = map((u,d,c)->mapcolor(dataset["ulsta"][u][:,:,d,c],dataset["minmaxcg"][c]),keys(dataset["ulsta"]),ud,uc)
ucmsta = map(i->alphamask(i,radius=0.5,sigma=10,masktype="Disk")[1],ucsta)


p=plotunitpositionimage(up,ucmsta,layer=layer)
save("UnitPosition_STA.svg",p)

p=plotunitlayerimage(ul,ucmsta,unitid=us)
save("UnitLayer_STA.svg",p)


uexs = mapreduce(u->map(j->j.ex,dataset["ucex"][u]),hcat,us)
uds = mapreduce(u->map(j->j.d,dataset["ucex"][u]),hcat,us)

function rftype(cws,cpd)
    a=cws[1];l=cws[2];m=cws[3];s=cws[4];c=1;d = cpd[argmax(abs.(exs))]
    t = join(filter(k->!isnothing(k),map((i,j)->i==0 ? nothing : "$j$(i > 0 ? "+" : "-")",cws,["A","L","M","S"])),"_")
    if l != 0
        if m != 0 # both significent l,m
            if l > 0 && m < 0
                t = "L+M-"
                c = 2
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l < 0 && m > 0
                t = "M+L-"
                c = 3
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l > 0 && m > 0
                if s < 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S-L+M+"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            elseif l < 0 && m < 0
                if s > 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S+L-M-"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            end
        else
            if a ==0
                t = "L+"
                c = 4
                d = cpd[4]
            else
            end
        end
    end

        if l>0 && m<0
            t = "L+M-"
            c = 2
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l<0 && m>0
            t="M+L-"
            c=3
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l>0 && m>0
            if s<0
                t="S-L+M+"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "+-"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        elseif l<0 && m<0
            if s>0
                t="S+L-M-"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "-+"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        end
    return (t=t,c=c,d=d)
end

tcd = [rftype(uexs[:,i],uds[:,i]) for i in 1:size(uexs,2)]
uc = map(i->i.c,tcd)
ud = map(i->i.d,tcd)

_,nsg,_,_=histrv(unitposition[unitgood,2],0,3500,binwidth=60)
_,ns,_,_=histrv(up[:,2],0,3500,binwidth=60)

bar(0:60:3500,[nsg ns],orientation=:h,bar_position=:stack,xlims=[-0.3,10])

plot([nsg ns],range(0,step=40,length=length(ns)),fillalpha=0.4,fill=true,labels=["SU" "SU with STA"],linealpha=0,xlims=[-0.3,7])
hline!([layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(-0.3,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

savefig("unit_layer_dist.svg")



import Base:sin
sin(x,y;fx=1,fy=1) = sin(2π*(fx*x+fy*y))




t = dataset["ulsta"][92][:,:,10,2]

Plots.heatmap(t)

ft=abs.(fft(t))
Plots.heatmap(ft)

freqs = fftfreq(25,30)





lmpi = [Tuple(argmax(clc))...]

lmi = clc[:,:,lmpi[3:end]...]

heatmap(lmi,aspect_ratio=:equal,frame=:all,color=:coolwarm,yflip=true)
scatter!([lmpi[2]],[lmpi[1]],markersize=10)
seg = seeded_region_growing(lmi,[(CartesianIndex(lmpi[1:2]...),2),(CartesianIndex(1,1),1)])

lr = labels_map(seg)
heatmap(lr,aspect_ratio=:equal,frame=:all,color=:grays)

ri = mapreduce(i->[i...],hcat, Tuple.(findall(lr.>1)))


t = extrema(ri,dims=2)
tr = map(i->i[2]-i[1],t)
rsize = floor(Int,maximum(tr)/2)

lp = floor.(Int,mean.(t))[:]

hspan!([lp[1]-rsize,lp[1]+rsize],alpha=0.2)
vspan!([lp[2]-rsize,lp[2]+rsize],alpha=0.2)


















@manipulate for u in uids
t1 = stas[1]["usta"][u]
t2=stas[2]["usta"][u]
t3=stas[3]["usta"][u]
t4=stas[4]["usta"][u]

sds = Array{Float64}(undef,imagesize...,41)
for d in 1:41
t = t1[d,:]
tt = fill(mean(t),imagesize)
tt[xi]=t
sds[:,:,d]= mapwindow(std,tt,(5,5))
end

mi = [Tuple(argmax(sds))...][1:2]
display(mi)
mir= map(i->filter(j->j>0,(-5:5) .+i),mi)
display(mir)
sds1 = [sum(sds[mir...,d]) for d in 1:41]
# bm1 = mean(t1[1,:])
# display(bm1)
# bsd1 = std(t1[1,:])
# display(bsd1)
# bm2 = mean(t2[1,:])
# bsd2 = std(t2[1,:])
# bm3 = mean(t3[1,:])
# bsd3 = std(t3[1,:])
# bm4 = mean(t4[1,:])
# bsd4 = std(t4[1,:])
#
#
# sds1=sum((t1.-bm1)./bsd1.^2,dims=2)[:]
# sds2=sum((t2.-bm2)./bsd2.^2,dims=2)[:]
# sds3=sum((t3.-bm3)./bsd3.^2,dims=2)[:]
# sds4=sum((t4.-bm4)./bsd4.^2,dims=2)[:]

plot(sds1)


# c = [:black :red :lightgreen :blue]
# plot([sds1 sds2 sds3 sds4],color=c,lw=2,grid=false)
# a=[sds1[[1,2,end,end-1]],sds2[[1,2,end,end-1]],sds3[[1,2,end,end-1]],sds4[[1,2,end,end-1]]]'
# hline!(mean.(a).+5*std.(a),color=c)
end













skm = pyimport("skimage.metrics")

@manipulate for u in uids
t1 = stas[1]["usta"][u][1,:]
t2 = stas[1]["usta"][u][25,:]
ti1=fill(mean(t1),imagesize)
ti1[xi]=t1
ti2=fill(mean(t2),imagesize)
ti2[xi]=t2
ssim,di=skm.structural_similarity(ti1,ti2,full=true)

heatmap(di,clims=(0.4,1))
end
