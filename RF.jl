using NeuroAnalysis,Statistics,StatsBase,FileIO,Images,StatsPlots,Interact,ImageSegmentation,LsqFit,ProgressMeter,
Distances,Random,Optim,RCall
import NeuroAnalysis: isresponsive

dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL5"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)
layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")


## RF Functions
"extrema value of max amplitude and it's delay index"
function exd(sta)
    exi = [argmin(sta),argmax(sta)]
    ex = sta[exi]
    absexi = argmax(abs.(ex))
    (ex=ex[absexi],d=exi[absexi][3])
end
"join color channel stas"
function joinsta(stas;siteid="",ccode = Dict("DKL_X"=>'A',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S'))
    sizepx = stas[1]["sizepx"]
    sizedeg = stas[1]["sizedeg"]
    delays = stas[1]["delays"]
    ppd = first(sizepx./sizedeg)
    xi = stas[1]["xi"]
    notxi = setdiff(1:prod(sizepx),xi)
    cii=CartesianIndices(sizepx)
    bdi = filter!(i->!isnothing(i),indexin([-100:0;200:300],delays))

    dataset = Dict("sizedeg"=>sizedeg,"sizepx"=>sizepx,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi,"siteid"=>siteid)
    colors = [i["color"] for i in stas]
    colorcodes = map(i->ccode[i],colors)
    corder = sortperm(colorcodes)

    dataset["log"] = [stas[i]["log"] for i in corder]
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    dataset["minmaxcolor"] = [(stas[i]["mincolor"],stas[i]["maxcolor"]) for i in corder]
    usta = Dict();ucexd=Dict();uresponsive=Dict()

    uids = mapreduce(c->keys(stas[c]["usta"]),intersect,corder)
    for u in uids
        csta = Array{Float64}(undef,sizepx...,length(delays),length(corder))
        cexd = []
        for c in corder
            zsta = stas[c]["usta"][u]
            for d in eachindex(delays)
                csta[cii[notxi],d,c].=mean(zsta[d,:])
                csta[cii[xi],d,c] = zsta[d,:]
            end
            push!(cexd,exd(csta[:,:,:,c]))
        end
        usta[u] = csta
        ucexd[u] = cexd
        uresponsive[u] = false
    end
    dataset["usta"] = usta
    dataset["ucexd"] = ucexd
    dataset["uresponsive"] = uresponsive
    return dataset
end
"local RMS contrast of each image, highlighting local structure regions"
function localcontrast(csta;ws=0.5,ppd=45)
    dims = size(csta)
    clc = Array{Float64}(undef,dims)
    w = round(Int,ws*ppd)
    w = iseven(w) ? w+1 : w
    for d in 1:dims[3], c in 1:dims[4]
        clc[:,:,d,c] = mapwindow(std,csta[:,:,d,c],(w,w))
    end
    return clc
end
function localcontrast(data::Matrix;ws=0.5,ppd=45)
    w = round(Int,ws*ppd)
    w = iseven(w) ? w+1 : w
    mapwindow(std,data,(w,w))
end
"peak ROI region and its delay"
function peakroi(clc)
    ds = size(clc)[1:2]
    pi = [Tuple(argmax(clc))...]
    plc = clc[:,:,pi[3:end]...]
    segs = seeded_region_growing(plc,[(CartesianIndex(1,1),1),(CartesianIndex(1,ds[2]),1),(CartesianIndex(ds[1],1),1),
                (CartesianIndex(ds...),1),(CartesianIndex(pi[1:2]...),2)])
    idx = findall(labels_map(segs).==2)
    idxlims = dropdims(extrema(mapreduce(i->[Tuple(i)...],hcat,idx),dims=2),dims=2)
    hw = map(i->i[2]-i[1],idxlims)
    center = round.(Int,mean.(idxlims))
    radius = round(Int,maximum(hw)/2)
    return (i=idx,center=center,radius=radius,pd=pi[3])
end
function peakroi(data::Matrix)
    ds = size(data)
    pi = [Tuple(argmax(data))...]
    segs = seeded_region_growing(data,[(CartesianIndex(1,1),1),(CartesianIndex(1,ds[2]),1),(CartesianIndex(ds[1],1),1),
                            (CartesianIndex(ds...),1),(CartesianIndex(pi...),2)])
    idx = findall(labels_map(segs).==2)
    idxlims = dropdims(extrema(mapreduce(i->[Tuple(i)...],hcat,idx),dims=2),dims=2)
    hw = map(i->i[2]-i[1],idxlims)
    center = round.(Int,mean.(idxlims))
    radius = round(Int,maximum(hw)/2)
    return (i=idx,center=center,radius=radius)
end
"check local mean and contrast in a region significently different from baseline"
function isresponsive(sta,idx,bdi;mfactor=3,cfactor=3)
    lc = [std(sta[idx,j]) for j in 1:size(sta,3)]
    lm = [mean(sta[idx,j]) for j in 1:size(sta,3)]
    lcmaxd = argmax(lc); lcmax = lc[lcmaxd]
    lmmaxd = argmax(lm); lmmax = lm[lmmaxd]
    lmmind = argmin(lm); lmmin = lm[lmmind]
    bmm = mean(lm[bdi]);bmsd=std(lm[bdi])
    bcm = mean(lc[bdi]);bcsd=std(lc[bdi])

    (!(lcmaxd in bdi) && lcmax > bcm+cfactor*bcsd) ||
    (!(lmmaxd in bdi) && lmmax > bmm+mfactor*bmsd) ||
    (!(lmmind in bdi) && lmmin < bmm-mfactor*bmsd)
end
"check unit sta responsive and cut local sta"
function responsivesta!(dataset;ws=0.5,mfactor=3.5,cfactor=3.5,roimargin=0.2,peakroirlim=1.5)
    sizepx = dataset["sizepx"]
    ppd = dataset["ppd"]
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    ucexd = dataset["ucexd"]
    uresponsive=dataset["uresponsive"]

    ulsta=Dict();ulcexd=Dict();ulcroi=Dict();ulroi=Dict();ucresponsive=Dict()
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        clc = localcontrast(usta[u],ws=ws,ppd=ppd)
        cproi = map(c->peakroi(clc[:,:,:,c]),1:size(clc,4))
        ucresponsive[u] = map(c->cproi[c].radius < peakroirlim*ppd ? isresponsive(usta[u][:,:,:,c],cproi[c].i,bdi,mfactor=mfactor,cfactor=cfactor) : false,1:size(clc,4))
        if any(ucresponsive[u])
            uresponsive[u] = true
            vroi = cproi[ucresponsive[u]]
            vcs = mapfoldl(r->r.center,hcat,vroi,init=zeros(Int,2,0))
            center = round.(Int,mean(vcs,dims=2)[:])
            cdev = maximum(Distances.colwise(Euclidean(),vcs,center))
            radius = round(Int,(maximum(map(r->r.radius,vroi))+cdev)*(1+roimargin))
            vr = map(i->intersect(center[i].+(-radius:radius),1:sizepx[i]),1:2)
            radius = (minimum ∘ mapreduce)((r,c)->abs.([r[begin],r[end]].-c),append!,vr,center)
            ulroi[u] = (center=center,diameterdeg=(2radius+1)/ppd)
            ulsta[u] = usta[u][map(i->i.+(-radius:radius),center)...,:,:]
            ulcexd[u] = map(i->exd(ulsta[u][:,:,:,i]),1:size(ulsta[u],4))
            ulcroi[u] = cproi
        else
            uresponsive[u] = false
        end
        next!(p)
    end
    dataset["ulsta"]=ulsta
    dataset["ulcroi"]=ulcroi
    dataset["ulroi"]=ulroi
    dataset["ulcexd"]=ulcexd
    dataset["ucresponsive"]=ucresponsive

    return dataset
end
"gabor RF model"
rfgabor(x,y,p...) = gaborf(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5],θ=p[6],f=p[7],phase=p[8])
# rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3]*p[5],θₑ=p[6],aᵢ=p[7],μᵢ₁=p[2]+p[8],σᵢ₁=p[3]*p[9],μᵢ₂=p[4]+p[10],σᵢ₂=p[3]*p[9]*p[5],θᵢ=p[6])
"dog RF model"
rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[3]*p[6],μᵢ₂=p[4],σᵢ₂=p[3]*p[6],θᵢ=0)
"fit a 2D model to an image"
function modelfit(data,ppd;model=:gabor)
    drpx = (size(data)[1]-1)/2
    radius = drpx/ppd
    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-(drpx+1) (drpx+1)])/ppd
    y = vec(data)

    # try estimate solution
    roi = peakroi(localcontrast(data,ws=0.5,ppd=ppd))
    alb,aub = abs.(extrema(data[roi.i]))
    ab = max(alb,aub)
    r = roi.radius/ppd
    c = [roi.center[2] - (drpx+1), -roi.center[1] + (drpx+1)]/ppd

    rlt = mfun = missing
    try
        if model == :dog
            if aub >= alb
                ai = 3.5alb
                ae = aub + ai
                es = 0.2r;esl=0.15r;esu=0.3r
                ier=2;ierl = 1.1;ieru = 3
            else
                ae = 3.5aub
                ai = alb + ae
                es = 0.4r;esl=0.16r;esu=0.6r
                ier=0.5;ierl = 0.3;ieru = 0.9
            end
            # lb=[0,          -0.4sr,    0.1sr,   -0.4sr,    0.5,    0,     0,       -0.1sr,     0.1,    -0.1sr]
            # ub=[10,         0.4sr,    0.5sr,    0.4sr,    2,      π,     Inf,      0.1sr,     10,       0.1sr]
            # p0=[0,       0,        0.3sr,    0,        1,      π/4,   aei[2],    0,         0.25,       0]
            ub=[1.5ae,    0.36r,    esu,    0.36r,     1.5ai,    ieru]
            lb=[0.5ae,   -0.36r,    esl,   -0.36r,     0.5ai,    ierl]
            p0=[ae,       0,        es,     0,          ai,      ier]
            mfun = (x,p) -> rfdog.(x[:,1],x[:,2],p...)
        elseif model == :gabor
            ori,sf = f1orisf(powerspectrum2(data,ppd)...)
            ub=[1.5ab,   0.5r+c[1],   0.6r,    0.5r+c[2],    0.6r,      π,     10,     1]
            lb=[0.5ab,  -0.5r+c[1],   0.1r,   -0.5r+c[2],    0.1r,      0,     0.1,    0]
            p0=[ab,      c[1],        0.3r,     c[2],        0.3r,      ori,   sf,     0.5]
            mfun = (x,p) -> rfgabor.(x[:,1],x[:,2],p...)
            ofun = (p;x=x,y=y) -> sum((y.-mfun(x,p)).^2)
        end
        if !ismissing(mfun)
            ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.9),Optim.Options(iterations=200000))
            param=ofit.minimizer; converged = Optim.converged(ofit); yy = mfun(x,param); resid = y .- yy
            rt = rcopy(R"cor.test($y,$yy,alternative='greater',method='pearson')")
            r = rt[:estimate]; pr = rt[:p_value]

            # mfit = curve_fit(mfun,x,y,param,lower=lb,upper=ub,maxIter=5000)
            # param=mfit.param; r = cor(y,mfun(x,param)); converged = mfit.converged; resid = mfit.resid

            rlt = (model=model,radius=radius,param=param,converged=converged,resid=resid,r=r,pr=pr)
        end
    catch exc
        display.(stacktrace(catch_backtrace()))
    end
    return rlt
end
"fit responsive sta to each type of models"
function rffit!(dataset;model=[:gabor,:dog])
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
        rs = dataset["ucresponsive"][u]
        ds = map(i->i.d,dataset["ulcexd"][u])
        for m in model
            urf[u][m] = [rs[i] ? modelfit(ulsta[u][:,:,ds[i],i],ppd,model=m) : missing for i in 1:length(rs)]
        end
        next!(p)
    end
    return dataset
end
function rfimage(rffit,x,y;yflip=false)
    if rffit.model == :dog
        rfi = [rfdog(i,j,rffit.param...) for j in y, i in x]
    else
        rfi = [rfgabor(i,j,rffit.param...) for j in y, i in x]
    end
    yflip && (rfi=reverse(rfi,dims=1))
    rfi
end
function rfimage(rffit,radiusdeg;ppd=45)
    x=y= -radiusdeg:1/ppd:radiusdeg
    rfimage(rffit,x,y,yflip=true)
end
rfimage(rffit;ppd=rffit.ppd) = rfimage(rffit,rffit.radius,ppd=ppd)
function rfcolorimage(img;color=:coolwarm)
    cg = cgrad(color)
    m = maximum(abs.(img))
    alphamask(get(cg,img,(-m,m)),radius=0.5,sigma=0.35,masktype="Gaussian")[1]
end

## process all stas
testids = ["$(siteid)_HartleySubspace_$i" for i in 1:4]
testn=length(testids)
stadir = joinpath(resultsitedir,"sta")
isdir(stadir) || mkpath(stadir)
lstadir = joinpath(resultsitedir,"lsta")
isdir(lstadir) || mkpath(lstadir)
rffitdir = joinpath(resultsitedir,"rffit")
isdir(rffitdir) || mkpath(rffitdir)

dataset = joinsta(load.(joinpath.(resultsitedir,testids,"sta.jld2")),siteid=siteid)
dataset = responsivesta!(dataset)
dataset = rffit!(dataset,model=[:gabor])

save(joinpath(resultsitedir,"stadataset.jld2"),"dataset",dataset)
dataset = load(joinpath(resultsitedir,"stadataset.jld2"),"dataset")

## stas
plotstas=(u,d;dir=nothing)->begin
    usta = dataset["usta"][u]
    delays = dataset["delays"]
    sizepx = dataset["sizepx"]
    sizedeg = dataset["sizedeg"]
    uresponsive = dataset["uresponsive"][u]
    clim = maximum(map(i->abs(i.ex),dataset["ucexd"][u]))
    xlims = [0,sizedeg[2]]
    ylims = [0,sizedeg[1]]
    x = range(xlims...,length=sizepx[2])
    y = range(ylims...,length=sizepx[1])

    title = "Unit_$(u)_STA_$(delays[d])"
    p = Plots.plot(layout=(1,testn),legend=false,size=(400testn,600),titlefontcolor=uresponsive ? :green : :match,title=title)
    foreach(c->Plots.heatmap!(p,subplot=c,x,y,usta[:,:,d,c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlims,ylims=ylims,xticks=xlims,yticks=ylims,yflip=true,xlabel=dataset["color"][c]),1:testn)
    isnothing(dir) ? p : savefig(joinpath(dir,"$title.png"))
end

@manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    plotstas(u,d)
end
for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    u != 112 && continue
    plotstas(u,d,dir=stadir)
end
## responsive local stas
plotlstas=(u;dir=nothing)->begin
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    diameterpx = size(ulsta)[1]
    diameterdeg = dataset["ulroi"][u].diameterdeg
    ucresponsive = dataset["ucresponsive"][u]
    ulcex = map(i->i.ex,dataset["ulcexd"][u])
    ulcd = map(i->i.d,dataset["ulcexd"][u])
    clim = maximum(abs.(ulcex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = Plots.plot(layout=(1,testn+1),legend=false,size=(350(testn+1),600))
    Plots.bar!(p,subplot=1,dataset["color"],ulcex,frame=:zerolines,ylabel="extrema")
    foreach(c->Plots.heatmap!(p,subplot=c+1,xy,xy,ulsta[:,:,ulcd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[c] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,xlabel=dataset["color"][c],title="Unit_$(u)_STA_$(delays[ulcd[c]])"),1:testn)
    isnothing(dir) ? p : savefig(joinpath(dir,"Unit_$u.png"))
end

@manipulate for u in sort(collect(keys(dataset["ulsta"])))
    plotlstas(u)
end
for u in sort(collect(keys(dataset["ulsta"])))
    plotlstas(u,dir=lstadir)
end
## rf fit
plotrffit=(u,m;dir=nothing)->begin
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    ppd = dataset["ppd"]
    diameterpx = size(ulsta)[1]
    diameterdeg = dataset["ulroi"][u].diameterdeg
    ucresponsive = dataset["ucresponsive"][u]
    ulcex = map(i->i.ex,dataset["ulcexd"][u])
    ulcd = map(i->i.d,dataset["ulcexd"][u])
    clim = maximum(abs.(ulcex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = Plots.plot(layout=(4,testn),legend=false,size=(400testn,4*250))
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,ulsta[:,:,ulcd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[c] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ulcd[c]])"),1:testn)

    umfit=dataset["urf"][u][m]
    rpx = (diameterpx-1)/2
    x=y=(-rpx:rpx)./ppd
    xylims=[round.(extrema(x),digits=2)...]
    rfs = map(f->ismissing(f) ? missing : rfimage(f,x,y),umfit)
    foreach(c->ismissing(rfs[c]) ? Plots.plot!(p,subplot=c+testn,frame=:none) :
    Plots.heatmap!(p,subplot=c+testn,x,y,rfs[c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylims,ylims=xylims,xticks=xylims,yticks=xylims,xlabel=dataset["color"][c],titlefontcolor=ucresponsive[c] ? :green : :match,title="Fit_$(m)"),1:testn)

    foreach(c->ismissing(rfs[c]) ? Plots.plot!(p,subplot=c+2testn,frame=:none) :
    Plots.histogram!(p,subplot=c+2testn,umfit[c].resid,frame=:semi,xlabel="Residual",grid=false),1:testn)

    foreach(c->ismissing(rfs[c]) ? Plots.plot!(p,subplot=c+3testn,frame=:none) :
    Plots.scatter!(p,subplot=c+3testn,vec(ulsta[:,:,ulcd[c],c]),vec(rfs[c]),frame=:semi,grid=false,
    xlabel="y",ylabel="predict y",titlefontcolor=umfit[c].pr < 0.05 ? :green : :match,title="r = $(round(umfit[c].r,digits=3))",markerstrokewidth=0,markersize=1),1:testn)
    isnothing(dir) ? p : savefig(joinpath(dir,"Unit_$(u)_Fit_$(m).png"))
end

@manipulate for u in sort(collect(keys(dataset["urf"]))),m in collect(keys(first(values(dataset["urf"]))))
    plotrffit(u,m)
end
for u in sort(collect(keys(dataset["urf"]))),m in collect(keys(first(values(dataset["urf"]))))
    plotrffit(u,m,dir=rffitdir)
end
## rf space of unit


crti = map(t->findfirst(t.==crtypes),values(ucrt))

spike = load(joinpath(resultsitedir,testids[1],"spike.jld2"),"spike")
unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

ui = indexin(keys(ucrt),unitid)
up=unitposition[ui,:]


crtcs[2.0,3.0]

plotunitposition(up,color=crtcs[crti])


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

spike = load(joinpath(resultsitedir,testids[1],"spike.jld2"),"spike")
unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

us = collect(keys(dataset["ulsta"]))
up = unitposition[indexin(us,unitid),:]
ul = unitlayer[indexin(us,unitid)]



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


tcd = [rftype(uexs[:,i],uds[:,i]) for i in 1:size(uexs,2)]
uc = map(i->i.c,tcd)
ud = map(i->i.d,tcd)

_,nsg,_,_=histrv(unitposition[unitgood,2],0,3500,binwidth=60)
_,ns,_,_=histrv(up[:,2],0,3500,binwidth=60)

bar(0:60:3500,[nsg ns],orientation=:h,bar_position=:stack,xlims=[-0.3,10])

plot([nsg ns],range(0,step=40,length=length(ns)),fillalpha=0.4,fill=true,labels=["SU" "SU with STA"],linealpha=0,xlims=[-0.3,7])
hline!([layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(-0.3,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

savefig("unit_layer_dist.svg")
