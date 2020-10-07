using NeuroAnalysis,StatsBase,FileIO,Images,StatsPlots,Interact,ProgressMeter

dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL1"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

## Functions
"join color channel stas"
function joinsta(stas;ccode = Dict("DKL_X"=>'A',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S'),btw=[-100:0;200:300])
    siteid = stas[1]["siteid"]
    sizepx = stas[1]["sizepx"]
    sizedeg = stas[1]["sizedeg"]
    delays = stas[1]["delays"]
    ppd = first(sizepx./sizedeg)
    xi = stas[1]["xi"]
    notxi = setdiff(1:prod(sizepx),xi)
    cii=CartesianIndices(sizepx)
    bdi = filter!(i->!isnothing(i),indexin(btw,delays))

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
    p = ProgressMeter.Progress(length(uids),desc="Join STAs ... ")
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
        next!(p)
    end
    dataset["usta"] = usta
    dataset["ucexd"] = ucexd
    dataset["uresponsive"] = uresponsive
    return dataset
end
"check unit sta responsive and cut local sta"
function responsivesta!(dataset;w=0.5,mfactor=3.5,sdfactor=3.5,roimargin=0.1,peakroirlim=1.5)
    sizepx = dataset["sizepx"]
    ppd = dataset["ppd"]
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    ucexd = dataset["ucexd"]
    uresponsive=dataset["uresponsive"]

    ulsta=Dict();ulcexd=Dict();ulcroi=Dict();ulroi=Dict();ucresponsive=Dict()
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        clc = localcontrast(usta[u],round(Int,w*ppd))
        cproi = map(c->peakroi(clc[:,:,:,c]),1:size(clc,4))
        ucresponsive[u] = map(c->cproi[c].radius < peakroirlim*ppd ? isresponsive(usta[u][cproi[c].i,:,c],bdi;mfactor,sdfactor) : false,1:size(clc,4))
        if any(ucresponsive[u])
            uresponsive[u] = true
            vroi = cproi[ucresponsive[u]]
            roi = mergeroi(vroi,sizepx;roimargin)
            ulroi[u] = (roi...,centerdeg=roi.center/ppd,radiusdeg=roi.radius/ppd)
            ulsta[u] = usta[u][map(i->i.+(-roi.radius:roi.radius),roi.center)...,:,:]
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
    siteroi = mergeroi(values(ulroi),sizepx;roimargin)
    dataset["siteroi"] = (siteroi...,centerdeg=siteroi.center/ppd,radiusdeg=siteroi.radius/ppd)

    return dataset
end
"fit responsive sta to each type of models"
function rffit!(dataset;model=[:gabor,:dog])
    ulsta = dataset["ulsta"]
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
            urf[u][m] = [rs[i] ? fitmodel2(m,ulsta[u][:,:,ds[i],i],ppd) : missing for i in 1:length(rs)]
        end
        next!(p)
    end
    return dataset
end
function rfimage(fit;ppd=30,tightrf=false)
    if tightrf
        param = deepcopy(fit.param)
        if fit.model == :dog
        else
            param[[2,4]].=0
            radius = 3*maximum(param[[3,5]])
            fit = (;fit.model,fit.fun,param,radius)
        end
    end
    x=y= -fit.radius:1/ppd:fit.radius
    predict(fit,x,y,yflip=true)
end

## process all stas
testids = ["$(siteid)_HartleySubspace_$i" for i in 1:4]
testn=length(testids)
stadir = joinpath(siteresultdir,"sta")
isdir(stadir) || mkpath(stadir)
lstadir = joinpath(siteresultdir,"lsta")
isdir(lstadir) || mkpath(lstadir)
rffitdir = joinpath(siteresultdir,"rffit")
isdir(rffitdir) || mkpath(rffitdir)

dataset = joinsta(load.(joinpath.(siteresultdir,testids,"sta.jld2")))
dataset = responsivesta!(dataset)
dataset = rffit!(dataset,model=[:gabor])

save(joinpath(siteresultdir,"stadataset.jld2"),"dataset",dataset)
save(joinpath(siteresultdir,"siteroi.jld2"),"siteid",dataset["siteid"],"siteroi",dataset["siteroi"])
dataset = load(joinpath(siteresultdir,"stadataset.jld2"),"dataset")

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
    p = plot(layout=(1,testn),legend=false,size=(400testn,600),titlefontcolor=uresponsive ? :green : :match,title=title)
    foreach(c->heatmap!(p[c],x,y,usta[:,:,d,c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlims,ylims=ylims,xticks=xlims,yticks=ylims,yflip=true,xlabel=dataset["color"][c]),1:testn)
    isnothing(dir) ? p : savefig(joinpath(dir,"$title.png"))
end

@manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    plotstas(u,d)
end
for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    u == 112 || continue
    plotstas(u,d,dir=stadir)
end
## responsive local stas
plotlstas=(u;dir=nothing)->begin
    ppd = dataset["ppd"]
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    diameterpx = size(ulsta)[1]
    diameterdeg = diameterpx/ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcex = map(i->i.ex,dataset["ulcexd"][u])
    ulcd = map(i->i.d,dataset["ulcexd"][u])
    clim = maximum(abs.(ulcex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = plot(layout=(1,testn+1),legend=false,size=(350(testn+1),600))
    bar!(p[1],dataset["color"],ulcex,frame=:zerolines,ylabel="extrema")
    foreach(c->heatmap!(p[c+1],xy,xy,ulsta[:,:,ulcd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[c] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,
    xlabel=dataset["color"][c],title="Unit_$(u)_STA_$(delays[ulcd[c]])"),1:testn)
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
    diameterdeg = diameterpx/ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcex = map(i->i.ex,dataset["ulcexd"][u])
    ulcd = map(i->i.d,dataset["ulcexd"][u])
    clim = maximum(abs.(ulcex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = plot(layout=(4,testn),legend=false,size=(400testn,4*250))
    foreach(c->heatmap!(p[1,c],xy,xy,ulsta[:,:,ulcd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[c] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ulcd[c]])"),1:testn)

    umfit=dataset["urf"][u][m]
    rpx = (diameterpx-1)/2
    x=y=(-rpx:rpx)./ppd
    xylims=[round.(extrema(x),digits=2)...]
    rfs = map(f->ismissing(f) ? missing : predict(f,x,y),umfit)
    foreach(c->ismissing(rfs[c]) ? plot!(p[2,c],frame=:none) :
    heatmap!(p[2,c],x,y,rfs[c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylims,ylims=xylims,xticks=xylims,yticks=xylims,xlabel=dataset["color"][c],titlefontcolor=ucresponsive[c] ? :green : :match,title="Fit_$(m)"),1:testn)

    foreach(c->ismissing(rfs[c]) ? plot!(p[3,c],frame=:none) :
    histogram!(p[3,c],umfit[c].resid,frame=:semi,xlabel="Residual",grid=false),1:testn)

    foreach(c->ismissing(rfs[c]) ? plot!(p[4,c],frame=:none) :
    scatter!(p[4,c],vec(ulsta[:,:,ulcd[c],c]),vec(rfs[c]),frame=:semi,grid=false,
    xlabel="y",ylabel="predicted y",titlefontcolor=ucresponsive[c] ? :green : :match,title="r = $(round(umfit[c].r,digits=3))",markerstrokewidth=0,markersize=2),1:testn)
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

spike = load(joinpath(siteresultdir,testids[1],"spike.jld2"),"spike")
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

spike = load(joinpath(siteresultdir,testids[1],"spike.jld2"),"spike")
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
