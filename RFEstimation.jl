using NeuroAnalysis,StatsBase,FileIO,Images,StatsPlots,Interact,ProgressMeter

ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym")

function joinsta(indir;ccode=ccode,btw=-100:0,datafile="sta.jld2")
    stas=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            push!(stas,load(joinpath(root,datafile)))
        end
    end
    siteid = stas[1]["siteid"]
    sizepx = stas[1]["sizepx"]
    sizedeg = stas[1]["sizedeg"]
    delays = stas[1]["delays"]
    ppd = (sizepx[1]-1)/sizedeg[1]
    xi = stas[1]["xi"]
    notxi = setdiff(1:prod(sizepx),xi)
    cii=CartesianIndices(sizepx)
    bdi = filter!(i->!isnothing(i),indexin(btw,delays))

    dataset = Dict("sizedeg"=>sizedeg,"sizepx"=>sizepx,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi,"siteid"=>siteid)
    colors = [i["color"] for i in stas]
    colorcodes = map(i->ccode[i],colors)
    corder = sortperm(colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    dataset["log"] = [stas[i]["log"] for i in corder]
    dataset["minmaxcolor"] = [(stas[i]["mincolor"],stas[i]["maxcolor"]) for i in corder]
    usta = Dict();uresponsive=Dict()

    uids = mapreduce(i->keys(stas[i]["usta"]),intersect,corder)
    p = ProgressMeter.Progress(length(uids),desc="Join STAs ... ")
    for u in uids
        csta = Array{Float64}(undef,sizepx...,length(delays),length(corder))
        for j in eachindex(corder)
            zsta = stas[corder[j]]["usta"][u]
            for d in eachindex(delays)
                csta[cii[notxi],d,j].=mean(zsta[d,:])
                csta[cii[xi],d,j] = zsta[d,:]
            end
        end
        usta[u] = csta
        uresponsive[u] = false
        next!(p)
    end
    dataset["usta"] = usta
    dataset["uresponsive"] = uresponsive
    return dataset
end
"check unit sta responsive and cut local sta"
function responsivesta!(dataset;w=0.5,sdfactor=3.5,roimargin=0.1,peakroirlim=1.5)
    sizepx = dataset["sizepx"]
    ppd = dataset["ppd"]
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    uresponsive=dataset["uresponsive"]

    ulsta=Dict();ulcsdd=Dict();ulcroi=Dict();ulroi=Dict();ucresponsive=Dict()
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        clc = localcontrast(usta[u],round(Int,w*ppd))
        cproi = map(i->peakroi(clc[:,:,:,i]),1:size(clc,4))
        crsdd = map(i->isresponsive(usta[u][cproi[i].i,:,i],bdi;sdfactor),1:size(clc,4))
        ucresponsive[u] = map(i->cproi[i].radius <= peakroirlim*ppd ? crsdd[i].r : false,1:size(clc,4))
        if any(ucresponsive[u])
            uresponsive[u] = true
            vroi = cproi[ucresponsive[u]]
            roi = mergeroi(vroi,sizepx;roimargin)
            ulroi[u] = (roi...,centerdeg=roi.center/ppd,radiusdeg=roi.radius/ppd)
            ulsta[u] = usta[u][map(i->i.+(-roi.radius:roi.radius),roi.center)...,:,:]
            ulcroi[u] = cproi
            ulcsdd[u] = crsdd
        else
            uresponsive[u] = false
        end
        next!(p)
    end
    dataset["ulsta"]=ulsta
    dataset["ulcroi"]=ulcroi
    dataset["ulroi"]=ulroi
    dataset["ulcsdd"]=ulcsdd
    dataset["ucresponsive"]=ucresponsive
    return dataset
end
"fit responsive sta to models"
function fitsta!(dataset;model=[:gabor,:dog])
    ulsta = dataset["ulsta"]
    ppd = dataset["ppd"]

    if !haskey(dataset,"ulfit")
        dataset["ulfit"]=Dict()
    end
    ulfit = dataset["ulfit"]
    p = ProgressMeter.Progress(length(ulsta),desc="Fit STAs ... ")
    for u in keys(ulsta)
        if !haskey(ulfit,u)
            ulfit[u] = Dict()
        end
        rs = dataset["ucresponsive"][u]
        ds = map(i->i.d,dataset["ulcsdd"][u])
        for m in model
            ulfit[u][m] = [rs[i] ? fitmodel2(m,ulsta[u][:,:,ds[i],i],ppd) : missing for i in 1:length(rs)]
        end
        next!(p)
    end
    return dataset
end

## process all stas
resultroot = "../Result"
subject = "AF8";recordsession = "HLV1";recordsite = "ODR6"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

dataset = joinsta(siteresultdir)
dataset = responsivesta!(dataset)
dataset = fitsta!(dataset,model=[:dog,:gabor])

save(joinpath(siteresultdir,"stadataset.jld2"),"dataset",dataset)
dataset = load(joinpath(siteresultdir,"stadataset.jld2"),"dataset")

## stas
plotstas=(u,d;dir=nothing)->begin
    cn = length(dataset["ccode"])
    usta = dataset["usta"][u]
    delays = dataset["delays"]
    sizepx = dataset["sizepx"]
    sizedeg = dataset["sizedeg"]
    uresponsive = dataset["uresponsive"][u]
    clim = maximum(abs.(extrema(usta)))
    xlims = [0,sizedeg[2]]
    ylims = [0,sizedeg[1]]
    x = range(xlims...,length=sizepx[2])
    y = range(ylims...,length=sizepx[1])

    title = "Unit_$(u)_STA_$(delays[d])"
    p = plot(layout=(1,cn),legend=false,size=(400cn,600),titlefontcolor=uresponsive ? :green : :match,title=title)
    foreach(i->heatmap!(p[i],x,y,usta[:,:,d,i],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlims,ylims=ylims,xticks=xlims,yticks=ylims,yflip=true,xlabel=dataset["color"][i]),1:cn)
    isnothing(dir) ? p : savefig(joinpath(dir,"$title.png"))
end

@manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    plotstas(u,d)
end

stadir = joinpath(siteresultdir,"sta")
isdir(stadir) || mkpath(stadir)
for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    u == 112 || continue
    plotstas(u,d,dir=stadir)
end
## responsive local stas
plotlstas=(u;dir=nothing)->begin
    cn = length(dataset["ccode"])
    ppd = dataset["ppd"]
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    diameterpx = size(ulsta)[1]
    diameterdeg = diameterpx/ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcd = map(i->i.d,dataset["ulcsdd"][u])
    ulcdex = map((d,i)->exd(ulsta[:,:,d,i]).ex,ulcd,1:cn)
    clim = maximum(abs.(ulcdex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = plot(layout=(1,cn+1),legend=false,size=(350(cn+1),600))
    bar!(p[1],dataset["color"],ulcdex,frame=:zerolines,ylabel="extrema")
    foreach(i->heatmap!(p[i+1],xy,xy,ulsta[:,:,ulcd[i],i],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[i] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,
    xlabel=dataset["color"][i],title="Unit_$(u)_STA_$(delays[ulcd[i]])"),1:cn)
    isnothing(dir) ? p : savefig(joinpath(dir,"Unit_$(u)_STA.png"))
end

@manipulate for u in sort(collect(keys(dataset["ulsta"])))
    plotlstas(u)
end

lstadir = joinpath(siteresultdir,"lsta")
isdir(lstadir) || mkpath(lstadir)
for u in sort(collect(keys(dataset["ulsta"])))
    plotlstas(u,dir=lstadir)
end
## fit stas
plotfitstas=(u,m;dir=nothing)->begin
    cn = length(dataset["ccode"])
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    ppd = dataset["ppd"]
    diameterpx = size(ulsta)[1]
    diameterdeg = diameterpx/ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcd = map(i->i.d,dataset["ulcsdd"][u])
    ulcdex = map((d,i)->exd(ulsta[:,:,d,i]).ex,ulcd,1:cn)
    clim = maximum(abs.(ulcdex))
    xylims = [0,round(diameterdeg,digits=1)]
    xy = range(xylims...,length=diameterpx)

    p = plot(layout=(4,cn),legend=false,size=(400cn,4*250))
    foreach(i->heatmap!(p[1,i],xy,xy,ulsta[:,:,ulcd[i],i],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[i] ? :green : :match,xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ulcd[i]])"),1:cn)

    umfit=dataset["ulfit"][u][m]
    rpx = (diameterpx-1)/2
    x=y=(-rpx:rpx)./ppd
    xylims=[round.(extrema(x),digits=2)...]
    fs = map(f->ismissing(f) ? missing : predict(f,x,y),umfit)
    foreach(i->ismissing(fs[i]) ? plot!(p[2,i],frame=:none) :
    heatmap!(p[2,i],x,y,fs[i],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylims,ylims=xylims,xticks=xylims,yticks=[],xlabel=dataset["color"][i],titlefontcolor=ucresponsive[i] ? :green : :match,title="Fit_$(m)"),1:cn)

    foreach(i->ismissing(fs[i]) ? plot!(p[3,i],frame=:none) :
    histogram!(p[3,i],umfit[i].resid,frame=:semi,xlabel="Residual",grid=false),1:cn)

    foreach(i->ismissing(fs[i]) ? plot!(p[4,i],frame=:none) :
    scatter!(p[4,i],vec(ulsta[:,:,ulcd[i],i]),vec(fs[i]),frame=:semi,grid=false,
    xlabel="y",ylabel="ŷ",titlefontcolor=ucresponsive[i] ? :green : :match,title="r = $(round(umfit[i].r,digits=3))",markerstrokewidth=0,markersize=2),1:cn)
    isnothing(dir) ? p : savefig(joinpath(dir,"Unit_$(u)_Fit_$(m).png"))
end

@manipulate for u in sort(collect(keys(dataset["ulfit"]))),m in collect(keys(first(values(dataset["ulfit"]))))
    plotfitstas(u,m)
end

lstafitdir = joinpath(siteresultdir,"lstafit")
isdir(lstafitdir) || mkpath(lstafitdir)
for u in sort(collect(keys(dataset["ulfit"]))),m in collect(keys(first(values(dataset["ulfit"]))))
    plotfitstas(u,m,dir=lstafitdir)
end
