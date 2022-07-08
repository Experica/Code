using NeuroAnalysis,FileIO,JLD2,Images,StatsBase,StatsPlots,Interact,ProgressMeter,DataFrames,XLSX

ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")

getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function joinsta(indir;ccode=ccode,btw=-100:0,rtw=10:130,datafile="sta.jld2")
    stas=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            push!(stas,load(joinpath(root,datafile)))
        end
    end
    ugood1 = stas[1]["ugood"]
    siteid = stas[1]["siteid"]
    sizepx = stas[1]["sizepx"]
    sizedeg = stas[1]["sizedeg"]
    delays = stas[1]["delays"]
    ppd = (sizepx[1]-1)/sizedeg[1]
    xi = stas[1]["xi"]
    notxi = setdiff(1:prod(sizepx),xi)
    cii=CartesianIndices(sizepx)
    bdi = filter!(i->!isnothing(i),indexin(btw,delays))
    rdi = filter!(i->!isnothing(i),indexin(rtw,delays))

    dataset = Dict("sizedeg"=>sizedeg,"sizepx"=>sizepx,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi,"rdi"=>rdi,"siteid"=>siteid)
    colors = [i["exenv"]["color"] for i in stas]
    colorcodes = map(i->getccode(i,ccode),colors)
    corder = sortperm(colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    dataset["eye"] = [stas[i]["exenv"]["eye"] for i in corder]
    dataset["minmaxcolor"] = [(stas[i]["exenv"]["mincolor"],stas[i]["exenv"]["maxcolor"]) for i in corder]
    ugood=Dict();usta = Dict();uresponsive=Dict()

    uns = map(i->length(keys(stas[i]["usta"])),corder)
    @info "Number of units for each RF test: $uns "
    uids = mapreduce(i->keys(stas[i]["usta"]),intersect,corder) # only for unit that spikes in all tests
    @showprogress "Join STAs ... " for u in uids
        csta = Array{Float64}(undef,sizepx...,length(delays),length(corder))
        for j in eachindex(corder)
            zsta = stas[corder[j]]["usta"][u]
            @views for d in eachindex(delays)
                m = median(zsta[d,:]);sd = mad(zsta[d,:],normalize=true)
                csta[cii[notxi],d,j] = m .+ sd*randn(length(notxi)) # white noise with the same distribution
                csta[cii[xi],d,j] = zsta[d,:]
            end
        end
        ugood[u] = ugood1[u]
        usta[u] = csta
        uresponsive[u] = false
    end
    dataset["ugood"] = ugood
    dataset["usta"] = usta
    dataset["uresponsive"] = uresponsive
    return dataset
end
"check unit sta responsive and cut local sta"
function responsivesta!(dataset;w=0.5,sdfactor=6,roimargin=0.1,peakroirlim=(0.15,1.5))
    sizepx = dataset["sizepx"]
    ppd = dataset["ppd"]
    bdi = dataset["bdi"]
    rdi = dataset["rdi"]
    usta = dataset["usta"]
    uresponsive=dataset["uresponsive"]

    ulsta=Dict();ulcsdd=Dict();ulcroi=Dict();ulroi=Dict();ucresponsive=Dict()
    @showprogress "Test STAs ... " for u in keys(usta)
        clc = localcontrast(usta[u],round(Int,w*ppd))
        @views cproi = map(i->peakroi(clc[:,:,:,i]),1:size(clc,4))
        @views crsdd = map(i->isresponsive(usta[u][cproi[i].i,:,i];bi=bdi,ri=rdi,sdfactor),1:size(clc,4))
        ucresponsive[u] = map(i->(peakroirlim[1] <= cproi[i].radius/ppd <= peakroirlim[2]) ? crsdd[i].r : false,1:size(clc,4))
        if any(ucresponsive[u])
            uresponsive[u] = true
            vroi = cproi[ucresponsive[u]]
            roi = mergeroi(vroi;roimargin,imgsize=sizepx,issquare=true)
            ulroi[u] = (roi...,centerdeg=roi.center./ppd,radiideg=roi.radii./ppd,radiusdeg=roi.radius/ppd)
            ulsta[u] = usta[u][map((c,r)->c.+(-r:r),roi.center,roi.radii)...,:,:]
            ulcroi[u] = cproi
            ulcsdd[u] = crsdd
        else
            uresponsive[u] = false
        end
    end
    dataset["ulsta"]=ulsta
    dataset["ulcroi"]=ulcroi
    dataset["ulroi"]=ulroi
    dataset["ulcsdd"]=ulcsdd
    dataset["ucresponsive"]=ucresponsive
    return dataset
end
"fit responsive sta to models"
function fitsta!(dataset;model=[:gabor,:edog])
    ulsta = dataset["ulsta"]
    ppd = dataset["ppd"]

    if !haskey(dataset,"ulfit")
        dataset["ulfit"]=Dict()
    end
    ulfit = dataset["ulfit"]
    @showprogress "Fit STAs ... " for u in keys(ulsta)
        if !haskey(ulfit,u)
            ulfit[u] = Dict()
        end
        rs = dataset["ucresponsive"][u]
        ds = map(i->i.d,dataset["ulcsdd"][u])
        for m in model
            ulfit[u][m] = [rs[i] ? fitmodel2(m,ulsta[u][:,:,ds[i],i],ppd) : missing for i in 1:length(rs)]
        end
    end
    return dataset
end

plotstas=(dataset,u,d;dir=nothing,figfmt=[".png"])->begin
    cn = length(dataset["ccode"])
    usta = dataset["usta"][u]
    ugood = dataset["ugood"][u]
    delays = dataset["delays"]
    sizepx = dataset["sizepx"]
    sizedeg = dataset["sizedeg"]
    uresponsive = dataset["uresponsive"][u]
    clim = maximum(abs.(extrema(usta)))
    xlims = [0,sizedeg[2]]
    ylims = [0,sizedeg[1]]
    x = range(xlims...,length=sizepx[2])
    y = range(ylims...,length=sizepx[1])

    title = "$(ugood ? "S" : "M")U$(u)_STA_$(delays[d])"
    p = plot(;layout=(1,cn),legend=false,size=(400cn,450),titlefontcolor=uresponsive ? :green : :match,title)
    @views foreach(i->heatmap!(p[i],x,y,usta[:,:,d,i];aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),
        xlims,ylims,xticks=xlims,yticks=ylims,yflip=true,xlabel=dataset["color"][i]),1:cn)
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotlstas=(dataset,u;dir=nothing,figfmt=[".png"])->begin
    cn = length(dataset["ccode"])
    ppd = dataset["ppd"]
    ulsta = dataset["ulsta"][u]
    ugood = dataset["ugood"][u]
    delays = dataset["delays"]
    diameterpx = size(ulsta)[1:2]
    diameterdeg = diameterpx./ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcd = map(i->i.d,dataset["ulcsdd"][u])
    @views ulcdex = map((d,i)->exd(ulsta[:,:,d,i]).ex,ulcd,1:cn)
    clim = maximum(abs.(ulcdex))
    xlims = [0,round(diameterdeg[2],digits=1)]
    ylims = [0,round(diameterdeg[1],digits=1)]
    x = range(xlims...,length=diameterpx[2])
    y = range(ylims...,length=diameterpx[1])

    p = plot(layout=(1,cn+1),grid=false,legend=false,size=(350(cn+1),450))
    bar!(p[1],dataset["color"],ulcdex,frame=:zerolines,ylabel="extrema",leftmargin=8mm)
    @views foreach(i->heatmap!(p[i+1],x,y,ulsta[:,:,ulcd[i],i];aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[i] ? :green : :match,xlims,ylims,xticks=xlims,yticks=false,yflip=true,
    xlabel=dataset["color"][i],title="$(ugood ? "S" : "M")U$(u)_STA_$(delays[ulcd[i]])"),1:cn)
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ugood ? "S" : "M")U$(u)_STA$ext")),figfmt)
end

plotfitstas=(dataset,u,m;dir=nothing,figfmt=[".png"])->begin
    cn = length(dataset["ccode"])
    ulsta = dataset["ulsta"][u]
    ugood = dataset["ugood"][u]
    delays = dataset["delays"]
    ppd = dataset["ppd"]
    diameterpx = size(ulsta)[1:2]
    diameterdeg = diameterpx./ppd
    ucresponsive = dataset["ucresponsive"][u]
    ulcd = map(i->i.d,dataset["ulcsdd"][u])
    @views ulcdex = map((d,i)->exd(ulsta[:,:,d,i]).ex,ulcd,1:cn)
    clim = maximum(abs.(ulcdex))
    xlims = [0,round(diameterdeg[2],digits=1)]
    ylims = [0,round(diameterdeg[1],digits=1)]
    x = range(xlims...,length=diameterpx[2])
    y = range(ylims...,length=diameterpx[1])

    p = plot(layout=(3,cn),legend=false,size=(300cn,3*320))
    @views foreach(i->heatmap!(p[1,i],x,y,ulsta[:,:,ulcd[i],i];aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),
    titlefontcolor=ucresponsive[i] ? :green : :match,xlims,ylims,xticks=xlims,yticks=false,yflip=true,
    title="$(ugood ? "S" : "M")U$(u)_STA_$(delays[ulcd[i]])"),1:cn)

    umfit=dataset["ulfit"][u][m]
    rspx = (diameterpx.-1)./2
    x = range(-rspx[2],rspx[2],length=diameterpx[2])./ppd
    y = range(-rspx[1],rspx[1],length=diameterpx[1])./ppd
    xlims=[round.(extrema(x),digits=2)...]
    ylims=[round.(extrema(y),digits=2)...]
    fs = map(f->ismissing(f) ? missing : predict(f,x,y),umfit)
    foreach(i->ismissing(fs[i]) ? plot!(p[2,i],frame=:none) :
    heatmap!(p[2,i],x,y,fs[i];aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),xlims,ylims,
    xticks=xlims,yticks=false,xlabel=dataset["color"][i],titlefontcolor=ucresponsive[i] ? :green : :match,title="Fit_$(m)"),1:cn)

    @views foreach(i->ismissing(fs[i]) ? plot!(p[3,i],frame=:none) :
    heatmap!(p[3,i],x,y,ulsta[end:-1:1,:,ulcd[i],i].-fs[i];aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),xlims,ylims,
    xticks=false,yticks=false,xlabel="Residual",titlefontcolor=ucresponsive[i] ? :green : :match,title="r = $(round(umfit[i].r,digits=3))"),1:cn)
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ugood ? "S" : "M")U$(u)_Fit_$(m)$ext")),figfmt)
end

stainfo = (siteresultdir) -> begin
dataset = joinsta(siteresultdir)
dataset = responsivesta!(dataset)
dataset = fitsta!(dataset,model=[:edog,:gabor])
jldsave(joinpath(siteresultdir,"stadataset.jld2");dataset)

lstadir = joinpath(siteresultdir,"lsta")
rm(lstadir,force=true,recursive=true)
mkpath(lstadir)
for u in keys(dataset["ulsta"])
    plotlstas(dataset,u,dir=lstadir)
end

lstafitdir = joinpath(siteresultdir,"lstafit")
rm(lstafitdir,force=true,recursive=true)
mkpath(lstafitdir)
for u in keys(dataset["ulfit"]),m in keys(first(values(dataset["ulfit"])))
    plotfitstas(dataset,u,m,dir=lstafitdir)
end
end


## process all stas
resultroot = "Z:/"
subject = "AG1";recordsession = "V1";recordsite = "ODL1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
dataset = load(joinpath(siteresultdir,"stadataset.jld2"),"dataset")

## stas
@manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    plotstas(dataset,u,d)
end

# stadir = joinpath(siteresultdir,"sta")
# isdir(stadir) || mkpath(stadir)
# for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
#     u == 265 || continue
#     plotstas(dataset,u,d,dir=stadir)
# end

## responsive local stas
@manipulate for u in sort(collect(keys(dataset["ulsta"])))
    plotlstas(dataset,u)
end

## fit local stas
@manipulate for u in sort(collect(keys(dataset["ulfit"]))),m in collect(keys(first(values(dataset["ulfit"]))))
    plotfitstas(dataset,u,m)
end



## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
@showprogress "Batch All STAs ... " for r in eachrow(penetration)
    stainfo(joinpath(resultroot,r.Subject_ID,r.siteid))
end
