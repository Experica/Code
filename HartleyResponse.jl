using NeuroAnalysis,FileIO,JLD2,Images,StatsBase,StatsPlots,ProgressMeter,DataFrames,XLSX

ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")

getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function cond2hartley(cond)
    kcond = select(cond,[:Ori,:SpatialFreq,:SpatialPhase]=>ByRow((o,s,p)->map(i->round(i,digits=5),sin2cas(deg2rad(o),s,p)))=>AsTable)
    kx = levels(kcond.kx)
    ky = levels(kcond.ky)
    phase = levels(kcond.phase)
    kx == ky || error("kx and ky values are not match, may not be a valid hartley subspace")

    hi = [(findfirst(i->i==r.ky,ky),findfirst(i->i==r.kx,kx),findfirst(i->i==r.phase,phase)) for r in eachrow(kcond)]
    (;ky,kx,phase,hi)
end

function joinhartley(indir;ccode=ccode,btw=-100:0,rtw=10:130,datafile="sta.jld2")
    stas=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            push!(stas,load(joinpath(root,datafile)))
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
    @showprogress "Join Responses ... " for u in uids
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
    return dataset
end



hartleyinfo = (siteresultdir;figfmt=[".png"]) -> begin
    dataset = joinhartley(siteresultdir)
    # dataset = responsivesta!(dataset)
    # dataset = fitsta!(dataset,model=[:edog,:gabor])
    jldsave(joinpath(siteresultdir,"hartleydataset.jld2");dataset)

    lstadir = joinpath(siteresultdir,"lsta")
    rm(lstadir,force=true,recursive=true)
    mkpath(lstadir)
    for u in keys(dataset["ulsta"])
        plotlstas(dataset,u;dir=lstadir,figfmt)
    end

    lstafitdir = joinpath(siteresultdir,"lstafit")
    rm(lstafitdir,force=true,recursive=true)
    mkpath(lstafitdir)
    for u in keys(dataset["ulfit"]),m in keys(first(values(dataset["ulfit"])))
        plotfitlstas(dataset,u,m;dir=lstafitdir,figfmt)
    end
end

resultroot = "Z:/"


## process all Hartley responses of a RecordSite
subject = "AG1";recordsession = "V1";recordsite = "ODL1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
dataset = load(joinpath(siteresultdir,"hartleydataset.jld2"),"dataset")
