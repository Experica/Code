using NeuroAnalysis,FileIO,JLD2,DataFrames,StatsPlots,StatsBase,ProgressMeter,VegaLite,Dierckx,Combinatorics,XLSX

ccode = Dict("DKL_X"=>"A","DKL_Y"=>"Y","DKL_Z"=>"S","LMS_Xmcc"=>"L","LMS_Ymcc"=>"M","LMS_Zmcc"=>"S",
             "LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")

getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function collectcondtest(indir;id="Flash2Color",datafile="factorresponse.jld2")
    rs=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            exenv = load(joinpath(root,datafile),"exenv")
            (haskey(exenv,"ID") && exenv["ID"] == id) || continue
            push!(rs,load(joinpath(root,datafile)))
        end
    end
    return rs
end

function coloronoff(x,p)
    if p=="A"
        s = diff(gray.(Gray.(x)))[1] > 0 
    elseif p=="Y"
        s = diff(comp1.(LMS.(x)))[1] > 0 
    elseif p=="S"
        s = diff(comp3.(LMS.(x)))[1] > 0 
    end
    ss = s ? ["-","+"] : ["+","-"]
    map(i->p*i,ss)
end

function joinflash2color_od(rs;ccode=ccode,cs=["A","A"])
    dataset = Dict{String,Any}("siteid"=>rs[1]["siteid"])
    ti = [r["testindex"] for r in rs]
    length(rs)<=3 && return nothing
    ei = indexin([0,1],ti)
    rs = rs[ei]
    colors = [r["exenv"]["color"] for r in rs]
    dataset["ccode"] = map(i->getccode(i,ccode),colors)
    dataset["color"] = colors
    dataset["eye"] = ["N", "D"]

    fa = [map(c->RGBA(c...),r["fa"].Color) for r in rs]
    fac = map(coloronoff,fa,cs)
    dataset["fa"] = fa
    dataset["fac"] = fac
    dataset["fac_od"] = map((c,e)->map(i->"$(e)_$i",c),fac,dataset["eye"])

    uns = map(r->length(r["unitid"]),rs)
    @info "Number of units for each test: $uns"
    uids = mapreduce(r->r["unitid"],intersect,rs) # only include units that spikes in all tests
    uis = map(r->indexin(uids,r["unitid"]),rs)
    dataset["ugood"] = Dict(uids.=>rs[1]["unitgood"][uis[1]])

    getru = (n;k=nothing) -> isnothing(k) ? map((r,i)->r[n][i],rs,uis) : map((r,i)->r[n][k][i],rs,uis)
    dataset["frs"] = Dict(uids.=>zip(getru("frs")...))
    dataset["responsive"] = Dict(uids.=>zip(getru("responsive")...))
    dataset["modulative"] = Dict(uids.=>zip(getru("modulative")...))
    dataset["enoughresponse"] = Dict(uids.=>zip(getru("enoughresponse")...))
    dataset["maxi"] = Dict(uids.=>zip(getru("maxi")...))
    dataset["fpsth"] = (psth=Dict(uids.=>zip(getru("fpsth";k=:psth)...)), x = rs[1]["fpsth"].x)

    return dataset
end

function joinflash2color(rs;ccode=ccode,cs=["A","Y","S"])
    dataset = Dict{String,Any}("siteid"=>rs[1]["siteid"])
    ti = [r["testindex"] for r in rs]
    ei = length(rs)>3 ? indexin([1,2,3],ti) : trues(3)
    rs = rs[ei]
    colors = [r["exenv"]["color"] for r in rs]
    colorcodes = map(i->getccode(i,ccode),colors)
    corder = indexin(cs,colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    rs = rs[corder]
    fa = [map(c->RGBA(c...),r["fa"].Color) for r in rs]
    fac = map(coloronoff,fa,cs)
    dataset["fa"] = fa
    dataset["fac"] = fac

    uns = map(r->length(r["unitid"]),rs)
    @info "Number of units for each test: $uns"
    uids = mapreduce(r->r["unitid"],intersect,rs) # only include units that spikes in all tests
    uis = map(r->indexin(uids,r["unitid"]),rs)
    dataset["ugood"] = Dict(uids.=>rs[1]["unitgood"][uis[1]])

    getru = (n;k=nothing) -> isnothing(k) ? map((r,i)->r[n][i],rs,uis) : map((r,i)->r[n][k][i],rs,uis)
    dataset["frs"] = Dict(uids.=>zip(getru("frs")...))
    dataset["responsive"] = Dict(uids.=>zip(getru("responsive")...))
    dataset["modulative"] = Dict(uids.=>zip(getru("modulative")...))
    dataset["enoughresponse"] = Dict(uids.=>zip(getru("enoughresponse")...))
    dataset["maxi"] = Dict(uids.=>zip(getru("maxi")...))
    dataset["fpsth"] = (psth=Dict(uids.=>zip(getru("fpsth";k=:psth)...)), x = rs[1]["fpsth"].x)

    return dataset
end

plotflash2color_od = (dataset,u;dir=nothing,figfmt=[".png"],csc=["N_A+","N_A-","D_A+","D_A-"])->begin
    fa = vcat(dataset["fa"]...)
    fac_od = vcat(dataset["fac_od"]...)
    um = vcat(map(i->i.m,dataset["frs"][u])...)
    use = vcat(map(i->i.se,dataset["frs"][u])...)
    ug = dataset["ugood"][u]
    ii= indexin(csc,fac_od)

    p = plot(;grid=false,legend=false,tickdir=:out,xlabel="Color",ylabel="Response (spike/s)",size=(600,450))
    bar!(p,fac_od[ii],um[ii];yerror=use[ii],color=fa[ii],msc=:gray,msw=2)
    
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ug ? "S" : "M")U$(u)_od$ext")),figfmt)
end

plotflash2color = (dataset,u;dir=nothing,figfmt=[".png"],csc=["A+","A-","Y+","Y-","S+","S-"])->begin
    fa = vcat(dataset["fa"]...)
    fac = vcat(dataset["fac"]...)
    um = vcat(map(i->i.m,dataset["frs"][u])...)
    use = vcat(map(i->i.se,dataset["frs"][u])...)
    ug = dataset["ugood"][u]
    ii= indexin(csc,fac)

    p = plot(;grid=false,legend=false,tickdir=:out,xlabel="Color",ylabel="Response (spike/s)",size=(600,450))
    bar!(p,fac[ii],um[ii];yerror=use[ii],color=fa[ii],msc=:gray,msw=2)
    
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ug ? "S" : "M")U$(u)$ext")),figfmt)
end

flash2colordataset_od = (siteresultdir;figfmt=nothing) -> begin
    rs = collectcondtest(siteresultdir,id="Flash2Color")
    dataset = joinflash2color_od(rs)
    isnothing(dataset) && return dataset
    jldsave(joinpath(siteresultdir,"flash2colordataset_od.jld2");dataset)
    if !isnothing(figfmt)
        dir = joinpath(siteresultdir,"flash2color");mkpath(dir)
        for u in sort(collect(keys(dataset["ugood"])))
            plotflash2color_od(dataset,u;dir,figfmt)
        end
    end
end

flash2colordataset = (siteresultdir;figfmt=nothing) -> begin
    rs = collectcondtest(siteresultdir,id="Flash2Color")
    dataset = joinflash2color(rs)
    jldsave(joinpath(siteresultdir,"flash2colordataset.jld2");dataset)
    if !isnothing(figfmt)
        dir = joinpath(siteresultdir,"flash2color");mkpath(dir)
        for u in sort(collect(keys(dataset["ugood"])))
            plotflash2color(dataset,u;dir,figfmt)
        end
    end
end

resultroot = "Z:/"

## process all Flash2Color of a RecordSite
subject = "AG1";recordsession = "V1";recordsite = "ODL1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

plotflash2color_od(dataset,750)
plotflash2color(dataset,750)

## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch All Flash2Colors ... " for r in eachrow(penetration)
    # flash2colordataset_od(joinpath(resultroot,r.Subject_ID,r.siteid);figfmt=[".png",".svg"])
    # flash2colordataset(joinpath(resultroot,r.Subject_ID,r.siteid);figfmt=[".png",".svg"])
    flash2colordataset(joinpath(resultroot,r.Subject_ID,r.siteid))
end



