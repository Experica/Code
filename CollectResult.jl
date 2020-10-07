using NeuroAnalysis,DataFrames,FileIO

# Collect all results and merge into database

"Collect Recording Sites"
function collectsite(indir;site=DataFrame(),datafile="siteroi.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            data = load(joinpath(root,datafile))
            siteroi = data["siteroi"]
            siteid = data["siteid"]
            df = DataFrame(site=siteid,roicenter=[siteroi.centerdeg],roiradius=siteroi.radiusdeg)
            append!(site,df)
        end
    end
    return unique!(site,:site)
end

"Collect Units"
function collectunit(indir;cell=DataFrame(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spikedata = load(joinpath(root,datafile))
            spike = spikedata["spike"]
            siteid = spikedata["siteid"]
            sui = findall(spike["unitgood"])
            df = DataFrame(site=siteid,id=["$(siteid)_SU$u" for u in spike["unitid"][sui]],
             position = [spike["unitposition"][i:i,:] for i in sui])
            append!(cell,df)
        end
    end
    return unique!(cell,:id)
end

"Collect Layers and Merge into Cells"
function collectlayer(indir,cell;datafile="layer.jld2")
    ls = Dict()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            l = load(joinpath(root,datafile))
            ls[l["siteid"]] = l["layer"]
        end
    end
    getdepth = (s,p) -> begin
        haskey(ls,s) || return missing
        return -p[2] + ls[s]["Out"][1]
    end
    getlayer = (s,p) -> begin
        haskey(ls,s) || return missing
        return assignlayer(p[2],ls[s])
    end
    return transform(cell,[[:site,:position] => ByRow(getdepth) => :depth, [:site,:position] => ByRow(getlayer) => :layer])
end

"Collect Circuits and Merge into Cells"
function collectcircuit(indir,cell;datafile="sitecircuit.jld2")
    cs = Dict()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile))
            cs[c["siteid"]] = c["circuit"]
        end
    end
    getoutput = (s,id) -> begin
        haskey(cs,s) || return missing
        ts=[]
        foreach(p->"$(s)_SU$(p[1])"==id && push!(ts,"$(s)_SU$(p[2])"),cs[s].projs)
        isempty(ts) ? missing : ts
    end
    getinput = (s,id) -> begin
        haskey(cs,s) || return missing
        ss=[]
        foreach(p->"$(s)_SU$(p[2])"==id && push!(ss,"$(s)_SU$(p[1])"),cs[s].projs)
        isempty(ss) ? missing : ss
    end
    getprojtype = (s,id) -> begin
        haskey(cs,s) || return missing
        for i in cs[s].eunits
            "$(s)_SU$i"==id && return 'E'
        end
        for i in cs[s].iunits
            "$(s)_SU$i"==id && return 'I'
        end
        return missing
    end
    return transform(cell,[[:site,:id] => ByRow(getoutput) => :output, [:site,:id] => ByRow(getinput) => :input,
                            [:site,:id] => ByRow(getprojtype) => :projtype])
end

"Collect Color and Merge into Cells"
function collectcolor(indir;cell=DataFrame(id=[]),datafile="colordataset.jld2")
    df=DataFrame()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile))
            siteid = dataset["siteid"]
            sdf=DataFrame(id=[])
            for i in eachindex(dataset["log"])
                id = ["$(siteid)_SU$u" for u in dataset["uid"][i]]
                tn = contains(dataset["log"][i],"DKL") ? "dkl_" : "hsl_"
                frf = dataset["frf"][i]
                sdf = outerjoin(sdf,DataFrame(:id=>id,Symbol(tn,"frf")=>frf),on=:id)
            end
            append!(df,sdf)
        end
    end
    return outerjoin(cell,df,on=:id)
end

"Collect OriSF and Merge into Cells"
function collectorisf(indir;cells=DataFrame(id=[]),model=:gabor,datafile="stadataset.jld2")
    rfs=DataFrame()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_SU$u" for u in keys(dataset["urf"])]
            rc = Any[map((c,r)->r ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) for u in keys(dataset["urf"])]
            rf = Any[map(f->ismissing(f) ? missing : (;f.model,f.radius,f.param,f.r),dataset["urf"][u][model]) for u in keys(dataset["urf"])]
            ui = map(i->!isempty(skipmissing(i)),rc)
            append!(rfs,DataFrame(id=id[ui],rc=rc[ui],rf=rf[ui]))
        end
    end
    return outerjoin(cells,unique!(rfs,:id),on=:id)
end

"Collect RFs and Merge into Cells"
function collectrf(indir;cell=DataFrame(id=[]),model=:gabor,datafile="stadataset.jld2")
    rfs=DataFrame()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_SU$u" for u in keys(dataset["ulroi"])]
            roic = [v.centerdeg for v in values(dataset["ulroi"])]
            roir = [v.radiusdeg for v in values(dataset["ulroi"])]
            rc = Any[map((c,r)->r ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) for u in keys(dataset["ulroi"])]
            rf=missing
            if haskey(dataset,"urf")
                rf = Any[map(f->ismissing(f) ? missing : (;f.model,f.fun,f.radius,f.param,f.r),dataset["urf"][u][model]) for u in keys(dataset["ulroi"])]
            end
            append!(rfs,DataFrame(id=id,roicenter=roic,roiradius=roir,rc=rc,rf=rf))
        end
    end
    return outerjoin(cell,unique!(rfs,:id),on=:id)
end

## Collect Results
resultroot = "../Result"

site = collectsite(joinpath(resultroot,"AF5"))
save(joinpath(resultroot,"site.jld2"),"site",site)


cell = collectunit(joinpath(resultroot,"AF5"))
cells = collectlayer(joinpath(resultroot,"AF5"),cells)
cells = collectcircuit(joinpath(resultroot,"AF5"),cells)

cells = collectcolor(resultroot,cells=cells)
cells = collectorisf(resultroot,cells=cells)
cell = collectrf(resultroot,cell=cell,model=:gabor)

save(joinpath(resultroot,"cell.jld2"),"cell",cell)
