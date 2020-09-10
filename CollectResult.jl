using NeuroAnalysis,DataFrames,FileIO

# Collect all results and merge into a database

"Collect Units"
function collectunit(indir;cells=DataFrame(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spikedata = load(joinpath(root,datafile))
            spike = spikedata["spike"]
            siteid = spikedata["siteid"]
            sui = findall(spike["unitgood"])
            cs = DataFrame(site=siteid,id=["$(siteid)_SU$u" for u in spike["unitid"][sui]],
             position = [spike["unitposition"][i:i,:] for i in sui])
            append!(cells,cs)
        end
    end
    return unique!(cells,:id)
end

"Collect Layers and Merge into Cells"
function collectlayer(indir,cells;datafile="layer.jld2")
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
    return transform(cells,[[:site,:position] => ByRow(getdepth) => :depth, [:site,:position] => ByRow(getlayer) => :layer])
end

"Collect Circuits and Merge into Cells"
function collectcircuit(indir,cells;datafile="sitecircuit.jld2")
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
    return transform(cells,[[:site,:id] => ByRow(getoutput) => :output, [:site,:id] => ByRow(getinput) => :input,
                            [:site,:id] => ByRow(getprojtype) => :projtype])
end

"Collect Color and Merge into Cells"
function collectcolor(indir;cells=DataFrame(id=[]),datafile="colordataset.jld2")
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
    return outerjoin(cells,df,on=:id)
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

"Collect RF and Merge into Cells"
function collectrf(indir;cells=DataFrame(id=[]),model=:gabor,datafile="stadataset.jld2")
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

## Collect Results
resultroot = "../Result"

cells = collectunit(joinpath(resultroot,"AF5"))
cells = collectlayer(joinpath(resultroot,"AF5"),cells)
cells = collectcircuit(joinpath(resultroot,"AF5"),cells)

cells = collectcolor(resultroot,cells=cells)
cells = collectorisf(resultroot,cells=cells)
cells = collectrf(resultroot,cells=cells,model=:gabor)

save(joinpath(resultroot,"cells.jld2"),"cells",cells)
