using NeuroAnalysis,DataFrames,FileIO,Statistics,LightGraphs,MetaGraphs

# Collect all results and merge into database
ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym")

"Collect Units"
function collectunit!(indir;cell=DataFrame(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spikedata = load(joinpath(root,datafile))
            spike = spikedata["spike"]
            siteid = spikedata["siteid"]
            sui = findall(spike["unitgood"])
            df = DataFrame(site=siteid,id=["$(siteid)_SU$u" for u in spike["unitid"][sui]],
             position = [spike["unitposition"][i:i,:] for i in sui])
            append!(cell,df,cols=:union)
        end
    end
    return unique!(cell,:id)
end

"Collect Layers and Merge into Cell"
function collectlayer(indir;cell=DataFrame(site=[],position=[]),datafile="layer.jld2")
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

"Collect CondTests and Merge into Cell"
function collectcondtest(indir;cell=DataFrame(id=[]),ccode=ccode,datafile="factorresponse.jld2")
    sdf = Dict()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            fr = load(joinpath(root,datafile))
            siteid = fr["siteid"]
            haskey(sdf,siteid) || (sdf[siteid]=DataFrame(id=[]))
            id = ["$(siteid)_SU$u" for u in fr["unitid"]]
            c = ccode[fr["color"]]
            tk = ["responsive","modulative","enoughresponse"]
            isempty(fr["f1f0"]) || push!(tk,"f1f0")
            df1 = ("$(k)!$c"=>fr[k] for k in tk)
            df2 = ("frf_$(k)!$c"=>fr["factorresponsefeature"][k] for k in keys(fr["factorresponsefeature"]))
            df3 = ("pzfr_$(k)!$c"=>map((i,p,m)->(m[i...].-mean(p[i...]))/std(p[i...]),fr["optfri"][k],fr["pfms"],fr["fms"]) for k in keys(fr["fa"]))
            df4 = ("f_$(k)!$c"=>fill(fr["fa"][k],length(id)) for k in keys(fr["fa"]))
            df = DataFrame(df1...,df2...,df3...,df4...)
            df.id = id
            sdf[siteid] = outerjoin(sdf[siteid],df,on=:id)
        end
    end
    return outerjoin(cell,reduce((i,j)->append!(i,j,cols=:union),values(sdf)),on=:id)
end

"Collect Circuits and Merge into CellGraph"
function collectcircuitgraph(indir;datafile="sitecircuit.jld2")
    cs = Dict()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile))
            cs[c["siteid"]] = c["circuit"]
        end
    end
    vs = DataFrame();es=DataFrame()
    for k in keys(cs)
        e=DataFrame(site=k,src=map(first,cs[k].projs),dst=map(last,cs[k].projs),w=cs[k].projweights)
        v=DataFrame(site=k,u=union(e.src,e.dst))
        v.type = map(i->i in cs[k].eunits ? 'E' : i in cs[k].iunits ? 'I' : missing,v.u)
        append!(vs,v);append!(es,e)
    end
    g = MetaDiGraph(nrow(vs))
    foreach(i->set_prop!(g,i,:site,vs.site[i]) && set_prop!(g,i,:id,"$(vs.site[i])_SU$(vs.u[i])") && set_prop!(g,i,:type,vs.type[i]),1:nv(g))
    foreach(r->add_edge!(g,findfirst((vs.site.==r.site) .& (vs.u.==r.src)),findfirst((vs.site.==r.site) .& (vs.u.==r.dst)),:weight,r.w),eachrow(es))
    return g
end

"Collect STAs and Merge into Cell"
function collectsta(indir;cell=DataFrame(id=[]),datafile="stadataset.jld2")
    dfs=DataFrame()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_SU$u" for u in keys(dataset["ulroi"])]
            roic = [v.centerdeg for v in values(dataset["ulroi"])]
            roir = [v.radiusdeg for v in values(dataset["ulroi"])]
            rc = Any[map((c,r)->r ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) for u in keys(dataset["ulroi"])]
            df = DataFrame(id=id,roicenter=roic,roiradius=roir,rc=rc)
            if haskey(dataset,"ulfit")
                foreach(m->df[!,"sta!$m"]=Any[map(f->ismissing(f) ? missing : (;(k=>f[k] for k in setdiff(keys(f),[:resid]))...),dataset["ulfit"][u][m]) for u in keys(dataset["ulroi"])],
                        keys(first(values(dataset["ulfit"]))))
            end
            append!(dfs,df,cols=:union)
        end
    end
    return outerjoin(cell,unique!(dfs,:id),on=:id)
end

## Collect Results
resultroot = "../Result"
indir = joinpath(resultroot,"AF5")

cell = collectunit!(indir)
cell = collectlayer(indir;cell)
cell = collectcondtest(indir;cell)
cell = collectsta(indir;cell)

save(joinpath(resultroot,"cell.jld2"),"cell",cell)
cell = load(joinpath(resultroot,"cell.jld2"),"cell")


cellgraph = collectcircuitgraph(indir)
save(joinpath(resultroot,"cellgraph.jld2"),"cellgraph",cellgraph)
cellgraph = load(joinpath(resultroot,"cellgraph.jld2"),"cellgraph")
