using NeuroAnalysis,DataFrames,FileIO,JLD2,StatsBase,ProgressMeter,XLSX,LinearAlgebra,Graphs,MetaGraphs,
    VegaLite,CairoMakie,GraphMakie
import StatsPlots

"Merge and check projections of one recording site"
function mergeprojection(indir;check=true,datafile="projection.jld2",debug=true,usez=true)
    projs=[];projeis=[];lagis=[];lagvs=[];lagzs=[];siteid=nothing;x=nothing
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile))
            siteid=c["siteid"];x=c["x"]
            append!(projs,c["projs"]);append!(projeis,c["projeis"]);
            append!(lagis,c["lagis"]);append!(lagvs,c["lagvs"]);append!(lagzs,c["lagzs"])
        end
    end
    if check
        projs,projeis,lagis,lagvs,lagzs = checkprojection!(projs,projeis,lagis,lagvs,lagzs;debug,usez)
    end
    return (;projs,projeis,lagis,lagvs,lagzs,siteid,x)
end

function sitecircuit(indir;figfmt = [".png"])
    sc = mergeprojection(indir)
    save(joinpath(indir,"sitecircuit.jld2"),"circuit",sc)

    unit = load(joinpath(indir,"unit.jld2"),"unit")
    unitid = unit["unitid"];unitgood=unit["unitgood"];unitposition=unit["unitposition"]
    layer = load(joinpath(indir,"layer.jld2"),"layer")
    plotcircuit(unitposition,unitid,sc.projs;unitgood,layer,projtypes=sc.projeis,projweights=sc.lagzs,showunit=:circuit)
    foreach(ext->savefig(joinpath(indir,"UnitPositionCircuit$ext")),figfmt)
end

resultroot = "Z:/"
figfmt = [".svg",".png"]

## Batch RecordSites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch Site Circuit ... " for r in eachrow(penetration)
    sitecircuit(joinpath(resultroot,r.Subject_ID,r.siteid);figfmt)
end



## Collect All Site Circuits
function collectcircuit(indir;datafile="sitecircuit.jld2")
    cs = []
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            push!(cs,load(joinpath(root,datafile),"circuit"))
        end
    end
    vs = DataFrame();es=DataFrame()
    getid = (siteid,id)->"$(siteid)_SU$id"
    for c in cs
        e=DataFrame(siteid=c.siteid,src=first.(c.projs),dst=last.(c.projs),projei=replace(c.projeis,true=>"E",false=>"I"),
                    lag=map(i->abs(c.x[i]),c.lagis),lagv=c.lagvs,lagz=c.lagzs)
        transform!(e,[:siteid,:src]=>ByRow(getid)=>last,[:siteid,:dst]=>ByRow(getid)=>last)
        v = unique!(outerjoin(select(e,[:siteid,:projei],:src=>:id),select(e,:siteid,:dst=>:id),on=[:siteid,:id]))
        append!(vs,v);append!(es,e)
    end
    return (;vs,es)
end

sucircuit = collectcircuit(resultroot)
jldsave(joinpath(resultroot,"sucircuit.jld2");sucircuit)
sucircuit = load(joinpath(resultroot,"sucircuit.jld2"),"sucircuit")

allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
layer,nlbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","nlbt")
layercolor = StatsPlots.palette(:tab10).colors.colors
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
allcsu = subset(allcu,:good)
vs = innerjoin(allcsu,sucircuit.vs,on=[:siteid,:id])
leftjoin!(vs,penetration[:,[:siteid,:od,:cofd,:pid]],on=:siteid)
es = filter(r->(r.src in vs.id) && (r.dst in vs.id),sucircuit.es)
transform!(es,:src=>ByRow(s->vs.layer[findfirst(vs.id.==s)])=>:srclayer,:dst=>ByRow(d->vs.layer[findfirst(vs.id.==d)])=>:dstlayer)
transform!(es,[:srclayer,:dstlayer]=>ByRow((s,d)->"$s ➡ $d")=>:layerproj)
jldsave(joinpath(resultroot,"csucircuit.jld2");vs,es)


# projection statistics
cdir = joinpath(resultroot,"Circuit");mkpath(cdir)

projpair=transform!(outerjoin(
transform!(combine(groupby(allcsu,:siteid),nrow=>:ncsu),:ncsu=>ByRow(i->binomial(i,2))=>:npair),
combine(groupby(es,:siteid),nrow=>:nproj),
on=:siteid),:nproj=>(i->replace!(i,missing=>0))=>:nproj,[:nproj,:npair]=>ByRow((i,j)->i/j*100)=>:projpair)

pl = projpair |> [@vlplot(:bar,x={"siteid"},y={"ncsu",title="Number of Cortical Single Unit"});
                  @vlplot(:bar,x={"siteid"},y={"projpair",title="Projections / Pairs (%)"})]
foreach(ext->save(joinpath(cdir,"csu_projpair$ext"),pl),figfmt)

pl = vs |> @vlplot(:bar,y={"spiketype"},x={"count()",title="Number of Cortical Single Unit of Circuit"},color={"projei"})
foreach(ext->save(joinpath(cdir,"csu_spiketype_projei$ext"),pl),figfmt)

pl = es |> @vlplot(:bar,y={"layerproj"},x={"count()",title="Number of Projections"})
foreach(ext->save(joinpath(cdir,"csu_layerproj$ext"),pl),figfmt)

pl = es |> [@vlplot(:bar,x={"lag"},y={"count()",title="Number of Projections"})
            @vlplot(:bar,x={"lagv",bin={maxbins=50}},y={"count()",title="Number of Projections"})
            @vlplot(:bar,x={"lagz",bin={maxbins=50,extent=[-30,70]}},y={"count()",title="Number of Projections"})]
foreach(ext->save(joinpath(cdir,"csu_projparam$ext"),pl),figfmt)

# layer projection graph
function layerprojgraph(es;ln=sort!(filter(l->l ∉ ["WM","Out"],collect(keys(layer)))))
    lps = combine(groupby(es,[:srclayer,:dstlayer]),nrow=>:nlp)
    g = MetaDiGraph(length(ln))
    foreach(i->set_prop!(g,i,:name,ln[i]),1:nv(g))
    foreach(r->add_edge!(g,findfirst(ln.==r.srclayer),findfirst(ln.==r.dstlayer),:np,r.nlp),eachrow(lps))
    g,ln
end

plotlayerprojgraph = (g;layout = _ -> map(i->(0,mean(nlbt[i:i+1])),1:nv(g)),node_color=layercolor,nlabels=ln) -> begin
    np = [get_prop(g,e,:np) for e in edges(g)]
    el = [abs(e.src-e.dst) for e in edges(g)]
    f, ax, p = graphplot(g;layout,
    node_color,node_size=30,node_strokewidth=0.5,
    nlabels,nlabels_align=(:center,:center),nlabels_fontsize=8,
    edge_width=0.02np,edge_color=RGBAf(0,0,0,0.9),arrow_size=1.5log2.(np),arrow_shift=0.75,
    selfedge_size=0.03,selfedge_direction=Point2f(0,-1),selfedge_width=1.2,
    curve_distance=0.036el)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    ax.yreversed=true
    f
end

lpg,ln = layerprojgraph(es)
f=plotlayerprojgraph(lpg)
foreach(ext->save(joinpath(cdir,"csu_layerprojgraph$ext"),f),figfmt)

foreach(l->begin
lpg,_ = layerprojgraph(filter(r->r.srclayer==l || r.dstlayer==l,es))
f=plotlayerprojgraph(lpg)
foreach(ext->save(joinpath(cdir,"csu_layerprojgraph_$(replace(l,'/'=>'-'))$ext"),f),figfmt)
end,ln)















function circuit2graph(vs,es)
    g = MetaDiGraph(nrow(vs))
    foreach(i->set_props!(g,i,Dict(k=>vs[i,k] for k in propertynames(vs))),1:nv(g))
    foreach(r->add_edge!(g,findfirst(vs.id.==r.src),findfirst(vs.id.==r.dst),Dict(:weight=>r.w,:l=>r.l)),eachrow(es))
    g
end



unitgraph = circuit2graph(vs,es)
jldsave(joinpath(resultroot,"unitgraph.jld2");unitgraph)
unitgraph = load(joinpath(resultroot,"unitgraph.jld2"),"unitgraph")



## projection of RFs
function plotprojsta(ps,stacell;src=nothing,dst=nothing,ptype="4Cb ⟶ 23",dir=nothing)
    mi(x) = isnothing(x) ? missing : x
    isnothing(src) || isnothing(dst) || (ptype = "$(src) ⟶ $(dst)")
    vpj=dropmissing!(transform!(filter(r->r.projtype == ptype,ps),:src=>ByRow(i->mi(findfirst(j->j==i,stacell.id)))=>:srci,
            :dst=>ByRow(i->mi(findfirst(j->j==i,stacell.id)))=>:dsti),[:srci,:dsti])
    n = nrow(vpj)
    n==0 && return nothing
    p=Plots.plot(layout=(n,2),frame=:none,size=(2*450,n*200),leg=false,titlefontsize=10)
    for i in 1:n
        Plots.plot!(p[i,1],stacell.srimg[vpj.srci[i]],title=vpj.src[i])
        Plots.plot!(p[i,2],stacell.srimg[vpj.dsti[i]],title=vpj.dst[i])
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"$ptype.png"))
end

pstadir = joinpath(resultroot,"psta")
isdir(pstadir) || mkpath(pstadir)
for pt in levels(projsummary.projtype)
    plotprojsta(projsummary,rswcell,ptype=pt,dir=pstadir)
end



function vergenttree(g,v;dir=:in)
    vs = dir==:in ? inneighbors(g,v) : outneighbors(g,v)
    length(vs) < 2 && return missing
    (vs=get_prop.([g],vs,:id), v=get_prop(g,v,:id),dir)
end

vergentsummary = [collect(skipmissing(vergenttree(cellgraph,i,dir=:in) for i in 1:nv(cellgraph)));
                collect(skipmissing(vergenttree(cellgraph,i,dir=:out) for i in 1:nv(cellgraph)))]

function plottreesta(tree,stacell;dir=nothing)
    mi(x) = isnothing(x) ? missing : x
    vi=findfirst(i->i==tree.v,stacell.id)
    isnothing(vi) && return nothing
    vsi = collect(skipmissing(mi(findfirst(i->i==v,stacell.id)) for v in tree.vs))
    n = length(vsi)
    n<2 && return nothing
    p=Plots.plot(layout=(n+1,1),frame=:none,size=(400,(n+1)*200),leg=false,titlefontsize=10)
    title = "$(tree.v) - $(stacell.layer[vi]) - $(tree.dir)"
    Plots.plot!(p[1,1],stacell.srimg[vi],title=title,titlefontsize=15)
    for i in 1:n
        Plots.plot!(p[i+1,1],stacell.srimg[vsi[i]],title="$(stacell.id[vsi[i]]) - $(stacell.layer[vsi[i]])")
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"$title.png"))
end

iostadir = joinpath(resultroot,"iosta")
isdir(iostadir) || mkpath(iostadir)
for t in vergentsummary
    plottreesta(t,rswcell,dir=iostadir)
end






save("UMAP_c_sta_layer.svg",plotunitlayerimage(rswcell.layer,rswcell.imgs,width=2000,markersize=30))


vi = filter(i->"4Cb"==get_prop(cellgraph,i,:layer),1:nv(cellgraph))

outneighbors(cellgraph,9)

outneighbors(cellgraph,14)



t=bfs_tree(cellgraph,9)
t1=dfs_tree(cellgraph,9)
t=union(t,t1)

t=diffusion(cellgraph,1,100,initial_infections=[9])

tt = reduce(union,t)

sg,vmap = induced_subgraph(cellgraph,tt)
sg,vmap = induced_subgraph(cellgraph,collect(edges(t)))
graphplot(sg,nodeshape=:circle,nodecolor=:lightgray,markersize=3,method=:spectrual,names=vmap,dim=3,linewidth=3)






