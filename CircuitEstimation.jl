using NeuroAnalysis,DataFrames,Statistics,FileIO,JLD2,StatsBase,StatsPlots,VegaLite,ProgressMeter,XLSX,LinearAlgebra,
    Graphs,MetaGraphs,GraphMakie,GraphRecipes
import GLMakie as mk

"Merge and check circuits of one recording site"
function mergecircuit(indir;check=true,datafile="circuit.jld2")
    projs=[];eunits=[];iunits=[];projlags=[];projweights=[];siteid=nothing
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile))
            siteid=c["siteid"]
            append!(projs,c["projs"]);append!(eunits,c["eunits"]);append!(iunits,c["iunits"])
            append!(projlags,c["projlags"]);append!(projweights,c["projweights"])
        end
    end
    if check
        projs,eunits,iunits,projlags,projweights = checkcircuit(projs,eunits,iunits,projlags,projweights,debug=true)
    end
    return (;projs,eunits,iunits,projlags,projweights,siteid)
end

function sitecircuit(indir;figfmt = [".png"])
    sc = mergecircuit(indir)
    save(joinpath(indir,"sitecircuit.jld2"),"circuit",sc)

    unit = load(joinpath(indir,"unit.jld2"),"unit")
    unitid = unit["unitid"];unitgood=unit["unitgood"];unitposition=unit["unitposition"]
    layer = load(joinpath(indir,"layer.jld2"),"layer")
    plotcircuit(unitposition,unitid,sc.projs;unitgood,sc.eunits,sc.iunits,layer)
    foreach(ext->savefig(joinpath(indir,"UnitPositionCircuit$ext")),figfmt)
end


resultroot = "Z:/"
figfmt = [".png"]

## Batch RecordSites
# penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
# @showprogress "Batch Site Circuit ... " for r in eachrow(penetration)
#     sitecircuit(joinpath(resultroot,r.Subject_ID,r.siteid))
# end


## Collect All Site Circuits
function collectcircuit(indir;datafile="sitecircuit.jld2")
    cs = []
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile),"circuit")
            push!(cs,c)
        end
    end
    vert = DataFrame();edge=DataFrame()
    getid = (siteid,id)->"$(siteid)_SU$id"
    for c in cs
        e=DataFrame(siteid=c.siteid,src=map(first,c.projs),dst=map(last,c.projs),l=c.projlags,w=c.projweights)
        v=DataFrame(siteid=c.siteid,id=union(e.src,e.dst))
        v.projtype = map(i->i in c.eunits ? 'E' : i in c.iunits ? 'I' : missing,v.id)
        transform!(e,[:siteid,:src]=>ByRow(getid)=>last,[:siteid,:dst]=>ByRow(getid)=>last)
        transform!(v,[:siteid,:id]=>ByRow(getid)=>last)
        append!(vert,v);append!(edge,e)
    end
    return (;vert,edge)
end

unitcircuit = collectcircuit(resultroot)
jldsave(joinpath(resultroot,"unitcircuit.jld2");unitcircuit)
unitcircuit = load(joinpath(resultroot,"unitcircuit.jld2"),"unitcircuit")


function circuit2graph(vs,es)
    g = MetaDiGraph(nrow(vs))
    foreach(i->set_props!(g,i,Dict(k=>vs[i,k] for k in propertynames(vs))),1:nv(g))
    foreach(r->add_edge!(g,findfirst(vs.id.==r.src),findfirst(vs.id.==r.dst),Dict(:weight=>r.w,:l=>r.l)),eachrow(es))
    g
end
function layergraph(es)
    ls = sort(union(es.srclayer,es.dstlayer))
    lps = combine(groupby(es,[:srclayer,:dstlayer]),nrow=>:nlp)
    g = MetaDiGraph(length(ls))
    foreach(i->set_prop!(g,i,:layer,ls[i]),1:nv(g))
    foreach(r->add_edge!(g,findfirst(ls.==r.srclayer),findfirst(ls.==r.dstlayer),:weight,r.nlp),eachrow(lps))
    g
end

allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
layer = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate")
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
penetration = select(penetration,[:siteid,:od,:cofd],:cofd=>ByRow(i->i∈["L/M","S/LM","L/M, S/LM"])=>:incofd,
                :cofd=>ByRow(i->i∈["B","W","B, W"])=>:inbw,:cofd=>ByRow(i->i=="None")=>:none)
csu = leftjoin(subset(allcu,:good),unitcircuit.vert,on=[:siteid,:id])
leftjoin!(csu,penetration,on=:siteid)
es = filter(r->(r.src in csu.id) && (r.dst in csu.id),unitcircuit.edge)
vs = subset(csu,:id=>ByRow(i->i in es.src || i in es.dst))
transform!(es,:src=>ByRow(s->vs.layer[findfirst(vs.id.==s)])=>:srclayer,:dst=>ByRow(d->vs.layer[findfirst(vs.id.==d)])=>:dstlayer)
transform!(es,[:srclayer,:dstlayer]=>ByRow((s,d)->"$s ➡ $d")=>:layerproj)


projpair=transform!(outerjoin(
    transform!(combine(groupby(csu,:siteid),nrow=>:ncsu),:ncsu=>ByRow(i->binomial(i,2))=>:np),
    combine(groupby(es,:siteid),nrow=>:ne),on=:siteid),
    [:ne,:np]=>ByRow((i,j)->i/j*100)=>:pp)

projpair |> [@vlplot(:bar,x={"siteid"},y={"ncsu",title="Number of Cortical Single Unit"});
             @vlplot(:bar,x={"siteid"},y={"pp",title="Projections / Pairs (%)"})]

es |> @vlplot(:bar,y={"layerproj"},x={"count()",title="Number of Projections"})


lg = layergraph(es)
GraphRecipes.graphplot(lg,nodesize=0.1,nodeweights=fill(2,9),nodeshape=:circle,names=map(i->get_prop(g,i,:layer),1:nv(g)),method=:shell,
    fontsize=10,ew=(s,d,w)->0.005*get_prop(g,s,d,:weight))

unitgraph = circuit2graph(vs,es)
jldsave(joinpath(resultroot,"unitgraph.jld2");unitgraph)
unitgraph = load(joinpath(resultroot,"unitgraph.jld2"),"unitgraph")



p=GraphRecipes.graphplot(g,ms=0.1,nodeshape=:rect,names=map(i->get_prop(g,i,:layer),1:nv(g)),method=:shell,
    fontsize=10,ew=(s,d,w)->0.005*get_prop(g,s,d,:weight),curvature=0.2)
p=graphplot(g,arrow_show=true,nlabels=map(i->get_prop(g,i,:layer),1:nv(g)))

nodecolor = [in(i,eunits) ? RGB(1,0.2,0.2) : in(i,iunits) ? RGB(0.2,0.2,1) : RGB(0.4,0.4,0.4) for i in 1:nn]

p=gplot(ug,unitposition[unitgood,1][1:nn],unitposition[unitgood,2][1:nn],nodelabel=1:nn,edgestrokec="gray30",nodefillc=map(c->coloralpha(c,0.5),nodecolor),
nodelabelc=nodecolor,nodesize=3,arrowangleoffset=10/180*pi,arrowlengthfrac=0.025)

display(p)


gg = MetaDiGraph(3)
add_edge!(gg,1,2,:weight,1)
add_edge!(gg,1,3,:weight,5)
add_edge!(gg,3,2,:weight,10)


p=GraphRecipes.graphplot(gg,x=[1,1,1],y=[1,2,3],method=:chorddaigram,curves=true,curvature=0.05,names=1:3,ew=(s,d,w)->get_prop(gg,s,d,:weight))

pyplot()

p=graphplot(ug,linewidth=1,framestyle=:axes,linecolor=:gray,marker=:circle,markersize=2,markerstrokewidth=0,arrow=arrow(:closed,:head,1,1),
x=unitposition[:,1][1:nn],y=unitposition[:,2][1:nn],names=unitid[1:nn],fontsize=2)
hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)






layercolor=Dict("Out"=>RGB(0.1,0.1,0.1),
                "1"=>RGB(0.2,0.9,0.9),
                "23"=>RGB(0.12,0.46,0.7),
                "2"=>[2550,1800],
                "3"=>[2550,1800],
                "4AB"=>RGB(1,0.5,0.05),
                "4A"=>[2575,1800],
                "4B"=>[2430,1800],
                "4C"=>[1800,1800],
                "4Ca"=>RGB(0.17,0.62,0.17),
                "4Cb"=>RGB(0.83,0.15,0.15),
                "56"=>RGB(0.58,0.4,0.74),
                "5"=>[1300,1800],
                "6"=>[1500,1800],
                "WM"=>RGB(0.9,0.9,0.9))

graphplot(cellgraph,nodeshape=:circle,edgecolor=:gray10,linewidth=0.3,curvature_scalar=0.01,arrow=0.13,method=:stress,
        markerstrokecolor=:gray10,markerstrokewidth=0.3,nodecolor=map(i->layercolor[get_prop(cellgraph,i,:layer)],1:nv(cellgraph)))
foreach(ext->savefig(joinpath(resultroot,"graph_layer$ext")),figfmt)









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






g = [0 1 1;
     0 0 1;
     0 1 0]

graphplot(g, names=1:3,edgecolor=:red,linewidth=3, curvature_scalar=0.1,arrow=0.4,dim=2,markersize=0.3)
