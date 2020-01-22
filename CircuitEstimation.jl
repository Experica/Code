using NeuroAnalysis,Statistics,FileIO,Plots,LightGraphs,MetaGraphs

# Combine all circuit estimations of one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)
layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")

# Merge Circuits
function mergecircuit(resultsitedir;ischeck=true)
    projs=[];eunits=[];iunits=[];projweights=[];unitids=[];unitpositions=[];unitgoods=[]
    for (root,dirs,files) in walkdir(resultsitedir)
        if "circuit.jld2" in files
            c = load(joinpath(root,"circuit.jld2"))
            if !isempty(c["projs"])
                append!(projs,c["projs"]);append!(eunits,c["eunits"]);append!(iunits,c["iunits"]);append!(projweights,c["projweights"])
                s = load(joinpath(root,"spike.jld2"),"spike")
                push!(unitids,s["unitid"]);push!(unitpositions,s["unitposition"]);push!(unitgoods,s["unitgood"])
            end
        end
    end
    unitid=reduce(intersect,unitids)
    uidi=indexin(unitid,unitids[1])
    unitposition = unitpositions[1][uidi,:]
    unitgood = unitgoods[1][uidi]

    intersect!(eunits,eunits,unitid)
    intersect!(iunits,iunits,unitid)
    ui = indexin(unique(projs),projs)
    projs=projs[ui];projweights=projweights[ui]
    ivi=map(p->!isempty(setdiff(p,unitid)),projs)
    deleteat!(projs,ivi);deleteat!(projweights,ivi)
    if ischeck
        projs,eunits,iunits,projweights=checkcircuit(projs,eunits,iunits,projweights)
    end
    return projs,eunits,iunits,projweights,unitid,unitgood,unitposition
end



projs,eunits,iunits,projweights,unitid,unitgood,unitposition = mergecircuit(resultsitedir)



# Layer Circuit Graph
plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,projweights=projweights,layer=layer,showmode=:circuit)
foreach(i->savefig(joinpath(resultsitedir,"Layer_Circuit$i")),[".png",".svg"])
save(joinpath(resultsitedir,"site_circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights,"unitid",unitid,"unitgood",unitgood,"unitposition",unitposition)



## Graph Analysis
G = SimpleDiGraphFromIterator(Edge(i) for i in projs)




neighbors(G,102)



















nodecolor = [in(i,eunits) ? RGB(1,0.2,0.2) : in(i,iunits) ? RGB(0.2,0.2,1) : RGB(0.4,0.4,0.4) for i in 1:nn]

p=gplot(ug,unitposition[unitgood,1][1:nn],unitposition[unitgood,2][1:nn],nodelabel=1:nn,edgestrokec="gray30",nodefillc=map(c->coloralpha(c,0.5),nodecolor),
nodelabelc=nodecolor,nodesize=3,arrowangleoffset=10/180*pi,arrowlengthfrac=0.025)




pyplot()

p=graphplot(ug,linewidth=1,framestyle=:axes,linecolor=:gray,marker=:circle,markersize=2,markerstrokewidth=0,arrow=arrow(:closed,:head,1,1),
x=unitposition[:,1][1:nn],y=unitposition[:,2][1:nn],names=unitid[1:nn],fontsize=2)
hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
