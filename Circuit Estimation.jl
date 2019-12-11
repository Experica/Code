using NeuroAnalysis,Statistics,FileIO,Plots,LightGraphs

# Combine all circuit estimations for one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
settimeunit(1)

subject = "AE9";recordsession = "";recordsite = "u003"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
sitedir = joinpath(resultroot,subject,siteid)
layer = load(joinpath(sitedir,"layer.jld2"),"layer")

# merge circuits
testid="$(siteid)_000"
spike = load(joinpath(sitedir,testid,"spike.jld2"),"spike")
projs,eunits,iunits = load(joinpath(sitedir,testid,"circuit.jld2"),"projs","eunits","iunits")
vprojs,veunits,viunits=checkcircuit(projs,eunits,iunits)
unitgood=trues(length(spike["unitspike"]))

# Unit Layer Circuit Graph
plotcircuit(spike["unitposition"],vprojs,findall(unitgood),unitid=spike["unitid"],eunits=veunits,iunits=viunits,layer=layer)
foreach(i->savefig(joinpath(sitedir,"$(testid)_layer_circuit$i")),[".png",".svg"])










ug = SimpleDiGraphFromIterator((Edge(i) for i in projs))
nn = nv(ug)
nodecolor = [in(i,eunits) ? RGB(1,0.2,0.2) : in(i,iunits) ? RGB(0.2,0.2,1) : RGB(0.4,0.4,0.4) for i in 1:nn]

p=gplot(ug,unitposition[unitgood,1][1:nn],unitposition[unitgood,2][1:nn],nodelabel=1:nn,edgestrokec="gray30",nodefillc=map(c->coloralpha(c,0.5),nodecolor),
nodelabelc=nodecolor,nodesize=3,arrowangleoffset=10/180*pi,arrowlengthfrac=0.025)





pyplot()

p=graphplot(ug,linewidth=1,framestyle=:axes,linecolor=:gray,marker=:circle,markersize=2,markerstrokewidth=0,arrow=arrow(:closed,:head,1,1),
x=unitposition[:,1][1:nn],y=unitposition[:,2][1:nn],names=unitid[1:nn],fontsize=2)
hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
