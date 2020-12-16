using NeuroAnalysis,Statistics,FileIO,Plots

"Merge and Check Circuits of one recording site"
function mergecircuit(indir;check=true,datafile="circuit.jld2")
    projs=[];eunits=[];iunits=[];projweights=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            c = load(joinpath(root,datafile))
            append!(projs,c["projs"]);append!(eunits,c["eunits"]);append!(iunits,c["iunits"]);append!(projweights,c["projweights"])
        end
    end
    if check
        projs,eunits,iunits,projweights=checkcircuit(projs,eunits,iunits,projweights)
    end
    return (;projs,eunits,iunits,projweights)
end

"Batch Site Circuits"
function batchsitecircuit(indir)
    for site in readdir(indir)
        circuit = mergecircuit(joinpath(indir,site))
        save(joinpath(indir,site,"sitecircuit.jld2"),"circuit",circuit,"siteid",site)
    end
end



resultroot = "../Result"
batchsitecircuit(joinpath(resultroot,"AF5"))

## Layer Circuit Graph
plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,projweights=projweights,layer=layer,showmode=:circuit)
foreach(i->savefig(joinpath(siteresultdir,"Layer_Circuit$i")),[".png",".svg"])









nodecolor = [in(i,eunits) ? RGB(1,0.2,0.2) : in(i,iunits) ? RGB(0.2,0.2,1) : RGB(0.4,0.4,0.4) for i in 1:nn]

p=gplot(ug,unitposition[unitgood,1][1:nn],unitposition[unitgood,2][1:nn],nodelabel=1:nn,edgestrokec="gray30",nodefillc=map(c->coloralpha(c,0.5),nodecolor),
nodelabelc=nodecolor,nodesize=3,arrowangleoffset=10/180*pi,arrowlengthfrac=0.025)




pyplot()

p=graphplot(ug,linewidth=1,framestyle=:axes,linecolor=:gray,marker=:circle,markersize=2,markerstrokewidth=0,arrow=arrow(:closed,:head,1,1),
x=unitposition[:,1][1:nn],y=unitposition[:,2][1:nn],names=unitid[1:nn],fontsize=2)
hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
