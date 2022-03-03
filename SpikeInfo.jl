using NeuroAnalysis,FileIO,Statistics,StatsPlots,Plots,StatsBase,ProgressMeter,XLSX

function mergespike!(indir;unit=Dict(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spike,siteid = load(joinpath(root,datafile),"spike","siteid")
            delete!.([spike],["unitspike","isspikesorted","unitsync"])
            if haskey(unit,"siteid")
                if siteid == unit["siteid"]
                    newids = setdiff(spike["unitid"],unit["unitid"])
                    if !isempty(newids)
                        newididx = indexin(newids,spike["unitid"])
                        append!(unit["unitid"],newids)
                        append!(unit["unitgood"],spike["unitgood"][newididx])
                        unit["unitposition"] = [unit["unitposition"];spike["unitposition"][newididx,:]]
                        unit["unitwaveform"] = [unit["unitwaveform"];spike["unitwaveform"][newididx,:]]
                        unit["unittemplatewaveform"] = [unit["unittemplatewaveform"];spike["unittemplatewaveform"][newididx,:]]
                        foreach(k->append!(unit["unitfeature"][k],spike["unitfeature"][k][newididx]),keys(spike["unitfeature"]))
                        foreach(k->append!(unit["unittemplatefeature"][k],spike["unittemplatefeature"][k][newididx]),keys(spike["unittemplatefeature"]))
                        append!(unit["unittemplateamplitude"],spike["unittemplateamplitude"][newididx])
                    end
                end
            else
                unit = spike
                unit["siteid"] = siteid
            end
        end
    end
    return unit
end


## Merge Spike Data of a RecordSite
resultroot = "Z:/"
indir = joinpath(resultroot,"AG2","AG2_V1_ODR1")


spikeinfo = (indir) -> begin

unit = mergespike!(indir)
save(joinpath(indir,"unit.jld2"),"unit",unit)

unitid = unit["unitid"];unitwave = unit["unitwaveform"]
unittempwave=unit["unittemplatewaveform"];unitgood=unit["unitgood"]
unitposition=unit["unitposition"];unitfeature = unit["unitfeature"]
unittempfeature=unit["unittemplatefeature"];unittempamp = unit["unittemplateamplitude"]
siteid = unit["siteid"];figfmt = [".png"]


## Unit Position
plotunitposition(unitposition;unitgood,chposition=unit["chposition"],size=(400,700))
foreach(ext->savefig(joinpath(indir,"UnitPosition$ext")),figfmt)

## Unit Density
muc=1.2
r=100
w = replace(unitgood,0=>muc)

n,y = unitdensity(unitposition[:,2];w,bw=80,step=40,r)
plot(n*1e9,y;xlabel="Density (unit/mm³)",ylabel="Depth (μm)",leg=false,size=(400,700),
    lw=2,left_margin=4Plots.mm,title="Count(MU)=$muc, r=$(r)μm")
foreach(ext->savefig(joinpath(indir,"UnitDensity$ext")),figfmt)

df=Dict("feature"=>["Density";;],"depth"=>y,"depthfeature"=>[n;;])

## Unit Feature
un,wn = size(unitwave)
wx = range(-wn/2,step=1,length=wn)

uwy = [1.2*unitwave[i,j]+unitposition[i,2] for i in 1:un,j in 1:wn]
uwx = [0.05*wx[j]+unitposition[i,1] for i in 1:un,j in 1:wn]

plot(uwx',uwy',leg=false,size=(400,700),xlabel="X (μm)",ylabel="Y (μm)",grid=false,lw=2,left_margin=4Plots.mm,
    color=permutedims(map(i->i ? :darkgreen : :gray30,unitgood)))
foreach(ext->savefig(joinpath(indir,"UnitWaveform$ext")),figfmt)


plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread" "leftspread" "rightspread"];xlabel="Spread (μm)",df)
foreach(ext->savefig(joinpath(indir,"UnitFeature_spread$ext")),figfmt)

plotdepthfeature(unittempfeature,unitposition,["uppvinv" "downpvinv"],kw=1000;xlabel="1/Propagation Speed (ms/μm)",df)
foreach(ext->savefig(joinpath(indir,"UnitFeature_invspeed$ext")),figfmt)


plotdepthfeature(unitfeature,unitposition,["amplitude" "peaktroughratio"];kw=[0.1 1],xlabel="A.U.",df)
foreach(ext->savefig(joinpath(indir,"UnitFeature_amp$ext")),figfmt)

plotdepthfeature(unitfeature,unitposition,["halftroughwidth" "halfpeakwidth" "duration"];kw=1000,xlabel="Time (ms)",df)
foreach(ext->savefig(joinpath(indir,"UnitFeature_width$ext")),figfmt)

plotdepthfeature(unitfeature,unitposition,["repolarrate" "recoverrate"];kw=1e-3,xlabel="A.U.",df)
foreach(ext->savefig(joinpath(indir,"UnitFeature_rate$ext")),figfmt)

save(joinpath(indir,"unitdepthfeature.jld2"),"df",df)
end


plotdepthfeature = (feature,position,ks;kw=ones(1,length(ks)),size=(350,700),xlabel="",df=nothing) -> begin
    kn = length(ks)
    vs = hcat(map(k->feature[k],ks)...).*kw
    cs = permutedims(palette(:default).colors.colors[1:kn])
    p=scatter(vs,position[:,2],label=ks,markerstrokewidth=0,alpha=0.6,markersize=3,size=size,
        xlabel=xlabel,ylabel="Depth (μm)",color=cs,left_margin=4Plots.mm)

    yms = [unitdensity(position[:,2],w=vs[:,i],wfun=mean,bw=80,step=40) for i in 1:kn]
    y = yms[1].y
    yms = replace!(hcat(map(i->i.n,yms)...),NaN=>0,Inf=>0,-Inf=>0)
    if df isa Dict
        if haskey(df,"feature")
            df["feature"] = [df["feature"] ks]
            df["depth"] = y
            df["depthfeature"] = [df["depthfeature"] yms]
        else
            df["feature"] = ks
            df["depth"] = y
            df["depthfeature"] = yms
        end
    end

    plot!(p,[zeros(1,kn);yms;zeros(1,kn)],[minimum(y);y;maximum(y)], st=:shape,lw=0,label=false,alpha=0.2,color=cs)
    plot!(p,yms,y;label=false,lw=2,color=cs)
    p
end


## Penetration Sites
# penetration = unique!(meta[:,[:Subject_ID,:RecordSession,:RecordSite]])
# penetration = transform!(penetration,All()=>ByRow((a,b,c)->join(filter!(!isempty,[a,b,c]),"_"))=>:siteid)
# XLSX.writetable(joinpath(resultroot,"penetration.xlsx"),collect(eachcol(penetration)),names(penetration))

## Batch RecordSites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
@showprogress "Batch Spike Info ... " for r in eachrow(penetration)
    spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid))
end




## Collect Units of All RecordSites
function collectunit!(indir;unit=Dict(),datafile="unit.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            u = load(joinpath(root,datafile),"unit")
            delete!(u,"chposition")
            u["siteid"] = fill(u["siteid"],length(u["unitid"]))
            u["id"] = map((s,g,i)->g ? "$(s)_SU$i" : "$(s)_MU$i",u["siteid"],u["unitgood"],u["unitid"])
            if haskey(unit,"id")
                append!(unit["siteid"],u["siteid"])
                append!(unit["id"],u["id"])
                append!(unit["unitid"],u["unitid"])
                append!(unit["unitgood"],u["unitgood"])
                unit["unitposition"] = [unit["unitposition"];u["unitposition"]]
                unit["unitwaveform"] = [unit["unitwaveform"];u["unitwaveform"]]
                unit["unittemplatewaveform"] = [unit["unittemplatewaveform"];u["unittemplatewaveform"]]
                foreach(k->append!(unit["unitfeature"][k],u["unitfeature"][k]),keys(u["unitfeature"]))
                foreach(k->append!(unit["unittemplatefeature"][k],u["unittemplatefeature"][k]),keys(u["unittemplatefeature"]))
                append!(unit["unittemplateamplitude"],u["unittemplateamplitude"])
            else
                unit = u
            end
        end
    end
    return unit
end

allunit = collectunit!(resultroot)
save(joinpath(resultroot,"allunit.jld2"),"allunit",allunit)





## Single Unit Type
allunit = load(joinpath(resultroot,"allunit.jld2"),"allunit")
unitid = allunit["unitid"];unitwave = allunit["unitwaveform"]
unittempwave=allunit["unittemplatewaveform"];unitgood=allunit["unitgood"]
unitposition=allunit["unitposition"];unitfeature = allunit["unitfeature"]
unittempfeature=allunit["unittemplatefeature"];unittempamp = allunit["unittemplateamplitude"]
siteid = allunit["siteid"];figfmt = [".png"]

f1 = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate","ttrough","amplitude"]
f2 = ["upspread","downspread","leftspread","rightspread","uppvinv","downpvinv"]
f = [f1;f2]
F = Float64[hcat(map(k->unitfeature[k],f1)...);;hcat(map(k->unittempfeature[k],f2)...)][unitgood,:]
replace!(F,NaN=>0,Inf=>0,-Inf=>0)

# Y = Float64[abs.(unitfeature[f1[1]]);;hcat(map(k->unitfeature[k],f1[2:end])...);;
#   hcat(map(k->replace(unittempfeature[k],NaN=>0,Inf=>0,-Inf=>0),f2)...)][unitgood,:]

foreach(i->F[:,i]=zscore(F[:,i]),1:size(F,2))
dotplot(permutedims(f),F,leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0)),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
foreach(ext->savefig(joinpath(resultroot,"SpikeFeature$ext")),figfmt)


using Interact,UMAP,Clustering,Distances,PyCall

ft = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate","upspread","downspread","uppvinv","downpvinv"]
ft = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate","upspread","downspread"]
ft = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","upspread","downspread"]
Ft = F[:,indexin(ft,f)]


scatter(Ft[:,1],Ft[:,3],marker=(1,stroke(0)))



D = pairwise(Euclidean(),Ft,dims=1)

cr = hclust(D,branchorder=:optimal)
cid = cutree(cr,k=4)

cr = kmeans(Ft',4)


cr = kmeans(Ft5,6)
cid = assignments(cr)


crs = [kmeans(Ft',i) for i in 1:12]

mi = [mutualinfo(crs[i],crs[end]) for i in 1:11]
ri = [randindex(crs[i],crs[end])[1] for i in 1:11]
vi = [vmeasure(crs[i],crs[end]) for i in 1:11]
si = [mean(silhouettes(crs[i],D)) for i in 2:11]


plot(mi,leg=false)
plot(ri,leg=false)
plot(vi,leg=false)
plot(si,leg=false)

D = pairwise(Euclidean(),Ft,dims=1)
cr = dbscan(D,0.4,20)
cid = cr.assignments


pc = pyimport("pyclustering")

op = pyimport("sklearn.cluster")

Ft5 = umap(Ft', 5, n_neighbors=20, min_dist=1e-6,n_epochs=150)
scatter(Ft5[2,:],Ft5[1,:],marker=(1,stroke(0)))
D = pairwise(Euclidean(),Ft5,dims=2)


cr = dbscan(D,0.75,150)
cid = cr.assignments


clust = op.OPTICS(min_samples=15,xi=0.05,max_eps=1)
cid = clust.fit_predict(Ft)




unique(cid)

plot(cr)






D = rand(10, 10)
D += D' # symmetric distance matrix (optional)
result = hclust(D, linkage=:single)

plot(result)









cm = cgrad(:turbo)

@manipulate for i in 100:300, n in 5:50, d in 0.001:0.001:3
    Ft2 = umap(Ft', 2, n_neighbors=n, min_dist=d, n_epochs=i)
    # scatter(Ft2[2,:], Ft2[1,:],leg=false,grid=false,marker=(1,stroke(0)),size=(600,600),
    #         xticks=[],yticks=[],ratio=:equal,xlabel="UMAP Dimension 2",ylabel="UMAP Dimension 1")
    scatter(Ft2[2,:], Ft2[1,:],group=cid,leg=true,grid=false,marker=(1,stroke(0)),size=(600,600),palette=:tab10,
            xticks=[],yticks=[],ratio=:equal,xlabel="UMAP Dimension 2",ylabel="UMAP Dimension 1")
end



Ft2 = umap(Ft', 2, n_neighbors=15, min_dist=0.05, n_epochs=150)
scatter(Ft2[2,:], Ft2[1,:],group=cid,leg=false,grid=false,marker=(2,stroke(0)),size=(500,450),palette=:tab10,
        xticks=[],yticks=[],ratio=:equal,xlabel="UMAP Dimension 2",ylabel="UMAP Dimension 1")
foreach(ext->savefig(joinpath(resultroot,"UnitCluster$ext")),figfmt)


groupedhist(Ft[:,3],group=cid,bar_position = :stack)

dotplot(permutedims(ft),Ft[cid.==1,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:blue),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),Ft[cid.==2,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:orange),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),Ft[cid.==3,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:green),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),Ft[cid.==4,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:purple),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))



clim=4;ms=0.1
dotplot(permutedims(ft),clamp!(Ft[cid.==1,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[1]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),clamp!(Ft[cid.==2,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[2]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),clamp!(Ft[cid.==3,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[3]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),clamp!(Ft[cid.==4,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[4]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),clamp!(Ft[cid.==5,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[5]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
dotplot!(permutedims(ft),clamp!(Ft[cid.==6,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[6]),
        xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
foreach(ext->savefig(joinpath(resultroot,"SpikeFeature_UnitCluster$ext")),figfmt)

cidx = [findall(cid.==i) for i in 1:length(unique(cid))]

cws = map(i->unitwave[unitgood,:][sample(i,200,replace=false),:],cidx)
ncws = map(i->0.8i./map(j->max(abs.(j)...),extrema(i,dims=2)),cws)

cmws = map(i->dropdims(mean(i,dims=1),dims=1),cws)
csews = map(i->dropdims(std(i,dims=1)/sqrt(200),dims=1),cws)
ncmws = map(i->0.8i./max(abs.(extrema(i))...),cmws)


cs = permutedims(coloralpha.(palette(:tab10).colors.colors[1:6],0.9))
p=plot(leg=false)
foreach(i->plot!(ncws[i]',color=cs[i],lw=3),[2,4])
p


p=plot(leg=false)
foreach(i->plot!(ncmws[i],color=cs[i],lw=3),1:6)
p
foreach(ext->savefig(joinpath(resultroot,"SpikeWave_UnitCluster$ext")),figfmt)



p=plot(cws,group=repeat(1:5,inner=100))
p=plot(cws,group=1:81,leg=false)
repeat([1 2 3 4 5],inner=(1,100))
