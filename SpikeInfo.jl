using NeuroAnalysis,FileIO,JLD2,Statistics,StatsPlots,StatsBase,Images,ProgressMeter,DataFrames,XLSX,Dierckx

function mergespikeinfo!(indir;unit=Dict(),datafile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            spike,siteid = load(joinpath(root,datafile),"spike","siteid")
            foreach(k->delete!(spike,k),["unitspike","isspikesorted","unitsync","t0"])
            if haskey(unit,"siteid")
                if siteid == unit["siteid"]
                    newids = setdiff(spike["unitid"],unit["unitid"])
                    if !isempty(newids)
                        newididx = indexin(newids,spike["unitid"])
                        append!(unit["unitid"],newids)
                        append!(unit["unitgood"],spike["unitgood"][newididx])
                        unit["unitposition"] = [unit["unitposition"];spike["unitposition"][newididx,:]]
                        unit["unitwaveforms"] = [unit["unitwaveforms"];spike["unitwaveforms"][newididx,:,:]]
                        unit["unitwaveform"] = [unit["unitwaveform"];spike["unitwaveform"][newididx,:]]
                        foreach(k->append!(unit["unitfeature"][k],spike["unitfeature"][k][newididx]),keys(spike["unitfeature"]))
                        append!(unit["unittemplatemeanamplitude"],spike["unittemplatemeanamplitude"][newididx])
                        unit["unittemplateposition"] = [unit["unittemplateposition"];spike["unittemplateposition"][newididx,:]]
                        unit["unittemplatewaveform"] = [unit["unittemplatewaveform"];spike["unittemplatewaveform"][newididx,:]]
                        foreach(k->append!(unit["unittemplatefeature"][k],spike["unittemplatefeature"][k][newididx]),keys(spike["unittemplatefeature"]))
                        foreach(k->append!(unit["qm"][k],spike["qm"][k][newididx]),keys(spike["qm"]))
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

resultroot = "Z:/"

# ## Get Penetration Sites from Metadata
# penetration = unique!(meta[:,[:Subject_ID,:RecordSession,:RecordSite]])
# penetration = transform!(penetration,All()=>ByRow((a,b,c)->join(filter!(!isempty,[a,b,c]),"_"))=>:siteid)
# XLSX.writetable(joinpath(resultroot,"penetration.xlsx"),collect(eachcol(penetration)),names(penetration))

## Merge ALL Spike Info of a RecordSite
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch Merging Spike Info ... " for r in eachrow(penetration)
    indir = joinpath(resultroot,r.Subject_ID,r.siteid)
    unit = mergespikeinfo!(indir)
    save(joinpath(indir,"unit.jld2"),"unit",unit)
end


## Spike Info of a RecordSite
figfmt = [".svg"]
indir = joinpath(resultroot,"AG2","AG2_V1_ODL3")
layer = load(joinpath(indir,"layer.jld2"),"layer")

plotdepthfeature = (feature,position,ks;good=trues(size(position,1)),kw=ones(1,length(ks)),size=(350,700),xlabel="",df=nothing,layer=nothing,grid=true,xticks=:auto,vl=[]) -> begin
    kn = length(ks)
    vs = hcat(map(k->feature[k][good],ks)...).*kw
    cs = permutedims(palette(:default).colors.colors[1:kn])
    p=plot(;size,grid,tickdir=:out)
    isempty(vl) || vline!(p,vl;linecolor=:gray10,legend=false,lw=1)
    scatter!(p,vs,position[good,2];label=ks,markerstrokewidth=0,alpha=0.6,markersize=3,xticks,
        xlabel,ylabel="Depth (μm)",color=cs,left_margin=4Plots.mm)

    yms = [unitdensity(position[good,2],w=vs[:,i],wfun=mean,bw=40,step=20,s=1.4) for i in 1:kn]
    y = yms[1].y
    yms = hcat(map(i->i.n,yms)...)
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
    if !isnothing(layer)
        xmin,xmax=extrema(vs)
        xm = 0.02(xmax-xmin)
        xmin-=7xm; xmax+=xm
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    p
end

spikeinfo = (indir;layer=nothing,figfmt=[".png"]) -> begin

    unit = load(joinpath(indir,"unit.jld2"),"unit")

    unitid = unit["unitid"];unitgood=unit["unitgood"]
    unittempposition=unit["unittemplateposition"];unittempamp = unit["unittemplatemeanamplitude"]
    unittempwave=unit["unittemplatewaveform"];unittempfeature=unit["unittemplatefeature"]
    unitposition=unit["unitposition"];unitwave = unit["unitwaveform"];unitfeature = unit["unitfeature"]
    unitqm=unit["qm"];siteid = unit["siteid"]

    if layer == :batch
        layerpath = joinpath(indir, "layer.jld2")
        layer = ispath(layerpath) ? load(layerpath, "layer") : nothing
    end

    # Unit Position
    plotunitposition(unitposition;unitgood,chposition=unit["chposition"],size=(400,700),layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitPosition$ext")), figfmt)

    # Unit Density
    muc=1.2
    r=100
    w = replace(unitgood,0=>muc)

    # n,y = unitdensity(unitposition[unitgood,2];bw=40,step=20,r,s=1.4)
    n,y = unitdensity(unitposition[:,2];w,bw=40,step=20,r,s=1.4)
    plot(1e9n,y;xlabel="Density (unit/mm³)",ylabel="Depth (μm)",leg=false,size=(350,700),grid=false,tickdir=:out,
        lw=2,left_margin=4Plots.mm,title="Count(MU)=$muc, r=$(r)μm")
    if !isnothing(layer)
        xmin,xmax = extrema(1e9n)
        xm = 0.02(xmax-xmin)
        xmin-=xm;xmax+=xm
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!([l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitDensity$ext")), figfmt)

    df=Dict("feature"=>["Density";;],"depth"=>y,"depthfeature"=>[1e9n;;],"siteid"=>siteid)

    # Unit Feature
    un,wn = size(unitwave)
    wx = range(-wn/2,step=1,length=wn)

    uwy = [1e-2*unitwave[i,j]+unitposition[i,2] for i in 1:un,j in 1:wn]
    uwx = [0.05*wx[j]+unitposition[i,1] for i in 1:un,j in 1:wn]

    plot(uwx',uwy',leg=false,size=(350,700),xlabel="X (μm)",ylabel="Y (μm)",grid=false,lw=1,left_margin=4Plots.mm,alpha=0.6,tickdir=:out,
        color=permutedims(map(i->i ? RGB(0,0.55,0) : RGB(0),unitgood)),title="Unit Waveform")
    if !isnothing(layer)
        xmin=1.5;xmax=43
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!([l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann,xlims=(xmin,xmax),
        yticks=0:200:3820)
    end
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitWaveform$ext")), figfmt)


    # n=findfirst(unitid.==432)
    # plot!(uwx[n,:],uwy[n,:],leg=false,lw=1,color=:red)
    # plot(0.1*(-40:40) .+ unit["chposition"][:,1]',  0.03*unit["unitwaveforms"][n,:,:] .+ unit["chposition"][:,2]';
    #     leg=false,ratio=0.1,size=(600,3840),color=:dodgerblue,lw=2,frame=:none)

    plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread"];xlabel="Spread (μm)",layer,grid=false,vl=[0])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")Spread$ext")), figfmt)
    plotdepthfeature(unitfeature,unitposition,["duration"];kw=1000,xlabel="Time (ms)",layer,grid=false,vl=[0])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")Duration$ext")), figfmt)
    plotdepthfeature(unitfeature,unitposition,["peaktroughratio"];xlabel="Peak/Trough",layer,grid=false,xticks=-1:-2:-15,vl=[-1])
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")PTRatio$ext")), figfmt)

    # plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread" "leftspread" "rightspread"];xlabel="Spread (μm)",df,layer,good=unitgood)
    plotdepthfeature(unittempfeature,unitposition,["upspread" "downspread" "leftspread" "rightspread"];xlabel="Spread (μm)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_spread$ext")), figfmt)

    # plotdepthfeature(unittempfeature,unitposition,["uppvinv" "downpvinv"],kw=1000;xlabel="1/Propagation Speed (ms/μm)",df,layer,good=unitgood)
    plotdepthfeature(unittempfeature,unitposition,["uppvinv" "downpvinv"],kw=1000;xlabel="1/Propagation Speed (ms/μm)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_invspeed$ext")), figfmt)


    # plotdepthfeature(unitfeature,unitposition,["amplitude" "peaktroughratio"];kw=[1e-3 1],xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["amplitude" "peaktroughratio"];kw=[1e-3 1],xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_amp$ext")), figfmt)

    # plotdepthfeature(unitfeature,unitposition,["halftroughwidth" "halfpeakwidth" "duration"];kw=1000,xlabel="Time (ms)",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["halftroughwidth" "halfpeakwidth" "duration"];kw=1000,xlabel="Time (ms)",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_width$ext")), figfmt)

    # plotdepthfeature(unitfeature,unitposition,["repolarrate" "recoverrate"];kw=1e-3,xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitfeature,unitposition,["repolarrate" "recoverrate"];kw=1e-3,xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_rate$ext")), figfmt)

    # plotdepthfeature(unitqm,unitposition,["fr" "fp" "pisi"];xlabel="A.U.",df,layer,good=unitgood)
    plotdepthfeature(unitqm,unitposition,["fr" "fp" "pisi"];xlabel="A.U.",df,layer)
    foreach(ext -> savefig(joinpath(indir, "$(isnothing(layer) ? "" : "Layer_")UnitFeature_qm$ext")), figfmt)

    save(joinpath(indir,"unitdepthfeature.jld2"),"df",df)
end

## Batch Spike Info
@showprogress "Batch Spike Info ... " for r in eachrow(penetration)
    # spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid))
    spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)
end



## Collect Units of All RecordSites
function collectunit!(indir;unit=Dict(),datafile="unit.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            u = load(joinpath(root,datafile),"unit")
            foreach(k->delete!(u,k),["chposition","unitwaveforms","unittemplatemeanamplitude"])
            u["siteid"] = fill(u["siteid"],length(u["unitid"]))
            u["id"] = map((s,g,i)->"$(s)_$(g ? "S" : "M")U$i",u["siteid"],u["unitgood"],u["unitid"])
            if haskey(unit,"id")
                append!(unit["siteid"],u["siteid"])
                append!(unit["id"],u["id"])
                append!(unit["unitid"],u["unitid"])
                append!(unit["unitgood"],u["unitgood"])
                unit["unitposition"] = [unit["unitposition"];u["unitposition"]]
                unit["unitwaveform"] = [unit["unitwaveform"];u["unitwaveform"]]
                foreach(k->append!(unit["unitfeature"][k],u["unitfeature"][k]),keys(u["unitfeature"]))
                unit["unittemplateposition"] = [unit["unittemplateposition"];u["unittemplateposition"]]
                unit["unittemplatewaveform"] = [unit["unittemplatewaveform"];u["unittemplatewaveform"]]
                foreach(k->append!(unit["unittemplatefeature"][k],u["unittemplatefeature"][k]),keys(u["unittemplatefeature"]))
                foreach(k->append!(unit["qm"][k],u["qm"][k]),keys(u["qm"]))
            else
                unit = u
            end
        end
    end
    return unit
end

allunit = collectunit!(resultroot)

# Add layer info for all units
alllayer = load(joinpath(resultroot,"alllayer.jld2"),"alllayer")
layertemplate,lbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","lbt")
transform!(alllayer,"1"=>ByRow(i->(x->last(i).-x))=>"e2cfun")
transform!(alllayer,vcat.("e2cfun",names(alllayer,Not([:siteid,:e2cfun,:GM]))) .=> ByRow((f,x)->f(x)) => last)
layerboundary = combine(alllayer,names(alllayer,Not([:siteid,:e2cfun,:GM])) .=> ByRow(last)=>identity)
alllayer.tcfun = [gettcfun(collect(r),lbt,true) for r in eachrow(layerboundary)] # gettcfun defined in "LayerEstimation.jl"

usi = map(i->findfirst(j->j==i,alllayer.siteid),allunit["siteid"])
allunit["unitaligndepth"] = [alllayer.e2cfun[usi[i]](allunit["unitposition"][i,2]) |> alllayer.tcfun[usi[i]] for i in eachindex(usi)]
allunit["unitlayer"] = assignlayer.(allunit["unitaligndepth"],[layertemplate])

# assign GM units that's been assigned as WM units
for i in eachindex(usi)
    allunit["unitlayer"][i] == "WM" || continue
    gmb = alllayer.GM[usi[i]]
    ismissing(gmb) && continue
    (gmb[begin]<= allunit["unitposition"][i,2] < gmb[end]) && (allunit["unitlayer"][i]="GM")
end

jldsave(joinpath(resultroot,"allunit.jld2");allunit)



## Unit Feature
allunit = load(joinpath(resultroot,"allunit.jld2"),"allunit")
unitid = allunit["id"];unitwave = allunit["unitwaveform"];unittempfeature = allunit["unittemplatefeature"]
unitgood=allunit["unitgood"];unitposition=allunit["unitposition"];unitfeature = allunit["unitfeature"]
unitaligndepth=allunit["unitaligndepth"];unitlayer=allunit["unitlayer"];unitqm=allunit["qm"]
figfmt = [".svg",".png"]
unitfeaturedir = joinpath(resultroot,"UnitFeature")
mkpath(unitfeaturedir)

f1 = ["duration","peaktroughratio","lefthalftroughwidth","righthalftroughwidth","halftroughwidth","lefthalfpeakwidth","righthalfpeakwidth","halfpeakwidth","repolarrate","recoverrate","ttrough","amplitude"]
f2 = ["upspread","downspread","leftspread","rightspread","uppvinv","downpvinv"]
fqm= ["fr","pisi","fp"]
f = [f1;f2;fqm]
F = Float64[hcat(map(k->unitfeature[k],f1)...);;hcat(map(k->unittempfeature[k],f2)...);;hcat(map(k->unitqm[k],fqm)...)]
replace!(F,NaN=>0,Inf=>0,-Inf=>0)

@views Fz = hcat(map(i->zscore(F[:,i]),1:size(F,2))...)
dotplot(permutedims(1:length(f)),Fz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="Unit Spike Feature",size=(750,600),ann=(2,15,"n=$(size(Fz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"unit_feature$ext")),figfmt)


cui = 0 .<= unitaligndepth .<= 1 # cortical units
csui = cui .& unitgood # cortical single units
gmui = unitlayer.=="GM" # GM units
gmsui = gmui .& unitgood # GM single units

cuF = F[cui,:]
@views cuFz = hcat(map(i->zscore(cuF[:,i]),1:size(cuF,2))...)
dotplot(permutedims(1:length(f)),cuFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="Cortical Unit Spike Feature",size=(750,600),ann=(2,15,"n=$(size(cuFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"cu_feature$ext")),figfmt)

csuF = F[csui,:]
@views csuFz = hcat(map(i->zscore(csuF[:,i]),1:size(csuF,2))...)
dotplot(permutedims(1:length(f)),csuFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="Cortical Single-Unit Spike Feature",size=(750,600),ann=(2,15,"n=$(size(csuFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_feature$ext")),figfmt)

gmuF = F[gmui,:]
@views gmuFz = hcat(map(i->zscore(gmuF[:,i]),axes(gmuF,2))...)
dotplot(permutedims(1:length(f)),gmuFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="GM Unit Spike Feature",size=(750,600),ann=(2,15,"n=$(size(gmuFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"gmu_feature$ext")),figfmt)

gmsuF = F[gmsui,:]
@views gmsuFz = hcat(map(i->zscore(gmsuF[:,i]),axes(gmsuF,2))...)
dotplot(permutedims(1:length(f)),gmsuFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="GM Single-Unit Spike Feature",size=(750,600),ann=(2,15,"n=$(size(gmsuFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_feature$ext")),figfmt)


## Single-Unit Classification
using UMAP,Clustering,Distances,VegaLite

wn = size(unitwave,2)
wx = collect(range(-wn/2,step=1,length=wn))
csuwy = permutedims(unitwave[csui,:])
csuwx = repeat(wx,outer=(1,size(csuwy,2)))
gmsuwy = permutedims(unitwave[gmsui,:])
gmsuwx = repeat(wx,outer=(1,size(gmsuwy,2)))

# cf = filter(i->i ∉ ["ttrough","fp","leftspread","rightspread"],f) # all possible features
# cf = filter(i->i ∉ ["ttrough","fp","leftspread","rightspread","fr","pisi","amplitude"],f) # all probable features
cf = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate","uppvinv","downpvinv","upspread","downspread"] # 1D+2D features
# cf = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate","uppvinv","downpvinv"] # 1D+2D features
# cf = ["duration","peaktroughratio","halftroughwidth","halfpeakwidth","repolarrate","recoverrate"] # 1D features

cFz = csuFz[:,indexin(cf,f)]
cFzD = pairwise(Euclidean(),cFz,dims=1)

dotplot(permutedims(1:length(cf)),cFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(cf),xformatter=i->cf[Int(i)],
        xlabel="Cortical Single-Unit Clustering Feature",size=(750,600),ann=(2,10,"n=$(size(cFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature$ext")),figfmt)

cFz2 = umap(cFz', 2, n_neighbors=25, min_dist=0.8, n_epochs=300,metric=Euclidean())
scatter(cFz2[2,:], cFz2[1,:],leg=false,msw=0,ms=2,ma=0.8,size=(600,600),frame=:none,ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_umap$ext")),figfmt)
plot(cFz2[2,:]' .+ 0.85e-2*csuwx, cFz2[1,:]' .+ 2.8e-5*csuwy,color=:darkgreen,alpha=0.8,leg=false,frame=:none,lw=0.5,size=(800,800),ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_umap_wave$ext")),figfmt)


gmFz = gmsuFz[:,indexin(cf,f)]
gmFzD = pairwise(Euclidean(),gmFz,dims=1)

dotplot(permutedims(1:length(cf)),gmFz,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(cf),xformatter=i->cf[Int(i)],
        xlabel="GM Single-Unit Clustering Feature",size=(750,600),ann=(2,10,"n=$(size(gmFz,1))"))
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature$ext")),figfmt)

gmFz2 = umap(gmFz', 2, n_neighbors=25, min_dist=0.8, n_epochs=300,metric=Euclidean())
scatter(gmFz2[2,:], gmFz2[1,:],leg=false,msw=0,ms=2,ma=0.8,size=(600,600),frame=:none,ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_umap$ext")),figfmt)
plot(gmFz2[2,:]' .+ 0.85e-2*gmsuwx, gmFz2[1,:]' .+ 2.8e-5*gmsuwy,color=:darkgreen,alpha=0.8,leg=false,frame=:none,lw=0.5,size=(800,800),ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_umap_wave$ext")),figfmt)

# try to get optimal number of clusers
function noclu(F,FD;ks=2:10,itr=100)
    ok=ones(itr)
    for i in 1:itr
        cs = [kmeans(F',k) for k in ks]
        ss = [mean(Clustering.silhouettes(c,FD)) for c in cs]
        ok[i] = ks[argmax(ss)]
    end
    ok
end

ks = 3:6
oks = noclu(cFz,cFzD;ks)
density(oks,xticks=ks,xlabel="k",ylabel="Silhouette",leg=false)

k=5
kc = kmeans(cFz',k)
clu = assignments(kc)

ks = 3:6
oks = noclu(gmFz,gmFzD;ks)
density(oks,xticks=ks,xlabel="k",ylabel="Silhouette",leg=false)

k=5
kc = kmeans(gmFz',k)
gmclu = assignments(kc)

# hc = hclust(cFzD,linkage=:ward,branchorder=:optimal)
# plot(hc,xticks=false)
# clu = cutree(hc;k)

scatter(cFz2[2,:], cFz2[1,:],group=clu,c=cgrad(:tab10)[clu],leg=:inline,frame=:none,msw=0,ms=2,ma=0.8,size=(600,600),ratio=1,legendfontsize=12)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_umap_clu$ext")),figfmt)
plot(cFz2[2,:]' .+ 0.85e-2*csuwx, cFz2[1,:]' .+ 2.8e-5*csuwy,c=cgrad(:tab10)[clu'],alpha=0.8,leg=false,frame=:none,lw=0.5,size=(800,800),ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_umap_wave_clu$ext")),figfmt)

scatter(gmFz2[2,:], gmFz2[1,:],group=gmclu,c=cgrad(:tab10)[gmclu],leg=:inline,frame=:none,msw=0,ms=2,ma=0.8,size=(600,600),ratio=1,legendfontsize=12)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_umap_clu$ext")),figfmt)
plot(gmFz2[2,:]' .+ 0.85e-2*gmsuwx, gmFz2[1,:]' .+ 2.8e-5*gmsuwy,c=cgrad(:tab10)[gmclu'],alpha=0.8,leg=false,frame=:none,lw=0.5,size=(800,800),ratio=1)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_umap_wave_clu$ext")),figfmt)


plotclufeature=(F,f,clu;lim=maximum(abs.(F)))->begin
    cluid = sort(unique(clu))
    p = plot(leg=false,grid=false,size=(750,600),xlabel="Spike Feature",ylabel="Z Score",xrotation=20,ann=(2,lim-1,"n=$(size(F,1))"))
    foreach(i->dotplot!(permutedims(1:length(f)),clamp!(F[clu.==i,:],-lim,lim),marker=(2,0.3, stroke(0),cgrad(:tab10)[i]),xticks=1:length(f),xformatter=i->f[Int(i)]),cluid)
    p
end

plotclufeature(cFz,cf,clu,lim=7)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_clu$ext")),figfmt)
plotclufeature(gmFz,cf,gmclu,lim=7)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_clu$ext")),figfmt)

pl = flatten(DataFrame(f=cf,F=[cFz[:,i] for i in eachindex(cf)],clu=fill(clu,length(cf))),[:F,:clu]) |> @vlplot(width=200,height=300,
    mark={:boxplot, extent="min-max"},
    x={"clu:n",axis={title="Cluster",titleFontSize=24,titleFontWeight="normal",labelAngle=0,labelFontSize=18}},
    y={:F, axis={title="Z Score",titleFontSize=20,titleFontWeight="normal",grid=false,labelFontSize=14}},
    color={"clu:n",scale={scheme=:category10}},
    column={"f",sort=cf,header={labelFontSize=26,labelFontWeight="bold",title=""}}
)
foreach(ext->save(joinpath(unitfeaturedir,"csu_clu_cfeature$ext"),pl),figfmt)


# @df DataFrame(f="duration",clu=clu,F=csuF[:,findfirst(i->i=="pisi",f)]) boxplot(:clu,:F,fillalpha=0.75, linewidth=2,ylabel="Duration (μm)",xlabel="Cluster")

# DataFrame(f="duration",clu=clu,F=csuF[:,findfirst(i->i=="duration",f)]) |> @vlplot(width=200,height=300,
# mark={:boxplot, extent="min-max"},
# x={"clu:n",axis={title="Cluster",titleFontSize=24,titleFontWeight="normal",labelAngle=0,labelFontSize=18}},
# y={:F, axis={title="Z Score",titleFontSize=20,titleFontWeight="normal",grid=false,labelFontSize=14}},
# color={"clu:n",scale={scheme=:category10}}
# )


plotcluwaveform=(wys,clu;fs=30e3,isalign=true,isnorm=true,ismean=true,n=50,color=:tab10,clucode=nothing)->begin
    wn,un=size(wys)
    cluid = sort(unique(clu))
    cluname = isnothing(clucode) ? cluid : map(i->clucode[i],cluid)
    wys = hlpass(wys',fs,low=8500)'
    if isalign
        wx = range(-round(Int,wn/2),step=1,length=wn)
        iwx = -wn:0.25:wn # double x range, upsampling 4 times
        fs *= 4
        @views iwys = hcat(map(i->Spline1D(wx,wys[:,i],k=3,bc="nearest")(iwx),1:un)...)
        r = 4wx[begin]:4wx[end] # original x range
        @views ws = hcat(map(i->iwys[r .+ argmax(abs.(iwys[:,i])),i],1:un)...)
    else
        ws = wys
    end
    isnorm && (ws ./= maximum(abs.(ws),dims=1))
    ci = [c.==clu for c in cluid]
    cws = map(i->ws[:,i],ci)
    if ismean
        cwm = hcat(map(i->mean(i,dims=2),cws)...)
        cwse = hcat(map(i->std(i,dims=2)/sqrt(size(i,2)),cws)...)
        plot(cwm;ribbon=cwse,label=permutedims(cluname),color=cgrad(color)[cluid'],lw=2,xformatter=i->round(i/fs*1000,digits=2),xlabel="Time (ms)",ylabel="a.u.")
    else
        colors=cgrad(color)[cluid]
        @views csw = map(i->i[:,sample(1:size(i,2),min(n,size(i,2)),replace=false)],cws)
        p=plot(leg=false,xlabel="Time (ms)",ylabel="a.u.",xformatter=i->round(i/fs*1000,digits=2))
        for i in eachindex(csw)
            plot!(p,csw[i],color=colors[i],alpha=0.3,lw=1)
        end
        for i in eachindex(csw)
            plot!(p,mean(csw[i],dims=2),color=colors[i],alpha=1,lw=3)
        end
        p
    end
end

clucode = Dict(1=>"tp",2=>"E",3=>"A",4=>"I",5=>"En")
plotcluwaveform(csuwy,clu;ismean=true,clucode)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_clu_mwave$ext")),figfmt)
plotcluwaveform(csuwy,clu;ismean=false)
foreach(ext->savefig(joinpath(unitfeaturedir,"csu_cfeature_clu_wave$ext")),figfmt)

gmclucode = Dict(1=>"E",2=>"I",3=>"En",4=>"tp",5=>"A")
plotcluwaveform(gmsuwy,gmclu;ismean=true,clucode=gmclucode)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_clu_mwave$ext")),figfmt)
plotcluwaveform(gmsuwy,gmclu;ismean=false)
foreach(ext->savefig(joinpath(unitfeaturedir,"gmsu_cfeature_clu_wave$ext")),figfmt)


# Add spike type info for cortical units
allcu = DataFrame(siteid = allunit["siteid"][cui],id = unitid[cui],good=unitgood[cui],aligndepth=unitaligndepth[cui],layer=unitlayer[cui])
leftjoin!(allcu,DataFrame(id=unitid[csui],spiketype=map(i->clucode[i],clu)),on=:id)
jldsave(joinpath(resultroot,"allcu.jld2");allcu)

# Add spike type info for GM units
allgmu = DataFrame(siteid = allunit["siteid"][gmui],id = unitid[gmui],good=unitgood[gmui],aligndepth=unitaligndepth[gmui],layer=unitlayer[gmui])
leftjoin!(allgmu,DataFrame(id=unitid[gmsui],spiketype=map(i->gmclucode[i],gmclu)),on=:id)
jldsave(joinpath(resultroot,"allgmu.jld2");allgmu)