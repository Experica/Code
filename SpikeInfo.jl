using NeuroAnalysis,FileIO,Statistics,StatsPlots,StatsBase,Images,ProgressMeter,DataFrames,XLSX,Dierckx

function mergespike!(indir;unit=Dict(),datafile="spike.jld2")
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


## Merge Spike Data of a RecordSite
resultroot = "Z:/"
figfmt = [".svg"]
indir = joinpath(resultroot,"AG1","AG1_V1_ODL2")
layer = load(joinpath(indir,"layer.jld2"),"layer")


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


## Penetration Sites
# penetration = unique!(meta[:,[:Subject_ID,:RecordSession,:RecordSite]])
# penetration = transform!(penetration,All()=>ByRow((a,b,c)->join(filter!(!isempty,[a,b,c]),"_"))=>:siteid)
# XLSX.writetable(joinpath(resultroot,"penetration.xlsx"),collect(eachcol(penetration)),names(penetration))

## Batch RecordSites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch Merge Spike ... " for r in eachrow(penetration)
    indir = joinpath(resultroot,r.Subject_ID,r.siteid)
    unit = mergespike!(indir)
    save(joinpath(indir,"unit.jld2"),"unit",unit)
end
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

## Layer info for all units
usi = map(i->findfirst(j->j==i,alllayer.siteid),allunit["siteid"])
allunit["unitaligndepth"] = [alllayer.e2cfun[usi[i]](allunit["unitposition"][i,2]) |> alllayer.tcfun[usi[i]] for i in eachindex(usi)]
allunit["unitlayer"] = assignlayer.(allunit["unitaligndepth"],[layertemplate])
save(joinpath(resultroot,"allunit.jld2"),"allunit",allunit)





## Single Unit Feature
allunit = load(joinpath(resultroot,"allunit.jld2"),"allunit")
unitid = allunit["id"];unitwave = allunit["unitwaveform"];unittempfeature = allunit["unittemplatefeature"]
unitgood=allunit["unitgood"];unitposition=allunit["unitposition"];unitfeature = allunit["unitfeature"]
unitaligndepth=allunit["unitaligndepth"];unitlayer=allunit["unitlayer"];unitqm=allunit["qm"]
fs=allunit["fs"];figfmt = [".png"]

f1 = ["duration","peaktroughratio","lefthalftroughwidth","righthalftroughwidth","halftroughwidth","lefthalfpeakwidth","righthalfpeakwidth","halfpeakwidth","repolarrate","recoverrate","ttrough","amplitude"]
f2 = ["upspread","downspread","leftspread","rightspread","uppvinv","downpvinv"]
fqm= ["fr","pisi","fp"]
f = [f1;f2;fqm]
F = Float64[hcat(map(k->unitfeature[k],f1)...);;hcat(map(k->unittempfeature[k],f2)...);;hcat(map(k->unitqm[k],fqm)...)]
replace!(F,NaN=>0,Inf=>0,-Inf=>0)

cui = 0 .<= unitaligndepth .<= 1 # cortical units
csui = cui .& unitgood # cortical single units
csuF = F[csui,:]

foreach(i->csuF[:,i]=zscore(csuF[:,i]),1:size(csuF,2))
dotplot(permutedims(1:length(f)),csuF,leg=false,grid=false,ylabel="Z Score",xrotation=25,marker=(1, stroke(0)),xticks=1:length(f),xformatter=i->f[Int(i)],
        xlabel="Spike Feature",size=(750,600),ann=(18,15,"n=$(size(csuF,1))"))
foreach(ext->savefig(joinpath(resultroot,"csu_feature$ext")),figfmt)

wn = size(unitwave,2)
wx = collect(range(-wn/2,step=1,length=wn))
csuwy = permutedims(unitwave[csui,:])
csuwx = repeat(wx,outer=(1,size(csuwy,2)))

## Single Unit Classification
using UMAP,Clustering,Distances

ft = ["duration","lefthalftroughwidth","righthalftroughwidth","halftroughwidth",
    "lefthalfpeakwidth","righthalfpeakwidth","halfpeakwidth","repolarrate","recoverrate",
    "uppvinv","downpvinv"]
ft = ["duration","lefthalftroughwidth","righthalftroughwidth","halftroughwidth",
    "lefthalfpeakwidth","righthalfpeakwidth","halfpeakwidth","repolarrate","recoverrate"]
ft = ["duration","halftroughwidth","halfpeakwidth","repolarrate","recoverrate"]
ft = ["duration","halftroughwidth","halfpeakwidth"]
Ft = csuF[:,indexin(ft,f)]

Ft2 = umap(Ft', 2, n_neighbors=25, min_dist=0.8, n_epochs=300,metric=Euclidean())
scatter(Ft2[1,:], Ft2[2,:],leg=false,msw=0,ms=3,ma=0.8,size=(600,600),frame=:none,ratio=1)
plot(Ft2[1,:]' .+ 1e-2*csuwx, Ft2[2,:]' .+ 3e-5*csuwy,color=:darkgreen,alpha=0.8,leg=false,frame=:none,lw=0.5,size=(800,800),ratio=1)
foreach(ext->savefig(joinpath(resultroot,"csu_umap$ext")),figfmt)

FtD = pairwise(Euclidean(),Ft,dims=1)
hc = hclust(FtD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=5)
clu = replace(clu,1=>"En",2=>"E",3=>"I",4=>"Ei",5=>"Ib")

scatter(Ft2[1,:], Ft2[2,:],group=clu,leg=:inline,frame=:none,msw=0,ms=3,ma=0.8,size=(600,600),ratio=1)
foreach(ext->savefig(joinpath(resultroot,"csu_umap_clu$ext")),figfmt)





# groupedhist(Ft[:,7],group=clu,bar_position = :stack)

# dotplot(permutedims(ft),Ft[cid.==1,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:blue),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# # dotplot!(permutedims(ft),Ft[cid.==2,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:orange),
# #         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),Ft[cid.==3,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:green),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),Ft[cid.==4,:],leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(1, stroke(0),:purple),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))



# clim=4;ms=0.1
# dotplot(permutedims(ft),clamp!(Ft[clu.==1,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[1]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),clamp!(Ft[clu.==2,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[2]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),clamp!(Ft[clu.==3,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[3]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),clamp!(Ft[clu.==4,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[4]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),clamp!(Ft[clu.==5,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[5]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# dotplot!(permutedims(ft),clamp!(Ft[clu.==6,:],-clim,clim),leg=false,grid=false,ylabel="Z Score",xrotation=20,marker=(ms, stroke(0),cs[6]),
#         xlabel="Spike Feature",size=(750,600),ann=(15,17,"n=$(size(F,1))"))
# foreach(ext->savefig(joinpath(resultroot,"SpikeFeature_UnitCluster$ext")),figfmt)





plotcluwaveform=(wys,clu;isalign=true,isnorm=true,ismean=true,n=30,color=:default)->begin
    wn,un=size(wys)
    cluid = unique(clu)
    ci = [c.==clu for c in cluid]
    if isalign
        wx = range(-round(Int,wn/2),step=1,length=wn)
        iwx = -wn:0.25:wn
        @views iwys = hcat(map(i->Spline1D(wx,wys[:,i],k=3,bc="nearest")(iwx),1:un)...)
        r = range(4wx[begin],4wx[end],step=1)
        @views ws = hcat(map(i->iwys[r .+ argmax(abs.(iwys[:,i])),i],1:un)...)
    else
        ws = wys
    end
    cws = map(i->ws[:,i],ci)
    isnorm && map!(i->i./maximum(abs.(i),dims=1),cws,cws)
    if ismean
        cwm = hcat(map(i->mean(i,dims=2),cws)...)
        cwse = hcat(map(i->std(i,dims=2)/sqrt(size(i,2)),cws)...)
        plot(cwm;ribbon=cwse,label=permutedims(cluid),palette=color,lw=2)
    else
        colors = palette(color,length(cluid)).colors.colors
        @views csw = map(i->i[:,sample(1:size(i,2),min(n,size(i,2)),replace=false)],cws)
        p=plot(leg=false)
        for i in eachindex(csw)
            plot!(p,csw[i],color=colors[i],alpha=0.8)
        end
        p
    end
end

plotcluwaveform(csuwy,clu,ismean=true)
foreach(ext->savefig(joinpath(resultroot,"csu_clu_spikewave$ext")),figfmt)


allcu = DataFrame(siteid = allunit["siteid"][cui],id = unitid[cui],good=unitgood[cui],aligndepth=unitaligndepth[cui],layer=unitlayer[cui])
leftjoin!(allcu,DataFrame(id=unitid[csui],type=clu),on=:id)
save(joinpath(resultroot,"allcu.jld2"),"allcu",allcu)

allcsu = subset(allcu,:good)
plotdepthhist(allcsu;g=:type,layer)