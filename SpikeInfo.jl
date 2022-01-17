using NeuroAnalysis,FileIO,Statistics,Plots,ProgressMeter,XLSX

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

## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
@showprogress "Batch Spike Info ... " for r in eachrow(penetration)
    spikeinfo(joinpath(resultroot,r.Subject_ID,r.siteid))
end

## Collect Units
function collectunit!(indir;unit=Dict(),datafile="unit.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            u = load(joinpath(root,datafile),"unit")
            delete!(u,"chposition")
            u["siteid"] = fill(u["siteid"],length(u["unitid"]))
            u["id"] = map((s,g,i)->g ? "$(s)_SU$i" : "$(s)_MU$i",u["siteid"],u["unitgood"],u["unitid"])
            if haskey(unit,"id")
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


f1 = ["duration" "peaktroughratio" "halftroughwidth" "halfpeakwidth" "repolarrate" "recoverrate"]
f2 = ["upspread" "downspread" "leftspread" "rightspread" "uppvinv" "downpvinv"]
Y = Float64[abs.(unitfeature[f1[1]]);;hcat(map(k->unitfeature[k],f1[2:end])...);;
     hcat(map(k->replace(unittempfeature[k],NaN=>0,Inf=>0,-Inf=>0),f2)...)][unitgood,:]
fs = [f1 f2]


using StatsPlots,StatsBase,Interact,UMAP

foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))
violin(fs,Y,leg=false,ylabel="Z Score",xrotation=20,xlabel="Spike Feature")
foreach(ext->savefig(joinpath(indir,"UnitFeature_dist$ext")),figfmt)


cm = cgrad(:fire)

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3, c in 1:length(fs)
    Y2 = umap(Y', 2, n_neighbors=n, min_dist=d, n_epochs=i)
    scatter(Y2[1,:], Y2[2,:],leg=false,color=cm[clampscale(Y[:,c],3)],markerstrokewidth=0)
end
