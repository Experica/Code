using NeuroAnalysis,Statistics,FileIO,Plots,Images,Interact

# Combine all color tests of one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL1"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)
figfmt = [".png",".svg"]

layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")

testids = ["$(siteid)_Color_$i" for i in 1:2]
testn=length(testids)

## All factor response
frs = load.(joinpath.(resultsitedir,testids,"factorresponse.jld2"))
testlogs = ["$(i["log"])_$(i["color"])" for i in frs]
cms = [occursin("DKL",i) ? cgrad(ColorMaps["dkl_mcchue_l0"].colors) : cgrad(ColorMaps["hsl_mshue_l0.4"].colors)  for i in testlogs]

vus = map(i->i["responsive"].&i["modulative"].&i["unitgood"].&i["enoughresponse"],frs)
ohs = map((i,v)->map(j->j.oh,i["factorresponsefeature"][:HueAngle][v]),frs,vus)
hcvs = map((i,v)->map(j->j.hcv,i["factorresponsefeature"][:HueAngle][v]),frs,vus)

## All spike
spks = load.(joinpath.(resultsitedir,testids,"spike.jld2"),"spike")
upos = map((s,v)->s["unitposition"][v,:],spks,vus)

## plot
plothuehist=(;w=700,h=700)->begin
    p=plot(layout=(1,testn),link=:none,legend=false,grid=true,size=(testn*w,h))
    for i in 1:testn
        ys,ns,ws,is = epochspiketrain(ohs[i],0:30:360)
        x = 2*π*cms[i].values
        y = fill(1.1maximum(ns),length(x))
        scatter!(p,subplot=i,x,y,projection=:polar,color=cms[i].colors.colors,markerstrokewidth=0,markersize=12)
        plot!(p,subplot=i,deg2rad.(mean.([ws;ws[1]])),[ns;ns[1]],title=testlogs[i],projection=:polar,linewidth=5,linecolor=:gray30)
    end
    p
end

plothuehist()
foreach(ext->savefig(joinpath(resultsitedir,"OptHueHist$ext")),figfmt)

plothueposition=(;w=800,h=700)->begin
    p=plot(layout=(1,testn),legend=false,grid=false,size=(testn*w,h))
    for i in 1:testn
        xlims = extrema(upos[i][:,1]).+[-2,1]
        if !isnothing(layer)
            lx = xlims[1]+1
            hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,
            annotations=[(lx,layer[k][1],text(k,7,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray70,legend=false)
        end
        scatter!(p,subplot=i,upos[i][:,1],upos[i][:,2],title=testlogs[i],markersize=(1 .-hcvs[i])*5 .+2, color=cms[i][ohs[i]/360],
        xlims=xlims,markerstrokewidth=0)
    end
    p
end

plothueposition()
foreach(ext->savefig(joinpath(resultsitedir,"OptHue_UnitPosition$ext")),figfmt)
