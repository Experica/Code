using Query,VegaLite,DataFrames,Plots

## Load Cells
resultroot = "../Result"
figfmt = [".png",".svg"]
cells = load(joinpath(resultroot,"cells.jld2"),"cells")

dklcg = cgrad(ColorMaps["dkl_mcchue_l0"].colors)
hslcg = cgrad(ColorMaps["hsl_mshue_l0.4"].colors)
dklcs = "#" .* hex.(dklcg.colors.colors,:rgb)
hslcs = "#" .* hex.(hslcg.colors.colors,:rgb)

## plot
dklcells = select!(dropmissing(cells,[:dkl_frf]),Not([:dkl_frf,:hsl_frf]),
                    :dkl_frf=>ByRow(x->x.oh)=>:dkl_oh,:dkl_frf=>ByRow(x->x.hcv)=>:dkl_hcv,
                    :dkl_frf=>ByRow(x->x.fit.ph)=>:dkl_ph,:dkl_frf=>ByRow(x->x.fit.hhw)=>:dkl_hhw,
                    :dkl_frf=>ByRow(x->x.fit.gvm.r)=>:dkl_gvm_r)
hslcells = select!(dropmissing(cells,[:hsl_frf]),Not([:dkl_frf,:hsl_frf]),
                    :hsl_frf=>ByRow(x->x.oh)=>:hsl_oh,:hsl_frf=>ByRow(x->x.hcv)=>:hsl_hcv,
                    :hsl_frf=>ByRow(x->x.fit.ph)=>:hsl_ph,:hsl_frf=>ByRow(x->x.fit.hhw)=>:hsl_hhw,
                    :hsl_frf=>ByRow(x->x.fit.gvm.r)=>:hsl_gvm_r)


dklcells |> @vlplot(:bar,x={"dkl_oh",bin={maxbins=20}},y={"count()"},row={"site"})
hslcells |> @vlplot(:bar,x={"hsl_oh",bin={maxbins=20}},y={"count()"},row={"site"})


plothuehist=(df;w=700,h=700,hue="oh",cs="dkl",cg=dklcg)->begin
    sites = sort(unique(df.site));n=length(sites);cshue="$(cs)_$hue"
    p=plot(layout=(n,1),legend=false,grid=true,size=(w,n*h))
    for i in 1:n
        ys,ns,ws,is = epochspiketrain(df[df.site.==sites[i],cshue],0:30:360)
        if !isnothing(cg)
            x = 2*π*cg.values
            y = fill(1.1maximum(ns),length(x))
            scatter!(p,subplot=i,x,y,projection=:polar,color=cg.colors.colors,markerstrokewidth=0,markersize=12)
        end
        plot!(p,subplot=i,deg2rad.(mean.([ws;ws[1]])),[ns;ns[1]],title="$(sites[i])--$cshue",projection=:polar,linewidth=5,linecolor=:gray30)
    end
    p
end

plothuehist(dklcells,hue="oh",cs="dkl",cg=dklcg)
foreach(ext->savefig(joinpath(resultroot,"DKLOptHueHist$ext")),figfmt)
plothuehist(hslcells,hue="oh",cs="hsl",cg=hslcg)
foreach(ext->savefig(joinpath(resultroot,"HSLOptHueHist$ext")),figfmt)


dklcells |> @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"},
                        color={"dkl_oh",scale={range=dklcs}},column={"site"})

dklcells |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"dkl_oh",scale={range=dklcs}},column={"site"})

hslcells |> @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"},
                        color={"hsl_oh",scale={range=hslcs}},column={"site"})

hslcells |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"hsl_oh",scale={range=hslcs}},column={"site"})


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
foreach(ext->savefig(joinpath(siteresultdir,"OptHue_UnitPosition$ext")),figfmt)
