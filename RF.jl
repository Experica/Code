using NeuroAnalysis,Statistics,StatsBase,FileIO,JLD2,Images,StatsPlots,Interact,Distances,Combinatorics,
    VegaLite,DataFrames,Clustering,UMAP,MultivariateStats,HypothesisTests


function collectsta(indir;unit=DataFrame(),datafile="stadataset.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            # STA was batched only for SU
            # SU was chosen that spikes in all tests
            # lroi only for STA responsive in at least one tests
            id = ["$(siteid)_SU$u" for u in keys(dataset["ulroi"])]
            roic = [v.centerdeg for v in values(dataset["ulroi"])]
            roir = [v.radiusdeg for v in values(dataset["ulroi"])]
            cr = Any[map((c,r)->r ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) for u in keys(dataset["ulroi"])]
            df = DataFrame(siteid=siteid,id=id,roicenter=roic,roiradius=roir,cr=cr)
            if haskey(dataset,"ulfit")
                foreach(m->df[!,"sta!$m"]=Any[dataset["ulfit"][u][m] for u in keys(dataset["ulroi"])], keys(first(values(dataset["ulfit"]))))
            end
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

## STA
resultroot = "Z:/"
figfmt = [".png"]

staunit = collectsta(resultroot)


sym(r) = r < 1 ? 1 - 1/r : r - 1
abssym(r) = abs(sym(r))
signsym(r) = sign(sym(r))
oddeven(p) = p < 0.5 ? 4abs(p-0.25) : 4abs(p-0.75)
onoff(p) = p < 0.5 ? 1 - 4abs(p-0.25) : 4abs(p-0.75) - 1
dogfun = (x,y,p) -> dogf.(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[6],μᵢ₂=p[4],σᵢ₂=p[6],θᵢ=0)
gaborfun = (x,y,p) -> gaborf.(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5],θ=p[6],f=p[7],phase=p[8])

function fitpredict(model,param;ppd=45,tight=false)
    if model==:dog
        fit=(;model,param,fun=dogfun)
    elseif model==:gabor
        fit=(;model,param,fun=gaborfun)
    end

   fitpredict(fit;ppd,tight)
end
function fitpredict(fit;ppd=45,tight=false)
   if tight
       param = deepcopy(fit.param)
       if fit.model == :dog
           param[[2,4]].=0
           radius = 3*maximum(param[[6,3]])
           fit = (;fit.model,fit.fun,param,radius)
       elseif fit.model == :gabor
           param[[2,4]].=0
           radius = 3*maximum(param[[3,5]])
           fit = (;fit.model,fit.fun,param,radius)
       end
   end
   x=y= -fit.radius:1/ppd:fit.radius
   predict(fit,x,y,yflip=true)
end
function fitnorm(fit;sdfactor=3)
    param = deepcopy(fit.param)
    if fit.model == :dog
        param[[2,4]].=0
        radius = sdfactor*maximum(param[[6,3]])
    elseif fit.model == :gabor
        param[[2,4]].=0
        radius = sdfactor*maximum(param[[3,5]])
    end
    sqrt(hcubature(x->fit.fun(x...,param)^2,[-radius,-radius],[radius,radius])[1])
end

function stadogunit(unit;model=:dog,imgsize=(48,48))
    df = select!(dropmissing!(flatten(dropmissing(unit,["cr","sta!$model"]),["cr","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->i.model)=>:model,
        "sta!$model"=>ByRow(i->i.param)=>:param,
        # "sta!$model"=>ByRow(i->imresize(fitpredict(i,tight=true),imgsize))=>"img","sta!$model"=>ByRow(i->fitnorm(i))=>:W,
        "sta!$model"=>ByRow(i->i.r)=>:r,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->i.aic)=>:aic,
        "sta!$model"=>ByRow(i->i.adjr2)=>:adjr2,
        "sta!$model"=>ByRow(i->i.param[1])=>:ae,
        "sta!$model"=>ByRow(i->i.param[2])=>:x,
        "sta!$model"=>ByRow(i->i.param[3])=>:se,
        "sta!$model"=>ByRow(i->i.param[4])=>:y,
        "sta!$model"=>ByRow(i->i.param[5])=>:ai,
        "sta!$model"=>ByRow(i->i.param[6])=>:si,
        "sta!$model"=>ByRow(i->4*maximum(i.param[[3,6]]))=>:diameter,
        "sta!$model"=>ByRow(i->sym(i.param[6]/i.param[3]))=>:cs,
        "sta!$model"=>ByRow(i->abs(sym(i.param[6]/i.param[3])))=>:op,
        "sta!$model"=>ByRow(i->sym(i.param[1]/i.param[5]))=>:onoff,
        "sta!$model"=>ByRow(i->abs(sym(i.param[1]/i.param[5])))=>:amp,
        "sta!$model"=>ByRow(i->signsym(i.param[1]/i.param[5]))=>:sign)
end

function stagaborunit(unit;model=:gabor,imgsize=(48,48))
    df = select!(dropmissing!(flatten(dropmissing(unit,["cr","sta!$model"]),["cr","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->i.model)=>:model,
        "sta!$model"=>ByRow(i->i.param)=>:param,
        # "sta!$model"=>ByRow(i->imresize(fitpredict(i,tight=true),imgsize))=>"img","sta!$model"=>ByRow(i->fitnorm(i))=>:W,
        "sta!$model"=>ByRow(i->i.r)=>:r,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->i.aic)=>:aic,
        "sta!$model"=>ByRow(i->i.adjr2)=>:adjr2,
        "sta!$model"=>ByRow(i->i.param[1])=>:a,
        "sta!$model"=>ByRow(i->i.param[2])=>:x,
        "sta!$model"=>ByRow(i->i.param[3])=>:sx,
        "sta!$model"=>ByRow(i->i.param[4])=>:y,
        "sta!$model"=>ByRow(i->i.param[5])=>:sy,
        "sta!$model"=>ByRow(i->rad2deg(i.param[6]))=>:ori,
        "sta!$model"=>ByRow(i->i.param[7])=>:sf,
        "sta!$model"=>ByRow(i->i.param[8])=>:phase,
        "sta!$model"=>ByRow(i->4*maximum(i.param[[3,5]]))=>:diameter,
        "sta!$model"=>ByRow(i->sym(i.param[3]/i.param[5]))=>:el,
        "sta!$model"=>ByRow(i->abs(sym(i.param[3]/i.param[5])))=>:rd,
        "sta!$model"=>ByRow(i->4*i.param[5]*i.param[7])=>:cyc,
        "sta!$model"=>ByRow(i->oddeven(i.param[8]))=>:oddeven,
        "sta!$model"=>ByRow(i->onoff(i.param[8]))=>:onoff,
        "sta!$model"=>ByRow(i->sign(onoff(i.param[8])))=>:sign)
end

dogunit = stadogunit(staunit)
gaborunit = stagaborunit(staunit)

gaborunit |> Voyager()






# Test goodness of fit
gdogi = dogunit.r .> 0.8
ggabori = gaborunit.r .> 0.7



# BIC to select model
scatter(dogunit.bic,gaborunit.bic,leg=false,marker=(3,0.5,stroke(0)),ratio=:equal,grid=false,xlims=(0,2e5),
    ylims=(0,2e5),xticks=[],yticks=[],xlabel="BIC (DoG)",ylabel="BIC (Gabor)")
plot!(x->x,0,2e5,color=:black)

dogi = dogunit.bic .< gaborunit.bic

# Spatial Pattern
gdogunit = dogunit[dogi .& gdogi,:]
ggaborunit = gaborunit[.!dogi .& ggabori,:]




imgsize=(48,48)
gdogunit.img = imresize.(fitpredict.(gdogunit.model,gdogunit.param,tight=true),[imgsize])
gdogimg = map(i->float32.(alphamask(i)),gdogunit.img)
gdogimg = map(i->float32.(i),gdogunit.img)
scatter(gdogunit.cs,gdogunit.onoff,yflip=true,xlims=(-0.4,0.4),ylims=(-1.5,1.5),leg=false,xticks=[0],yticks=[0],
    marker=(0.8,stroke(0)),xlabel="CenterSurround",ylabel="OnOff")

scatter(ggaborunit.phase,ggaborunit.cyc,leg=false,yflip=true,ylims=(0,2),
    marker=(3,0.8,stroke(0)),xlabel="Phase",ylabel="Cycle")

# Spatial Pattern Clustering
features = [:cs,:onoff,:op,:amp]
F = mapfoldl(i->gdogunit[:,i],hcat,features)
foreach(i->F[:,i]=zscore(F[:,i]),1:size(F,2))

cr = kmeans(F',4)
cid = assignments(cr)


using GLMakie



@manipulate for i in 100:300, n in 5:50, d in 0.001:0.001:3
    F2 = umap(F', 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Plots.scatter(F2[2,:], F2[1,:],leg=false,grid=false,marker=(3,stroke(0)),size=(600,600),
            xticks=[],yticks=[],ratio=:equal,xlabel="UMAP Dimension 2",ylabel="UMAP Dimension 1")
    # scatter(Ft2[2,:], Ft2[1,:],group=cid,leg=true,grid=false,marker=(1,stroke(0)),size=(600,600),palette=:tab10,
    #         xticks=[],yticks=[],ratio=:equal,xlabel="UMAP Dimension 2",ylabel="UMAP Dimension 1")
end



@manipulate for i in 100:300, n in 5:50, d in 0.001:0.001:3
   F2 = umap(F', 2, n_neighbors=n, min_dist=d, n_epochs=i)
   Makie.scatter(F2[2,:], F2[1,:], marker=gdogimg,markersize=30)
end






scatter(gdogunit.cs,gdogunit.onoff,yflip=true,xlims=(-0.5,0.5),ylims=(-2.5,2.5),leg=true,xticks=[-0.2,0,0.2],#group=cid,
    marker=(0.8,stroke(0)),palette=:tab10,xlabel="CenterSurround",ylabel="OnOff")

scatter(gdogunit.cs,gdogunit.onoff,yflip=true,leg=true,xticks=[0],yticks=[0],group=cid,
    marker=(0.8,stroke(0)),palette=:tab10,xlabel="CenterSurround",ylabel="OnOff")








vdogcell |> [@vlplot(:bar,x={"r",bin={step=0.05,extent=[0,1]},title="r"},y={"count()",title="Number of Linear Filter"});
                @vlplot(mark={:line,size=2},x=:r,transform=[{sort=[{field=:r}],window=[{op="count",as="cum"}],frame=[nothing,0]}],
                y={"cum",title="Cumulated Number of Linear Filter"})]

doggaborgoodness = DataFrame(dogr=dogcell.r,gaborr=gaborcell.r,dogbic=dogcell.bic,gaborbic=gaborcell.bic,dogfvu=dogcell.fvu,gaborfvu=gaborcell.fvu)
doggaborgoodness |>
    [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.5e5}]},x=:x,y=:x}, {:point,x={:gaborbic,title="BIC (Gabor)"},y={:dogbic,title="BIC (DoG)"}}]);
    @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:gaborfvu,title="FVU (Gabor)"},y={:dogfvu,title="FVU (DoG)"}}]);
    @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:gaborr,title="R (Gabor)"},y={:dogr,title="R (DoG)"}}])]



srfcells |> [@vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Spatial RF"},color="rc");
       @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Spatial RF"},color="rc")]

srfcells |> [@vlplot(:point,x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
          @vlplot(:point,x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
          @vlplot(:point,x={"odd",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"})]

gaborcell |> [@vlplot(mark={:point,size=20},x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
       @vlplot(mark={:point,size=20},x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
       @vlplot(mark={:point,size=20},x={"oddeven",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"})]

gaborcell |> [@vlplot(:point,x={"diameter",title="Diameter (Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"el",title="Ellipticalness"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"ori",title="Ori (Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"sf",title="SpatialFreq (Cycle/Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"phase",title="Phase"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"oddeven",title="OddEven"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"onoff",title="OnOff"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}})]

dogunit |> [@vlplot(mark={:bar,binSpacing=0},x={"diameter",bin={step=0.03},title="Diameter (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ae",bin={maxbin=30},title="Ae"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ai",bin={maxbin=30},title="Ai"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"centersurround",bin={step=0.03},title="CenterSurround"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"opponency",bin={step=0.03},title="Opponency"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"onoff",bin={step=0.2},title="OnOff"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"amplitude",bin={step=0.2},title="OnOff Amplitude"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"sign",bin={maxbin=3},title="OnOff Sign"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}})]

gaborcell |> [@vlplot(mark={:bar,binSpacing=0},x={"diameter",bin={step=0.03},title="Diameter (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"el",bin={step=0.1},title="Ellipticalness"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"rd",bin={step=0.1},title="Roundness"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ori",bin={step=3.6},title="Ori (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"sf",bin={step=0.1},title="SpatialFreq (Cycle/Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"cyc",bin={step=0.1},title="Cycle"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"phase",bin={step=0.02},title="Phase"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"oddeven",bin={step=0.02},title="OddEven"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"onoff",bin={step=0.04},title="OnOff"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"sign",bin={maxbin=3},title="OnOff Sign"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}})]

gaborcell |> [@vlplot(:bar,x={"ar",bin={maxbins=20},title="Aspect Ratio"},y={"count()",title="Number of Spatial RF"},color={"rc"});
       @vlplot(:bar,x={"nc",bin={maxbins=20},title="Number of Cycle"},y={"count()",title="Number of Spatial RF"},color={"rc"});
       @vlplot(:bar,x={"odd",bin={maxbins=15},title="Oddness"},y={"count()",title="Number of Spatial RF"},color={"rc"})]








## UMAP of STA
function d2clu(x;r=1)
   cs = dbscan(x,r)
   cid = zeros(Int,size(x,2))
   foreach(i->cid[cs[i].core_indices] .= i, 1:length(cs))
   cid
end
function srimg(c,i)
   mapreduce((p,q)->float32.(alphamask(q,color=spectralcm[p])),hcat,c,i)
end
dogimg = map(i->float32.(alphamask(i)),dogcell.img)
gaborimg = map(i->float32.(alphamask(i)),gaborcell.img)
spectralcm = Dict('A'=>ColorMaps["dkl_mcclumiso"].colors,'L'=>ColorMaps["lms_mccliso"].colors,
                       'M'=>ColorMaps["lms_mccmiso"].colors,'S'=>ColorMaps["lms_mccsiso"].colors)
dogimgc = map((i,c)->float32.(alphamask(i,color=spectralcm[c])),dogcell.img,dogcell.rc)
gaborimgc = map((i,c)->float32.(alphamask(i,color=spectralcm[c])),gaborcell.img,gaborcell.rc)

@manipulate for i in 1:length(dogimg)
   dogimg[i]
end

@manipulate for i in 1:length(gaborimg)
   gaborimg[i]
end

# DoG
features = [:cs,:op,:amp]
Y = mapfoldl(i->dogcell[:,i],hcat,features)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
   Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
   Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(permutedims(Y), 2, n_neighbors=30, min_dist=0.4, n_epochs=300,metric=Euclidean(),init=:spectral)
p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=2.5),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_dog_shape_umap_clu.png"))
dogcell.umapclu = d2clu(Y2,r=2.5)



save("UMAP_sta.svg",plotunitpositionimage(Y2',dogimg))


dogcell |> [@vlplot(mark={:bar,binSpacing=0},x={"op",bin={step=0.03},title="Opponency"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"});
      @vlplot(mark={:bar,binSpacing=0},x={"amp",bin={step=0.2},title="OnOff Amplitude"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"})]
dogcell.sop = dogcell.umapclu .< 3
dogcell.spatial = ["dog-"*(i ? "op" : "nop") for i in dogcell.sop]



# Gabor
features = [:rd,:cyc,:oddeven]
Y = mapfoldl(i->gaborcell[:,i],hcat,features)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
   Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
   Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(permutedims(Y), 2, n_neighbors=16, min_dist=0.05, n_epochs=300,metric=Euclidean(),init=:spectral)
p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=1.2),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_gabor_shape_umap_clu.png"))
gaborcell.umapclu = d2clu(Y2,r=1.2)

save("UMAP_sta.svg",plotunitpositionimage(Y2',gaborimg))


gaborcell |> [@vlplot(mark={:bar,binSpacing=0},x={"cyc",bin={step=0.1},title="Cycle"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"});
      @vlplot(mark={:bar,binSpacing=0},x={"oddeven",bin={step=0.02},title="OddEven"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"})]
gaborcell.sop = gaborcell.umapclu .< 3
gaborcell.spatial = ["gabor-"*(i == 3 ? "nop" : i == 2 ? "op_odd" : "op_even") for i in gaborcell.umapclu]



save(joinpath(resultroot,"dogcell.jld2"),"dogcell",dogcell)
save(joinpath(resultroot,"gaborcell.jld2"),"gaborcell",gaborcell)










## t-SNE of responsive sta
lsta = mapreduce(u->mapreduce((c,r,exd)->r ? [dataset["ulsta"][u][:,:,exd.d,c]] : [],
append!,1:testn,dataset["ucresponsive"][u],dataset["ulcexd"][u]),
append!,keys(dataset["ulsta"]))

maxsizepx=mapreduce(i->maximum(size(i)),max,lsta)
map!(i->imresize(i,maxsizepx,maxsizepx),lsta,lsta)

@manipulate for i in eachindex(lsta)
    heatmap(lsta[i],yflip=true,aspect_ratio=:equal,frame=:none,color=:coolwarm)
end

# PCA reduction
X = mapreduce(vec,hcat,lsta)
foreach(i->X[i,:]=zscore(X[i,:]),1:size(X,1))
M = fit(PCA,X,maxoutdim=60,pratio=0.85,mean=0)
Y = MultivariateStats.transform(M,X)

@manipulate for i in 1000:6000, p in 2:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

# plot t-SNE sta
mlsta = map(i->alphamask(get(cgrad(:coolwarm),i,:extrema),radius=0.5,sigma=0.35,masktype="Gaussian")[1],lsta)

Y2 = tsne(permutedims(Y), 2, 0, 5000, 8)
p = plotunitpositionimage(Y2,mlsta)
save("tSNE_rlSTA.svg",p)

## t-SNE of responsive rf
m=:gabor
rfs = mapreduce(u->mapreduce(f->f.r > 0.6 ? [(;(n=>f[n] for n in (:model,:radius,:param))...)] : [],
append!,skipmissing(dataset["urf"][u][m])),
append!,keys(dataset["urf"]))

X = mapreduce(f->f.param,hcat,rfs)
Y = X
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 1000:6000, p in 2:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

mrf = map(f->alphamask(get(cgrad(:coolwarm),rfimage(f,f.radius),:extrema),radius=0.5,sigma=0.35,masktype="Gaussian")[1],rfs)

Y2 = tsne(permutedims(Y), 2, 0, 5000, 10)
p = plotunitpositionimage(Y2,mrf)
save("tSNE_rf_$m.svg",p)

## type of responsiveness
CTYPES = sort(join.(sort.(combinations(unique(values(ccode))))))
CTCOLORS = cgrad(:rainbow,length(CTYPES),categorical = true).colors
crtype = Dict(u=>join(sort(dataset["ccode"][dataset["ucresponsive"][u]])) for u in keys(dataset["ulsta"]))

p = DataFrame(ResponseType=collect(values(crtype))) |> @vlplot(:bar,y=:ResponseType,x="count()",width=600,height=400)

## type of model responsiveness
m=:gabor
cmrtype = filter!(p->!isempty(p.second),Dict(u=>join(sort(dataset["ccode"][map(f->!ismissing(f) && f.r > 0.6,dataset["urf"][u][m])])) for u in keys(dataset["urf"])))

p = DataFrame(ModelResponseType=collect(values(cmrtype))) |> @vlplot(:bar,y=:ModelResponseType,x="count()",width=600,height=400)

## rf space of unit
modelnull=Dict(:gabor=>fill(missing,8),:dog=>[])
function rf2unitspace(ps;model=:gabor)
    if model == :gabor
        as = ps[1:8:end]
        maxamp = maximum(skipmissing(as))
        foreach(i->ps[i]/=maxamp,1:8:length(ps))
        ps=replace(ps,missing=>0)
    elseif model == :dog
    end
    ps
end
"map RF models to unit feature space"
function rfunitspace(dataset;model=:gabor)
    corder = sortperm(dataset["ccode"])
    ufs=Dict()
    for u in keys(dataset["urf"])
        ps = vcat(map(f->ismissing(f) || f.r < 0.6 ? modelnull[model] : f.param,dataset["urf"][u][model])[corder]...)
        if !isempty(skipmissing(ps))
            ufs[u] = rf2unitspace(ps,model=model)
        end
    end
    ufs
end


rfspace = rfunitspace(dataset,model=:gabor)
Xg = [cmrtype[u] for u in keys(rfspace)]
X = hcat(values(rfspace)...)
Y = X
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 1000:6000, p in 2:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:circle,stroke(0)),group=Xg,frame=:none,leg=:inline)
    # savefig("test.svg")
end







XD=pairwise(Euclidean(),X,dims=2)
hc = hclust(XD,branchorder=:optimal)
p=plot(hc,xticks=false)

uids = collect(keys(rfspace))
cids = cutree(hc,h=2000)




scatter(X2D[:,1], X2D[:,2], marker=(6,6,:auto,stroke(0)),marker_z=cids,
color=cgrad(:seaborn_bright,length(unique(cids)),categorical=true),
series_annotations=text.(uids,5,:gray10,:center))




n=length(urfs)
rs = [kmeans(X,k,display=:final,maxiter=2000,init=:kmpp) for k in 1:15]
XD = Array{Float64}(undef,n,n)


## rf space of unit










crti = map(t->findfirst(t.==crtypes),values(ucrt))

spike = load(joinpath(resultsitedir,testids[1],"spike.jld2"),"spike")
unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

ui = indexin(keys(ucrt),unitid)
up=unitposition[ui,:]


crtcs[2.0,3.0]

plotunitposition(up,color=crtcs[crti])








## t

m = :dog
ucrf = Dict()
for u in keys(dataset["urf"])
    umfit = dataset["urf"][u][m]
    for i in eachindex(umfit)
        isnothing(umfit[i]) && continue
        k = (uid=u,d=dataset["ulcexd"][u][i].d,c=i)
        ucrf[k] = umfit[i].r
    end
end

gaborr = collect(values(ucrf))
dogr = collect(values(ucrf))

rs = [kmeans([gaborr dogr]',k,display=:final,maxiter=1000,init=:rand) for k in 1:10]

tloss = map(r->r.totalcost,rs)
vinfo = [Clustering.varinfo(rs[i],rs[i+1]) for i in 1:length(rs)-1]
minfo = [Clustering.mutualinfo(rs[i],rs[i+1]) for i in 1:length(rs)-1]
vm = [Clustering.vmeasure(rs[i],rs[i+1],β=2) for i in 1:length(rs)-1]

scatter(gaborr,dogr,xlims=[-0.5, 1],ylims=[-0.5,1],aspect_ratio=:equal,marker_z=assignments(r),color=:rainbow,
markerstrokewidth=0,markersize=10,leg=false,xlabel="Gabor Correlation",ylabel="Dog Correlation")



plot(tloss[6:end])
plot(vinfo)
plot(minfo)
plot(vm)






m = :gabor
ucfittedgabor = Dict()
for u in keys(dataset["urf"])
    umfit = dataset["urf"][u][m]
    for i in eachindex(umfit)
        ucfittedgabor[u]=map(f->isnothing(f) || f.r<0.6 ? false : true,umfit)
    end
end

ucfitted = Dict(u=>ucfitteddog[u] .| ucfittedgabor[u]  for u in keys(ucfitteddog))




fittedtype = hcat(values(ucfitted)...)
fittedtypecs = [kmeans(fittedtype,k,display=:final,maxiter=1000,init=:rand) for k in 1:10]
fminfo = [Clustering.mutualinfo(fittedtypecs[i],fittedtypecs[i+1]) for i in 1:length(fittedtypecs)-1]
plot(fminfo)









##
function getbestrf(dataset;rt=0.65)
    urf = dataset["urf"]
    nc = length(dataset["color"])
    ubrf=Dict()
    for u in keys(urf)
        # u!=96 && continue
        mfs=[]
        for m in keys(urf[u])
            push!(mfs, map(f->begin
             if isnothing(f)
                 nothing
             else
                 (param=f.param,r=f.r,m=m)
             end
         end,urf[u][m]))
         end
         # display(mfs)
         crf = []
         for i in 1:nc
             push!(crf,mapreduce(mf->mf[i],(m1,m2)->begin
             if isnothing(m1)
                 if isnothing(m2)
                     nothing
                 else
                     m2.r > rt ? m2 : nothing
                 end
             else
                 if isnothing(m2)
                     m1.r > rt ? m1 : nothing
                 else
                     t = m1.r > m2.r ? m1 : m2
                     t.r > rt ? t : nothing
                 end
             end
         end,mfs))
         end
         ubrf[u]=crf
    end
    return ubrf
end


usrf = getbestrf(dataset)


urf = map(i->i[:gabor][2],values(dataset["urf"]))
urft = map(i->i.converged,urf)
vurf = urf[urft]

vups = mapreduce(i->i.param,hcat,vurf)

using Clustering

r = kmeans(vups[[3,5],:],4,display=:final)

r = kmeans([abs.(vups[3:3,:]./vups[5:5,:].-1);maximum(vups[[3,5],:],dims=1).*vups[7:7,:]],3,display=:final)

r = kmeans(maximum(vups[[3,5],:],dims=1).*vups[7:7,:],2,display=:final)


ci = assignments(r)


cu = map(i->us[findall(i.==ci)],1:3)


sort(cu[1])




Plots.scatter(vups[3,:],vups[7,:],group=ci)






spike = load(joinpath(resultsitedir,testids[1],"spike.jld2"),"spike")
unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

us = collect(keys(dataset["ulsta"]))
up = unitposition[indexin(us,unitid),:]
ul = unitlayer[indexin(us,unitid)]



uc = [argmax(map(j->abs(j.ex),dataset["ucex"][u])) for u in keys(dataset["ulsta"])]
ud = map((u,c)->dataset["ucex"][u][c].d,keys(dataset["ulsta"]),uc)

ucsta = map((u,d,c)->mapcolor(dataset["ulsta"][u][:,:,d,c],dataset["minmaxcg"][c]),keys(dataset["ulsta"]),ud,uc)
ucmsta = map(i->alphamask(i,radius=0.5,sigma=10,masktype="Disk")[1],ucsta)


p=plotunitpositionimage(up,ucmsta,layer=layer)
save("UnitPosition_STA.svg",p)

p=plotunitlayerimage(ul,ucmsta,unitid=us)
save("UnitLayer_STA.svg",p)


uexs = mapreduce(u->map(j->j.ex,dataset["ucex"][u]),hcat,us)
uds = mapreduce(u->map(j->j.d,dataset["ucex"][u]),hcat,us)

function rftype(cws,cpd)
    a=cws[1];l=cws[2];m=cws[3];s=cws[4];c=1;d = cpd[argmax(abs.(exs))]
    t = join(filter(k->!isnothing(k),map((i,j)->i==0 ? nothing : "$j$(i > 0 ? "+" : "-")",cws,["A","L","M","S"])),"_")
    if l != 0
        if m != 0 # both significent l,m
            if l > 0 && m < 0
                t = "L+M-"
                c = 2
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l < 0 && m > 0
                t = "M+L-"
                c = 3
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l > 0 && m > 0
                if s < 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S-L+M+"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            elseif l < 0 && m < 0
                if s > 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S+L-M-"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            end
        else
            if a ==0
                t = "L+"
                c = 4
                d = cpd[4]
            else
            end
        end
    end

        if l>0 && m<0
            t = "L+M-"
            c = 2
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l<0 && m>0
            t="M+L-"
            c=3
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l>0 && m>0
            if s<0
                t="S-L+M+"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "+-"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        elseif l<0 && m<0
            if s>0
                t="S+L-M-"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "-+"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        end
    return (t=t,c=c,d=d)
end

tcd = [rftype(uexs[:,i],uds[:,i]) for i in 1:size(uexs,2)]
uc = map(i->i.c,tcd)
ud = map(i->i.d,tcd)

_,nsg,_,_=histrv(unitposition[unitgood,2],0,3500,binwidth=60)
_,ns,_,_=histrv(up[:,2],0,3500,binwidth=60)

bar(0:60:3500,[nsg ns],orientation=:h,bar_position=:stack,xlims=[-0.3,10])

plot([nsg ns],range(0,step=40,length=length(ns)),fillalpha=0.4,fill=true,labels=["SU" "SU with STA"],linealpha=0,xlims=[-0.3,7])
hline!([layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(-0.3,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

savefig("unit_layer_dist.svg")







t = dataset["ulsta"][92][:,:,10,2]

Plots.heatmap(t)

ft=abs.(fft(t))
Plots.heatmap(ft)

freqs = fftfreq(25,30)





lmpi = [Tuple(argmax(clc))...]

lmi = clc[:,:,lmpi[3:end]...]

heatmap(lmi,aspect_ratio=:equal,frame=:all,color=:coolwarm,yflip=true)
scatter!([lmpi[2]],[lmpi[1]],markersize=10)
seg = seeded_region_growing(lmi,[(CartesianIndex(lmpi[1:2]...),2),(CartesianIndex(1,1),1)])

lr = labels_map(seg)
heatmap(lr,aspect_ratio=:equal,frame=:all,color=:grays)

ri = mapreduce(i->[i...],hcat, Tuple.(findall(lr.>1)))


t = extrema(ri,dims=2)
tr = map(i->i[2]-i[1],t)
rsize = floor(Int,maximum(tr)/2)

lp = floor.(Int,mean.(t))[:]

hspan!([lp[1]-rsize,lp[1]+rsize],alpha=0.2)
vspan!([lp[2]-rsize,lp[2]+rsize],alpha=0.2)


















@manipulate for u in uids
t1 = stas[1]["usta"][u]
t2=stas[2]["usta"][u]
t3=stas[3]["usta"][u]
t4=stas[4]["usta"][u]

sds = Array{Float64}(undef,imagesize...,41)
for d in 1:41
t = t1[d,:]
tt = fill(mean(t),imagesize)
tt[xi]=t
sds[:,:,d]= mapwindow(std,tt,(5,5))
end

mi = [Tuple(argmax(sds))...][1:2]
display(mi)
mir= map(i->filter(j->j>0,(-5:5) .+i),mi)
display(mir)
sds1 = [sum(sds[mir...,d]) for d in 1:41]
# bm1 = mean(t1[1,:])
# display(bm1)
# bsd1 = std(t1[1,:])
# display(bsd1)
# bm2 = mean(t2[1,:])
# bsd2 = std(t2[1,:])
# bm3 = mean(t3[1,:])
# bsd3 = std(t3[1,:])
# bm4 = mean(t4[1,:])
# bsd4 = std(t4[1,:])
#
#
# sds1=sum((t1.-bm1)./bsd1.^2,dims=2)[:]
# sds2=sum((t2.-bm2)./bsd2.^2,dims=2)[:]
# sds3=sum((t3.-bm3)./bsd3.^2,dims=2)[:]
# sds4=sum((t4.-bm4)./bsd4.^2,dims=2)[:]

plot(sds1)


# c = [:black :red :lightgreen :blue]
# plot([sds1 sds2 sds3 sds4],color=c,lw=2,grid=false)
# a=[sds1[[1,2,end,end-1]],sds2[[1,2,end,end-1]],sds3[[1,2,end,end-1]],sds4[[1,2,end,end-1]]]'
# hline!(mean.(a).+5*std.(a),color=c)
end













skm = pyimport("skimage.metrics")

@manipulate for u in uids
t1 = stas[1]["usta"][u][1,:]
t2 = stas[1]["usta"][u][25,:]
ti1=fill(mean(t1),imagesize)
ti1[xi]=t1
ti2=fill(mean(t2),imagesize)
ti2[xi]=t2
ssim,di=skm.structural_similarity(ti1,ti2,full=true)

heatmap(di,clims=(0.4,1))
end
