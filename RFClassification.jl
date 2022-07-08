using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,LinearAlgebra,Interact,Images,UMAP,
    Clustering,Distances,VegaLite,Combinatorics,XLSX,DataStructures
import GLMakie,Contour

## Collect All STAs
function collectsta(indir;unit=DataFrame(),datafile="stadataset.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            delays = dataset["delays"]
            unitid = collect(keys(dataset["usta"]))
            unitgood = map(u->dataset["ugood"][u],unitid)
            unitresponsive = map(u->dataset["uresponsive"][u],unitid)
            id = map((u,g)->"$(siteid)_$(g ? "S" : "M")U$u",unitid,unitgood)
            roi = map((u,r)->r ? (;dataset["ulroi"][u].centerdeg,dataset["ulroi"][u].radiideg,dataset["ulroi"][u].radiusdeg) : missing,unitid,unitresponsive)
            rc = map((u,r)->r ? map((c,cr)->cr ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) : missing,unitid,unitresponsive)
            dc = map((u,r)->r ? map((d,cr)->cr ? delays[d.d] : missing,dataset["ulcsdd"][u],dataset["ucresponsive"][u]) : missing,unitid,unitresponsive)
            df = DataFrame(siteid=siteid,id=id,responsive=unitresponsive,roi=roi,rc=rc,dc=dc)
            if haskey(dataset,"ulfit")
                foreach(m->df[!,"sta!$m"]=map((u,r)->r ? map(f->ismissing(f) ? missing : (;(k=>f[k] for k in keys(f))...),dataset["ulfit"][u][m]) : missing,
                        unitid,unitresponsive), keys(first(values(dataset["ulfit"]))))
            end
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
figfmt = [".png"]

staunit = collectsta(resultroot)
# save(joinpath(resultroot,"staunit.jld2"),"staunit",staunit)

staunit = load(joinpath(resultroot,"staunit.jld2"),"staunit")
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
layer = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate")
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
penetration = select(penetration,[:siteid,:od,:cofd],:cofd=>ByRow(i->i∈["L/M","S/LM","L/M, S/LM"])=>:incofd,
                :cofd=>ByRow(i->i∈["B","W","B, W"])=>:inbw,:cofd=>ByRow(i->i=="None")=>:none)

## STA responsive
stacu = innerjoin(allcu,staunit,on=[:siteid,:id])
leftjoin!(stacu,penetration,on=:siteid)
transform!(stacu,:rc=>ByRow(i->ismissing(i) ? missing : join(skipmissing(i)))=>:RTYPE)
stacu.CTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([L,M,S]+)"=>s"\1"),stacu.RTYPE)
stacsu = subset(stacu,:good)
starcsu = dropmissing(stacsu,:rc)
stanrcsu = subset(stacsu,:responsive=>.!)
jldsave(joinpath(resultroot,"stanrcsu.jld2");stanrcsu)

plotlstaroi = (unit;siteid="AG1_V1_ODL1",cg=cgrad(:matter),dir=nothing,figfmt=[".png"],title="$(siteid)_lsta_roi")->begin
    cunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true)
    x = map(i->i.centerdeg[2],cunit.roi)
    y = map(i->i.centerdeg[1],cunit.roi)
    rx = map(i->i.radiideg[2],cunit.roi)
    ry = map(i->i.radiideg[1],cunit.roi)
    roix = [-1,1,1,-1,-1] .* rx' .+ x'
    roiy = [1,1,-1,-1,1] .* ry' .+ y'
    p=plot(roix,roiy;frame=:zerolines,tickdir=:out,xmirror=true,lw=1,ratio=1,la=0.8,
    xlabel="Stimulus_X (deg)",ylabel="Stimulus_Y (deg)",yflip=true,leg=false,color=cg[cunit.aligndepth]')
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotlstaroi(starcsu;siteid="AG2_V1_ODR13")
# foreach(siteid->plotlstaroi(starcsu;siteid,dir=joinpath(resultroot,"RF")),levels(starcsu.siteid))


plotdepthhist = (unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,title="CorticalDepth_$g",ylabel="Number of Units",palette=:tab20,leg=:best) -> begin
p = groupedhist(unit.aligndepth;group=unit[!,g],barpositions=:stack,bin=0:0.02:1,permute=(:y,:x),xflip=true,grid=false,legendfontsize=6,
    xlims=(0,1),lw=0,size=(350,500),tickor=:out,xlabel="Cortical Depth",ylabel,palette,leg)
if !isnothing(layer)
    ann = [(1,mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
    hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
end
isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

# foreach(siteid->plotdepthhist(filter(r->r.siteid==siteid,stacsu);dir=resultroot,g=:responsive,title="$(siteid)_CorticalDepth_$g",layer),levels(stacsu.siteid))
# foreach(siteid->plotdepthhist(filter(r->r.siteid==siteid,stacsu);dir=resultroot,g=:RTYPE,title="$(siteid)_CorticalDepth_$g",layer),levels(stacsu.siteid))

plotdepthhist(stacsu;g=:responsive,layer)
plotdepthhist(starcsu;g=:RTYPE,layer)
plotdepthhist(starcsu;g=:CTYPE,layer)

p=select(starcsu,[:RTYPE,:siteid,:layer]) |>
[@vlplot(:bar,y={"RTYPE"},x={"count()"});
@vlplot(:bar,y={"siteid"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}});
@vlplot(:bar,y={"layer"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}})]

# COFD, Achromatic and None
plotdepthhist(subset(starcsu,:none);g=:RTYPE,layer)

p=select(subset(starcsu,:none),[:RTYPE,:layer]) |>
[@vlplot(:bar,x={"RTYPE"},y={"count()"});
@vlplot(:bar,y={"layer"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}})]


## Spatial Pattern
staedogunit = (unit;model=:edog) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:cor,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:ae,
        "sta!$model"=>ByRow(i->i.param[7])=>:ai,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σxe,
        "sta!$model"=>ByRow(i->i.param[5])=>:r,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[8])=>:ρ)
end
stagaborunit = (unit;model=:gabor) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:cor,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:a,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σx,
        "sta!$model"=>ByRow(i->i.param[5])=>:r,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[7])=>:f,
        "sta!$model"=>ByRow(i->i.param[8])=>:p)
end
# brs(r) = r > 1 ? r - 1 : 1 - 1/r
# sinoddeven(p) = p < 0.5 ? 4abs(p-0.25) : 4abs(p-0.75)
# sinonoff(p) = p < 0.5 ? 1 - 4abs(p-0.25) : 4abs(p-0.75) - 1
function fitpredict(fit;ppd=45,rp=150,tight=false,metric=false)
   if tight
       param = deepcopy(fit.param)
       ppd = max(ppd,rp/min(fit.radii...)) # variable ppd to have rp pixels in local roi radius for all colors
       if fit.model == :edog
           param[[2,4]].=0
           σxe = param[3]
           r = param[5]
           ρ = param[8]
           radius = 3*max(σxe,r*σxe,ρ*σxe,ρ*r*σxe) # predict img radius
           radii = (radius,radius)
           fit = (;fit.model,fit.fun,fit.mfun,param,radii)
       elseif fit.model == :gabor
           param[[2,4]].=0
           σx = param[3]
           r = param[5]
           radius = 3*max(σx,r*σx)
           radii = (radius,radius)
           fit = (;fit.model,fit.fun,fit.mfun,param,radii)
       end
   end
   y = -fit.radii[1]:1/ppd:fit.radii[1]
   x = -fit.radii[2]:1/ppd:fit.radii[2]
   img = predict(fit,x,y,xygrid=true,yflip=true)
   if metric
       mask = tight ? predict((;fun=fit.mfun,param),x,y,xygrid=true,yflip=true) : trues(size(img))
       mimg = img[mask]
       onidx = mimg.>=0
       onap = count(onidx)/length(onidx)
       W = sum(mimg.^2)
       onwp = sum(mimg[onidx].^2)/W
       return (;img,onap,onwp,W)
   else
       return  img
   end
end
function fitpredict(fit,μx,μy,r;rp=170,ppd=max(50,rp/r))
    param = deepcopy(fit.param)
    param[2]=param[2]-μx
    param[4]=param[4]-μy
    fit = (;fit.model,fit.fun,fit.mfun,param,radii=(r,r))

    y = -fit.radii[1]:1/ppd:fit.radii[1]
    x = -fit.radii[2]:1/ppd:fit.radii[2]
    img = predict(fit,x,y,xygrid=true,yflip=true)
end

# cortical single units
edogcsu = load(joinpath(resultroot,"edogcsu.jld2"),"edogcsu")
gaborcsu = load(joinpath(resultroot,"gaborcsu.jld2"),"gaborcsu")

edogcsu = staedogunit(stacsu)
transform!(edogcsu,
    [:σxe,:r,:ρ]=>ByRow((σxe,r,ρ)->5*max(σxe,r*σxe,ρ*σxe,ρ*r*σxe))=>:diameter,
    [:ae,:ai]=>ByRow(/)=>:ρa,
    # :r=>ByRow(brs∘inv)=>:el,
    :r=>ByRow(log2∘inv)=>:el,
    :r=>ByRow(abs∘log2∘inv)=>:rd,
    # :ρ=>ByRow(brs)=>:cs,
    :ρ=>ByRow(log2)=>:cs,
    # :ρ=>ByRow(abs∘brs)=>:op,
    # [:ae,:ai]=>ByRow(brs∘/)=>:onoff,
    [:ae,:ai]=>ByRow(log2∘/)=>:onoff,
    [:ae,:ai]=>ByRow(sign∘log2∘/)=>:onoffsign,
    # [:ae,:ai]=>ByRow(abs∘brs∘/)=>:amp,
    :fit=>ByRow(i->fitpredict(i,tight=true,metric=true))=>[:img,:onap,:onwp,:W])

gaborcsu = stagaborunit(stacsu)
transform!(gaborcsu,
    [:σx,:r]=>ByRow((σx,r)->5*max(σx,r*σx))=>:diameter,
    [:σx,:r,:f]=>ByRow((σx,r,f)->5*r*σx*f)=>:cyc,
    # :r=>ByRow(brs∘inv)=>:el,
    :r=>ByRow(log2∘inv)=>:el,
    :r=>ByRow(abs∘log2∘inv)=>:rd,
    :p=>ByRow(p->sin(2π*p))=>:onoff,
    :p=>ByRow(p->(sign∘sin)(2π*p))=>:onoffsign,
    :p=>ByRow(p->(abs∘sin)(2π*p))=>:oddeven,
    # :p=>ByRow(sinoddeven)=>:oddeven,
    # :p=>ByRow(sinonoff)=>:onoff,
    :fit=>ByRow(i->fitpredict(i,tight=true,metric=true))=>[:img,:onap,:onwp,:W])

gn = [:id,:rc,:cor,:bic,:fvu]
csumg = innerjoin(edogcsu[:,gn],gaborcsu[:,gn],on=[:id,:rc],renamecols=:edog=>:gabor)


@df edogcsu scatter(:ρ,:ρa,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xticks=0:0.5:4,xlim=(-0.1,4.1),xlabel="ρ",ylabel="ρa")
@df edogcsu scatter(:cs,:onoff,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-10,10),
    xlabel="CenterSurround",ylabel="OnOff")
@df edogcsu corrplot(cols([:diameter,:cs,:onap,:onwp]),bin=50,leg=false,size=(850,650),grid=false)


@df gaborcsu scatter(:p,:cyc,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xlabel="p",ylabel="cyc")
@df gaborcsu scatter(:onoff,:cyc,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xlabel="OnOff",ylabel="cyc")
@df gaborcsu corrplot(cols([:diameter,:onoff,:onap,:onwp]),bin=50,leg=false,size=(850,650),grid=false)


@df csumg scatter(:bicgabor,:bicedog,ratio=:equal,leg=false,ms=3,msw=0,ma=0.7,grid=false)
lim=1.5e5
plot!([0,lim],[0,lim],color=:black,lw=1,xlim=[0,lim],ylim=[0,lim],xlabel="BIC (Gabor)",ylabel="BIC (eDoG)",xticks=false,yticks=false)
foreach(ext->savefig(joinpath(resultroot,"csu_bic$ext")),figfmt)

@df csumg scatter(:fvugabor,:fvuedog,ratio=:equal,leg=false,ms=3,msw=0,ma=0.7,grid=false)
plot!([0,1],[0,1],color=:black,lw=1,xlim=[0,1],xlabel="FVU (Gabor)",ylabel="FVU (eDoG)")
foreach(ext->savefig(joinpath(resultroot,"csu_fvu$ext")),figfmt)

@df csumg scatter(:corgabor,:coredog,ratio=:equal,leg=false,ms=3,msw=0,ma=0.7,grid=false)
plot!([0,1],[0,1],color=:black,lw=1,xlim=[0,1],xlabel="R (Gabor)",ylabel="R (eDoG)")
foreach(ext->savefig(joinpath(resultroot,"csu_cor$ext")),figfmt)



## Spatial Pattern Classification
function d2clu(x;r=1)
    cs = dbscan(x,r)
    cid = zeros(Int,size(x,2))
    foreach(i->cid[cs[i].core_indices] .= i, 1:length(cs))
    cid
end

ccodecm = Dict('A'=>ColorMaps["dkl_mcclumiso"],'L'=>ColorMaps["lms_mccliso"],'M'=>ColorMaps["lms_mccmiso"],'S'=>ColorMaps["lms_mccsiso"])
maskimg(i;imgsize=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,imgsize);masktype))
maskimg(i,c;imgsize=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,imgsize);color=ccodecm[c].colors,masktype))

edogimg = map(maskimg,edogcsu.img)
gaborimg = map(maskimg,gaborcsu.img)
edogimgc = map((i,c)->maskimg(i,c),edogcsu.img,edogcsu.rc)
gaborimgc = map((i,c)->maskimg(i,c),gaborcsu.img,gaborcsu.rc)

@manipulate for i in 1:length(edogimg)
    edogimg[i]
    # edogimgc[i]
end
@manipulate for i in 1:length(gaborimg)
    gaborimg[i]
    # gaborimgc[i]
end

# eDoG
features = [:onoffsign,:onap,:onwp]
Y = hcat(map(i->edogcsu[:,i],features)...)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[2,:], Y2[1,:],markersize=5)
    # Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(permutedims(Y), 2, n_neighbors=25, min_dist=1.8, n_epochs=300,metric=Euclidean())
scatter(Y2[1,:], Y2[2,:],ms=3,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

p = GLMakie.scatter(Y2[1,:], Y2[2,:], marker=edogimg, ms=30)
save(joinpath(resultroot,"sta_dog_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=2.5),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
dogcell.umapclu = d2clu(Y2,r=2.5)


save("UMAP_sta.svg",plotunitpositionimage(Y2',edogimg))



YD = pairwise(Euclidean(),Y,dims=1)
hc = hclust(YD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=4)
clu = replace(clu,2=>"eDoG-Off",3=>"eDoG-On",4=>"Off",1=>"On")

scatter(Y2[1,:], Y2[2,:],msw=0,ms=3,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600))
foreach(ext->savefig(joinpath(resultroot,"sta_edogcsu_pattern_umap_clu$ext")),figfmt)

edogcsu.pattern = clu
jldsave(joinpath(resultroot,"edogcsu.jld2");edogcsu)

dogcell |> [@vlplot(mark={:bar,binSpacing=0},x={"op",bin={step=0.03},title="Opponency"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"});
       @vlplot(mark={:bar,binSpacing=0},x={"amp",bin={step=0.2},title="OnOff Amplitude"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"})]



# Gabor
features = [:onoff,:onap,:onwp]
Y = hcat(map(i->gaborcsu[:,i],features)...)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(permutedims(Y), 2, n_neighbors=30, min_dist=1.5, n_epochs=300,metric=Euclidean())
scatter(Y2[1,:], Y2[2,:],ms=3,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=1.2),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_gabor_shape_umap_clu.png"))
gaborcell.umapclu = d2clu(Y2,r=1.2)

save("UMAP_sta.svg",plotunitpositionimage(Y2',gaborimg))




cr = kmeans(Y[:,1:2]',5)
clu = assignments(cr)

YD = pairwise(Euclidean(),Y[:,1:2],dims=1)
hc = hclust(YD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=5)
clu = replace(clu,2=>"Gabor-Even-Off",3=>"Gabor-Even-On",1=>"Gabor-Odd",4=>"On",5=>"Off")

scatter(Y2[1,:], Y2[2,:],msw=0,ms=3,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600))
foreach(ext->savefig(joinpath(resultroot,"sta_gaborcsu_pattern_umap_clu$ext")),figfmt)

gaborcsu.pattern = clu
jldsave(joinpath(resultroot,"gaborcsu.jld2");gaborcsu)


gaborcell |> [@vlplot(mark={:bar,binSpacing=0},x={"cyc",bin={step=0.1},title="Cycle"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"});
       @vlplot(mark={:bar,binSpacing=0},x={"oddeven",bin={step=0.02},title="OddEven"},y={"count()",title="Number of Spatial RF"},color={"umapclu:n"})]



names(gaborcsu)

## Model Selection and Goodness
gaborbetter = csumg.bicgabor .< csumg.bicedog
smcsu = sort!(append!(gaborcsu[gaborbetter,:],edogcsu[.!gaborbetter,:],cols=:union),:rc)
transform!(smcsu,:pattern=>ByRow(i->!startswith(i,'O'))=>:so)
srmcsu = filter(r->r.cor > 0.6,smcsu)

plotdepthhist(smcsu;g=:pattern,layer,ylabel="Number of Spatial Pattern")
plotdepthhist(srmcsu;g=:pattern,layer,ylabel="Number of Spatial Pattern")
plotdepthhist(srmcsu;g=:so,layer,ylabel="Number of Spatial Pattern")

foreach(siteid->plotdepthhist(srmcsu;dir=resultroot,siteid,g=:pattern,layer,ylabel="Number of Spatial Pattern"),levels(srmcsu.siteid))
foreach(siteid->plotdepthhist(srmcsu;dir=resultroot,siteid,g=:sop,layer,ylabel="Number of Spatial Pattern"),levels(srmcsu.siteid))

p=select(srmcsu,[:pattern,:siteid,:layer]) |>
[@vlplot(:bar,y={"pattern"},x={"count()"});
@vlplot(:bar,y={"siteid"},x={"count()"},color={"pattern:n",scale={scheme=:category20}});
@vlplot(:bar,y={"layer"},x={"count()"},color={"pattern:n",scale={scheme=:category20}})]


# vdogcell |> [@vlplot(:bar,x={"r",bin={step=0.05,extent=[0,1]},title="r"},y={"count()",title="Number of Linear Filter"});
#                 @vlplot(mark={:line,size=2},x=:r,transform=[{sort=[{field=:r}],window=[{op="count",as="cum"}],frame=[nothing,0]}],
#                 y={"cum",title="Cumulated Number of Linear Filter"})]



plotmodelcontour = (unit;siteid="AG1_V1_ODL1",cg=cgrad(:matter),dir=nothing,figfmt=[".png"],title="$(siteid)_model_contour")->begin
    cunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true)
    t = 0:0.02:2π
    p=plot(frame=:zerolines,tickdir=:out,xmirror=true,lw=1,ratio=1,la=0.8,
        xlabel="Stimulus_X (deg)",ylabel="Stimulus_Y (deg)",leg=false)
    for r in eachrow(cunit)
        xys = map(i->r.fit.cfun(i,r.fit.param).+Tuple(reverse(r.roi.centerdeg)),t)
        plot!(p,xys,color=cg[r.aligndepth])
    end
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotmodelcontour(srmcsu,siteid="AG1_V1_ODL1")
# foreach(siteid->plotmodelcontour(srmcsu;siteid,dir=joinpath(resultroot,"RF")),levels(srmcsu.siteid))


plotdepthscatter = (unit;x=:diameter,g=:pattern,dir=nothing,siteid="all",figfmt=[".png"],layer=nothing,xlabel="",
    wfun=x->isempty(x) ? 0 : median(x),palette=:tab20,leg=:best) -> begin
cunit = siteid=="all" ? unit : filter(r->r.siteid==siteid,unit)
title="$(siteid)_CorticalDepth_$g"
xmin,xmax = extrema(cunit[!,x])
p = Plots.scatter(cunit[!,x],cunit.aligndepth;group=cunit[!,g],yflip=true,grid=false,legendfontsize=6,
    ms=3,msw=0,ma=0.7,ylims=(0,1),size=(350,500),tickor=:out,ylabel="Cortical Depth",xlabel,palette,leg)

n,y = unitdensity(cunit.aligndepth;w=cunit[!,x],wfun,spacerange=(0,1),bw=0.05,step=0.05)
Plots.plot!(p,n,y,color=:gray40,label="Average")
if !isnothing(layer)
    ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
    hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
end
isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotdepthscatter(srmcsu;x=:diameter,g=:pattern,layer,xlabel="Diameter")
plotdepthscatter(srmcsu;x=:rd,g=:pattern,layer,xlabel="Roundness")


## Relation of Spatial Pattern between Spectral Channels
function plotssmpair(ssm,mp;w=350,h=350,dlim=nothing,mulim=nothing)
    pssm = filter(r->contains(r.MRTYPE,mp.first) && contains(r.MRTYPE,mp.second),ssm)

    p = Plots.plot(layout=(3,1),leg=false,size=(w,3h))
    pn = :diameter
    x = getindex.(pssm[!,pn],mp.first)
    y = getindex.(pssm[!,pn],mp.second)
    lim = isnothing(dlim) ? max(maximum(x),maximum(y)) : dlim
    Plots.scatter!(p[1],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Diameter",
    ms=3,msw=0,ma=0.7,ratio=1)
    Plots.plot!(p[1],[0,lim],[0,lim],color=:gray30)

    x = rad2deg.(getindex.(pssm[!,:θ],mp.first))
    y = rad2deg.(getindex.(pssm[!,:θ],mp.second))
    xticks=yticks=0:45:180
    Plots.scatter!(p[2],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Orientation",
    xticks,yticks,ms=3,msw=0,ma=0.7,ratio=1)
    Plots.plot!(p[2],[0,180],[0,180],color=:gray30)

    x = getindex.(pssm[!,:μx],mp.second) .- getindex.(pssm[!,:μx],mp.first)
    y = getindex.(pssm[!,:μy],mp.second) .- getindex.(pssm[!,:μy],mp.first)
    lim = isnothing(mulim) ? maximum(abs.([x;y])) : mulim
    xlims=ylims=[-lim,lim]
    Plots.scatter!(p[3],x,y;xlabel="X$(mp.second) - X$(mp.first) (deg)",ylabel="Y$(mp.second) - Y$(mp.first) (deg)",
    title="Center Displacement",ms=3,msw=0,ma=0.7,ratio=1,xlims,ylims)
    Plots.vline!(p[3],[0],color=:gray30)
    Plots.hline!(p[3],[0],color=:gray30)
end


ssmcsu = rightjoin(combine(groupby(srmcsu,:id),:rc=>join=>:MRTYPE,
        [:rc,:diameter]=>((c,v)->OrderedDict((c.=>v)...))=>:diameter,
        [:rc,:θ]=>((c,v)->OrderedDict((c.=>v)...))=>:θ,
        [:rc,:μx]=>((c,v)->OrderedDict((c.=>v)...))=>:μx,
        [:rc,:μy]=>((c,v)->OrderedDict((c.=>v)...))=>:μy,
        [:rc,:W]=>((c,v)->OrderedDict((c.=>v)...))=>:W,
        :μx=>mean=>:ssmμx,
        :μy=>mean=>:ssmμy,
        [:μx,:μy,:diameter]=>((x,y,d)->maximum([abs.(x.-mean(x)).+0.6d;abs.(y.-mean(y)).+0.6d]))=>:ssmradius,
        [:rc,:fit]=>((c,v)->OrderedDict((c.=>v)...))=>:ssm,
        [:rc,:pattern]=>((c,p)->join(c.*'_'.*p,", "))=>:ss,
        [:rc,:so]=>((c,o)->join(c.*'_'.*replace(o,true=>"O",false=>"N"),", "))=>:sso),
        stacsu,on=:id)
transform!(ssmcsu,:MRTYPE=>ByRow(i->!ismissing(i))=>:mresponsive)
ssrmcsu = dropmissing(ssmcsu,:MRTYPE)
transform!(ssrmcsu,:MRTYPE=>ByRow(i->replace(i,r"A([L,M,S]+)"=>s"\1"))=>:MCTYPE,
        :ss=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:css,
        :sso=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:csso,
        :W=>ByRow(d->begin
            ub = reduce(max,values(d))
            OrderedDict(k=>v/ub for (k,v) in d)
        end)=>:nW,
        :W=>ByRow(d->begin
            L = haskey(d,'L') ? d['L'] : 0
            M = haskey(d,'M') ? d['M'] : 0
            S = haskey(d,'S') ? d['S'] : 0
            L==M==S==0 && return [missing,missing]
            Σ = L+M+S
            [L/Σ,M/Σ]
        end)=>[:wl,:wm])

plotdepthhist(ssmcsu;g=:mresponsive,layer)

plotssmpair(ssrmcsu,'A'=>'L')
plot((plotssmpair(ssrmcsu,i=>j,dlim=6,mulim=0.6) for (i,j) in combinations(['A','L','M','S'],2))...,
            layout=(1,6),size=(6*350,3*350))

# Cone Weight
csrmcsu = dropmissing(ssrmcsu,[:wl,:wm])
plot([0,1],[1,0],ratio=1,leg=false,xlims=(-0.01,1.01),ylims=(-0.01,1.01),color=:black,tickor=:out)
scatter!(csrmcsu.wl,csrmcsu.wm,ms=4,msw=0,color=:dodgerblue,xlabel="L Cone Weight",ylabel="M Cone Weight")

# typecell |>
#     @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0,y=1},{x=1,y=0}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=1}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=-1}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=0,y=-1},{x=1,y=0}]},x=:x,y=:y},
#     {:point,x={:wl,title="L Cone Weight"},y={:wm,title="M Cone Weight"},color={"layer",scale={scheme=:category10}}}])

## Cell Spectral-Spatial Classification
function cellsop(s)
    us = unique(s)
    length(us) == 1 ? us[1] : -1
end
function celltype(so,co)
    so == -1 && return -1
    so == 1 ? (co ? 1 : 3) : (co ? 2 : 0)
end

function plotssmimg(ssm;rand=false,type=:img,ni=nothing,maxn=50,col=5,h=130,w=type==:img ? 450 : 130,path=nothing,bg=:gray,
    imgsize=(55,55),titlefontsize=type==:img ? 14 : 6,hidetitle=false)
    nu = nrow(ssm); n = min(maxn,nu); n==0 && return
    nc = n < col ? n : col
    nr = n <= col ? 1 : ceil(Int,n/col)
    ni = rand ? sample(1:nu,n,replace=false) : isnothing(ni) ? (1:n) : ni
    p=plot(;layout=(nr,nc),size=(nc*w,nr*h),frame=:none,leg=false,titlefontsize)
    for i in eachindex(ni)
        title = hidetitle ? "" : "$(ssm.id[ni[i]])"
        ssmimg = ssm.ssmimg[ni[i]]
        if type == :img
            t = hcat((maskimg(v,c;imgsize,masktype="None") for (c,v) in ssmimg)...)
            plot!(p[i],t;title,bginside=bg,frame=:none)
        elseif type == :contour
            x = 1:imgsize[2]; y = 1:imgsize[1]; xlims=(x[begin],x[end]);ylims=(y[begin],y[end])
            for (c,v) in ssmimg
                t = imresize(v,imgsize)
                lb,ub = extrema(t)
                lim = max(abs(lb),abs(ub))
                levels=[0.1lb,0.1ub]
                colors = get(cgrad(ccodecm[c].colors),[lb,ub],(-lim,lim))
                for j in 1:2
                    for l in Contour.lines(Contour.contour(y,x,t,levels[j]))
                        ys,xs = Contour.coordinates(l)
                        plot!(p[i],xs,ys;color=colors[j],yflip=true,ratio=1,
                            xlims,ylims,lw=2,title,bginside=bg,frame=:none)
                    end
                end
            end
        end
    end
    isnothing(path) && return p
    mkpath(dirname(path));savefig(path)
end


transform!(ssrmcsu,[:ssm,:ssmμx,:ssmμy,:ssmradius]=>ByRow((m,μx,μy,r)->OrderedDict(k=>fitpredict(v,μx,μy,r) for (k,v) in m))=>:ssmimg)
p=select(ssrmcsu,[:ss,:sso,:css,:csso]) |> @vlplot(:bar,y={"count()"},x={"ss"})

foreach(t->plotssmimg(filter(r->r.sso == t,ssrmcsu);maxn=nrow(ssrmcsu),col=5,hidetitle=true,type=:contour,path=joinpath(resultroot,"RF","$(t).contour.png")),
        levels(ssrmcsu.sso))

plotssmimg(filter(r->r.sso == "S_O",ssrmcsu);maxn=1,col=1,hidetitle=true,rand=true)
plotssmimg(filter(r->r.sso == "A_O, L_N",ssrmcsu);maxn=3,col=1,hidetitle=true,type=:contour)
plotssmimg(filter(r->r.sso == "A_O, L_N" || r.sso == "A_N, L_O",ssrmcsu);maxn=20,col=1,hidetitle=true,type=:contour)
plotssmimg(filter(r->r.sso == "L_O, M_O, S_N",ssrmcsu);maxn=5,col=1,hidetitle=true,type=:img)
plotssmimg(filter(r->r.sso == "A_O, L_O, M_N, S_N",ssrmcsu);maxn=1,col=1,hidetitle=true,type=:img)



ni = [2,4]



savefig("test.png")
save("test.png",p)

tt = dropmissing(ssmcsu)
plotdepthhist(t;g=:sops,layer,leg=:outerright)

t=filter(r->r.MRTYPE == "L" && r.layer == "4Cβ",select(ssrmcsu,Not(r"\w*!\w*")))


t |> @vlplot(:bar,x={"count()"},y={"ss"})



A4cβ=filter(r->r.MRTYPE == "A" && r.layer == "4Cβ",select(ssmrcsu,Not(r"\w*!\w*"))) 



pyplot()



typecell = select!(leftjoin(combine(groupby(doggaborcell,:id),:sop=>(i->cellsop(i))=>:sop,:sign=>(i->[1,-1] ⊆ i)=>:cop,
                [:rc,:W,:sign]=>((c,w,s)->rcw(c,w,s,cone='L'))=>:wl,
                [:rc,:W,:sign]=>((c,w,s)->rcw(c,w,s,cone='M'))=>:wm,
                [:rc,:W,:sign]=>((c,w,s)->rcw(c,w,s,cone='S'))=>:ws,
                [:rc,:sign]=>((i,j)->join(map((p,q)->p*(q==1 ? "" : "̲"),i,j)))=>:sr,
                [:rc,:img]=>((c,i)->[srimg(c,i)])=>:srimg), cell,on=:id),Not(r"\w*!\w*"),[:sop,:cop]=>ByRow((i,j)->celltype(i,j))=>:celltype)


typecell |> [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"depth",bin={step=200},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"celltype:n"})]









save("UMAP_c_sta.svg",plotunitpositionimage(Y2',rswcell.srimg))

save("UMAP_c_sta_layer.svg",plotunitlayerimage(rswcell.layer,rswcell.imgs,width=2000,markersize=30))

## non-sta but orisf responsive






