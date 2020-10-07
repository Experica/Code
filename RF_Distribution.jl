using Combinatorics,Query,VegaLite,DataFrames,Clustering,MultivariateStats,UMAP,Makie

## UMAP of responsive sta
lsta = mapreduce(u->mapreduce((c,r,exd)->r ? [dataset["ulsta"][u][:,:,exd.d,c]] : [],
append!,1:testn,dataset["ucresponsive"][u],dataset["ulcexd"][u]),
append!,keys(dataset["ulsta"]))

maxsizepx=mapreduce(i->maximum(size(i)),max,lsta)
map!(i->imresize(i,maxsizepx,maxsizepx),lsta,lsta)
@manipulate for i in eachindex(lsta)
    heatmap(lsta[i],yflip=true,aspect_ratio=:equal,frame=:none,color=:coolwarm)
end

mlsta = map(i->begin
                   m = maximum(abs.(extrema(i)))
                   alphamask(get(cgrad(:coolwarm),i,(-m,m)),radius=0.5,sigma=0.35,masktype="Gaussian")[1]
               end,lsta)
@manipulate for i in eachindex(mlsta)
    mlsta[i]
end

# PCA reduction
X = mapreduce(vec,hcat,lsta)
foreach(i->X[i,:]=zscore(X[i,:]),1:size(X,1))
M = fit(PCA,X,maxoutdim=60,pratio=0.9,mean=0)
Y = MultivariateStats.transform(M,X)

# plot t-SNE sta
@manipulate for i in 1000:2000, p in 5:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = tsne(permutedims(Y), 2, 0, 2000, 10)
p = plotunitpositionimage(Y2,mlsta)
save("tSNE_rlSTA.svg",p)

# plot UMAP sta
@manipulate for i in 200:2000, n in 5:60, d in 0.001:0.001:10
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = umap(Y, 2, n_neighbors=10, min_dist=0.01, n_epochs=500)'
p = plotunitpositionimage(Y2,mlsta)
save("UMAP_rlSTA.svg",p)



## Collect RF and Merge into Cells
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

cells = load(joinpath(resultroot,"cells.jld2"),"cells")


cells |> [@vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"});
          @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"})]

## RF Spatial Pattern
srfcells = transform(dropmissing!(flatten(dropmissing(cells,:rc),[:rf,:rc]),:rf),:rf=>ByRow(x->x.r)=>:r)

srfcells |> @vlplot(height=400,width=550,:bar,x={"r",bin={step=0.05,extent=[0,1]},title="r"},y={"count()",title="Number of Spatial RF"})
srfcells |> @vlplot(height=400,width=550,mark={:line,size=3},x=:r,transform=[{
                    sort=[{field=:r}],
                    window=[{op="count",as="cum"}],
                    frame=[nothing,0]
                    }],y={"cum",title="Cumulated Number of Spatial RF"})

odd(p) = p < 0.5 ? 4abs(p-0.25) : 4abs(p-0.75)
onoff(p) = p < 0.5 ? 1 : -1
srfcells = transform(filter!(r->r.r>=0.6,srfcells),[:rf=>ByRow(x->x.param[5]/x.param[3])=>:ar,
                                                    :rf=>ByRow(x->4*x.param[5]*x.param[7])=>:nc,
                                                    :rf=>ByRow(x->odd(x.param[8]))=>:odd,
                                                    :rf=>ByRow(x->onoff(x.param[8]))=>:onoff,
                                                    :rf=>ByRow(x->rad2deg(x.param[6]))=>:ori,
                                                    :rf=>ByRow(x->4*x.param[3])=>:diameter,
                                                    :rf=>ByRow(x->x.param[2])=>:cx,
                                                    :rf=>ByRow(x->x.param[4])=>:cy,
                                                    :rf=>ByRow(x->x.param[1])=>:amp])


srfcells |> [@vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Spatial RF"},color="rc");
          @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Spatial RF"},color="rc")]

srfcells |> [@vlplot(:point,x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
             @vlplot(:point,x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
             @vlplot(:point,x={"odd",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"})]

srfcells |> [@vlplot(mark={:point,size=20},x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
          @vlplot(mark={:point,size=20},x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
          @vlplot(mark={:point,size=20},x={"odd",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"})]

srfcells |> [@vlplot(:bar,x={"ar",bin={maxbins=20},title="Aspect Ratio"},y={"count()",title="Number of Spatial RF"},color={"layer"});
             @vlplot(:bar,x={"nc",bin={maxbins=20},title="Number of Cycle"},y={"count()",title="Number of Spatial RF"},color={"layer"});
             @vlplot(:bar,x={"odd",bin={maxbins=15},title="Oddness"},y={"count()",title="Number of Spatial RF"},color={"layer"})]

srfcells |> [@vlplot(:bar,x={"ar",bin={maxbins=20},title="Aspect Ratio"},y={"count()",title="Number of Spatial RF"},color={"rc"});
          @vlplot(:bar,x={"nc",bin={maxbins=20},title="Number of Cycle"},y={"count()",title="Number of Spatial RF"},color={"rc"});
          @vlplot(:bar,x={"odd",bin={maxbins=15},title="Oddness"},y={"count()",title="Number of Spatial RF"},color={"rc"})]

srfimage = rfcolorimage.(rfimage.(srfcells.rf,ppd=40,onlyrf=true))
@manipulate for i in 1:length(srfimage)
    imresize(srfimage[i],50,50)
end


Y = [srfcells.ar srfcells.nc srfcells.odd]
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))
cluid = x->begin
    cs = dbscan(x,1,min_cluster_size=3)
    cid = zeros(Int,size(x,2))
    foreach(i->cid[cs[i].core_indices] .= i, 1:length(cs))
    cid
end
@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(Y', 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:auto,stroke(0)),group=cluid(Y2),frame=:none,leg=:inline)
end


Y2 = umap(Y', 2, n_neighbors=10, min_dist=0.02, n_epochs=200)
srfcells.sclu = cluid(Y2)
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:auto,stroke(0)),group=srfcells.sclu,frame=:none,leg=:inline)
save("UMAP_srf.svg",plotunitpositionimage(Y2',srfimage))

srfcells |> [@vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"},color={"sclu:n"});
          @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"sclu:n"})]

## Relationship between Cone Spatial RF
lm = innerjoin(filter(r->r.rc=='L',srfcells),filter(r->r.rc=='M',srfcells),on=:id,makeunique=true)
lm = transform(lm,[:cx,:cx_1]=>ByRow((i,j)->i-j)=>:xdisplace,[:cy,:cy_1]=>ByRow((i,j)->i-j)=>:ydisplace)
lm |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=4}]},x=:x,y=:x}, {:point,x={:ar,title="L"},y={:ar_1,title="M"}}],title="Aspect Ratio");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=5}]},x=:x,y=:x}, {:point,x={:nc,title="L"},y={:nc_1,title="M"}}],title="Number of Cycle");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:odd,title="L"},y={:odd_1,title="M"}}],title="Oddness")]

lm |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=180}]},x=:x,y=:x}, {:point,x={:ori,title="L"},y={:ori_1,title="M"}}],title="Orientation");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-0.3,y=0},{x=0.3,y=0}]},x=:x,y=:y},
        {mark={:line,size=1,color="black"},data={values=[{x=0,y=-0.3},{x=0,y=0.3}]},x=:x,y=:y},{:point,x={:xdisplace,title="X(L-M)"},y={:ydisplace,title="Y(L-M)"}}],title="RF Displacement");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.5}]},x=:x,y=:x}, {:point,x={:diameter,title="L"},y={:diameter_1,title="M"}}],title="Diameter")]

ls = innerjoin(filter(r->r.rc=='L',srfcells),filter(r->r.rc=='S',srfcells),on=:id,makeunique=true)
ls = transform(ls,[:cx,:cx_1]=>ByRow((i,j)->i-j)=>:xdisplace,[:cy,:cy_1]=>ByRow((i,j)->i-j)=>:ydisplace)
ls |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=4}]},x=:x,y=:x}, {:point,x={:ar,title="L"},y={:ar_1,title="S"}}],title="Aspect Ratio");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=5}]},x=:x,y=:x}, {:point,x={:nc,title="L"},y={:nc_1,title="S"}}],title="Number of Cycle");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:odd,title="L"},y={:odd_1,title="S"}}],title="Oddness")]

ls |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=180}]},x=:x,y=:x}, {:point,x={:ori,title="L"},y={:ori_1,title="S"}}],title="Orientation");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-0.3,y=0},{x=0.3,y=0}]},x=:x,y=:y},
        {mark={:line,size=1,color="black"},data={values=[{x=0,y=-0.3},{x=0,y=0.3}]},x=:x,y=:y},{:point,x={:xdisplace,title="X(L-S)"},y={:ydisplace,title="Y(L-S)"}}],title="RF Displacement");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.5}]},x=:x,y=:x}, {:point,x={:diameter,title="L"},y={:diameter_1,title="S"}}],title="Diameter")]

ms = innerjoin(filter(r->r.rc=='M',srfcells),filter(r->r.rc=='S',srfcells),on=:id,makeunique=true)
ms = transform(ms,[:cx,:cx_1]=>ByRow((i,j)->i-j)=>:xdisplace,[:cy,:cy_1]=>ByRow((i,j)->i-j)=>:ydisplace)
ms |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=4}]},x=:x,y=:x}, {:point,x={:ar,title="M"},y={:ar_1,title="S"}}],title="Aspect Ratio");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=5}]},x=:x,y=:x}, {:point,x={:nc,title="M"},y={:nc_1,title="S"}}],title="Number of Cycle");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:odd,title="M"},y={:odd_1,title="S"}}],title="Oddness")]

ms |> [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=180}]},x=:x,y=:x}, {:point,x={:ori,title="M"},y={:ori_1,title="S"}}],title="Orientation");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-0.3,y=0},{x=0.3,y=0}]},x=:x,y=:y},
        {mark={:line,size=1,color="black"},data={values=[{x=0,y=-0.3},{x=0,y=0.3}]},x=:x,y=:y},{:point,x={:xdisplace,title="X(M-S)"},y={:ydisplace,title="Y(M-S)"}}],title="RF Displacement");
        @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.5}]},x=:x,y=:x}, {:point,x={:diameter,title="M"},y={:diameter_1,title="S"}}],title="Diameter")]


## Cone Weights
rcw = (c,a,s;cone='L')-> begin
    ci = findfirst(c.==cone)
    isnothing(ci) ? 0 : s[ci]*a[ci]/sum(a)
end
rcwcells = combine(groupby(filter(r->r.rc!='A',srfcells),:id),[:rc,:amp,:onoff]=>((c,a,s)->rcw(c,a,s,cone='L'))=>:wl,
                                                            [:rc,:amp,:onoff]=>((c,a,s)->rcw(c,a,s,cone='M'))=>:wm,
                                                            [:rc,:amp,:onoff]=>((c,a,s)->rcw(c,a,s,cone='S'))=>:ws)
rcwcells = leftjoin(rcwcells,cells,on=:id)
rcwcells |> @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0,y=1},{x=1,y=0}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=1}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=-1}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=0,y=-1},{x=1,y=0}]},x=:x,y=:y},
                            {:point,x={:wl,title="L Cone Weight"},y={:wm,title="M Cone Weight"},color=:layer}])


## type of responsiveness
CTYPES = sort(join.(sort.(combinations(['A','L','M','S']))))
CTCOLORS = cgrad(:rainbow,length(CTYPES),categorical = true).colors.colors

rdf = DataFrame(ResponseType=map(ct->join(skipmissing(ct)),crtype))
p = rdf |> @vlplot(:bar,y=:ResponseType,x="count()",width=600,height=400)

mrdf = DataFrame(ModelResponseType=filter!(i->!isempty(i),map(ct->join(skipmissing(ct)),cmrtype)))
p = mrdf |> @vlplot(:bar,y=:ModelResponseType,x="count()",width=600,height=400)

crdf = DataFrame(ConeResponseType=filter(i->!isempty(i),map(t->match(r"A?([L,M,S]*)",t).captures[1],mrdf.ModelResponseType)))
p = crdf |> @vlplot(:bar,y=:ConeResponseType,x="count()",width=600,height=400)
## ALMS space
X = map(fs->vcat(map(f->ismissing(f) || f.r <= 0.6 ? fill(missing,8) : f.param,fs)...),rfs)

Y = replace(hcat(filter(i->!isempty(skipmissing(i)),X)...),missing=>0)
Yg = mrdf.ModelResponseType
foreach(j->begin
    ai = 1:8:size(Y,1)
    maxamp = maximum(Y[ai,j])
    foreach(i->Y[i,j]/=maxamp,ai)
end,1:size(Y,2))

foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 200:1000, n in 5:60, d in 0.001:0.001:5
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    Plots.scatter(Y2[:,1], Y2[:,2], marker=(3,3,:circle,stroke(0)),group=Yg,frame=:none,leg=:inline)
end

## LMS space
LMSs = hcat(filter(i->!isempty(skipmissing(i)),map(fs->vcat(map(f->ismissing(f) || f.r <= 0.6 ? fill(missing,8) : f.param,fs[2:end])...),rfs))...)

σx = LMSs[3:8:end,:]
p = Makie.scatter(σx[1,:],σx[2,:],σx[3,:],markersize=0.01,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx=(σx;title="4σₓ") -> begin
    lmi = [!ismissing(σx[1,i]) && !ismissing(σx[2,i]) for i in 1:size(σx,2)]
    lm = σx[[1,2],lmi]

    lsi = [!ismissing(σx[1,i]) && !ismissing(σx[3,i]) for i in 1:size(σx,2)]
    ls = σx[[1,3],lsi]

    msi = [!ismissing(σx[2,i]) && !ismissing(σx[3,i]) for i in 1:size(σx,2)]
    ms = σx[[2,3],msi]

    lim = maximum(skipmissing(σx))
    p = Plots.plot(layout=(1,3),leg=false,xlims=[0,lim],ylims=[0,lim],title=title)
    Plots.plot!(p,subplot=1,x->x,0,lim)
    Plots.scatter!(p,subplot=1,lm[1,:],lm[2,:],aspect_ratio=1,xlabel="L",ylabel="M")
    Plots.plot!(p,subplot=2,x->x,0,lim)
    Plots.scatter!(p,subplot=2,ls[1,:],ls[2,:],aspect_ratio=1,xlabel="L",ylabel="S")
    Plots.plot!(p,subplot=3,x->x,0,lim)
    Plots.scatter!(p,subplot=3,ms[1,:],ms[2,:],aspect_ratio=1,xlabel="M",ylabel="S")
    p
end

plotσx(4σx)

σyoverσx = mapreduce(i->map(j->LMSs[8j+5,i]/LMSs[8j+3,i],0:2),hcat,1:size(LMSs,2))
p = Makie.scatter(σyoverσx[1,:],σyoverσx[2,:],σyoverσx[3,:],markersize=0.1,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx(σyoverσx,title="σy/σx")

f4σy = mapreduce(i->map(j->4LMSs[8j+5,i]*LMSs[8j+7,i],0:2),hcat,1:size(LMSs,2))
p = Makie.scatter(f4σy[1,:],f4σy[2,:],f4σy[3,:],markersize=0.1,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx(f4σy,title="f4σy")


ω = mapreduce(i->map(j->Ω(LMSs[8j+8,i]),0:2),hcat,1:size(LMSs,2))


plotσx(ω,title="ω")


θ = LMSs[6:8:end,:]

plotσx(180θ/π,title="θ")


dlm = map(i->sqrt(sum((LMSs[[2,4],i] .- LMSs[[2,4].+8,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dlm)),nbins=30,xlabel="LM Distance",leg=false)

dls = map(i->sqrt(sum((LMSs[[2,4],i] .- LMSs[[2,4].+16,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dls)),nbins=30,xlabel="LS Distance",leg=false)

dms = map(i->sqrt(sum((LMSs[[2,4].+8,i] .- LMSs[[2,4].+16,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dms)),nbins=30,xlabel="MS Distance",leg=false)

savefig("test.svg")








Y = replace(hcat(filter(i->!isempty(skipmissing(i[9:end])),X)...)[9:end,:],missing=>0)
Yg = crdf.ConeResponseType

Z = Y[8:8:end,:]

Plots.scatter(360(Zlm[1,:].-Zlm[2,:]),ones(20),nbins=30,projection=:polar)
