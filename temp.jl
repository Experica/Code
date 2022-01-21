using NeuroAnalysis,Makie,Colors,Plots,GeometryTypes

x = 1:10
y=1:10
u=fill(0.5,10)
v = fill(0.5,10)
arrows(x,y,u,v)





OriMarker = Shape([(0.5,0.2),(-0.5,0.2),(-0.5,-0.2),(0.5,-0.2)])



a=[0.5 -0.5 -0.5 0.5;
    0.2 0.2 -0.2 -0.2]
Shape(a)





unitposition.+[-1 0]

LineSegment((-1,1))

oris=rand(0:360,20)
unitposition = [ rand(20:70,20) rand(20:70,20)]
unitgood=rand(Bool,20)


Vec2f0[(0.5,0.2),(-0.5,0.2),(-0.5,-0.2),(0.5,-0.2)]



m = Point2f0[(0,0), (0, 1), (0.5, 0.5), (1, 1), (1, 0)]
meshscatter(rand(1:7,20),rand(0:3,20), marker=GLPlainMesh(m),rotations=rand(0:2π,20))



AbstractPlotting.set_theme!(
    plot = (show_axis = false, scale_plot = false),
    color = :turquoise1
)
poly(HyperRectangle(Vec2f0(0), Vec2f0(1,100)))



t = load("C:\\Users\\fff00\\Command\\Recordings\\0040.png")

tt=t[200:600,200:600]

gt=Gray.(tt)

ggt=gray.(gt)

plot(ggt[:,200])

using Makie,Plots,CairoMakie,GLMakie

AbstractPlotting.inline!(false)

p=Makie.scatter(1:2,1:2,marker=rand(RGB{Float32},10,10),markersize=(50,10))

Makie.scatter(1:2,1:2,marker=[rand(RGB{Float32},10,10), rand(RGB{Float32},20,10)],markersize=[(50,10),(10,50)])



imgs = [rand(RGB{Float32},20,10), rand(RGB{Float32},10,20)]
Makie.scatter(1:2,1:2,marker=imgs,markersize=size.(imgs))

CairoMakie.activate!()
GLMakie.activate!()

dataset = prepare("Y:\\AF5\\AF5_HLV1_ODL1_HartleySubspace_1.mat")
dataset = prepare("Y:\\AF5\\AF5_HLV1_ODL3_Flash2Color_2.mat")




using Interact,Plots,VegaLite,DataFrames,Distances,LinearAlgebra


using Interact,Plots
@manipulate for i in 1:10
    heatmap(rand(10,10),color=:temperaturemap)
end



Makie.surface(x,y,z,colormap=:coolwarm,shading=false)

sufstr=["_A","_L","_M","_S"];tcdir="E:\\Code\\temp"
for u in sort(collect(keys(dataset["ulsta"])))
    for t in [1,2,3,4]
    tdir = "Z:\\AF5\\AF5_HLV1_ODL1\\AF5_HLV1_ODL1_OriSF_$t"
    nm = "Single-Unit_$(u)_SpatialFreqTuning"
    fs = matchfile(Regex("$nm.png"),dir=tdir,adddir=true)
    if !isempty(fs)
   cp(fs[1],joinpath(tcdir,"$nm$(sufstr[t]).png"),force=true)
end
end
end

ct = Dict()
for u in sort(collect(keys(dataset["ulsta"])))
    # push!(ct,u=>"DO")
    ct[u]="DO"
end

ct

YAML.write_file("temp/cell type1.yaml",ct)

using YAML

us = collect(keys(ct))
up = up[si,:]
us = us[si]

si = sortperm(up[:,1])
jit = [cumsum(range(0.0000001,length=58,step=0.0000002)) .+ rand(0:0.00002:1,58) ones(58)]

plotunitposition(up.*jit,unitid=us,layer=layer,unitgood=trues(58),markersize=0,unitidsize=7)


savefig("unitposition.svg")



ct = YAML.load_file("cell type.yaml")
for u in keys(ct)
    tdir = "Z:\\AF5\\AF5_HLV1_ODL3\\AF5_HLV1_ODL3_Color_1"
    nm = "Single-Unit_$(u)_HueAngleTuning"
    fs = matchfile(Regex("$nm.png"),dir=tdir,adddir=true)
    if !isempty(fs)
       cp(fs[1],joinpath(tcdir,"$nm.png"),force=true)
    end
end









using NeuroAnalysis,Test,Plots,Interact,Statistics, BenchmarkTools,StatsPlots,DataFrames

t = prepare("Y:\\AF5\\AF5_HLV1_ODL3_Color_1.mat")


twf = t["spike_kilosort"]["templateswaveformfeature"]
twff = t["spike_kilosort"]["clusterwaveformfeature"]

twff = t["spike"]["unitfeature"]

sui = t["spike"]["unitgood"]

fns = collect(keys(twff))

@manipulate for x in fns, y in fns
    scatter(twff[x][sui],twff[y][sui],markersize=3,leg=false)
end

tt=t["lf"]["meta"]
tt["from"]

tt["roapgain"]
tt["savedchans"]
tt["snsSaveChanSubset"]

@manipulate for i = 1:size(twf, 1)
    ps, freq = powerspectrum(twf[i:i, :], 30000, freqrange = [3000, 12000])
    p = plot(layout=(4,1),size=(800,200*4))
    w = twf[i, :]
    plot!(subplot=1,w)
    plot!(subplot=2,freq, ps[1,:],ylims=[0,1e-5])
    fw = hlpass(twf[i:i,:],30000,high=150,low=12000)[1,:]
    plot!(subplot=3,fw)
    plot!(subplot=4,w.-fw,ylims=[-1e0,1e0])
end


rt = matread("ttt.mat")

save("ttt.jld2",t)


load("ttt.jld2","cg")

a=[1,2,3]

a[begin:begin]
a=rand(3,4)


a[begin:begin,:]

b = a[:,argmax(a,dims=2)]

[i[1]*i[2] for i in CartesianIndices(a)]




imagepatch(img,10,(0.5,0.5))
## Tuning Properties in layers
layer = load(joinpath(siteresultdir,"layer.jld2"),"layer")
# testids = ["$(siteid)_$(lpad(i,3,'0'))" for i in [8,12,13,14]]

testids = ["$(siteid)_OriSF_$i" for i in 1:5]
testn=length(testids)
ds = load.(joinpath.(siteresultdir,testids,"factorresponse.jld2"))
testtitles = ["$(i["color"])" for i in ds]
spikes = load.(joinpath.(siteresultdir,testids,"spike.jld2"),"spike")

f = :Ori
p=plot(layout=(1,testn),link=:all,legend=false,grid=false,xlims=(10,60))
for i in 1:testn
    vi = spikes[i]["unitgood"].&ds[i]["responsive"].&ds[i]["modulative"]
    if f in [:Ori,:Ori_Final]
        color = map((i,j)->j ? HSV(i.oo,1-i.ocv/1.5,1) : RGBA(0.5,0.5,0.5,0.1),ds[i]["factorstats"][:Ori_Final],vi)
    elseif f==:Dir
        color = map((i,j)->j ? HSV(i.od,1-i.dcv/1.5,1) : RGBA(0.5,0.5,0.5,0.1),ds[i]["factorstats"][:Ori_Final],vi)
    end
    scatter!(p,subplot=i,spikes[i]["unitposition"][:,1],spikes[i]["unitposition"][:,2],color=color,markerstrokewidth=0,markersize=2,title=testtitles[i],titlefontsize=8)
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(17,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
p
foreach(i->savefig(joinpath(siteresultdir,"Layer_UnitPosition_$(f)_Tuning$i")),[".png",".svg"])



if plot
    vi = uresponsive.&umodulative.&unitgood.&uenoughresponse
    pfactor=intersect([:Ori,:Ori_Final],keys(ufs))
    if length(pfactor)>0
        pfactor=pfactor[1]
        colorspace = ex["Param"]["ColorSpace"]
        hues = ex["Param"]["Color"]
        title = "UnitPosition_$(colorspace)_$(hues)_OptOriDirSF"
        p=plotunitpositionproperty(unitposition[vi,:],title=title,ori=map(i->i.oo,ufs[pfactor][vi]),os=map(i->1-i.ocv,ufs[pfactor][vi]),
        dir = map(i->i.od,ufs[pfactor][vi]),ds=map(i->1-i.dcv,ufs[pfactor][vi]),sf=map(i->i.osf,ufs[:SpatialFreq][vi]),layer=layer)
        foreach(i->save(joinpath(resultdir,"$title$i"),p),[".png",".svg"])

        oos = map(i->i.oo,ufs[pfactor][vi])
        oos = [oos;oos.+180]
        ys,ns,ws,is = histrv(oos,0,360,nbins=20)
        oa = deg2rad.(mean.([ws;ws[1]]))
        oh = [ns;ns[1]]./2

        ys,ns,ws,is = histrv(map(i->i.od,ufs[pfactor][vi]),0,360,nbins=20)
        da = deg2rad.(mean.([ws;ws[1]]))
        dh = [ns;ns[1]]

        title = "$(colorspace)_$(hues)_OptOriDirHist"
        Plots.plot([oa da],[oh dh],title=title,projection=:polar,linewidth=2,label=["Ori" "Dir"])
        foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
    end
    if haskey(ufs,:ColorID)
        plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.oh,1-i.hcv,1) : HSVA(0,0,0.5,0.2),ufs[:ColorID],umodulative))
        foreach(i->savefig(joinpath(resultdir,"UnitPosition_Hue_Tuning$i")),[".png",".svg"])
    end
end


## Test
ys = reshape2mask(epochsamplenp(mmlf,fs,[condon[1] condon[end]+1000],1:nch,meta=lfmeta,bandpass=[]),exchmask)

dy = 1e6*ys[:,1,:]
dy=csd(1e6*ys[:,1,:],h=hy)


ps,freq=powerspectrum(1e6dy,fs,freqrange=[0,10])

t=periodogram(dy[100,:];fs)
t=welch_pgram(dy[100,:];fs)
t=mt_pgram(dy[100,:];fs,nw=10)

plot(t.freq,t.power)

plotanalog(log10.(amp),x=freq1,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=:auto,color=:vik)
plotanalog(pha,x=freq1,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=:auto,color=:vik)


t=rfft(dy,2)
freq = rfftfreq(size(dy,2),fs)


t1=t[:,freq.<2]
freq1=freq[freq.<2]
amp=abs.(t1)
pha = angle.(t1)

f1i=argmin(abs.(freq.-0.33))
f2i=argmin(abs.(freq.-0.66))
f3i=argmin(abs.(freq.-0.99))

plot(depths,log10.(amp[:,f1i]))
plot(depths,log10.(amp[:,f2i]))
plot(depths,log10.(amp[:,f3i]))

plot(depths,pha[:,f1i])
plot(depths,pha[:,f2i])
plot(depths,pha[:,f3i])



## getnerate DKL plane
using MAT,ColorLab

fn = "StaticAdelson"
bodata = matread("C:\\Users\\fff00\\DownLoads\\$fn.mat")

imgs = map(i->(bodata["data"][:,:,i].+1)./2,1:size(bodata["data"],3))

saveunityrawtexture("C:\\Users\\fff00\\DownLoads\\$fn",imgs)

heatmap(bodata["data"][:,:,1])
heatmap(imgs[1])


## cell response
function tresponseness(mat,bindex)
    b = vec(mat[:,bindex])
    ht = [@views UnequalVarianceTTest(b,mat[:,i]) for i = 1:size(mat,2)] # Welch's t-test
    t = map(i->i.t,ht)
    replace!(t,NaN=>0)
    return (;m=abs.(t))
end

pnrow*hy

mfun = x->abs.(x.-mean(x[baseindex]))
uw = replace(unitgood,0=>1.5)
acmdepthpsth = Dict(condstring(r)=>spacepsth(map(pm->(;vmeanse(pm.mat[r.i,:];mfun).m,pm.x),unitepochpsth),unitposition,w=uw,lim=(0,3840),bw=100,step=50) for r in eachrow(cond))


acmdepthpsth = Dict(condstring(r)=>spacepsth(map(pm->(;tresponseness(pm.mat[r.i,:],baseindex).m,pm.x),unitepochpsth),unitposition,w=uw,lim=(0,3840),bw=100,step=50) for r in eachrow(cond))


dr = first(values(acmdepthpsth))
plotanalog(dr.psth)



ks = collect(keys(cmdepthpsth))
k=ks[2]


## lfp plot
using StatsBase,Plots
clfp = ys[:,1,:,4]

clfp = first(values(cmcys))
cmcysse = Dict(condstring(r)=>dropdims(std(ys[:,c,:,r.i],dims=3)/sqrt(length(r.i)),dims=3) for r in eachrow(cond))
clfpse = first(values(cmcysse))


offset = 0.00001
y = [clfp[i,j]+offset*(i-1) for i=1:192,j=1:374]

plot(y',leg=false,frame=:grid,grid=false,color=:black,size=(600,750),fillrange=range(0,length=192,step=offset)',fillalpha=0.03,
 ylabel="Depth (μm)",xlabel="Time (ms)",xticks=[],yticks=[],left_margin=4mm,
 annotations=[(1,-1e-4,Plots.text("0",10,:gray20,:center)),
            (374,-1e-4,Plots.text("150",10,:gray20,:center)),
            (-5,-0.5e-4,Plots.text("0",10,:gray20,:right)),
            (-5,1.9e-3,Plots.text("3840",10,:gray20,:right))])



plot(y',ribbon=clfpse',leg=false,frame=:grid,grid=false,color=:black,size=(600,750),fillalpha=0.2,
 ylabel="Depth (μm)",xlabel="Time (ms)",xticks=[],yticks=[],left_margin=4mm,
 annotations=[(1,-1e-4,Plots.text("0",10,:gray20,:center)),
            (374,-1e-4,Plots.text("150",10,:gray20,:center)),
            (-5,-0.5e-4,Plots.text("0",10,:gray20,:right)),
            (-5,1.9e-3,Plots.text("3840",10,:gray20,:right))])









savefig("test.svg")


cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],:CSD,h=hy),dims=3),dims=3) for r in eachrow(cond))

lfp = first(values(cmcys))
csdcsd = first(values(cmccsd))
csdlfp = csd(lfp,h=hy)
iclll = csd(lll,h=hy)

nlfp = stfilter(lfp,temporaltype=:sub,ti=baseindex)
ncsdlfp = stfilter(csdlfp,temporaltype=:sub,ti=baseindex)
ncsdlfp = csd(nlfp,h=hy)

plotanalog(lfp)
plotanalog(nlfp)
plotanalog(csdcsd)
plotanalog(csdlfp)
plotanalog(iclll)


plotanalog((iclll.-clll)[2:end-1,:])

extrema((iclll.-clll)[2:end-1,:])


## kCSD
using PyCall

d = permutedims(permutedims(collect(depths)))
t =  permutedims(permutedims((collect((0:373)/fs))))

kcsd = pyimport("kcsd")


k = kcsd.KCSD1D(d,lfp,
    h=hy,sigma=0.3,n_src_init=20000,
    gdx=20,
    src_type="gauss",R_init=1,lambd=0)

k.cross_validate(Rs=np.linspace(0.01, 0.15, 15))
rs = permutedims(permutedims(collect(0.01:0.01:0.15)))
rs = collect(0.01:0.01:0.15)
k.cross_validate(Rs=rs)



plotanalog(    k.values("CSD"))



ktt = zeros(192,374)
for i=1:374
k = kcsd.KCSD1D(d,lfp[:,i:i],
    h=hy,sigma=0.3,n_src_init=1000,
    src_type="gauss",R_init=1,lambd=0)
    k.values("CSD")
ktt[:,i] = k.values("CSD")
end
## gpcsd
gcsd = pyimport("gpcsd.gpcsd1d")
g = gcsd.GPCSD1D(ys[:,1,:,1:2:end],d,t)

g.fit(n_restarts=5)
gtt = g.sample_prior(100)
gtt = g.predict(d,t)

plotanalog(dropdims(mean(gtt,dims=3),dims=3))

##



Gray.(clampscale(responses[:,:,4],1.5))
maps = complexmap(responses,angles,filter=nothing,presdfactor=1.5,sufsdfactor=nothing)

## od

lod = load("Z:\\AF9\\AF9_V1V2_Full\\AF9_V1V2_Full_ISICycleOri_4\\isi.jld2")
rod = load("Z:\\AF9\\AF9_V1V2_Full\\AF9_V1V2_Full_ISICycleOri_3\\isi.jld2")
coneiso = load("Z:\\AF9\\AF9_V1V2_Full\\AF9_V1V2_Full_ISICycle2Color_0\\isi.jld2")

f1 = coneiso["F1"]
f1p = angle.(f1)

Gray.(clampscale(f1p,-π,π))
Gray.(clampscale(abs.(f1p),0,π))
Gray.(clampscale(-abs.(f1p),3))
isoimg = adjust_histogram(f1p, AdaptiveEqualization(nbins = 256,
minval=-π,maxval=π, rblocks = 8, cblocks = 8, clip = 0.1))

isoimg = adjust_histogram(abs.(f1p), AdaptiveEqualization(nbins = 64,
minval=0,maxval=π, rblocks = 4, cblocks = 4, clip = 0.5))

isoimg = adjust_histogram(clampscale(f1p,1.5), AdaptiveEqualization(nbins = 64,
minval=0,maxval=1, rblocks = 4, cblocks = 4, clip = 0.04))

Gray.(clampscale(isoimg))
Gray.(f1p./(2π))

f1pp = copy(f1p)
f1pp[f1pp.>π] = 2π .- f1pp[f1pp.>π]
Gray.(f1pp./π)

f1ppimg = adjust_histogram(f1pp, AdaptiveEqualization(nbins = 64,
minval=0,maxval=π, rblocks = 4, cblocks = 4, clip = 0.7))

Gray.(clampscale(f1ppimg))

od1 = lod["F1mag"] .- rod["F1mag"]
od2 = lod["F2mag"] .- rod["F2mag"]


od1 = clampscale(lod["F1mag"],1.5) .- clampscale(rod["F1mag"],1.5)
od2 = clampscale(lod["F2mag"]) .- clampscale(rod["F2mag"])

histogram(vec(od1))

Gray.(clampscale(od1,1.5))

Gray.(clampscale(od2,1.5))

img = adjust_histogram(od1, AdaptiveEqualization(nbins = 256, rblocks = 4, cblocks = 4, clip = 0.2))
Gray.(clampscale(img))

## delay
ncycle=9
responsedelay = 500





using  Combinatorics


a=[1,2,3]

filter!(i->issorted(i)&&all(j->j==1,diff(i)),collect(permutations(1:3,3)))



function segmentpermutations(r,l::Integer)
    n = length(r)
    l<1 && error("Segment length smaller than 1.")
    l>n && error("Segment length larger than whole data length.")
    s = 1:l
    [s.+i for i in 0:(n-l)]
end


segmentpermutations(1:10,1)

function osp(ns,l)
    ci = filter!(i->issorted(i)&&all(j->j==1,diff(i)),collect(permutations(ns,l)))
end

osp(1:10,9)

fis = [epoch2sampleindex([c[1]-1 1000*c[end]/modulatefreq].+(preicidur+responsedelay),framerate,maxsampleindex=nframe) for c in osp(1:10,1)]

fis=[]
append!(fis, [epoch2sampleindex([c[begin]-1 c[end]].*1000 ./modulatefreq.+(preicidur+responsedelay),framerate,maxsampleindex=nframe) for c in segmentpermutations(1:10,1)])

frameindex=fis[5]

Gray.(F1polarityce)


ifs = [dft_imager(imagefile[fi],w,h,framerate,baseimg,modulatefreq) for fi in fis]

ips = map(i->angle.(i[1]),ifs)
ims = map(i->abs.(i[1]),ifs)

rightips = ips
rightims = ims

leftips = ips
leftims = ims

F1ps = foreach(i->begin
    F1phase = angle.(ifs[i][1])
    F1polarity = (π .- abs.(F1phase)) ./ π
    # filter
    F1polarityce = adjust_histogram(F1polarity, AdaptiveEqualization(nbins = 256, rblocks = 8, cblocks = 8, clip = 0.1))
    F1polarityce = clampscale(F1polarityce,0,1)
    foreach(ext->save(joinpath(resultdir,"F1Polarity_ContrastEnhanced_$i$ext"),F1polarityce),figfmt)
end,eachindex(ifs))



tt = [circmean(map(p->p[i,j],ips)) for i=1:h,j=1:w]

tt = [circmean(map(p->p[i,j],ips),map(p->p[i,j],ims)) for i=1:h,j=1:w]


tt = sum(first.(ifs))

F1phase = angle.(tt)





t = [UnequalVarianceTTest(map(p->p[i,j],leftims),map(p->p[i,j],rightims)).t for i = 1:h,j=1:w]
t[isnan.(t)].=0

Gray.(clampscale(t,1.5))


##

ori = clampscale(F2phase01,2)

ori = imfilter(F2phase01,Kernel.gaussian(5))
ori = adjust_histogram(F2phase01, LinearStretching())
clamp01!(ori)

Gray.(ori)

map(a->HSV(360a,1,1),ori)

map((a,m)->HSV(360a,1,0.5+m/2),F2phase01,clampscale( F2mag01))


t = F2phase01[1200:1300,1200:1300]
tt = adjust_histogram(t, LinearStretching())
Gray.(tt)


extrema(t)
extrema(tt)




## ap plot
yy = pys[:,:,2]
yn,xn = size(yy)
offset = 7e-5
y = yy .+ range(start=0,step=offset,length=yn)


plot(y',leg=false,frame=:grid,grid=false,color=:black,size=(600,750),
 ylabel="Depth (μm)",xlabel="Time (ms)",xticks=[],yticks=[],left_margin=4mm,
 annotations=[(1,-4offset,Plots.text("0",10,:gray20,:center)),
            (xn,-4offset,Plots.text("150",10,:gray20,:center)),
            (-20,0,Plots.text("0",10,:gray20,:right)),
            (-20,(yn-1)*offset,Plots.text("3840",10,:gray20,:right))])

savefig("test.svg")


a=unitepochpsth[1]

typeof(a.x)

map(p->(;vmeanse(p.mat[1:2:end,:])...,p.x),unitepochpsth)


##
