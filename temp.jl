

ct = YAML.load_file("cell type.yaml")
for u in keys(ct)
    tdir = "Z:\\AF5\\AF5_HLV1_ODL3\\AF5_HLV1_ODL3_Color_1"
    nm = "Single-Unit_$(u)_HueAngleTuning"
    fs = matchfile(Regex("$nm.png"),dir=tdir,adddir=true)
    if !isempty(fs)
       cp(fs[1],joinpath(tcdir,"$nm.png"),force=true)
    end
end





imagepatch(img,10,(0.5,0.5))



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

gr()

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


using NeuroAnalysis


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


batchon = range(start=0,step=1500,length=800)
epochs = ref2sync(batchon.+epoch,dataset,ii)

ys = fill2mask(epochsamplenp(mmbinfile,fs,epochs.+t0,1:nch;meta=[],bandpass=[300,3000],whiten=nothing),exchmask,chmap=chmapraw)
ys = resample(ys,1//3,dims=3)
fs /= 3

ys = fill2mask(epochsamplenp(mmbinfile,fs,epochs,1:nch;meta,bandpass=[300,3000],whiten=nothing),exchmask,chmap=chmapraw)
ys = resample(ys,1//3,dims=3)
fs /= 3



cc = mapreduce(i->cor(pys[:,:,i],dims=2),(i,j)->cat(i,j,dims=3),1:size(pys,3))
@manipulate for i in 1:size(cc,3)
    # plotanalog(c[:,:,f],x=depths,y=depths,color=:turbo)
    heatmap(cc[:,:,i],color=:turbo)
end

lcc = localcoherence(cc,s=1.4)
plotanalog(lcc,y=depths,color=:turbo)

plot!(1e3mean(lcc,dims=2),depths)



pys = ys[:,1,:,1:2:end]
ps,freqs = powerspectrum(pys,fs,freqrange=[300,3000])

using Interact

@manipulate for i in 1:size(pys,3)
    @views plotanalog(ps[:,:,i],x=freqs,y=depths,color=:turbo,n=mean(ps[:,:,i],dims=2))
end



pss = dropdims(mean(ps,dims=2),dims=2)
plotanalog(pss;hy,color=:heat,n=mean(pss,dims=2),xlabel="",cunit=:fr)







pys = ys[:,1,:,1:2:end]

cs,freqs = coherencespectrum(pys,fs,freqrange=[300,3000])

pss=wpc[4][:,1:2:200]
pbss=aps[4]["wpbc"][:,1:2:200]
ds = depths[1]
pdc = abs.(pss.-pbss)

plotanalog(pss,x=1:100,y=ds,xlabel="Trials",color=:heat,cunit=:fr,n=mean(pss,dims=2))
plotanalog(pbss,x=1:100,y=ds,xlabel="Trials",color=:heat,cunit=:fr,n=mean(pbss,dims=2))

plotanalog(pdc,x=1:100,y=ds,xlabel="",color=:heat,cunit=:fr,n=mean(pdc,dims=2))
plotanalog(pdc,x=1:100,y=ds,xlabel="",n=mean(pdc,dims=2))





@manipulate for f in 1:length(freqs)
    plotanalog(cs[:,:,f,1],x=depths,y=depths,xlabel="Depth",xunit="μm",color=:heat,cunit=:fr,aspectratio=:equal)
    # heatmap(c[:,:,f],color=:turbo)
end

plotanalog(cs[:,:,43,1],x=depths,y=depths,xlabel="Depth",xunit="μm",color=:heat,cunit=:fr,aspectratio=:equal)
savefig("test.png")
freqs[43]


css = dropdims(mean(cs,dims=3),dims=3)
@manipulate for f in 1:size(css,3)
    plotanalog(css[:,:,f],x=depths,y=depths,xlabel="Depth",xunit="μm",color=:heat,cunit=:fr,aspectratio=:equal)
    # heatmap(c[:,:,f],color=:turbo)
end
plotanalog(css[:,:,43],x=depths,y=depths,xlabel="Depth",xunit="μm",color=:heat,cunit=:fr,aspectratio=:equal)


lcss = mapreduce(i->bandmean(css[:,:,i]),hcat,1:size(css,3))

plotanalog(lcss,x=freqs,y=depths,xlabel="Frequency",xunit="Hz",color=:heat,cunit=:fr,n=mean(lcss,dims=2))

plotanalog(lcss,y=depths,xlabel="",color=:heat,cunit=:fr,n=mean(lcss,dims=2))





##

x=-4:0.1:4
tv = 100*dogf.(x)
ob = 100*dogf.(x) .+ 10*randn(length(x))

chi2count(x)= round.(Int,x.-minimum(x))
function chi2p(x)
    t = x.-minimum(x)
    t ./= sum(t)
end

pvalue(ChisqTest(chi2count(ob),chi2p(tv)))
pvalue(ChisqTest(chi2count(ob)))


Plots.plot(chi2count(ob))
Plots.plot(chi2p(tv))

Plots.plot(ob)
Plots.plot(tv)


pvalue(ApproximateTwoSampleKSTest(tv,ob))


savefig("test.png")
save("test.png",p)

savefig("test.svg")
##
using MultivariateStats

tt = fit(ICA,cat(values(cdrms)...,dims=2),7,maxiter=200)

plot(tt.mean)

heatmap(first(values(cdrms)),color=:coolwarm)

plot!(abs.(tt.W*1e4),0:191,leg=false)

pp = cat(values(cdrms)...,dims=2)

heatmap(pp)

pyplot()
gr()

F2phase01e = adjust_histogram(F2phase01, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.8))
        clamp01!(ori)


F2phase01e = adjust_histogram(F2phase01, LinearStretching())

bglp = imfilter(F2phase01,Kernel.gaussian(50))

tt = clampscale(F2phase01.-bglp)
tte = adjust_histogram(tt, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.8))
tte = adjust_histogram(tt, LinearStretching())

tt = als(F2phase01,rb=15,cb=15)
Gray.(F2mag01e)
Gray.(F2phase01e)


function takeroi(imgsize,ppd;roi=missing,roimaxresize=32,issquare=false)
    if ismissing(roi)
        idxrange = map(i->1:i,imgsize)
        roisize = imgsize
    else
        cs = round.(Int,roi.centerdeg*ppd)
        # rs = round.(Int,roi.radiideg*ppd)
        r = round(Int,roi.radiusdeg*ppd)
        rs = (r,r)
        cs,rs,r = clamproi(cs,rs,imgsize;issquare)
        idxrange = map((c,r)->(-r:r).+c,cs,rs)
        roisize = length.(idxrange)
        roi = (centerdeg=cs./ppd,radiideg=rs./ppd,radiusdeg=r/ppd)
    end
    d = max(roisize...)
    f = min(roimaxresize,d)/d
    roiresize = round.(Int,roisize.*f)
    (;idxrange,roiresize,roi)
end





## adaptive linear straching

function als(img;rb=20,cb=20,pp=(1,99))
    h,w = size(img)
    rn = round(Int,h/rb)
    cn = round(Int,w/cb)
    rs = [intersect((1:rn).+i*rn,1:h) for i in 0:rb-1]
    cs = [intersect((1:cn).+i*cn,1:w) for i in 0:cb-1]
    m = similar(img)
    for i in rs, j in cs
        b = @views img[i,j]
        vb = vec(b)
        m[i,j]=adjust_histogram(b, LinearStretching(src_minval=percentile(vb,pp[1]),src_maxval=percentile(vb,pp[2])))
    end
    m
end



ff = dogfilter(F1phase01)
t = als(ff;rb=20,cb=20)

map(a->HSV(360a,1,1),t)

t=mapwindow(x->clampscale(x[50,50],extrema(x)...),F2phase01,(101,101))
Gray.(t)
map(i->cm[i],t)


tt = ahe(ff)


Gray.(tt)

mtt = clampscale(dogfilter(F2mag01),3)
Gray.(mtt)

Gray.(ff)


Gray.(F1phase01_dog)


ColorMaps["dkl_mcchue_l0"].colors

cm.[F1phase01_dog]

##  
F1mag01 = clampscale(F1mag,(3,97))
F1mag01 = clampscale(dogfilter(F1mag),3)
histogram(vec(F1mag))


histogram(vec(F2phase01))
histogram(vec(dogfilter(F2phase01)))
histogram(vec(pr))





ct = map(i->exp(im*2π*i),F2phase01)

knl = Kernel.DoG((3,3),(25,25),(151,151))
k=collect(knl)
ctt = imfilter(ct,knl)
ctt = mapwindow(r->sum(r.*k),ct,(151,151))

ft = mod2pi.(angle.(ctt)) ./ (2π)

Gray.(F2phase01)
Gray.(F1phase01)
Gray.(ft)
lft = als(F1phase01)
Gray.(lft)


t = clampscale(dogfilter(F1phase01),3)
lt = als(dogfilter(F1phase01))
histogram(vec(lt))
Gray.(lt)

map(a->HSV(360a,1,1),F2phase01)
map(a->HSV(360a,1,1),lft)
map(a->HSV(360a,1,1),t)


Gray.(F2mag01)

leftF2mag = F2mag
leftF1mag = F1mag

t=log2.(leftF1mag./F1mag)

tt=clampscale(t)

Gray.(tt)

Gray.(F2mag01)



