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



function plotunitposition(unitposition;unitgood=[],title="Unit_Position",color=nothing,alpha=0.4,marker=:hline,rot=range(0,π,length=20))
    if isnothing(color)
        if !isempty(unitgood)
            color = map(i->i ? :darkgreen : :gray30,unitgood)
        else
            color = :gray30
        end
        if !isnothing(alpha)
            color = coloralpha.(parse.(RGB,color),alpha)
        end
    end

    markshape = HyperRectangle(Vec2f0(0),Vec2f0(1,60))
    markshape = Vec2f0[(0.5,0.2),(-0.5,0.2),(-0.5,-0.2),(0.5,-0.2)]
    markshape = GLPlainMesh(Point2f0[(0,0), (0, 1), (0.5, 0.5), (1, 1), (1, 0)])
    p = Scene()
    Makie.meshscatter!(p,unitposition[:,1],unitposition[:,2],scale_plot=true,color=color,marker=markshape,rotations=rot,markersize=10,strokewidth=0,linewidth=0)

    t1= unitposition.+[-1 0]
    t2 = unitposition.+[1 0]
    pos = zeros(2*size(t1,1),size(t1,2))
    pos[1:2:end,:]=t1
    pos[2:2:end,:]=t2


    # Makie.linesegments!(p,pos[:,1],pos[:,2],color=:black,rotations=rot,linewidth=2)
    # axis = p.Axis
    # axis.names.axisnames=("Position_X (μm)","Position_Y (μm)")
    return Makie.title(p,title)
end

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





using JLD2



using Interact,Plots,VegaLite,DataFrames,Distances,LinearAlgebra
@vlplot(:point,rand(10),rand(10))

Dict(Pair.(1:3,2:4))

using Interact,Plots
@manipulate for i in 1:10
    heatmap(rand(10,10),color=:temperaturemap)
end


x = DataFrame(a=[missing,missing,['A','L',missing]])
save("test.jld2","x",x)

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



t = minmaxcolormap(:coolwarm,0,1)
t = minmaxcolorgradient(RGBA(0,0,0,1),RGBA(1,1,1,1))

save("ttt.yaml",t)

fieldnames(ColorGradient)

propertynames(t)

rt = YAML.load_file("ttt.yaml")


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





# Variables for testing
x0 = "foo"
x1 = rand(3000)
x2 = rand(3000,1)
x3 = rand(1,3000)
x4 = rand(1,1)
x5 = rand(3000,1,30,1,3)
x6 = [x0, x1, x2, x3, x4, x5]
x7 = fill(x4,200,20,30)
x8 = Dict(:a=>x4,:b=>x4)
x9 = Any[x4,x4]
x = Dict(:x0=>x0, :x1=>x1, :x2=>x2, :x3=>x3, :x4=>x4, :x5=>x5) # Dict of arrays and strings
y = reshape([x, x0, x1, x6], (1,4))                            # Array of dicts and arrays and strings
z = Dict(:y=>y, :x=>x, :x6=>x6, :x7=>x7, :x8=>x8, :x9=>x9)                                # Dict of dicts and arrays and strings

# Some tests to check that all the arrays are squeezed recursively
["$k: size conversion: $(size(v)) → $(size(f(v))))" for (k,v) in x if v isa Array]
f(z)     # is a vector
f(z)[:y][1][:x5] # is a (3,3,3) array


@btime f($(z))
@btime f2($(z))
@btime f3($(z))
@btime dropmatdim!($(z))

f2(x7)


any((<:).(eltype(x7),[Any,Array,Dict]))

x7[1] :: Array

x7 isa Array

dropmatdim!(z)






using MAT
d=matread("c:\\users\\fff00\\mattest.mat")




d=readmat("c:\\users\\fff00\\mattest.mat")



using Unitful

a=3.6u"ms"

ustrip(datatimeunit,edur)

a=[2,4.5,7,9.3]

b=a*u"s"


edur = 150u"ms"

datatimeunit = u"ms"

uconvert.(u"ms",b)


a=rand(10)

e=[0u"ms" 150u"ms"]

a=[[1.0,2.0],[3.0,4.0]]

typeof(a) <: Vector{Vector}

b=[1.0,2.0]

typeof(b) <: Vector



## test
scatter(rand(10),group=rand(1:2,10),color=[:coolwarm :reds])


heatmap(randn(30,30),color=:grays,frame=:none,leg=false,aspect_ratio=1)

savefig("test.png")



x=collect(-10:0.1:10)
mfun(x,p) = p[1]*sin.(p[2]*x .+ p[3])
y=mfun(x,[5,2,1]) .+ randn(length(x))

mfun(x,p) = vmf.(x,β=p[1],μ=p[2],κ=p[3],n=2)
y=mfun(x,[5,0,1]) .+ randn(length(x))

plot(x,[y mfun(x,mfit.param)])

mfit = curve_fit(mfun,x,y,[2,1.0,0],x_tol=1e-11)



mfit = curve_fit(mfun,x,y,1 ./ mfit.resid.^2,[2,1.0,0])



## temp
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


## optim


data = dataset["ulsta"][89][:,:,20,2]
function dataxy(data,ppd)
    rpx = (size(data)[1]-1)/2
    r = rpx/ppd
    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-(rpx+1) rpx+1])/ppd
    y = vec(data)
    return x,y
end

x,y = dataxy(data,45)

function objfun(p;x=x,y=y)
    model = rfgabor.(x[:,1],x[:,2],p...)
    sum((y.-model).^2)
end
lb = [1000.0, -0.4,  0.1,  -0.4,  0.1,  0,  0.1,    0]
ub = [2000.0,  0.4,  0.4,   0.4,  0.4,  π,   8,  1]

llb = [1.0,     -1,  0.1,  -1,  0.1,  0,  0.2,    0]
uub = [10000.0,  1,  0.6,   1,  0.6,  π,   8,  1]


res = optimize(objfun,llb,uub,[1500.0,0.0,0.2,0,0.2,1.4,1,0],SAMIN(),Optim.Options(iterations=10^5))

pyfit = (model=:gabor,param=res.minimizer)
heatmap(rfimage(pyfit,0.82),aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)


yp = rfgabor.(x[:,1],x[:,2],res.minimizer...)



rcopy(R"cor.test($y,$yp,alternative='greater',method='pearson')")


## rffit


lc = localcontrast(data,ws=0.5,ppd=45)
roi=peakroi(lc)
rg = zeros(size(lc))
rg[roi.i].=1
heatmap(data,aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)
heatmap(lc,aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)
heatmap(rg,aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)



frlt = modelfit(data,45)
heatmap(rfimage(frlt,frlt.radius),aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)


using PyCall
lmfit = pyimport("lmfit")


pyfit = (model=:gabor,param=py"fparams")
heatmap(rfimage(pyfit,0.82),aspect_ratio=1,frame=:none,color=:coolwarm,yflip=true)


py"""
x = $(x[:,1])
y = $(x[:,2])
data = $y

from numpy import power,exp, sin,cos,pi

from lmfit import minimize, Parameters


def residual(params, x, y, data):
    a = params['a']
    miux = params['miux']
    sigmax = params['sigmax']
    miuy = params['miuy']
    sigmay = params['sigmay']
    ori = params['ori']
    sf = params['sf']
    phase = params['phase']

    sinv = sin(ori)
    cosv = cos(ori)
    xx = cosv * x + sinv * y
    yy = cosv * y - sinv * x
    model = a*exp(-0.5*(power((xx-miux)/sigmax,2) + power((yy-miuy)/sigmay,2))) * sin(2*pi*(sf * yy + phase))

    return data-model

r=0.8
params = Parameters()
params.add('a', value=1500.0,min=1000)
params.add('miux', value=0.0,min=-0.3*r,max=0.3*r)
params.add('sigmax', value=0.2*r,min=0.1*r,max=0.3*r)
params.add('miuy', value=0.0,min=-0.3*r,max=0.3*r)
params.add('sigmay', value=0.2*r,min=0.1*r,max=0.3*r)
params.add('ori', value=1.4,min=0,max=pi)
params.add('sf', value=1.0,min=1,max=3)
params.add('phase', value=0.0,min=0,max=1)

out = minimize(residual, params, args=(x, y, data),method='ampgo')


t = out.params.valuesdict()
fparams = [t['a'],t['miux'],t['sigmax'],t['miuy'],t['sigmay'],t['ori'],t['sf'],t['phase']]
"""



## test

a=DataFrame(id=1:3,p=2:4)

b=DataFrame(id=1:2,r=4:5)
c=DataFrame(id=[])
outerjoin(a,b,on=:id)
crossjoin(a,b)

