using NeuroAnalysis,Makie,Colors,Plots,GeometryTypes

x = 1:10
y=1:10
u=fill(0.5,10)
v = fill(0.5,10)
arrows(x,y,u,v)


m = rand(RGB,10,10)
ms = [rand(RGBf0,10,10) for _ in 1:10]
mms =fill(m,10)


m=AbstractPlotting.logo()
fill(AbstractPlotting.logo(),10)



Makie.scatter(x,y,marker=ms,markersize=1)
Plots.scatter(x,y,marker=m,markersize=10)



function plotunitpositionold(unitposition;unitgood=[],chposition=[],unitid=[],layer=nothing,color=nothing,alpha=0.4,title="",marker=:circle)
    nunit = size(unitposition,1);ngoodunit = isempty(unitgood) ? nunit : count(unitgood);us = "$ngoodunit/$nunit"
    xlim = isempty(chposition) ? (minimum(unitposition[:,1])-5,maximum(unitposition[:,1])+5) : (minimum(chposition[:,1])-5,maximum(chposition[:,1])+5)
    p = Plots.plot(legend=:topright,xlabel="Position_X (um)",ylabel="Position_Y (um)",xlims=xlim)
    if !isempty(chposition)
        Plots.scatter!(p,chposition[:,1],chposition[:,2],markershape=:rect,markerstrokewidth=0,markersize=2,color=:grey60,label="Electrode")
    end
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
    if !isempty(unitid)
        Plots.scatter!(p,unitposition[:,1],unitposition[:,2],label=us,color=color,markerstrokewidth=0,markersize=6,series_annotations=text.(unitid,3,:gray10,:center),title=title)
    else
        Plots.scatter!(p,unitposition[:,1],unitposition[:,2],label=us,marker=marker,color=color,markerstrokewidth=0,markersize=5,title=title)
    end
    if !isnothing(layer)
        lx = xlim[1]+2
        Plots.hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lx,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    return p
end



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

plotunitposition(unitposition,unitgood=unitgood)

plotunitpositionold(unitposition,unitgood=unitgood,marker=OriMarker)


m = Point2f0[(0,0), (0, 1), (0.5, 0.5), (1, 1), (1, 0)]
meshscatter(rand(1:7,20),rand(0:3,20), marker=GLPlainMesh(m),rotations=rand(0:2π,20))



AbstractPlotting.set_theme!(
    plot = (show_axis = false, scale_plot = false),
    color = :turquoise1
)
poly(HyperRectangle(Vec2f0(0), Vec2f0(1,100)))

using NeuroAnalysis,DataFrames,YAML

t=mat2julia!(matread("hartleycond.mat")["c"])

YAML.write_file("Hartley_k[30,30].yaml",t)

c= DataFrame(t)
condin(c)



t = load("C:\\Users\\fff00\\Command\\Recordings\\0040.png")

tt=t[200:600,200:600]

gt=Gray.(tt)

ggt=gray.(gt)

plot(ggt[:,200])




g = grating(ori=0,sf=0,phase=0.5)

a=Gray.(g)


c=DataFrame(ex["Cond"])

imageset = map(i->GrayA.(grating(deg2rad(i.Ori),i.SpatialFreq,i.SpatialPhase)),eachrow(c))



unmaskindex

r=sta(x,y)

Plots.imageHack()

using Plots

plot(a,seriestype=:heatmap,color=:coolwarm,ratio=:equal,framestyle=:none,yflip=true)

plot(a,seriestype=:heatmap,color=:coolwarm,ratio=:equal,framestyle=:none,title="sdfsdf")

coloralpha(RGB(),0)

cgrad(:coolwarm)

C(g::ColorGradient) = RGB[g[z] for z=range(0,1,length=30)]

cgrad(:coolwarm) |> C




using VegaLite,DataFrames,Images

function plotunitposition(unitposition;ori=nothing,os=nothing,dir=nothing,ds=nothing,sf=nothing,vi=trues(size(unitposition,1)),width=500,height=400,title="")
    df = DataFrame(x=unitposition[:,1],y=unitposition[:,2],m=Any[:circle for _ in 1:size(unitposition,1)],a=0.0,sw=0,sf=0.0,s=10.0)
    if !isnothing(sf)
        df[vi,:sf]=sf[vi]
    end
    if !isnothing(ori)
        df[vi,:m] = :stroke
        df[vi,:a] = -ori[vi]
        df[vi,:sw] = 1
        if !isnothing(os)
            df[vi,:s]=os[vi]
        else
            df[vi,:s]=160
        end
    end
    if !isnothing(dir)
        t = copy(df)
        arrowpath = "M 0 -0.1 H 1 V -0.3 L 1.6 0 L 1 0.3 V 0.1 H 0 Z"
        t[vi,:m] = arrowpath
        t[vi,:a] = -dir[vi]
        t[vi,:sw] = 0
        if !isnothing(ds)
            t[vi,:s]=ds[vi]
        else
            t[vi,:s]=200
        end
        df = [df;t]
    end
    @vgplot(height=height,width=width,padding=5,data=[:df=>df],
    marks=[
    {
    type="symbol",
    from={data="df"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        shape={field="m"},
        angle={field="a"},
        size={field="s",scale="s"},
        strokeWidth={field="sw"},
        stroke={field="sf",scale="c"},
        fill={field="sf",scale="c"}
        }}
    }
    ],
    scales=[
    {
        name="x",
        nice=true,
        zero=false,
        range="width",
        domain={data="df",field="x"},
        type="linear",
        round=true
    },
    {
        name="y",
        nice=true,
        zero=false,
        range="height",
        domain={data="df",field="y"},
        type="linear",
        round=true
    },
    {
        name="c",
        nice=true,
        zero=false,
        round=false,
        type="linear",
        range={scheme = "yellowgreenblue"},
        domain={data="df",field="sf"}
    },
    {
        name="s",
        nice=true,
        zero=true,
        domain={data="df",field="s"},
        type="linear",
        round=true,
        range=[0, 200]
    }
    ],
    axes=[
    {
        domain=true,
        tickCount=5,
        grid=false,
        title="Position_X (μm)",
        scale="x",
        orient="bottom"
    },
    {
        domain=true,
        tickCount=5,
        grid=false,
        titlePadding=5,
        title="Position_Y (μm)",
        scale="y",
        orient="left"
    }
    ],
    title={
    text = title
    })
end


function foo(pos;ori=0,dir=90)
    data1 = DataFrame(x=pos[:,1],y=pos[:,2],shape=1,angle=ori)
    data2 = DataFrame(x=pos[:,1],y=pos[:,2],shape=2,angle=dir)
    data = [data1;data2]
    # data |> @vlplot(
    # mark={:point,filled=true,strokeWidth=4},
    # encoding={
    # x=:x,y=:y,shape={"shape:n",scale={domain=[1,2],range=[:stroke,:arrow]}},size={"shape:n",scale={domain=[1,2],range=[200,400]}},
    # angle="angle:q"
    # },
    # width=400,height=400 )

    arrowpath = "M 0 -0.04 H 1 V -0.09 L 1.2 0 L 1 0.09 V 0.04 H 0 Z"
    # @vlplot(data=data1,x=:x,y=:y,width=400,height=400)+
    # @vlplot(mark={:point,shape=cs,size=400,filled=true,angle=45})+
    # @vlplot(mark={:point,shape=:stroke,strokeWidth=1,size=200,angle=0})
    # #+@vlplot(mark={:text,text="→",size=100})

    @vgplot(
    height=400,
    width=400,
    padding=5,
    data=[:source=>data1,:source2=>data2],
    marks=[{
        encode={
            update={
                angle={field="angle"},
                shape={value="stroke"},
                x={field="x", scale="x"},
                y={field="y", scale="y"},
                size={value=200},
                strokeWidth={value=1},
                stroke={value="red"}
            }
        },
        from={data="source"},
        type="symbol"
    },
    {
        encode={
            update={
                angle={field="angle"},
                shape={value=arrowpath},
                x={field="x", scale="x"},
                y={field="y", scale="y"},
                size={value=400},
                strokeWidth={value=0},
                fill={value="blue"},
                stroke={value="blue"}
            }
        },
        from={data="source2"},
        type="symbol"
        }],
    axes=[
        {
            domain=false,
            tickCount=5,
            grid=true,
            title="Horsepower",
            scale="x",
            orient="bottom"
        },
        {
            domain=false,
            grid=true,
            titlePadding=5,
            title="Miles_per_Gallon",
            scale="y",
            orient="left"
        }
    ],
    scales=[
        {
            name="x",
            nice=true,
            zero=true,
            range="width",
            domain={data="source",field="x"},
            type="linear",
            round=true
        },
        {
            name="y",
            nice=true,
            zero=true,
            range="height",
            domain={data="source",field="y"},
            type="linear",
            round=true
        }
    ]
)
end

up = [rand(20:70,20) rand(20:3500,20)]
p=plotunitposition(up)
p=plotunitposition(up,ori=rand(0:180,20))
p=plotunitposition(up,ori=rand(0:180,20),dir=rand(0:360,20),sf=rand(20),os=rand(0:10,20),ds=rand(0:10,20),title="test sample")

p=foo(up,ori=rand(0:180,10),dir=rand(0:360,10))

save("testup.svg",p)
save("testup.png",p)


t=DataFrame(a=1:10,m=:sdc)

tt=copy(t)



t[:m]=:ttt

t


tt

a=[1,2,3]

-a
