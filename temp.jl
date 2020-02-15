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




(a= true ? 0 : 4, b=3)

a=[true,false,true]

b = any(a)
a



dogf(x,βₑ=1,βᵢ=1,μₑ=0,μᵢ=0,σₑ=1,σᵢ=1) = βₑ*exp(-(x-μₑ)^2 / (2σₑ*σₑ)) - βᵢ*exp(-(x-μᵢ)^2 / (2σᵢ*σᵢ))

plot(x->dogf(x),-10,10)











using Makie,Colors
Makie.scatter(1:2,1:2,marker=[rand(RGB,10,10) rand(RGB,10,10)],markersize=0.2)





Plots.plot(gaussianf,-3,3)
Plots.plot(gaborf,-3,3)
Plots.plot(dogf,-3,3)



x=y=-3:0.01:3
z = [gaussianf(i,j,σ₁=0.8,σ₂=0.5,θ=0.5) for j in reverse(y),i in x]
z = [gratingf(i,j,θ=0.5,f=0.5) for j in reverse(y),i in x]
z = [gaborf(i,j,σ₁=0.8,σ₂=0.5,θ=0.5,f=0.5) for j in reverse(y),i in x]
z = [dogf(i,j,σₑ₁=0.8,σₑ₂=0.5,σᵢ₁=1,σᵢ₂=0.7,θₑ=0.5,θᵢ=1.5) for j in reverse(y),i in x]


Plots.heatmap(z,aspect_ratio=:equal,frame=:none,yflip=true,leg=false,color=:plasma)



Makie.surface(x,y,z,colormap=:coolwarm,shading=false)
