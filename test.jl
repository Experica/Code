using Makie,Plots

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



function plotunitpositionold(unitposition;unitgood=[],chposition=[],unitid=[],layer=nothing,color=nothing,alpha=0.4,title="")
    nunit = size(unitposition,1);ngoodunit = isempty(unitgood) ? nunit : count(unitgood);us = "$ngoodunit/$nunit"
    xlim = isempty(chposition) ? (minimum(unitposition[:,1])-5,maximum(unitposition[:,1])+5) : (minimum(chposition[:,1])-5,maximum(chposition[:,1])+5)
    p = plot(legend=:topright,xlabel="Position_X (um)",ylabel="Position_Y (um)",xlims=xlim)
    if !isempty(chposition)
        scatter!(p,chposition[:,1],chposition[:,2],markershape=:rect,markerstrokewidth=0,markersize=2,color=:grey60,label="Electrode")
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
        scatter!(p,unitposition[:,1],unitposition[:,2],label=us,color=color,markerstrokewidth=0,markersize=6,series_annotations=text.(unitid,3,:gray10,:center),title=title)
    else
        scatter!(p,unitposition[:,1],unitposition[:,2],label=us,color=color,markerstrokewidth=0,markersize=5,title=title)
    end
    if !isnothing(layer)
        lx = xlim[1]+2
        hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lx,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    return p
end


function plotunitposition(unitposition;oris=nothing)
    p = Scene(xlabel="Position_X (um)",ylabel="Position_Y (um)")
    scatter!(p,unitposition[:,1],unitposition[:,2])
    return p
end

oris=rand(0:360,20)
unitposition = rand(20,2)

plotunitposition(unitposition)
