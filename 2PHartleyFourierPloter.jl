using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase
using Plots,StatsPlots,VegaLite,ImageView,Interact,CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random
using SchwarzChristoffel
using Luxor

# Expt info
disk = "K:"
subject = "AE6"  # Animal
recordPlane = ["001", "000"]  # I hard-code the plane Id here. The Id here is matched with the naming rule of multiple plane scanning


# pValue = 0.05

# oriaucThres = 0.7
# diraucThres = 0.7
# hueaucThres = 0.9

## Make folders and path
mainpath = joinpath(disk,subject, "2P_analysis")
dataExportFolder1 = joinpath(mainpath, "Summary", "DataExport")
roibkgFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_roibkg.jld2"),dir=dataExportFolder1,adddir=true)[1]
roi=load(roibkgFile,"roi")[1]
bkg=load(roibkgFile,"bkg")[1]

meta = readmeta(joinpath(dataExportFolder1,"metadata.mat"))
for t in 1:nrow(meta)
    meta.files[t]=meta.files[t][3:end]
end

## Query Tests
tests = @from i in meta begin
        @where i.Subject_ID == subject
        @where i.sourceformat == "Scanbox"
        @select {i.ID,i.RecordSite,i.filename}
        @collect DataFrame
        end
sort!(tests, [:RecordSite, :filename])

hartleytests = @from i in tests begin
        @where i.RecordSite != "u001"
        @where i.RecordSite != "u005"
        @where i.ID == "Hartley"
        @select {i.ID,i.RecordSite,i.filename}
        @collect DataFrame
        end

recordSession = unique!(hartleytests.RecordSite) # Unit
recordSession = map(i->i[2:4],recordSession)

hartleyExptId = hartleytests.filename    # Unit Ids have L, M, S, and Achromatic hartleys.
exptHartNum = length(hartleyExptId)

planeNum = Int64(length(hartleyExptId) / length(recordSession))
sessNum = length(recordSession)*planeNum
if length(recordPlane) != planeNum
        @warn "Warning: More plane acquired than you think!!!!"
        exit()
end

# paraName = string("oa",string(oriaucThres),"da",string(diraucThres),"ha",string(hueaucThres),"oc",string(oricvThres),
#            "dc",string(dsiThres),"hc",string(huecvThres),"osi",string(osiThres),"dsi",string(dsiThres),"hsi",string(huesiThres),"cpi",string(cpiThres))

# for i = 1:exptHartNum

display("Processing Hartley expt No: $i")
j=map(j->isodd(i) ? 1 : 2, i)
unitId = hartleyExptId[i][5:7]
siteId = join([unitId, "_",recordPlane[j]])
dataFolder = joinpath(mainpath, join(["U", unitId]), join([hartleyExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_tuning_result.jld2"),dir=dataFolder,adddir=true)[1]
result = load(dataFile)["result"]



for u in sort(collect(keys(dataset["ulsta"])))
        # u=11
        ulsta = dataset["ulsta"][u]
        delays = dataset["delays"]
        imagesize = size(ulsta)[1]
        stisize = dataset["stisize"]
        maskradius = dataset["maskradius"]
        truestimsz = stisize*maskradius*2
        ucex = map(i->i.ex,dataset["ulcex"][u])
        ucd = map(i->i.pd,dataset["ulcex"][u])
        clim = maximum(abs.(ucex))
        xylim = [0,round(truestimsz,digits=1)]
        xy = range(xylim...,length=imagesize)

        p = Plots.plot(layout=(1,5),legend=false,size=(1650,600))
        Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
        foreach(c->Plots.heatmap!(p,subplot=c+1,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:grid,
        color=:coolwarm,clims=(-clim,clim),xlims=xylim,ylims=xylim,xticks=[],yticks=[],
        yflip=true,xlabel=string(dataset["color"][c]),title="Cell_$(u)_STA_$(delays[ucd[c]])"),1:4)
        foreach(c->Plots.plot!(p,subplot=c+1,[0.3,0.6,0.9,1.2,1.5], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.3,0.6,0.9,1.2,1.5,1.8],["-0.9","-0.6","-0.3","0","0.3","0.6","0.9"]), label=""),1:4)
        foreach(c->Plots.plot!(p,subplot=c+1,[0.3,0.6,0.9,1.2,1.5], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
        foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.3,0.6,0.9,1.2,1.5,1.8],["-0.9","-0.6","-0.3","0","0.3","0.6","0.9"])),1)

        # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]), label=""),1:4)
        # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
        # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"])),1)

        # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"]), label=""),1:4)
        # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
        # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"])),1)

        p

        savefig(joinpath(resultFolderPlot,join([subject,"_U",recordSession,"_Plane",recordPlane, "_Cell",u,".png"])))
end








# end


# imshow(bkg)
color_hex_keys = ["#FF8080", "#FF80BF", "#FF80FF", "#BF80FF", "#8080FF", "#80BFFF", "#80FFFF", "#80FFBF", "#80FF80", "#BFFF80", "#FFFF80", "#FFBF80", "#808080"]   # DKL hues # DKL hues
poly=Any[]
filrange = []
for i = 1: length(roi)
        # i=1
        rr=roi[i]
        push!(filrange, [min(rr[:,2]...),max(rr[:,2]...)])
        ra=rr[1:2:end,:]
        rb=rr[2:2:end,:]
        rcb=reverse(rb,dims=1)
        # r=vcat(ra,rcb,[rr[1,1] rr[1,2]])
        r=vcat(ra,rcb,ra[1:2,:])
        x=r[:,1]  # actually this is column number, I revese it in matlab, PL
        y=r[:,2]  # actually this is row number
        p = Polygon(x,y)
        push!(poly,p)
end
fillrange=filrange[1]
# plot!(poly[i] for i in 1:size(poly,1))
# end@manipulate margin=0
for i in 1:size(poly,1)
        global  img3
        img3=plot!(poly[i], linecolor = parse(Colorant, color_hex_keys[2]),fillcolor=parse(Colorant, color_hex_keys[2]), fillalpha=1, fillrange=[filrange[i][1] filrange[i][2]],
        size=(size(bkg,2),size(bkg,1)), fmt=:svg, showaxis=false, dpi=300)
end

x = [-0,-0.2,0.3,0.2,-0.0]; y = [-0,0.1,0.2,0.5,0.1];
p = Polygon(x,y)
plot(poly[1], linecolor = :red, linewidth=0.1, fillcolor=:blue, fillalpha=1,fillrange=[filrange[1][1] filrange[1][2]], size=(size(bkg,2),size(bkg,1)),fmt=:svg, dpi=300)

savefig("test000.svg")
aaa=savefig(img3,"test11.svg")
savefig(img3,"test2.png")
imag=load("test11.svg")
plot(imag, size=(size(bkg,2),size(bkg,1)),fmt=:svg, showaxis=false,dpi=300)
savefig(aaa,"test11.png")

@png begin
        img3
end

@png begin
        tiles = Tiler(size(bkg,2),size(bkg,1), 1, 1, margin=0)
        # tile1, tile2 = collect(tiles)
        Luxor.translate(tile1[1])
        # randompoints = [Point(rand(-100:100), rand(-100:100)) for i in 1:10]

        # for i in 1:length(roi)
        temproi=roi[1]
        pl=Point[]
        for j in 1:size(temproi,1)
                pts = [Point(temproi[j,1], temproi[j,2])]
                append!(pl, pts)
        end
        Luxor.poly(pl, :fill)
        # end
        # gsave()

        Luxor.poly(randompoints, :stroke)

end



randompoints = [Point(rand(-100:100), rand(-100:100)) for i in 1:10]


ori_lut = ListedColormap(sio.loadmat(os.path.join(colormappath, 'ori_lut_alpha0.mat'), squeeze_me=True, struct_as_record=False)['lut'])
# cmap_patch = plt.cm.get_cmap('hsv')
pylt.cm.get_cmap("hot")
cmap_patch_cpi = pylt.cm.get_cmap('jet')
cmap_patch_osi = plt.cm.get_cmap('jet')
cmap_patch_hue = plt.cm.get_cmap('hsv_r')
cmap_patch_hue1 = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,'dkllut.mat'), squeeze_me=True, struct_as_record=False)['lut'])
# if colorSpace  == 'DKL':

color_hex_keys = ["#FF8080", "#FF80BF", "#FF80FF", "#BF80FF", "#8080FF", "#80BFFF", "#80FFFF", "#80FFBF", "#80FF80", "#BFFF80", "#FFFF80", "#FFBF80", "#808080"]   # DKL hues # DKL hues
# elif colorSpace == 'CIE':
color_hex_keys = plcolor.ListedColormap(['#af1600', '#8a4600', '#5a5d01', '#2a6600', '#006a00', '#006931', '#006464',
             '#0058b6', '#002DFF', '#6a2ade', '#97209b', '#aa1c50', '#808080'])   # HSL hues
# gray = '#B4B4B4'    # '#808080'


# boxplot is defined in StatsPlots
using StatsPlots
gr(leg=false, bg=:lightgrey)

# Create a filled contour and boxplot side by side.
plot(contourf(randn(10,20)), boxplot(rand(1:4,1000),randn(1000)))

# Add a histogram inset on the heatmap.
# We set the (optional) position relative to bottom-right of the 1st subplot.
# The call is `bbox(x, y, width, height, origin...)`, where numbers are treated as "percent of parent"
histogram!(randn(1000), inset = (1, bbox(0.05,0.05,0.5,0.25,:bottom,:right)), ticks=nothing, subplot=3, bg_inside=nothing)

# Add sticks floating in the window (inset relative to the window, as opposed to being relative to a subplot)
sticks!(randn(100), inset = bbox(0,-0.2,200px,100px,:center), ticks=nothing, subplot=4)
contourf(randn(10,20))

plot(Shape([25, 425, 425, 25], [-1, -5.24, 5.24, 1]))

plot(0:5,0:5)

x = LinRange(-2, 2, 40)
y = 2 .* x .+ 4
plot(x, y)


C(g::ColorGradient) = RGB[g[z] for z=range(0,stop=1,length=30)]
g = :inferno
cgrad(g) |> C

cgrad(g, scale=:log) |> C
cgrad(g, scale=:exp) |> C

cgrad(g, [0.01, 0.99]) |> C

using RDatasets
iris = dataset("datasets", "iris");

# load the StatsPlots recipes (for DataFrames) available via:
# Pkg.add("StatsPlots")
using StatsPlots

# Scatter plot with some custom settings
@df iris scatter(:SepalLength, :SepalWidth, group=:Species,
        title = "My awesome plot",
        xlabel = "Length", ylabel = "Width",
        m=(0.2, [:cross :hex :star5], 10),
        bg=RGB(.2,.2,.2))

# save a png
png("iris")


plot(rand(10))

scatter(rand(10,4), markershape = [:circle, :rect])
scatter(rand(10,4), markershape = [:circle :rect])

# 10 data points in 4 series
xs = 0 : 2π/10 : 2π
data = [sin.(xs) cos.(xs) 2sin.(xs) 2cos.(xs)]

# We put labels in a row vector: applies to each series
labels = ["Apples" "Oranges" "Hats" "Shoes"]

# Marker shapes in a column vector: applies to data points
markershapes = [:circle, :star5]

# Marker colors in a matrix: applies to series and data points
markercolors = [:green :orange :black :purple
                :red   :yellow :brown :white]

plot(xs, data, label = labels, shape = markershapes, color = markercolors,
     markersize = 10)



pyplot(leg = false, grid = false, xticks = nothing, yticks = nothing, size=(500,500))

function make_batman()
 p = P2[(0,0), (0.5, 0.2), (1, 0), (1,2),  (0.3,1.2), (0.2,2), (0,1.7)]
 m = P2[(p[i]+p[i+1])/2 for i=1:length(p)-1]
 m += P2[(0.2, 1), (0.4, 1), (2, 0), (0.5, -0.6), (0,0), (0,-0.15)]

 pts = P2[]
 for (i,mi) in enumerate(m)
     append!(pts, curve_points(BezierCurve(P2[p[i], m[i], p[i+1]])))
 end
 x, y = Plots.unzip(pts)
 Shape(vcat(x, -reverse(x)), vcat(y, reverse(y)))
end

x = range(0; stop=2*pi, length=1000); y = sin.(3 * x + 4 * cos.(2 * x));
plot(x, y, color="red", linewidth=2.0, linestyle="--")
title("A sinusoidally modulated sinusoid")

randompoints = [Point(rand(-100:100), rand(-100:100)) for i in 1:10]
randompoints[1]
@png begin
        tiles = Tiler(600, 250, 1, 2, margin=20)
        tile1, tile2 = collect(tiles)

        randompoints = [Point(rand(-100:100), rand(-100:100)) for i in 1:10]

        gsave()

        Luxor.translate(tile1[1])
        Luxor.poly(randompoints, :stroke)

        grestore()
end
gsave()
Luxor.translate(tile2[1])
Luxor.poly(randompoints, :fill)
grestore()

@png begin
    Luxor.text("Hello world")
    circle(Point(0, 0), 200, :stroke)
end
