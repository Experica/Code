using NeuroAnalysis,FileIO,JLD2,HypothesisTests,Images,ImageBinarization,ImageEdgeDetection,ImageMorphology,
StatsBase,StatsPlots,Contour,ImageDraw

function drawcontour!(img,ct,color=oneunit(eltype(img));fill=false)
    for line in lines(ct)
        xs, ys = coordinates(line)
        vert = map((i,j)->CartesianIndex(round(Int,i),round(Int,j)),xs,ys)
        draw!(img, Polygon(vert; fill),color)
    end
    return img
end

dataroot = "I:/"
resultroot = "Y:/"
subject = "Test";recordsession = "";recordsite = "Full"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)


## ocular dominance
function odmap(siteresultdir,lrtests::NTuple{2,String};datafile="isi.jld2")
    ds = load.(joinpath.(siteresultdir,lrtests,datafile))
    eyes = [i["exenv"]["eye"] for i in ds];eyeis = sortperm(eyes);eyesort=eyes[eyeis]
    eyesort == ["Left","Right"] || @warn "Tests aren't from Left and Right eyes, got $eyes instead."

    ler = ds[eyeis[1]]["epochresponse"];rer = ds[eyeis[2]]["epochresponse"]
    pairtest(ler,rer).stat
end

t = odmap(siteresultdir,("AG1_V1V2_Full_ISIEpochOri8_0","AG1_V1V2_Full_ISIEpochOri8_1"))
od = clampscale(dogfilter(t,lσ=50),3)
h,w = size(od)
save(joinpath(siteresultdir,"OD.png"),od)


## blood vessel and mask
bvdir = joinpath(dataroot,subject,"blood_vessel")
bvg = clampscale(readrawim_Mono12Packed(joinpath(bvdir,"blood_vessel-Frame0.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel_green.png"),bvg)
bvr = clampscale(readrawim_Mono12Packed(joinpath(bvdir,"blood_vessel-Frame1.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel_red.png"),bvr)

bvge = ahe(bvg,clip=0.9)
save(joinpath(siteresultdir,"bloodvessel_green_enhance.png"),bvge)
bvre = ahe(bvr,clip=0.9)
save(joinpath(siteresultdir,"bloodvessel_red_enhance.png"),bvre)

bvge = clampscale(dogfilter(bvg),3)
bvgm = binarize(Bool,bvge,AdaptiveThreshold(bvge))
save(joinpath(siteresultdir,"bloodvessel_green_mask.png"),bvgm)
bvre = clampscale(dogfilter(bvr),3)
bvrm = binarize(Bool,bvre,AdaptiveThreshold(bvre))

r = 4
foreach(_->dilate!(bvrm),1:r)
foreach(_->erode!(bvrm),1:r+2)
save(joinpath(siteresultdir,"bloodvessel_red_mask.png"),bvrm)

bvsdir = joinpath(dataroot,subject,"blood_vessel_after")
bvs = clampscale(readrawim_Mono12Packed(joinpath(bvsdir,"blood_vessel_after-Frame1.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel_scale.png"),bvs)


## OD border and contour
odmask = gray.(load(joinpath(siteresultdir,"mask_od.png"))) .==  1
ode = gaussianfilter(od,σ=7)
ode[.!odmask] .= 0.5
ode = ahe(ode,nblock=1,clip=0.9)

odb = detect_edges(ode, Canny(spatial_scale=10, high=ImageEdgeDetection.Percentile(70),low=ImageEdgeDetection.Percentile(30)))
save(joinpath(siteresultdir,"OD_border.png"),odb)
odbm = map(i->GrayA(i,i),odb)
save(joinpath(siteresultdir,"OD_border_mask.png"),odbm)

cl=[0.25,0.75] # contour levels for left(dark) and right(bright) domain
odct = contours(1:h,1:w,ode,cl)
odctl = Contour.levels(odct)[1]
odcth = Contour.levels(odct)[2]
odcl = drawcontour!(similar(od,Gray{Bool}),odctl)
odcr = drawcontour!(similar(od,Gray{Bool}),odcth)
save(joinpath(siteresultdir,"OD_contour_left.png"),odcl)
save(joinpath(siteresultdir,"OD_contour_right.png"),odcr)

odcml = drawcontour!(similar(od,GrayA{N0f8}),odctl)
odcmr = drawcontour!(similar(od,GrayA{N0f8}),odcth)
save(joinpath(siteresultdir,"OD_contour_mask_left.png"),odcml)
save(joinpath(siteresultdir,"OD_contour_mask_right.png"),odcmr)

odcfl = drawcontour!(similar(od,Gray{Bool}),odctl,fill=true)
odcfr = drawcontour!(similar(od,Gray{Bool}),odcth,fill=true)
save(joinpath(siteresultdir,"OD_contourfill_left.png"),odcfl)
save(joinpath(siteresultdir,"OD_contourfill_right.png"),odcfr)

odcfml = drawcontour!(similar(od,RGBA{N0f8}),odctl,RGBA{N0f8}(1,0,0,0.02),fill=true)
odcfmr = drawcontour!(similar(od,RGBA{N0f8}),odcth,RGBA{N0f8}(0,1,0,0.02),fill=true)
save(joinpath(siteresultdir,"OD_contourfill_mask_left.png"),odcfml)
save(joinpath(siteresultdir,"OD_contourfill_mask_right.png"),odcfmr)


## Cone-Opponent Functional Domain
function cofdmap(siteresultdir,test;datafile="isi.jld2")
    d = load(joinpath(siteresultdir,test,datafile))
    color = d["exenv"]["color"];minmaxcolor = d["exenv"]["minmaxcolor"]
    if haskey(d,"F1polarity")
        cofd = d["F1polarity"]
        color*="_p"
    elseif haskey(d,"t")
        cofd = d["t"]
        color*="_t"
    end
    (;cofd,color,minmaxcolor)
end

cofd,color,minmaxcolor = cofdmap(siteresultdir,"Test_Full_ISIEpoch2Color_0")
cofd = clampscale(dogfilter(cofd),3)
save(joinpath(siteresultdir,"COFD_$(color).png"),cofd)


## COFD contour
cofdmask = gray.(load(joinpath(siteresultdir,"mask_cofd.png"))) .== 1
cofdmask .&= bvrm
cofde = gaussianfilter(cofd)
cofde = ahe(cofde,nblock=30,clip=1)
cofde[.!cofdmask] .= 0.5

cl = [0.25,0.75] # contour levels for mincolor(dark) and maxcolor(bright) domain
cofdct = contours(1:h,1:w,cofde,cl)
cofdctl = Contour.levels(cofdct)[1]
cofdcth = Contour.levels(cofdct)[2]
cofdcl = drawcontour!(similar(cofd,Gray{Bool}),cofdctl)
cofdch = drawcontour!(similar(cofd,Gray{Bool}),cofdcth)
save(joinpath(siteresultdir,"COFD_$(color)_contour_low.png"),cofdcl)
save(joinpath(siteresultdir,"COFD_$(color)_contour_high.png"),cofdch)

cofdcml = drawcontour!(similar(cofd,RGBA{Float64}),cofdctl,minmaxcolor.mincolor)
cofdcmh = drawcontour!(similar(cofd,RGBA{Float64}),cofdcth,minmaxcolor.maxcolor)
save(joinpath(siteresultdir,"COFD_$(color)_contour_mask_low.png"),cofdcml)
save(joinpath(siteresultdir,"COFD_$(color)_contour_mask_high.png"),cofdcmh)

cofdcfl = drawcontour!(similar(cofd,Gray{Bool}),cofdctl,fill=true)
cofdcfh = drawcontour!(similar(cofd,Gray{Bool}),cofdcth,fill=true)
save(joinpath(siteresultdir,"COFD_$(color)_contourfill_low.png"),cofdcfl)
save(joinpath(siteresultdir,"COFD_$(color)_contourfill_high.png"),cofdcfh)

cofdcfml = drawcontour!(similar(cofd,RGBA{Float64}),cofdctl,coloralpha(minmaxcolor.mincolor,0.1),fill=true)
cofdcfmh = drawcontour!(similar(cofd,RGBA{Float64}),cofdcth,coloralpha(minmaxcolor.maxcolor,0.1),fill=true)
save(joinpath(siteresultdir,"COFD_$(color)_contourfill_mask_low.png"),cofdcfml)
save(joinpath(siteresultdir,"COFD_$(color)_contourfill_mask_high.png"),cofdcfmh)


