using NeuroAnalysis,FileIO,JLD2,HypothesisTests,Images,ImageBinarization,ImageEdgeDetection,
StatsBase,StatsPlots,PyCall,Contour

dataroot = "I:/"
dataexportroot = "Y:/"
resultroot = "Y:/"
subject = "AG1";recordsession = "V1V2";recordsite = "Full"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
skid = pyimport("skimage.draw")

function drawcontour!(img,ct;color=1)
    for line in lines(ct)
        xs, ys = coordinates(line)
        foreach((i,j)->img[round(Int,i),round(Int,j)] = color,xs,ys)
    end
    return img
end
function fillcontour!(img,ct;color=1)
    for line in lines(ct)
        xs, ys = coordinates(line)
        vert = map((i,j)->CartesianIndex(round(Int,i),round(Int,j)),xs,ys)
        draw!(img,vert,BoundaryFill(round(Int,mean(extrema(xs))),round(Int,mean(extrema(ys))),fill_value=color,boundary_value=color))
    end
    return img
end
function fillcontour!(img,ct,skid;color=1)
    px=Int[];py=Int[];d=size(img)
    for line in lines(ct)
        xs, ys = coordinates(line)
        pp = skid.polygon(xs.-1,ys.-1,d)
        append!(px,pp[1].+1);append!(py,pp[2].+1)
    end
    foreach((i,j)->img[i,j]=color,px,py)
    img
end

## blood vessel
bvdatadir = joinpath(dataroot,subject,"blood_vessel")
h = 2080
w = 2080
bv1 = clampscale(readrawim_Mono12Packed(joinpath(bvdatadir,"blood_vessel2-Epoch0-Frame4.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel1.png"),bv1)
bv2 = clampscale(readrawim_Mono12Packed(joinpath(bvdatadir,"blood_vessel2-Epoch0-Frame5.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel2.png"),bv2)

bv1e = adjust_histogram(bv1, AdaptiveEqualization(nbins = 256, rblocks = 12, cblocks = 12, clip = 0.85))
clamp01!(bv1e)
save(joinpath(siteresultdir,"bloodvessel1_Enhanced.png"),bv1e)

bv2e = adjust_histogram(bv2, AdaptiveEqualization(nbins = 256, rblocks = 12, cblocks = 12, clip = 0.85))
clamp01!(bv2e)
save(joinpath(siteresultdir,"bloodvessel2_Enhanced.png"),bv2e)

bv1m = binarize(-bv1e,AdaptiveThreshold())
save(joinpath(siteresultdir,"bloodvessel1_mask.png"),bv1m)

bv2m = binarize(-bv2e,AdaptiveThreshold())
save(joinpath(siteresultdir,"bloodvessel2_mask.png"),bv2m)


bvsdatadir = joinpath(dataroot,subject,"blood_vessel_after")
bv1s = clampscale(readrawim_Mono12Packed(joinpath(bvsdatadir,"blood_vessel_after-Epoch0-Frame2.Raw"),w,h))
save(joinpath(siteresultdir,"bloodvessel1_scale.png"),bv1s)



## ocular dominance
function odmap(siteresultdir,lrtests::NTuple{2,String};datafile="isi.jld2")
    ds = load.(joinpath.(siteresultdir,lrtests,datafile))
    eyes = [i["exenv"]["eye"] for i in ds];eyeis = sortperm(eyes);eyesort=eyes[eyeis]
    eyesort == ["Left","Right"] || @warn "Tests aren't from Left and Right eyes, got $eyes instead."

    ler = ds[eyeis[1]]["epochresponse"];rer = ds[eyeis[2]]["epochresponse"]
    pairtest(ler,rer)
end


t = odmap(siteresultdir,("AG1_V1V2_Full_ISIEpochOri8_0","AG1_V1V2_Full_ISIEpochOri8_1")).stat
od = clampscale(dogfilter(t,lσ=50),3)
save(joinpath(siteresultdir,"OD_LR.png"),od)



## OD inpainting
od = Gray.(load(joinpath(siteresultdir,"OD_inpaint1.png")))
od = imresize(gray.(od),h,w)


## OD border

Gray.(oda)
Gray.(clampscale(oda,0.5))



odt = gaussianfilter(od)
oda = ahe(odt,nblock=2,clip=0.9)
odb = detect_edges(od, Canny(spatial_scale=10, high=ImageEdgeDetection.Percentile(70),low=ImageEdgeDetection.Percentile(30)))
save(joinpath(siteresultdir,"OD_border.png"),odb)
odbm = map(i->GrayA(i,i),odb)
save(joinpath(siteresultdir,"OD_border_mask.png"),odbm)


cl=[0.25,0.75]
odct = contours(1:h,1:w,od,cl)
odc0 = drawcontour!(similar(od),Contour.levels(odct)[1])
odc1 = drawcontour!(similar(od),Contour.levels(odct)[2])
save(joinpath(siteresultdir,"OD_contour_$(eye0).png"),odc0)
save(joinpath(siteresultdir,"OD_contour_$(eye1).png"),odc1)

odcm0 = drawcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[1],color=RGBA{N0f8}(1,1,1,1))
odcm1 = drawcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[2],color=RGBA{N0f8}(1,1,1,1))
save(joinpath(siteresultdir,"OD_contour_mask_$(eye0).png"),odcm0)
save(joinpath(siteresultdir,"OD_contour_mask_$(eye1).png"),odcm1)

odcf0 = fillcontour!(similar(od),Contour.levels(odct)[1],skid)
odcf1 = fillcontour!(similar(od),Contour.levels(odct)[2],skid)
save(joinpath(siteresultdir,"OD_contourfill_$(eye0).png"),odcf0)
save(joinpath(siteresultdir,"OD_contourfill_$(eye1).png"),odcf1)

odcfm0 = fillcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[1],skid,color=RGBA{N0f8}(1,0,0,0.02))
odcfm1 = fillcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[2],skid,color=RGBA{N0f8}(0,1,0,0.02))
save(joinpath(siteresultdir,"OD_contourfill_mask_$(eye0).png"),odcfm0)
save(joinpath(siteresultdir,"OD_contourfill_mask_$(eye1).png"),odcfm1)



## Cone Domain
isi = load(joinpath(siteresultdir,"AG1_V1V2_Full_ISICycle2Color_28","isi.jld2"))
color = isi["exenv"]["color"];minmaxcolor = isi["exenv"]["minmaxcolor"]
cone = isi["F1polarity"]
cone = adjust_histogram(cone, AdaptiveEqualization(nbins = 256, rblocks = 30, cblocks = 30, clip = 0.1))
clamp01!(cone)
save(joinpath(siteresultdir,"$(color)_COFD.png"),cone)


## Cone inpainting
cone = Gray.(load(joinpath(siteresultdir,"$(color)_COFD_inpaint1.png")))
cone = imresize(gray.(cone),h,w)


## Cone Domain contour
cone = imfilter(cone,Kernel.gaussian(5))
cone = adjust_histogram(cone, AdaptiveEqualization(nbins = 256, rblocks = 2, cblocks = 2, clip = 0.9))

cl = [0.25,0.75]
cdct = contours(1:h,1:w,cone,cl)
cdc0 = drawcontour!(similar(cone),Contour.levels(cdct)[1])
cdc1 = drawcontour!(similar(cone),Contour.levels(cdct)[2])
save(joinpath(siteresultdir,"$(color)_COFD_contour_low.png"),cdc0)
save(joinpath(siteresultdir,"$(color)_COFD_contour_high.png"),cdc1)

cdcm0 = drawcontour!(similar(cone,RGBA{Float64}),Contour.levels(cdct)[1],color=minmaxcolor.mincolor)
cdcm1 = drawcontour!(similar(cone,RGBA{Float64}),Contour.levels(cdct)[2],color=minmaxcolor.maxcolor)
save(joinpath(siteresultdir,"$(color)_COFD_contour_mask_low.png"),cdcm0)
save(joinpath(siteresultdir,"$(color)_COFD_contour_mask_high.png"),cdcm1)

cdcf0 = fillcontour!(similar(cone),Contour.levels(cdct)[1],skid)
cdcf1 = fillcontour!(similar(cone),Contour.levels(cdct)[2],skid)
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_low.png"),cdcf0)
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_high.png"),cdcf1)

cdcfm0 = fillcontour!(similar(cone,RGBA{Float64}),Contour.levels(cdct)[1],skid,color=coloralpha(minmaxcolor.mincolor,0.1))
cdcfm1 = fillcontour!(similar(cone,RGBA{Float64}),Contour.levels(cdct)[2],skid,color=coloralpha(minmaxcolor.maxcolor,0.1))
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_mask_low.png"),cdcfm0)
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_mask_high.png"),cdcfm1)


