using NeuroAnalysis,FileIO,HypothesisTests,Images,ImageEdgeDetection,ImageBinarization,Contour,PyCall

resultroot = "Z:/"
subject = "AG2";recordsession = "";recordsite = ""
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

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
h = 2080
w = 2080
bv1 = clampscale(readrawim_Mono12Packed("I:\\AG2\\blood_vessel_before\\blood_vessel_before-Epoch0-Frame1.Raw",w,h))
save(joinpath(siteresultdir,"bloodvessel1.png"),bv1)
bv2 = clampscale(readrawim_Mono12Packed("I:\\AG1\\blood_vessel\\blood_vessel-Epoch0-Frame4.Raw",w,h))
save(joinpath(siteresultdir,"bloodvessel2.png"),bv2)

bv1ce = adjust_histogram(bv1, AdaptiveEqualization(nbins = 256, rblocks = 12, cblocks = 12, clip = 0.85))
clamp01!(bv1ce)
save(joinpath(siteresultdir,"bloodvessel1_Enhanced.png"),bv1ce)

bv2ce = adjust_histogram(bv2, AdaptiveEqualization(nbins = 256, rblocks = 12, cblocks = 12, clip = 0.85))
clamp01!(bv2ce)
save(joinpath(siteresultdir,"bloodvessel2_Enhanced.png"),bv2ce)

## blood vessel mask
bv1m = imfilter(-bv1ce,Kernel.DoG((20,20)))
bv1m = binarize(bv1m,Otsu())
save(joinpath(siteresultdir,"bloodvessel1_mask.png"),bv1m)

bv2m = imfilter(-bv2ce,Kernel.DoG((20,20)))
bv2m = binarize(bv2m,Otsu())
save(joinpath(siteresultdir,"bloodvessel2_mask.png"),bv2m)



## ocular dominence
eye1,epochresponse1 = load(joinpath(siteresultdir,"AG2_ISIEpochOri8_3","isi.jld2"),"eye","epochresponse")
eye2,epochresponse2 = load(joinpath(siteresultdir,"AG2_ISIEpochOri8_4","isi.jld2"),"eye","epochresponse")
h,w = size(epochresponse1)

ht = [@views UnequalVarianceTTest(epochresponse1[i,j,:],epochresponse2[i,j,:]) for i = 1:h,j=1:w] # Welch's t-test
t = map(i->i.t,ht)
replace!(t,NaN=>0)
p1 = map(i->isnan(i.t) ? 0.5 : pvalue(i,tail=:left),ht)
p2 = map(i->isnan(i.t) ? 0.5 : pvalue(i,tail=:right),ht)

od = adjust_histogram(clampscale(t), AdaptiveEqualization(nbins = 256, rblocks = 12, cblocks = 12, clip = 0.1))
clamp01!(od)
save(joinpath(siteresultdir,"OD.png"),od)

## OD inpainting
cv = pyimport("cv2")
od = cv.inpaint(round.(UInt8,od.*255),UInt8.(bv1m),5,cv.INPAINT_TELEA)
save(joinpath(siteresultdir,"OD_inpaint.png"),od)

od = Gray.(load(joinpath(siteresultdir, "OD_nvidia_inpaint.png")))
od = imresize(gray.(od),h,w)


Gray.(odc1)
## OD border
od = imfilter(od,Kernel.gaussian(5))
od = adjust_histogram(od, AdaptiveEqualization(nbins = 256, rblocks = 2, cblocks = 2, clip = 0.5))
odb = detect_edges(od, Canny(spatial_scale=15, high=ImageEdgeDetection.Percentile(80),low=ImageEdgeDetection.Percentile(70)))
save(joinpath(siteresultdir,"OD_border.png"),odb)
odbm = map(i->GrayA(i,i),odb)
save(joinpath(siteresultdir,"OD_border_mask.png"),odbm)

odct = contours(1:h,1:w,od,[0.25,0.75])
odc1 = drawcontour!(similar(od),Contour.levels(odct)[1])
odc2 = drawcontour!(similar(od),Contour.levels(odct)[2])
save(joinpath(siteresultdir,"OD_contour_left.png"),odc1)
save(joinpath(siteresultdir,"OD_contour_right.png"),odc2)

odcm1 = drawcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[1],color=RGBA{N0f8}(1,1,1,1))
odcm2 = drawcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[2],color=RGBA{N0f8}(1,1,1,1))
save(joinpath(siteresultdir,"OD_contour_mask_left.png"),odcm1)
save(joinpath(siteresultdir,"OD_contour_mask_right.png"),odcm2)

skid = pyimport("skimage.draw")
odcf1 = fillcontour!(similar(od),Contour.levels(odct)[1],skid)
odcf2 = fillcontour!(similar(od),Contour.levels(odct)[2],skid)
save(joinpath(siteresultdir,"OD_contourfill_left.png"),odcf1)
save(joinpath(siteresultdir,"OD_contourfill_right.png"),odcf2)

odcfm1 = fillcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[1],skid,color=RGBA{N0f8}(1,0,0,0.02))
odcfm2 = fillcontour!(similar(od,RGBA{N0f8}),Contour.levels(odct)[2],skid,color=RGBA{N0f8}(0,1,0,0.02))
save(joinpath(siteresultdir,"OD_contourfill_mask_left.png"),odcfm1)
save(joinpath(siteresultdir,"OD_contourfill_mask_right.png"),odcfm2)



## Cone Domain
eye,color,mincolor,maxcolor,cone = load(joinpath(siteresultdir,"AG1_V1V2_Full_ISICycle2Color_3","isi.jld2"),"eye","color","mincolor","maxcolor","F1polarity")
cone = adjust_histogram(cone, AdaptiveEqualization(nbins = 256, rblocks = 30, cblocks = 30, clip = 0.1))
clamp01!(cone)
save(joinpath(siteresultdir,"$(color)_COFD.png"),cone)

## Cone Domain contour
cone = imfilter(cone,Kernel.gaussian(5))
cone = adjust_histogram(cone, AdaptiveEqualization(nbins = 256, rblocks = 2, cblocks = 2, clip = 0.85))
nlevel = 40
cdct = contours(1:h,1:w,cone,nlevel)

ilevel = 12
cdc = drawcontour!(similar(cone),levels(cdct)[ilevel])
save(joinpath(siteresultdir,"$(color)_COFD_contour_$(round(ilevel/nlevel,digits=2)).png"),cdc)

cdcm = drawcontour!(similar(cone,RGBA{Float64}),levels(cdct)[ilevel],color=mincolor)
save(joinpath(siteresultdir,"$(color)_COFD_contour_$(round(ilevel/nlevel,digits=2))_mask.png"),cdcm)

cdcf = fillcontour!(similar(cone),levels(cdct)[ilevel],skid)
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_$(round(ilevel/nlevel,digits=2)).png"),cdcf)

cdcfm = fillcontour!(similar(cone,RGBA{Float64}),levels(cdct)[ilevel],skid,color=coloralpha(mincolor,0.1))
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_$(round(ilevel/nlevel,digits=2))_mask.png"),cdcfm)

ilevel = 28
cdc = drawcontour!(similar(cone),levels(cdct)[ilevel])
save(joinpath(siteresultdir,"$(color)_COFD_contour_$(round(ilevel/nlevel,digits=2)).png"),cdc)

cdcm = drawcontour!(similar(cone,RGBA{Float64}),levels(cdct)[ilevel],color=maxcolor)
save(joinpath(siteresultdir,"$(color)_COFD_contour_$(round(ilevel/nlevel,digits=2))_mask.png"),cdcm)

cdcf = fillcontour!(similar(cone),levels(cdct)[ilevel],skid)
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_$(round(ilevel/nlevel,digits=2)).png"),cdcf)

cdcfm = fillcontour!(similar(cone,RGBA{Float64}),levels(cdct)[ilevel],skid,color=coloralpha(maxcolor,0.1))
save(joinpath(siteresultdir,"$(color)_COFD_contourfill_$(round(ilevel/nlevel,digits=2))_mask.png"),cdcfm)
