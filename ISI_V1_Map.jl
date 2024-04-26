using NeuroAnalysis,FileIO,JLD2,HypothesisTests,Images,ImageBinarization,ImageEdgeDetection,ImageMorphology,
StatsBase,StatsPlots,Contour,ImageDraw,DataFrames,XLSX

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
subject = "AG6";recordsession = "";recordsite = "Full"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)


## orientation map
function orimap(siteresultdir,test;datafile="isi.jld2",filter=x->dogfilter(x,hσ=3,lσ=25))
    opt,op = load(joinpath(siteresultdir,test,datafile),"opt","op")
    os = vec(stack(op)')
    complexmap([opt;opt],2deg2rad.(os);rsign=repeat([-1,1],inner=length(opt)),filter)
end

cmap,amap,mmap = orimap(siteresultdir,"$(siteid)_ISIEpochOri8_0",filter=x->dogfilter(x,hσ=5,lσ=25))
orianglemap = map(a->HSV(rad2deg(a),1,1),amap)
oripolarmap = map((a,m)->HSV(rad2deg(a),1,m),amap,mmap)

save(joinpath(siteresultdir,"ori_anglemap.png"),orianglemap)
save(joinpath(siteresultdir,"ori_polarmap.png"),oripolarmap)


## orientation local homogeneity index of penetrations
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
pidx = findall(r->r.Subject_ID==subject && !ismissing(r.coord),eachrow(penetration))
pcoord = map(i->parse.(Int,split(i,", ")),penetration.coord[pidx])
ppum = penetration.ppmm[pidx[1]]/1000 # pixels/micrometer
rs = map(i->localhomoindex(amap,i;σ=120ppum),pcoord)

XLSX.openxlsx(joinpath(resultroot,"penetration.xlsx"),mode="rw") do f
    s = f["Sheet1"]
    ci = findfirst(i->i=="ori_lhi",names(penetration))
    for (i,r) in zip(pidx,rs)
        s[i+1,ci] = first(r)
    end
end

pmapdir = joinpath(resultroot,"PenetrationMap");mkpath(pmapdir)
poamap = [alphamask(map(a->HSV(rad2deg(a),1,1),last(r)))[1] for r in rs]
foreach(i->i[map(c->c:c+1,floor.(Int,size(i)./2))...].=HSVA(0.0,0,0.1,1),poamap)
foreach((m,id)->save(joinpath(pmapdir,"$(id)_ori.png"),m),poamap,penetration.siteid[pidx])


## ocular dominance
function odmap(siteresultdir,lrtests::NTuple{2,String};datafile="isi.jld2")
    ds = load.(joinpath.(siteresultdir,lrtests,datafile))
    eyes = [i["exenv"]["eye"] for i in ds];eyeis = sortperm(eyes);eyesort=eyes[eyeis]
    eyesort == ["Left","Right"] || @warn "Tests aren't from Left and Right eyes, got $eyes instead."

    ler = ds[eyeis[1]]["epochresponse"];rer = ds[eyeis[2]]["epochresponse"]
    pairtest(ler,rer).stat
end

t = odmap(siteresultdir,("$(siteid)_ISIEpochOri8_0","$(siteid)_ISIEpochOri8_1"))
od = clampscale(dogfilter(t,lσ=50),3)
h,w = size(od)
save(joinpath(siteresultdir,"OD.png"),od)



## blood vessel and mask
function readframe(dir,name;w=2080,h=2080,siteresultdir=nothing,whichframe=last)
    file = matchfile(Regex("$name-Frame\\d*.Raw");dir,join=true) |> whichframe
    f = clampscale(readrawim_Mono12Packed(file,w,h))
    isnothing(siteresultdir) || save(joinpath(siteresultdir,"$name.png"),f)
    f
end
function bvmask(bvdir,bvname;h=2080,w=2080,siteresultdir=nothing,r=nothing)
    bv = readframe(bvdir,bvname;w,h,siteresultdir)
    bve = clampscale(dogfilter(bv),3)
    bvm = binarize(Bool,bve,AdaptiveThreshold(bve))
    if !isnothing(r)
        foreach(_->dilate!(bvm),1:r)
        foreach(_->erode!(bvm),1:r+2)
    end
    if !isnothing(siteresultdir)
        save(joinpath(siteresultdir,"$(bvname)_ahe.png"),ahe(bv,clip=0.9))
        save(joinpath(siteresultdir,"$(bvname)_mask.png"),bvm)
    end
    bvm
end

bvdir = joinpath(dataroot,subject,"BloodVessel")

bvgm = bvmask(bvdir,"BloodVessel_Green";siteresultdir)
bvrm = bvmask(bvdir,"BloodVessel_Red";siteresultdir,r=4)

bvgm = bvmask(bvdir,"BloodVessel_Green_After";siteresultdir,r=6)
bvrm = bvmask(bvdir,"BloodVessel_Red_After";siteresultdir,r=4)

readframe(bvdir,"BloodVessel_Green_After_Scale";siteresultdir)



## OD border and contour
# od mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
odmask = gray.(load(joinpath(siteresultdir,"mask_OD.png"))) .==  1
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

cofd,color,minmaxcolor = cofdmap(siteresultdir,"$(siteid)_ISIEpochFlash2Color_0")
cofd = clampscale(dogfilter(cofd),3)
save(joinpath(siteresultdir,"COFD_$(color).png"),cofd)


## COFD contour
# cofd mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
cofdmask = gray.(load(joinpath(siteresultdir,"mask_COFD.png"))) .== 1
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

