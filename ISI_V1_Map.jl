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

function contourmap(img;w=2080,h=2080,cl=[0.25,0.75],name="OD",clname=["left","right"],clcolor=[RGBA{N0f8}(1,0,0,0.1),RGBA{N0f8}(0,1,0,0.1)],dir=nothing)
    ct = contours(1:h,1:w,img,cl)
    maps=Dict()
    for (i,ctl) in enumerate(Contour.levels(ct))
        title = join([name,"contour",clname[i]],'_')
        maps[title]=drawcontour!(similar(img,Gray{Bool}),ctl)
        title = join([name,"contourfill",clname[i]],'_')
        maps[title] = drawcontour!(similar(img,Gray{Bool}),ctl,fill=true)

        title = join([name,"contour_draw",clname[i]],'_')
        maps[title] = drawcontour!(similar(img,RGBA{N0f8}),ctl,convert(RGBA{N0f8}, clcolor[i]))
        title = join([name,"contourfill_draw",clname[i]],'_')
        maps[title] = drawcontour!(similar(img,RGBA{N0f8}),ctl,convert(RGBA{N0f8}, clcolor[i]),fill=true)
    end

    if !isnothing(dir)
        foreach(k->save(joinpath(dir,"$k.png"),maps[k]),keys(maps))
    end
    maps
end

dataroot = "D:/"
resultroot = "Z:/"
subject = "AG5";recordsession = "";recordsite = "Full"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)


## orientation map
function orimap(siteresultdir,test;datafile="isi.jld2",filter=x->dogfilter(x,hσ=3,lσ=25),mnorm=:r,mask=nothing,nsd=3)
    if eltype(test)==String
        opt1,op1 = load(joinpath(siteresultdir,test[1],datafile),"opt","op")
        os1 = vec(stack(op1)')
        opt2,op2 = load(joinpath(siteresultdir,test[2],datafile),"opt","op")
        os2 = vec(stack(op2)')
        return complexmap([opt1;opt1;opt2;opt2],deg2rad.([os1;os2]);n=2,nsd,rsign=repeat([-1,1],inner=length(opt1),outer=2),filter,mnorm,mask)
    end
    opt,op = load(joinpath(siteresultdir,test,datafile),"opt","op")
    os = vec(stack(op)')
    complexmap([opt;opt],deg2rad.(os);n=2,nsd,rsign=repeat([-1,1],inner=length(opt)),filter,mnorm,mask)
end

# ori mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
orimask = gray.(load(joinpath(siteresultdir,"mask_Ori.png"))) .==  1
# cmap,amap,mmap = orimap(siteresultdir,["$(siteid)_ISIEpochOri8_3","$(siteid)_ISIEpochOri8_4"],filter=x->dogfilter(x,hσ=5,lσ=25),mnorm=:r,mask=orimask,nsd=1.5)
cmap,amap,mmap = orimap(siteresultdir,"$(siteid)_ISIEpochOri8_0",filter=x->dogfilter(x,hσ=5,lσ=25),mnorm=:r,mask=orimask,nsd=1.5)
orianglemap = map(a->HSV(rad2deg(2a),1,1),amap)
oripolarmap = map((a,m)->HSV(rad2deg(2a),1,m),amap,mmap)
orimagmap = Gray.(mmap)

save(joinpath(siteresultdir,"Ori_anglemap.png"),orianglemap)
save(joinpath(siteresultdir,"Ori_polarmap.png"),oripolarmap)
save(joinpath(siteresultdir,"Ori_magmap.png"),orimagmap)
jldsave(joinpath(siteresultdir,"ori_map.jld2");cmap,amap,mmap)



## blood vessel and mask
function readframe(dir,name;w=2080,h=2080,siteresultdir=nothing,whichframe=last)
    file = matchfile(Regex("$name-Frame\\d*.Raw");dir,join=true) |> whichframe
    f = clampscale(readrawim_Mono12Packed(file,w,h))
    isnothing(siteresultdir) || save(joinpath(siteresultdir,"$name.png"),f)
    f
end
function bvmap(bvdir,bvname;h=2080,w=2080,siteresultdir=nothing,r=nothing,mask=nothing,nsd=3)
    bv = readframe(bvdir,bvname;w,h,siteresultdir)
    if !isnothing(mask)
        nbv = bv[.!mask]
        bv[.!mask].= randn(length(nbv))*std(nbv) .+ mean(nbv)
    end
    bve = clampscale(dogfilter(bv),nsd;mask)
    ms = ""
    bvm = binarize(Bool,bve,AdaptiveThreshold(bve))
    if !isnothing(r)
        foreach(_->dilate!(bvm),1:r)
        foreach(_->erode!(bvm),1:r+2)
        ms = "_coarse"
    end
    if !isnothing(siteresultdir)
        save(joinpath(siteresultdir,"$(bvname)_dog.png"),bve)
        save(joinpath(siteresultdir,"$(bvname)_ahe.png"),ahe(bv,clip=0.9))
        save(joinpath(siteresultdir,"$(bvname)_mask$ms.png"),bvm)
    end
    bvm
end

# bv mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
bvmask = gray.(load(joinpath(siteresultdir,"mask_BV.png"))) .==  1
bvdir = joinpath(dataroot,subject,"BloodVessel")

bvgm = bvmap(bvdir,"BloodVessel_Green";siteresultdir,mask=bvmask)
bvgm = bvmap(bvdir,"BloodVessel_Green";siteresultdir,r=6,mask=bvmask)
bvrm = bvmap(bvdir,"BloodVessel_Red";siteresultdir,r=5,mask=bvmask)

bvrm = bvmap(bvdir,"BloodVessel_Red_After";siteresultdir,r=5,mask=bvmask)
bvgm = bvmap(bvdir,"BloodVessel_Green_After";siteresultdir,mask=bvmask)
bvgm = bvmap(bvdir,"BloodVessel_Green_After";siteresultdir,r=6,mask=bvmask)

bvs = readframe(bvdir,"BloodVessel_Green_After_Scale";siteresultdir)



## ocular dominance of left(<0) and right(>0) eye
function odmap(siteresultdir,lrtests::NTuple{2,String};datafile="isi.jld2")
    ds = load.(joinpath.(siteresultdir,lrtests,datafile))
    eyes = [i["exenv"]["eye"] for i in ds];eyeis = sortperm(eyes);eyesort=eyes[eyeis]
    eyesort == ["Left","Right"] || @warn "Tests aren't from Left and Right eyes, got $eyes instead."

    ler = ds[eyeis[1]]["epochresponse"];rer = ds[eyeis[2]]["epochresponse"]
    pairtest(ler,rer).stat
end

# od mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
odmask = gray.(load(joinpath(siteresultdir,"mask_OD.png"))) .==  1
t = odmap(siteresultdir,("$(siteid)_ISIEpochOri8_2","$(siteid)_ISIEpochOri8_1"))
od = clampscale(dogfilter(t,lσ=50),2,mask=odmask)
odmagmap = clampscale(gaussianfilter(od,σ=7),1,mask=odmask)
save(joinpath(siteresultdir,"OD.png"),od)
jldsave(joinpath(siteresultdir,"od_map.jld2");od,odmagmap)

# OD border
odmagmap[.!odmask] .= 0.5
ode = ahe(odmagmap,nblock=1,clip=0.9)
odb = detect_edges(ode, Canny(spatial_scale=10, high=ImageEdgeDetection.Percentile(70),low=ImageEdgeDetection.Percentile(30)))
save(joinpath(siteresultdir,"OD_border.png"),odb)
odbd = map(i->GrayA(i,i),odb)
save(joinpath(siteresultdir,"OD_border_draw.png"),odbd)
# OD contour
odmaps = contourmap(ode;dir=siteresultdir)



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

# cofd mask drawn in GIMP/Photoshop to exclude non-relevant region(black)
cofdmask = gray.(load(joinpath(siteresultdir,"mask_COFD.png"))) .== 1
cofd,color,minmaxcolor = cofdmap(siteresultdir,"$(siteid)_ISIEpochFlash2Color_5")
cofd = clampscale(dogfilter(cofd),2,mask=cofdmask)
cofdmagmap = clampscale(gaussianfilter(cofd),2,mask=cofdmask)
save(joinpath(siteresultdir,"COFD_$(color).png"),cofd)
jldsave(joinpath(siteresultdir,"cofd_map_$(color).jld2");cofd,cofdmagmap,minmaxcolor)

# COFD contour for mincolor(dark) and maxcolor(bright) domain
cofde = ahe(cofdmagmap,nblock=30,clip=1)
cofde[.!(cofdmask.&bvgm)] .= 0.5
cofdmaps = contourmap(cofde;dir=siteresultdir,name="COFD_$color",clname=["low","high"],clcolor=minmaxcolor)







