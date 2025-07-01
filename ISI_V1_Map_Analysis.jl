using NeuroAnalysis,FileIO,JLD2,HypothesisTests,Images,ImageBinarization,ImageEdgeDetection,ImageMorphology,
StatsBase,CairoMakie,Contour,ImageDraw,DataFrames,XLSX,PairPlots

resultroot = "Z:/"
isidir = joinpath(resultroot,"ISI_V1_Map");mkpath(isidir)
alignment = DataFrame(XLSX.readtable(joinpath(resultroot,"ISIMapAlignment.xlsx"),"Sheet1"))
select!(alignment,Not([:Subject_ID,:ppmm]),:Subject_ID=>:subject,:ppmm=>(i->i/1000)=>:ppum)
transform!(alignment,Not([:subject,:siteid,:ppum,:bvname]).=>ByRow(x->ismissing(x) ? [0.,0] : parse.(Float64,split(strip(x,['(',')']),", "))).=>identity)


cs_cofd = ["DKL_X", "LMS_X", "LMS_Y", "LMS_Z", "DKL_RG", "DKL_BY"]
function loadcofd(dir;cs=cs_cofd)
    files = [c=>matchfile(Regex("cofd_map_$(c)\\w*.jld2");dir,join=true) for c in cs]
    ds = Dict(k=>isempty(v) ? nothing : load(v[1]) for (k,v) in files)
end

odmaps=Dict(r.subject=>load(joinpath(resultroot,r.subject,r.siteid,"od_map.jld2")) for r in eachrow(alignment))
orimaps = Dict(r.subject=>load(joinpath(resultroot,r.subject,r.siteid,"ori_map.jld2")) for r in eachrow(alignment))
cofdmaps=Dict(r.subject=>loadcofd(joinpath(resultroot,r.subject,r.siteid)) for r in eachrow(alignment))

cofdcode = Dict("DKL_X"=>"A","DKL_Z"=>"S","LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_RG"=>"L","DKL_BY"=>"S")
function colorminmax(x,c)
    ismissing(x) && return x
    if c=="A"
        s = gray.(Gray.(values(x)))
    elseif c=="L"
        s = comp1.(LMS.(values(x)))
    elseif c=="M"
        s = comp2.(LMS.(values(x))) 
    elseif c=="S"
        s = comp3.(LMS.(values(x))) 
    end
    s[2] > s[1]
end
getmap = (subject,type;c="DKL_X",roi=nothing,sigma_lhi=150) -> begin
    ppum = alignment.ppum[findfirst(alignment.subject.==subject)]
    if type=="od"
        m = odmaps[subject][type]
        a = alignment[findfirst(alignment.subject.==subject),type]
        isnothing(roi) || (m = m[map((i,j)->round.(Int,i.-j),roi,a)...])
    elseif type=="ori"
        amap = orimaps[subject]["amap"]
        mmap = orimaps[subject]["mmap"]
        # cmap = orimaps[subject]["cmap"]
        # mmap = clampscale(abs.(cmap),3)
        a = alignment[findfirst(alignment.subject.==subject),type]
        if isnothing(roi)
            m = (amap,mmap,nothing)
        else
            idx = map((i,j)->round.(Int,i.-j),roi,a)
            amap = amap[idx...];mmap = mmap[idx...]
            lhimap = [localhomoindex(amap,[i,j],σ=sigma_lhi*ppum,n=2).lhi for i in axes(amap,1),j in axes(amap,2)]
            m = (amap,mmap,lhimap)
        end
    elseif type=="cofd"
        t = cofdmaps[subject][c]
        m = t[type]
        colorminmax(t["minmaxcolor"],cofdcode[c]) || (m = 1 .- m)
        a = alignment[findfirst(alignment.subject.==subject),c]
        isnothing(roi) || (m = m[map((i,j)->round.(Int,i.-j),roi,a)...])
    end
    m
end

subject = "AG5"
# roi=nothing
roi=[700:1200,1200:1700]
map_od= getmap(subject,"od";roi)
map_ori,map_mori,map_lhiori = getmap(subject,"ori";roi)
map_cofd_A = getmap(subject,"cofd",c="DKL_X";roi)
map_cofd_L = getmap(subject,"cofd",c="LMS_X";roi)
map_cofd_M = getmap(subject,"cofd",c="LMS_Y";roi)
map_cofd_S = getmap(subject,"cofd",c="LMS_Z";roi)

roistr(roi) = replace(join(string.(roi),'_'),':'=>'-')
displaymap = (;size=(1000,800),dir=nothing,figfmt=[".png"]) -> begin
    f = Figure(;size)
    a,p =image(f[1,1],map_od',axis=(aspect=DataAspect(),title="OD",yreversed=true))
    hidedecorations!(a)
    a,p=image(f[1,2],map(a->HSV(rad2deg(2a),1,1),map_ori)',axis=(aspect=DataAspect(),title="ORI",yreversed=true))
    hidedecorations!(a)
    a,p=image(f[1,3],map_mori',axis=(aspect=DataAspect(),title="ORI_R",yreversed=true))
    hidedecorations!(a)
    if !isnothing(map_lhiori)
        a,p=image(f[1,4],map_lhiori',axis=(aspect=DataAspect(),title="ORI_LHI",yreversed=true))
        hidedecorations!(a)
    end

    a,p=image(f[2,1],map_cofd_A',axis=(aspect=DataAspect(),title="COFD_A",yreversed=true))
    hidedecorations!(a)
    a,p=image(f[2,2],map_cofd_L',axis=(aspect=DataAspect(),title="COFD_L",yreversed=true))
    hidedecorations!(a)
    a,p=image(f[2,3],map_cofd_M',axis=(aspect=DataAspect(),title="COFD_M",yreversed=true))
    hidedecorations!(a)
    a,p=image(f[2,4],map_cofd_S',axis=(aspect=DataAspect(),title="COFD_S",yreversed=true))
    hidedecorations!(a)
    isnothing(dir) ? f : foreach(ext->save(joinpath(dir,"$(subject)_$(roistr(roi))$ext"),f),figfmt)
end

displaymap()
displaymap(dir=isidir)


df = DataFrame(od=vec(map_od),ori=rad2deg.(vec(map_ori)),ori_r=vec(map_mori),ori_lhi=vec(map_lhiori),
        cofd_A=vec(map_cofd_A),cofd_L=vec(map_cofd_L),cofd_M=vec(map_cofd_M),cofd_S=vec(map_cofd_S))
plt = pairplot(df=>(PairPlots.Hist(colormap=:turbo),
            PairPlots.MarginHist()))
save(joinpath(isidir,"$(subject)_$(roistr(roi))_pp.png"),plt)















