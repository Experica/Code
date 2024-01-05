using NeuroAnalysis,FileIO,JLD2,DataFrames,StatsPlots,StatsBase,ProgressMeter,VegaLite,Dierckx,Combinatorics,XLSX

ccode = Dict("DKL_X"=>"A","DKL_Y"=>"Y","DKL_Z"=>"S","LMS_Xmcc"=>"L","LMS_Ymcc"=>"M","LMS_Zmcc"=>"S",
             "LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")

getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function collectcondtest(indir;id="OriSF",datafile="factorresponse.jld2")
    rs=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            exenv = load(joinpath(root,datafile),"exenv")
            (haskey(exenv,"ID") && exenv["ID"] == id) || continue
            push!(rs,load(joinpath(root,datafile)))
        end
    end
    return rs
end

function joinorisf(rs;ccode=ccode)
    siteid = rs[1]["siteid"]
    fa = rs[1]["fa"]

    dataset = Dict("siteid"=>siteid,"fa"=>fa)
    colors = [r["exenv"]["color"] for r in rs]
    colorcodes = map(i->getccode(i,ccode),colors)
    corder = sortperm(colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    rs = rs[corder]
    dataset["minmaxcolor"] = map(r->r["exenv"]["minmaxcolor"],rs)
    tfs = map(r->r["exenv"]["TemporalFreq"],rs)
    if allequal(tfs)
        tfs = tfs[1]
    else
        @warn "Different TemporalFreq: $tfs for Blocks of OriSF."
    end
    dataset["tf"]=tfs
    gts = map(r->r["exenv"]["GratingType"],rs)
    if allequal(gts)
        gts = gts[1]
    else
        @warn "Different GratingType: $gts for Blocks of OriSF."
    end
    dataset["GratingType"]=gts
    cds = map(r->r["exenv"]["conddur"],rs)
    if allequal(cds)
        cds = cds[1]
    else
        @warn "Different CondDur: $cds for Blocks of OriSF."
    end
    dataset["conddur"]=cds

    uns = map(r->length(r["unitid"]),rs)
    @info "Number of units for each OriSF test: $uns"
    uids = mapreduce(r->r["unitid"],intersect,rs) # only include units that spikes in all tests
    uis = map(r->indexin(uids,r["unitid"]),rs)
    dataset["ugood"] = Dict(uids.=>rs[1]["unitgood"][uis[1]])

    getru = (n;k=nothing) -> isnothing(k) ? map((r,i)->r[n][i],rs,uis) : map((r,i)->r[n][k][i],rs,uis)
    dataset["frs"] = Dict(uids.=>zip(getru("frs")...))
    dataset["responsive"] = Dict(uids.=>zip(getru("responsive")...))
    dataset["modulative"] = Dict(uids.=>zip(getru("modulative")...))
    dataset["enoughresponse"] = Dict(uids.=>zip(getru("enoughresponse")...))
    dataset["maxi"] = Dict(uids.=>zip(getru("maxi")...))
    dataset["fpsth"] = (psth=Dict(uids.=>zip(getru("fpsth";k=:psth)...)), x = rs[1]["fpsth"].x)
    dataset["frf"] = NamedTuple(k=>Dict(uids.=>zip(getru("frf";k)...)) for k in keys(fa))
    
    return dataset
end

function f1f0orisf!(dataset;btw=(-50,0))
    conddur = dataset["conddur"]
    tf = dataset["tf"]
    fpsth,x = dataset["fpsth"]
    bw = x[2]-x[1]
    fs = 1/(bw/1000)
    bi = findall(i->btw[begin]<=i<btw[end],x)
    ci = findall(i->0<=i<conddur,x)
    getF1F0 = (psth)->begin
        y = psth[ci].-mean(psth[bi])
        F1F0 = dft(y,fs,tf,0)
        m0 = abs(F1F0[2])
        m1 = abs(F1F0[1])
        a1 = mod2pi(angle(F1F0[1]))
        (;m0,m1,a1,f10=m1/m0)
    end
    f1f0 = Dict(u=>map(c->getF1F0.(c.m),v) for (u,v) in fpsth)
    f10 = Dict(u=>map(c->getproperty.(c,:f10),v) for (u,v) in f1f0)
    f1phase = Dict(u=>map(c->getproperty.(c,:a1),v) for (u,v) in f1f0)
    f1mag = Dict(u=>map(c->getproperty.(c,:m1),v) for (u,v) in f1f0)
    f1maxi = Dict(u=>map(c->[Tuple(argmax(getproperty.(c,:m1)))...],v) for (u,v) in f1f0)

    dataset["f1f0"] = f1f0
    dataset["f10"] = f10
    dataset["f1phase"] = f1phase
    dataset["f1mag"] = f1mag
    dataset["f1maxi"] = f1maxi
    return dataset
end

plotpsth = (dataset,u;isse=true,showfactor=true,dir=nothing,figfmt=[".png"],cs=["A","L","M","S"])->begin
    ccode = dataset["ccode"]
    cs = intersect(cs,ccode)
    ci = indexin(cs,ccode)
    cn = length(ci)

    conddur = dataset["conddur"]
    fa = dataset["fa"]
    fan1,fan2 = length.(values(fa))
    fpsth,x = dataset["fpsth"]
    psth = fpsth[u][ci]
    maxi = dataset["maxi"][u][ci]
    f1maxi = dataset["f1maxi"][u][ci]
    f1phase = dataset["f1phase"][u][ci]
    ug = dataset["ugood"][u]

    colors = map(c->c.maxcolor,dataset["minmaxcolor"][ci])
    "A" in cs && (colors[findfirst("A".==cs)]=RGBA(0.3,0.3,0.3,1))

    p = plot(;layout=(fan1,fan2),grid=false,legend=false,size=(200fan2,100fan1),
            link=:all,tickdir=:out,xlabel="Time (ms)",ylabel="Response (spike/s)")
    for j in 1:fan1, i in 1:fan2
        leftmargin = i == 1 ? 3Plots.mm : -18Plots.mm
        bottommargin = j==fan1 ? 3Plots.mm : -16Plots.mm
        topmargin = j==1 ? 6Plots.mm : :match
        rightmargin = i==fan2 ? 8Plots.mm : :match
        frame = (i ==1 && j==fan1) ? :auto : :none
        
        vspan!(p[j,i],[0,conddur];color=:gray92,leftmargin,bottommargin,topmargin,rightmargin,frame)
        for c in 1:cn
            ribbon = isse ? psth[c].se[j,i] : nothing
            y = psth[c].m[j,i]
            ann = ((0.1,0.9),("-",20,colors[c],:left,:vcenter,rad2deg(f1phase[c][j,i])))
            plot!(p[j,i],x,y;ribbon,color=colors[c],ann)
        end
    end
    
    ax = 0.06.*(0:cn-1)
    ax = ax .- mean(ax) .+ 0.5
    foreach((i,c,x)->annotate!(p[i...],((x,0.96),("∘",16,c,:center))),maxi,colors,ax)
    foreach((i,c,x)->annotate!(p[i...],((x,0.99),("▵",16,c,:center))),f1maxi,colors,ax)
    if showfactor
        foreach(i->annotate!(p[i,fan2],((1.1,0.5),("-ꜛ-",20,:center,fa[1][i]))),eachindex(fa[1]))
        foreach(i->annotate!(p[1,i],((0.5,1.2),("SF=$(round(fa[2][i],digits=1))",12,:center))),eachindex(fa[2]))
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ug ? "S" : "M")U$(u)_PSTH$ext")),figfmt)
end

orisfinfo = (siteresultdir;figfmt=nothing) -> begin
    rs = collectcondtest(siteresultdir,id="OriSF")
    dataset = joinorisf(rs)
    dataset = f1f0orisf!(dataset)
    jldsave(joinpath(siteresultdir,"orisfdataset.jld2");dataset)
    if !isnothing(figfmt)
        dir = joinpath(siteresultdir,"orisf");mkpath(dir)
        for u in sort(collect(keys(dataset["ugood"])))
            plotpsth(dataset,u;dir,figfmt)
        end
    end
end

resultroot = "Z:/"

## process all orisfs of a RecordSite
subject = "AG1";recordsession = "V1";recordsite = "ODL1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
# dataset = load(joinpath(siteresultdir,"orisfdataset.jld2"),"dataset")


plotpsth(dataset,748)
plotpsth(dataset,748,isse=false,cs=["L","M"])

## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch All OriSFs ... " for r in eachrow(penetration)
    # orisfinfo(joinpath(resultroot,r.Subject_ID,r.siteid);figfmt=[".png",".svg"])
    orisfinfo(joinpath(resultroot,r.Subject_ID,r.siteid))
end



