includet("Batch.jl")

## Prepare Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :stimuliroot => "../NaturalStimuli")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))

## Query Tests
tests = select!(filter(r->begin
                    startswith(r.Subject_ID, "AF5") &&
                    r.RecordSite == "ODL5" &&
                    r.ID == "Color" &&
                    r.sourceformat == "SpikeGLX"
                    end,meta),
                [:ID,:UUID,:files,:sourceformat])

## Batch Condition Tests
batchtests(tests,param,plot=true)

## Batch HartleySubspace Tests
param[:model]=[:STA]
param[:ppd] = 45
param[:blank] = (:Ori_Final,NaN)

param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100

batchtests(tests,param,plot=false)

## Collect units
function collectunit(indir;cells=DataFrame(),ufile="spike.jld2")
    for (root,dirs,files) in walkdir(indir)
        if ufile in files
            spike = load(joinpath(root,ufile),"spike")
            siteid = spike["siteid"]
            ui = findall(spike["unitgood"])
            cs = DataFrame(id=["$(siteid)_SU$u" for u in spike["unitid"][ui]],site=siteid,
             position = [spike["unitposition"][r:r,:] for r in ui])
            append!(cells,cs)
        end
    end
    return unique!(cells,:id)
end

cells = collectunit(joinpath(param[:resultroot],"AF5"))
save(joinpath(param[:resultroot],"cells.jld2"),"cells",cells)
