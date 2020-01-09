## Prepare Metadata and Batch Param
includet("batch.jl")
param = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :imageroot => "../NaturalStimuli")
meta = readmeta("$(param[:dataexportroot])/metadata.mat")

param[:layer] = layer



## Query Tests
tests = @from i in meta begin
        @where startswith(get(i.Subject_ID), "AF5")
        @where i.RecordSite == "ODL3"
        @where i.ID == "HartleySubspace"
        @where i.sourceformat == "SpikeGLX"
        @select {i.ID,i.UUID,i.files}
        @collect DataFrame
        end



## Flash of Two Colors
param = copy(env)
r,c = batchtests(tests,param,isplot=true)

## Hartley Subspace Parametric and Image Response
param[:model]=[:STA]

param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100
r,c = batchtests(tests,param,isplot=true)

## Drifting Grating with Ori and SpatialFreq
param = copy(env)
r,c = batchtests(tests,param,isplot=true)

# Drifting Grating with Ori, SpatialFreq and Color
param = copy(env)
param[:responsedelay] = 0.015
param[:blank] = (:ColorID, 36)
param[:corrcond] = :allcond
r,c = batchtests(tests,param,isplot=true)

## Color Test
param = copy(env)
r,c = batchtests(tests,param,isplot=true)



## Save Result
resultdir = joinpath(rroot,testtype)
!isdir(resultdir) && mkpath(resultdir)
save(joinpath(resultdir,"batchresult.jld2"),"ur",ur,"uc",uc)
