## Prepare Env and Metadata
includet("batch.jl")
env = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :imageroot => "../NaturalStimuli")
meta = readmeta("$(env[:dataexportroot])/metadata.mat")



env[:layer] = layer
## Query Test
tests = @from i in meta begin
        @where startswith(get(i.Subject_ID),"AF5")
        @where i.RecordSite=="ODL3"
        @where i.ID=="OriSF"
        @where i.sourceformat=="SpikeGLX"
        @select {i.ID,i.UUID,i.files}
        @collect DataFrame
        end



## Flash of Two Colors
param = copy(env)
r,c = batchtests(tests,param,isplot=true)

# Hartley Subspace Parametric and Image Response
param = copy(env)
param[:model]=[:ePPR]
param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100
ur,uc = batchtests(tests[4:4,:],param,isplot=true);

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
