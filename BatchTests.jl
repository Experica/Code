# Prepare Env and Metadata
includet("batch.jl")
env = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :imageroot => "../NaturalStimuli")
meta = readmeta("$(env[:dataexportroot])/metadata.mat")

# Query Test
test = "Flash"
tests=querymeta(meta,test=test,sourceformat="SpikeGLX",subject="AE9")

# Flash of Two Colors
param = copy(env)
r,c = batchtests(tests,param,isplot=true)


# Hartley Subspace Parametric and Image Response
param = copy(env)
param[:model]=[:ePPR]
param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100
ur,uc = batchtests(tests[4:4,:],param,isplot=true);

# Drifting Grating with Ori and SpatialFreq
param = copy(env)
param[:responsedelay] = 0.02
param[:blank] = (:Ori, "blank")
param[:corrcond] = :allcond
ur,uc = batchtests(tests[1:3,:],param,isplot=true);

# Drifting Grating with Ori, SpatialFreq and Color
param = copy(env)
param[:responsedelay] = 0.02
param[:blank] = (:ColorID, 36)
param[:corrcond] = :allcond
ur,uc = batchtests(tests[1:1,:],param,isplot=true);

# Save Result
using FileIO
resultdir = joinpath(rroot,testtype)
!isdir(resultdir) && mkpath(resultdir)
save(joinpath(resultdir,"batchresult.jld2"),"ur",ur,"uc",uc)
