## Prepare Param and Metadata
includet("Batch.jl")

param = Dict{Any,Any}(
    :dataexportroot => "O:\\AF4\\2P_analysis\\Summary\\DataExport",
    :interpolatedData => true,   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    :preOffset => 0.1,
    :responseOffset => 0.05,  # in sec
    :α => 0.05,   # p value
    :sampnum => 100,   # random sampling 100 times
    :fitThres => 0.5,
    :hueSpace => "DKL",   # Color space used? DKL or HSL
    :diraucThres => 0.8,   # if passed, calculate hue direction, otherwise calculate hue axis
    :oriaucThres => 0.5,
    :Respthres => 0.1,  # Set a response threshold to filter out low response cells?
    :blankId => 36,  # Blank Id; AF2AF3AF4=36, AE6AE7=34
    :excId => [27,28])  # Exclude some condition?

meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))
for t in 1:nrow(meta)
        meta.files[t]=meta.files[t][3:end]
end

## Query Tests
tests = @from i in meta begin
        # @where startswith(get(i.Subject_ID), "AF4")
        @where i.Subject_ID == "AF4"
        # @where (i.RecordSite == "u003" || i.RecordSite == "u004")
        # @where i.RecordSite == "u002"
        # @where i.RecordSite != "u010"
        @where i.ID == "Hartley"
        @where i.sourceformat == "Scanbox"
        @select {i.sourceformat,i.ID,i.files,i.RecordSite,i.filename,i.UUID}
        @collect DataFrame
        end

# tests.files = replace.(tests.files, Ref("H" => string(param[:dataexportroot][1])))   # Somtime the disk name changes becasue whatever reason
sort!(tests, [:RecordSite, :filename])
show(tests, allrows = true, allcols = true)
## Condition Tests
batchtests(tests,param,plot=false)

## HartleySubspace Parametric and Image Response
param[:model]=[:STA]
param[:stanorm] = nothing
param[:stawhiten] = nothing
param[:downsample] = 1  # for down-sampling stimuli image, 1 is no down-sampling
param[:hartleynorm] = false
param[:hartleyscale] = 1
param[:hartelyBlkId]=5641
param[:delayLB] = -0.066  # in sec; Usually do not need to change it
param[:delayUB] = 0.4   # in sec; Usually do not need to change it
param[:delayStep] = 0.033  # Bidirectional = 0.033 sec, Unidirectional = 0.066 sec
param[:delays] = param[:delayLB]:param[:delayStep]:param[:delayUB]
# param[:maskradius] = 1 #AE6=0.75 AE7=0.24 AF2=1 AF3=0.6 AF4=0.65

param[:model]=[:ePPR]
param[:downsample] = 2  # for down-sampling stimuli image, 1 is no down-sampling
param[:hartleynorm] = true
param[:hartleyscale] = 1
param[:eppr_ndelay]=1
param[:eppr_nft]=[3]
param[:eppr_lambda]=32  #16

batchtests(tests,param,plot=false)

tests
