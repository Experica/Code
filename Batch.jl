using NeuroAnalysis,FileIO,JLD2,Mmap,DataFrames,Statistics,ProgressMeter,Logging
import Base: close

# includet("Batch_Ripple.jl")
includet("Batch_Imager.jl")
includet("Batch_SpikeGLX.jl")
# includet("Batch_Scanbox.jl")

function batchtests(tests::DataFrame,param::Dict{Any,Any}=Dict{Any,Any}();log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),plot::Bool=true)
    p = ProgressMeter.Progress(size(tests,1),desc="Batch ... ",start=-1)
    for t in eachrow(tests)
        next!(p,showvalues = [(:Test, t.files)])
        try
            if t.ID=="OriGrating" && t.sourceformat=="Ripple"
                u,c=processori(dataset,resultroot,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t.ID=="Laser" && t.sourceformat=="Ripple"
                u,c=processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t.ID=="Image" && t.sourceformat=="Ripple"
                u,c=processimage(dataset,resultroot,env,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t.ID=="LaserImage" && t.sourceformat=="Ripple"
                u,c=processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t.ID in ["Flash","Flash2Color"] && t.sourceformat=="SpikeGLX"
                process_flash_spikeglx(t.files,param;uuid=t.UUID,log,plot)
            elseif t.ID in ["Hartley","HartleySubspace","Image"] && t.sourceformat=="SpikeGLX"
                process_hartley_spikeglx(t.files,param;uuid=t.UUID,log,plot)
            elseif t.ID in ["OriSF","OriSFColor","Color"] && t.sourceformat=="SpikeGLX"
                process_condtest_spikeglx(t.files,param;uuid=t.UUID,log,plot)
            elseif t.ID in ["CycleColorPlane"] && t.sourceformat=="SpikeGLX"
                process_cycle_spikeglx(t.files,param;uuid=t.UUID,log,plot)
            elseif t.ID in ["DirSF"] && t.sourceformat=="Scanbox"
                process_2P_dirsf(t.files,param,uuid=t.UUID,log=log,plot=plot)
            elseif t.ID in ["DirSFColor"] && t.sourceformat=="Scanbox"
                process_2P_dirsfcolor(t.files,param,uuid=t.UUID,log=log,plot=plot)
            elseif t.ID in ["Hartley"] && t.sourceformat=="Scanbox"
                # process_2P_hartleySTA(t.files,param,uuid=t.UUID,log=log,plot=plot)
                process_2P_hartleyFourier(t.files,param,uuid=t.UUID,log=log,plot=plot)
            elseif t.ID in ["ISICycle2Color","ISICycleOri"] && t.sourceformat=="Imager"
                process_cycle_imager(t.files,param;uuid=t.UUID,log,plot)
            elseif t.ID in ["ISIEpochOri8"] && t.sourceformat=="Imager"
                process_epoch_imager(t.files,param;uuid=t.UUID,log,plot)
            else
                display("Skip `$(t.ID)` from `$(t.sourceformat)`, no corresonding process function.")
            end
        catch exc
            display("============================================")
            display("Error In Processing: $(t.files)")
            display("============================================")
            display.(stacktrace(catch_backtrace()))
        end
    end
    finish!(p)
    close(log)
end

close(log::Dict{Any,AbstractLogger}) = foreach(l->close(l.stream),values(log))

function (log::Dict{Any,AbstractLogger})(key,f,logfile)
    if !haskey(log,key)
        log[key] = SimpleLogger(open(logfile,"w"))
    end
    logger = log[key]
    if f isa String
        println(logger.stream,f)
    else
        with_logger(f,logger)
    end
    flush(logger.stream)
end
