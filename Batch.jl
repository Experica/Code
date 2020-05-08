using NeuroAnalysis,Query,FileIO,ProgressMeter,Logging,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase#,ePPR
import Base: close

includet("Batch_Ripple.jl")
includet("Batch_SpikeGLX.jl")
includet("Batch_Scanbox.jl")

function batchtests(tests::DataFrame,param::Dict{Any,Any}=Dict{Any,Any}();log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),plot::Bool=true)
    p = ProgressMeter.Progress(size(tests,1),desc="Batch Tests ... ")
    for t in eachrow(tests)
        try
            if t[:ID]=="OriGrating"
                u,c=processori(dataset,resultroot,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="Laser"
                u,c=processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="Image"
                u,c=processimage(dataset,resultroot,env,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="LaserImage"
                u,c=processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID] in ["Flash","Flash2Color"] && t[:sourceformat]=="SpikeGLX"
                process_flash_spikeglx(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["Hartley","HartleySubspace"] && t[:sourceformat]=="SpikeGLX"
                process_hartley_spikeglx(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["OriSF","OriSFColor","Color"] && t[:sourceformat]=="SpikeGLX"
                process_condtest_spikeglx(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["DirSF"] && t[:sourceformat]=="Scanbox"
                process_2P_dirsf(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["DirSFColor"] && t[:sourceformat]=="Scanbox"
                process_2P_dirsfcolor(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["Hartley"] && t[:sourceformat]=="Scanbox"
                process_2P_hartleySTA(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
                # process_2P_hartleyFourier(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            end
        catch exc
            display("============================================")
            display("Error In Processing: $(t[:files])")
            display("============================================")
            display.(stacktrace(catch_backtrace()))
        end
        next!(p)
    end
    close(log)
end

function close(log::Dict{Any,AbstractLogger})
    for l in values(log)
        close(l.stream)
    end
end

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
