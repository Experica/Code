includet("Online/Online_Imager.jl")

## Manual online processing
resultroot = "Y:/"
testroot = "I:/Test/Test_ISI_Full/Test_ISI_Full_ISIEpochOri8_4"
testroot = "I:/Test/Test_ISI_Full/Test_ISI_Full_ISICycle2Color_10"
testroot = "I:/AG1/AG1_V1V2_Full/AG1_V1V2_Full_ISICycle2Color_2"
testroot = "I:/AG1/AG1_V1V2_Full/AG1_V1V2_Full_ISIEpochOri8_1"

online_epoch_imager(testroot,resultroot)
online_cycle_imager(testroot,resultroot)


## Automatic watching and online processing
function serve(rootdir;dt=5,resultroot="Z:/")

    dispatch = (t)->begin
        ex,s = watch_folder(rootdir,0)
        if s.renamed && !s.changed
            expath = joinpath(rootdir,ex)
            isdir(expath) || return
            if contains(ex,"Cycle")
                Threads.@spawn online_cycle_imager($expath,$resultroot)
            elseif contains(ex,"Epoch")
                Threads.@spawn online_epoch_imager($expath,$resultroot)
            end
        end
    end

    mkpath(rootdir)
    printstyled("Start Watching Experiment In Directory: $rootdir\n",color=:magenta,reverse=true)
    Timer(dispatch,0;interval=dt)
end

function unserve(rootdir,ST)
    printstyled("Stop Watching Experiment In Directory: $rootdir\n",color=:magenta,reverse=true)
    close(ST)
    sleep(0.2)
    unwatch_folder(rootdir)
end

resultroot = "Y:/"
rootdir = "I:/Test/Test_Full/"

const ST = serve(rootdir;resultroot);
unserve(rootdir,ST)

