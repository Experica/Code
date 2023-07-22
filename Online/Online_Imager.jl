
using FileWatching

function online_epoch_imager(testroot,resultroot;dt=5,tout=15)
    test = splitdir(testroot)[end]
    printstyled("Start Online Analysis: $test\n",color=:yellow,reverse=true)
    exfile = joinpath(testroot,"$test.yaml")
    metafile = joinpath(testroot,"$test.meta")

    while true
        timedwait(()->isfile(metafile),tout;pollint=dt) == :ok && break
        printstyled("No experiment and meta data detected, continue checking? (y|n)> ",color=:blue)
        if readline() != "y"
            printstyled("Abort Online Analysis: $test\n",color=:red,reverse=true)
            return
        end
    end

    printstyled("Preparing experiment and meta data ...\n",color=:cyan)

    ei = 0
    while true
        timedwait(()->isfile(joinpath(testroot,".Epoch$ei")),tout;pollint=dt) == :ok || break
        printstyled("Processing Epoch$ei ...\n",color=:red)
        epochdir = joinpath(testroot,"Epoch$ei")

        ei+=1
    end

    printstyled("Finish Online Analysis: $test\n",color=:green,reverse=true)
end

function online_cycle_imager(testroot,resultroot;dt=5,tout=15)
    test = splitdir(testroot)[end]
    printstyled("Start Online Analysis: $test\n",color=:yellow,reverse=true)
    exfile = joinpath(testroot,"$test.yaml")
    metafile = joinpath(testroot,"$test.meta")

    while true
        timedwait(()->isfile(metafile),tout;pollint=dt) == :ok && break
        printstyled("No experiment and meta data detected, continue checking? (y|n)> ",color=:blue)
        if readline() != "y"
            printstyled("Abort Online Analysis: $test\n",color=:red,reverse=true)
            return
        end
    end

    printstyled("Preparing experiment and meta data ...\n",color=:cyan)
    epochdir = joinpath(testroot,"Epoch0")
    fps = 10
    fi = 0
    df = floor(dt*fps)
    fmt = "Raw"

    while true
        if timedwait(()->isfile(joinpath(epochdir,"$test-Frame$(fi+df).$fmt")),tout;pollint=dt) == :ok
            fr = fi:fi+df
        else
            lf = length(readdir(epochdir))-1
            if lf>=fi
                fr = fi:lf
            else
                break
            end
        end
        printstyled("Processing Frames $fr ...\n",color=:red)

        fi=fr[end]+1
    end

    printstyled("Finish Online Analysis: $test\n",color=:green,reverse=true)
end



