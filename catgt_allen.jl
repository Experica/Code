using IniFile

## CatGT to pre-process flash segment of raw spike_band.dat
catgt = raw"C:\Users\fff00\CatGT-3.4\CatGT.exe"

# meta file template for allen raw spike_band.dat
metatempfile = joinpath(rawdatadir,"flash.imec0.ap.meta")
metatemp = read(Inifile(),metatempfile)

# flash segment in [1200, 1600] sec
@showprogress "CatGT Flash Segment ... " for r in eachrow(sps)
    indir = joinpath(rawdatadir,"$(r.session_id)","$(r.probe_id)")
    rawfile = joinpath(indir,"spike_band.dat")
    apfile = joinpath(indir,"flash_g0_t0.imec0.ap.bin")
    symlink(rawfile,apfile)
    apmetafile = joinpath(indir,"flash_g0_t0.imec0.ap.meta")
    
    # adjust meta info for each raw spike_band.dat
    apfilesize = stat(apfile).size
    set(metatemp,"fileName",apfile)
    set(metatemp,"fileSizeBytes",apfilesize)
    set(metatemp,"imSampRate",r.apfs)
    set(metatemp,"fileTimeSecs",apfilesize/2/384/r.apfs)
    open(apmetafile, "w+") do io
        write(io, metatemp)
    end

    catgtargs =
    """
    -dir=$indir -run=flash -no_run_fld -prb=0 \
    -g=0 -t=0 -zerofillmax=0 \
    -no_auto_sync \
    -startsecs=1200 -maxsecs=400 \
    -ap -apfilter=butter,12,500,0 -loccar=4,32 \
    """

    r.probe_id in [848037568,841431754,841435557] && (catgtargs *= "-chnexcl={0;265}") # bad channel for these probes

    run(Cmd(`$catgt $catgtargs`,windows_verbatim=true))
    rm(apfile)
    rm(apmetafile)
end

