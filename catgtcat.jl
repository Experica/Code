## concat binary files using CatGT
catgt = raw"C:\Users\fff00\CatGT-win\CatGT.exe"

indir = raw"E:\SpikeGLXData\AG1\AG1_V1_ODL1"
files = readdir(indir)
apfiles = filter(f->occursin(r"AG.*ap.bin",f),files)
apmetafiles = filter(f->occursin(r"AG.*ap.meta",f),files)
order = sortperm(mtime.(joinpath.(indir,apfiles)))
apfiles = apfiles[order]
apmetafiles = apmetafiles[order]
catfile = string(map(f->join(first.(split(f,'_')[end-1:end])),apfiles)...,".imec0.ap.bin")

n = length(apfiles)
catgtapfiles = [replace(apfiles[i],r".*imec0"=>"catgtrun_g$(i-1)_t0.imec0") for i in 1:n]
catgtapmetafiles = [replace(apmetafiles[i],r".*imec0"=>"catgtrun_g$(i-1)_t0.imec0") for i in 1:n]

appath = joinpath.(indir,apfiles)
apmetapath = joinpath.(indir,apmetafiles)
catgtappath = joinpath.(indir,catgtapfiles)
catgtapmetapath = joinpath.(indir,catgtapmetafiles)
catpath = joinpath(indir,catfile)
catmetapath = replace(catpath,"imec0.ap.bin"=>"imec0.ap.meta")
isfile(catpath) && rm(catpath)
isfile(catmetapath) && rm(catmetapath)


# foreach((t,l)->run(`cmd /c mklink $l $t`),appath,catgtappath)
# foreach((t,l)->run(`cmd /c mklink $l $t`),apmetapath,catgtapmetapath)
symlink.(appath,catgtappath)
symlink.(apmetapath,catgtapmetapath)


catgtargs =
"""
-dir=$indir -run=catgtrun -no_run_fld -prb=0 \
-g=0,$(n-1) -t=0 -zerofillmax=0 \
-ap -apfilter=butter,12,300,0 -gblcar \
"""

run(Cmd(`$catgt $catgtargs`,windows_verbatim=true))


rm.(catgtappath)
rm.(catgtapmetapath)
mv(joinpath(indir,"catgtrun_g0_tcat.imec0.ap.bin"),catpath)
mv(joinpath(indir,"catgtrun_g0_tcat.imec0.ap.meta"),catmetapath)
