using NeuroAnalysis,DataFrames,FileIO,Statistics,LightGraphs,MetaGraphs

# Collect all results and merge into database
ccode = Dict("DKL_X"=>'A',"DKL_Y"=>'Y',"DKL_Z"=>'S',"LMS_Xmcc"=>'L',"LMS_Ymcc"=>'M',"LMS_Zmcc"=>'S',
             "LMS_X"=>'L',"LMS_Y"=>'M',"LMS_Z"=>'S',"DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym")



"Collect CondTests and Merge into Cell"
function collectcondtest(indir;cell=DataFrame(id=[]),ccode=ccode,datafile="factorresponse.jld2")
    sdf = Dict()
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            fr = load(joinpath(root,datafile))
            siteid = fr["siteid"]
            haskey(sdf,siteid) || (sdf[siteid]=DataFrame(id=[]))
            id = ["$(siteid)_SU$u" for u in fr["unitid"]]
            c = ccode[fr["color"]]
            tk = ["responsive","modulative","enoughresponse"]
            isempty(fr["f1f0"]) || push!(tk,"f1f0")
            df1 = ("$(k)!$c"=>fr[k] for k in tk)
            df2 = ("frf_$(k)!$c"=>fr["factorresponsefeature"][k] for k in keys(fr["factorresponsefeature"]))
            df3 = ("pzfr_$(k)!$c"=>map((i,p,m)->(m[i...].-mean(p[i...]))/std(p[i...]),fr["optfri"][k],fr["pfms"],fr["fms"]) for k in keys(fr["fa"]))
            df4 = ("f_$(k)!$c"=>fill(fr["fa"][k],length(id)) for k in keys(fr["fa"]))
            df = DataFrame(df1...,df2...,df3...,df4...)
            df.id = id
            sdf[siteid] = outerjoin(sdf[siteid],df,on=:id)
        end
    end
    return outerjoin(cell,reduce((i,j)->append!(i,j,cols=:union),values(sdf)),on=:id)
end

