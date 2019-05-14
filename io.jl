using MAT, NeuroAnalysis.NAIO, DataFrames

"read `Matlab` exported data"
function readmat(f)
  d = matread(f)["dataset"]
  return d["data"],d["param"]
end

"convert `Matlab` struct of array to `DataFrame`"
matdictarray2df(d) = DataFrame(Any[squeeze(v,2) for v in values(d)],[Symbol(k) for k in keys(d)])

"Regular expression to match `VLab` data file names"
function vlabregex(testtype;subject="[0-9]",maxch=1,cell="[a-z]",maxrepeat=3,format="mat")
  mr = lpad(maxrepeat,2,0)
  Regex("^$subject+[A-Za-z]+[1-$maxch]$cell$testtype[0-$(mr[1])][0-$(mr[2])]\\.$format")
end

"Get matched file names in path"
function getdatafile(testtype;subject="[0-9]",path="./data")
  matchfile(vlabregex(testtype,subject=subject),path=path)
end
