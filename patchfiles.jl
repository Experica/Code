dir = "."
ext = ".yaml"
for (root, dirs, files) in walkdir(dir)
    fs = filter(f->endswith(f,ext),files)
    isempty(fs) && continue
    for f in fs
        fp = joinpath(root,f)
        display("Patching $fp ...")

        fc = read(fp,String)
        i = findfirst("      ID: ROGPG279Q",fc)
        ii = findnext("      Latency: 80",fc,i[end])
        pfc=fc[1:ii[end]-2] * "33" * fc[ii[end]+1:end]
        write(fp,pfc)
    end
end



fc = read("./AF5_HLV1_ODL1/AF5_HLV1_ODL1_OriSF_5.yaml",String)

i = findfirst("TimerDriftSpeed: ",fc)
ii = findnext("\r\n",fc,i[end])
mfc=fc[1:i[end]] * "4.89E-05\r\n" * fc[ii[end]+1:end]

write("test.yaml",mfc)
