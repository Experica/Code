using Statistics,FileIO,JLD2,Images,Plots

# Manual Layer Assignment
dataroot = "X:/"
dataexportroot = "Y:/"
resultroot = "Z:/"

subject = "AG2";recordsession = "V1";recordsite = "ODR1"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
figfmt = [".png"]


ii = '0'
test = "Flash2Color"
testids = ["$(siteid)_$(test)_$i" for i in 0:3]
absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
## AP dRMS
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = repeat(["$(i["eye"])_$(i["color"])" for i in aps],inner=2)
colors = mapreduce(i->[RGBA(parse.(Float64,split(match(r".*Color=\[(.*)\]",k).captures[1],", "))...) for k in keys(i["cdrms"])],append!,aps)
drms = mapreduce(i->[1e6v for v in values(i["cdrms"])],append!,aps)
rmstime = mapreduce(i->[i["time"],i["time"]],append!,aps)
rmsdepths = mapreduce(i->[i["depths"],i["depths"]],append!,aps)
rmslim=absmax(drms)

## Unit dPSTH
dpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
dpsth = mapreduce(i->[v.psth for v in values(i["cdpsth"])],append!,dpsths)
psthtime = mapreduce(i->[v.x for v in values(i["cdpsth"])],append!,dpsths)
psthdepths = mapreduce(i->[v.y for v in values(i["cdpsth"])],append!,dpsths)
psthlim=absmax(dpsth)

## LFP and CSD
lfps = load.(joinpath.(siteresultdir,testids,"lfp$ii.jld2"))
lfp = mapreduce(i->[1e6v for v in values(i["cmlfp"])],append!,lfps)
dcsd = mapreduce(i->[imfilter(v,Kernel.gaussian((1,0)),Fill(0)) for v in values(i["cdcsd"])],append!,lfps)
lfptime = mapreduce(i->[i["time"],i["time"]],append!,lfps)
lfpdepths = mapreduce(i->[i["depths"],i["depths"]],append!,lfps)
lfplim=absmax(lfp)
csdlim=absmax(dcsd)

## Unit Feature
uf = load(joinpath(siteresultdir,"unitdepthfeature.jld2"),"df")
ufy = uf["depth"]
udf = (;(Symbol(k)=>uf["depthfeature"][:,k.==uf["feature"][:]] for k in uf["feature"])...)

plotflashdepthinfo=(;w=150,h=500,r=4,layer=nothing,clim=:absmax)->begin
    n = length(drms)
	cs = permutedims(palette(:tab10).colors.colors[1:4])
    p=plot(layout=(r,n),link=:y,legend=false,grid=false,size=(n*w,r*h))
    for i in 1:n, j in 1:r
		if j == 1
			resp = drms;time = rmstime;depths = rmsdepths;lim = rmslim;color=:coolwarm;title=titles
		elseif j == 2
			resp = dpsth;time = psthtime;depths = psthdepths;lim = psthlim;color=:coolwarm;title=fill("",n)
		elseif j == 3
			resp = dcsd;time = lfptime;depths = lfpdepths;lim = csdlim;color=cgrad(:coolwarm,rev=true);title=fill("",n)
		end
		if j==4
			lmargin = i==1 ? 12Plots.mm : -4Plots.mm
			if i == 1
				ux = 1e9*udf.Density
				uxmin = minimum(ux)
				plot!(p[j,i],ux,ufy;xlabel="Density",ylabel="Depth (μm)",leg=false,color=cs[i],
				    lw=1.5,yticks=0:250:ufy[end],xticks=[ceil(uxmin,sigdigits=1),floor(maximum(ux),sigdigits=1)],
					left_margin=lmargin,bottom_margin=3Plots.mm,tickor=:out)
			elseif i == 2
				ux = [udf.upspread;reverse(udf.downspread)]
				uy = [ufy;reverse(ufy)]
				plot!(p[j,i],ux,uy,st=:shape,xlabel="Spread",lw=0,yticks=false,alpha=0.8,xticks=[0],color=cs[i],
				left_margin=lmargin,bottom_margin=3Plots.mm,tickor=:out)
			elseif i == 3
				ux = [0;udf.peaktroughratio;0]
				uy = [ufy[1];ufy;ufy[end]]
				plot!(p[j,i],ux,uy,st=:shape,xlabel="PTRatio",lw=0,yticks=false,alpha=0.8,xticks=[0],color=cs[i],
				left_margin=lmargin,bottom_margin=3Plots.mm,tickor=:out)
			elseif i == 4
				ux = [0;udf.duration;0]
				uy = [ufy[1];ufy;ufy[end]]
				plot!(p[j,i],ux,uy,st=:shape,xlabel="Duration",lw=0,yticks=false,alpha=0.8,xticks=[0],color=cs[i],
				left_margin=lmargin,bottom_margin=3Plots.mm,tickor=:out)
			else
				plot!(p[j,i],frame=:none,xticks=false,yticks=false,left_margin=lmargin)
			end
			if i <= 4
				if !isnothing(layer)
					anno = i==1 ? [(uxmin,mean(layer[k]),text(k,6,:gray10,:center,:left)) for k in keys(layer)] : []
		            hline!(p[j,i],[l[2] for l in values(layer)],linecolor=:gray25,legend=false,annotations=anno)
		        end
			end
		else
	        yticks = i==1 ? (0:250:depths[i][end]) : false
	        xticks = 0:50:time[i][end]
	        lmargin = i==1 ? 12Plots.mm : -4Plots.mm
	        xlabel = (j==3 && i==1) ? "Time (ms)" : ""
	        ylabel = i==1 ? "Depth (μm)" : ""
			canno = [(5,25,Plots.text("■",15,colors[i],:left,:bottom))]
			if clim == :link
				clims = (-lim,lim)
			elseif clim == :absmax
				lim = maximum(abs.(resp[i]))
				clims = (-lim,lim)
			else
				clims = :auto
			end
	        heatmap!(p[j,i],time[i],depths[i],resp[i];color,clims,
	        title=title[i],titlefontsize=10,yticks,xticks,tickor=:out,xlabel,ylabel,
	        annotation=canno,left_margin=lmargin)
			if !isnothing(layer)
				# anno = i==1 ? [(time[i][1]+3,mean(layer[k]),text(k,6,:gray10,:center,:left)) for k in keys(layer)] : []
				anno = [(time[i][1]+3,mean(layer[k]),text(k,6,:gray10,:center,:left)) for k in keys(layer)]
	            hline!(p[j,i],[l[2] for l in values(layer)],linecolor=:gray25,legend=false,lw=0.5,annotations=anno)
	        end
		end
    end
    p
end

## Assign Layer
layer = load(joinpath(siteresultdir,"layer.jld2"),"layer")
layer = Dict()

layer["1"] = [1800,1900]
layer["2"] = [1700,1800]
layer["3"] = [1350,1700]
layer["4A"] = [1300,2550]
layer["4B"] = [1150,1300]
layer["4Cα"] = [2150,2280]
layer["4Cβ"] = [950,2150]
layer["5"] = [750,2030]
layer["6"] = [500,750]
layer["WM"] = [0,500]

plotflashdepthinfo(;layer,clim=:absmax)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_LayerAssignment$ext")),figfmt)

## Finalize Layer
layer = checklayer!(layer)
jldsave(joinpath(siteresultdir,"layer.jld2");layer,siteid)
