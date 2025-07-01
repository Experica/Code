### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b6e9ff40-2e73-11ef-3453-d55ec196af6b
begin
	import Pkg
	Pkg.activate(Base.current_project())
	using PlutoUI,NeuroAnalysis,FileIO,JLD2,Images,DataFrames,XLSX,StatsBase,AlgebraOfGraphics,Plots
	import CairoMakie as mk
end

# ╔═╡ c006ddce-951b-4395-9b0e-706c1ef2c753
begin
	resultroot = "Z:/"
	alignment = DataFrame(XLSX.readtable(joinpath(resultroot,"ISIMapAlignment.xlsx"),"Sheet1"))
	rename!(alignment,:Subject_ID=>:subject,:siteid=>:isisiteid)
	select!(alignment,Not(:ppmm),:ppmm=>(i->i/1000)=>:ppum)
	transform!(alignment,Not([:subject,:isisiteid,:ppum,:bvname]).=>ByRow(x->ismissing(x) ? [0.,0] : parse.(Float64,split(strip(x,['(',')']),", "))).=>identity)
	
	penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
	transform!(penetration,:coord=>ByRow(c->ismissing(c) ? c : parse.(Float64,split(c,", ")))=>identity,:pid=>ByRow(i->split(i)[1])=>identity)
	rename!(penetration,:Subject_ID=>:subject)

	pdf = dropmissing(leftjoin(penetration[:,[:subject,:siteid,:pid,:coord]],alignment,on=[:subject]),[:subject,:coord])
	pmapdir = joinpath(resultroot,"PenetrationMap");mkpath(pmapdir)
end;

# ╔═╡ 7729558e-2c19-4891-9527-4367146b4e0e
begin
	saveresult=@bind saving CheckBox(default=false)
	md"#### Save All Results? $saveresult"
end

# ╔═╡ fdc2399c-f601-49d0-af34-18a34e54ccb1
md"# Orientation LHI and LR"

# ╔═╡ 84e0ab89-59db-4ca1-8a2f-5baa4b0f6528
orimaps = Dict(r.subject=>load(joinpath(resultroot,r.subject,r.isisiteid,"ori_map.jld2")) for r in eachrow(alignment))

# ╔═╡ c6c5f833-fe69-4d23-b4b0-78af456e780e
function plotlocalmap(li,lm;pid=pdf.pid,liname="LHI",dir=nothing,figfmt=[".png"])
	f=mk.Figure()
	a1=mk.Axis(f[1,1])
	is = map(i->ismissing(i) ? NaN : first(i),li)
	ms = map(i->ismissing(i) ? [] : i,lm)
	mk.scatter!(a1,is,marker=ms,markersize=30)
	a1.xticks=(1:4:length(pid),pid[1:4:end])
	a1.xlabel="Penetration";a1.ylabel=liname
	a2=mk.Axis(f[1,2])
	mk.hist!(a2,filter(!isnan,is),direction=:x)
	mk.hideydecorations!(a2,grid=false)
	mk.colsize!(f.layout,1,mk.Relative(6/7))
	isnothing(dir) || foreach(ext->save(joinpath(dir,"$liname$ext"),f),figfmt)
	f
end

# ╔═╡ f24ba5ea-a486-4065-b79e-cd6c1aba7233
begin
	sigma_lhi = @bind hs Slider(5:5:300,show_value=true,default=150)
	md"#### σ for LHI: $sigma_lhi μm"
end

# ╔═╡ 1b62a606-28a1-4c07-aa18-6382543c2b99
begin
	sigma_lr = @bind rs Slider(5:5:300,show_value=true,default=40)
	md"#### σ for LR: $sigma_lr μm"
end

# ╔═╡ 4a009372-619f-4986-b488-5df05e6aace6
begin
	po = [orimaps[r.subject]["amap"][round.(Int,r.coord.-r.ori)...] for r in eachrow(pdf)] # local orientation
	
	lhis = [localhomoindex(orimaps[r.subject]["amap"],r.coord.-r.ori;n=2,σ=hs*r.ppum) for r in eachrow(pdf)] # local homogeneity index

	lrs = [localaverage(orimaps[r.subject]["mmap"],r.coord.-r.ori;σ=rs*r.ppum) for r in eachrow(pdf)] # local r(1-circ_var)
	
	amap_ori = [alphamask(map(a->HSV(rad2deg(2a),1,1),orimaps[pdf.subject[i]]["amap"][lhis[i].roi...]))[1] for i in eachindex(lhis)]
	foreach(i->i[map(c->c:c+1,floor.(Int,size(i)./2))...].=HSVA(0,0,0.1,1),amap_ori)

	pmap_ori = [alphamask(map((a,m)->HSV(rad2deg(2a),1,m),orimaps[pdf.subject[i]]["amap"][lrs[i].roi...],orimaps[pdf.subject[i]]["mmap"][lrs[i].roi...]))[1] for i in eachindex(lrs)]
	foreach(i->i[map(c->c:c+1,floor.(Int,size(i)./2))...].=HSVA(0,0,0.9,1),pmap_ori)
end

# ╔═╡ 5b8b5618-cec1-47dd-bcff-f3c8ce845522
lhis[1]

# ╔═╡ d9395ae1-23f0-4805-a3ab-8002211954a5
let
	plt = data(DataFrame(LHI=first.(lhis),LR=first.(lrs)))*mapping(:LR,:LHI)*(linear()+visual(mk.Scatter,color=:gray40))
	f=draw(plt,figure=(size=(600,500),))
end

# ╔═╡ a4235ccb-8b31-4102-8b13-ce1c10177f10
begin
	srange = 5:5:400
	slhis = stack(s->[localhomoindex(orimaps[r.subject]["amap"],r.coord.-r.ori;n=2,σ=s*r.ppum)[1] for r in eachrow(pdf)], srange) # local homogeneity index
	slrs = stack(s->[localaverage(orimaps[r.subject]["mmap"],r.coord.-r.ori;σ=s*r.ppum)[1] for r in eachrow(pdf)], srange) # local r(1-circ_var)
	cors = [cor(slrs[:,j],slhis[:,i]) for i in eachindex(srange),j in eachindex(srange)]
	pmdir = saving ? pmapdir : nothing
	if saving
		foreach((id,m)->save(joinpath(pmdir,"$(id)_ori_lhi.png"),m) ,pdf.siteid,amap_ori)
		foreach((id,m)->save(joinpath(pmdir,"$(id)_ori_lr.png"),m) ,pdf.siteid,pmap_ori)
		jldsave(joinpath(resultroot,"ori_lhi_lr.jld2");po,srange,slhis,slrs,siteid=pdf.siteid)
	end
end

# ╔═╡ 06f1add9-32d8-4c2b-815c-835e0f254318
plotlocalmap(lhis,amap_ori,dir=pmdir)

# ╔═╡ 5c9333e3-6543-4c16-b7a7-a8d5f79342f8
plotlocalmap(lrs,pmap_ori,liname="LR",dir=pmdir)

# ╔═╡ e692b13f-4bea-4d05-a043-dfdbc10eece1
let
	f=mk.Figure()
	mk.Axis(f[1,1],xlabel="σ_LR (μm)",ylabel="σ_LHI (μm)",title="Correlation between LR and LHI",xticks=srange[1:10:end],yticks=srange[1:10:end])
	p = mk.contourf!(srange,srange,cors,levels=80,colormap=:turbo)
	mk.Colorbar(f[1,2],p)
	f
end

# ╔═╡ bc3f7c8a-22a1-41cb-a0d9-8e27ddf2d060
md"# OD Intensity"

# ╔═╡ 53d7e60b-b541-46dd-8b02-0d4e649b4efb
odmaps=Dict(r.subject=>load(joinpath(resultroot,r.subject,r.isisiteid,"od_map.jld2")) for r in eachrow(alignment))

# ╔═╡ 48889a7f-12fe-4084-99da-37eb44a2f0b6
begin
	sigma_od = @bind s_od Slider(5:5:300,show_value=true,default=50)
	md"#### σ for ODI: $sigma_od μm"
end

# ╔═╡ d8b4c153-dfc2-41a7-82d7-0b32eadd13b1
begin
	odis = [localaverage(odmaps[r.subject]["odmagmap"],r.coord.-r.od;σ=s_od*r.ppum) for r in eachrow(pdf)] # local od intensity
	
	map_od = [alphamask(Gray.(odmaps[pdf.subject[i]]["odmagmap"][odis[i].roi...]))[1] for i in eachindex(odis)]
	foreach(i->i[map(c->c:c+1,floor.(Int,size(i)./2))...].=GrayA(0.5,1),map_od)
end

# ╔═╡ 43b9c22b-165f-4407-80d8-c9bef1b4bdbc
plotlocalmap(odis,map_od;liname="OD Intensity",dir=pmdir)

# ╔═╡ fe4d198b-0ddc-4fd4-acb2-73d369592e2a
begin
	if saving
		foreach((id,m)->save(joinpath(pmdir,"$(id)_od.png"),m) ,pdf.siteid,map_od)
		
		sodis = stack(s->[localaverage(odmaps[r.subject]["odmagmap"],r.coord.-r.od;σ=s*r.ppum)[1] for r in eachrow(pdf)], srange)
		
		jldsave(joinpath(resultroot,"odi.jld2");srange,sodis,siteid=pdf.siteid)
	end
end

# ╔═╡ 2429d158-3dcf-45bc-b25a-9c2e4dddebab
md"# COFD Intensity"

# ╔═╡ 3a0cdbb0-d14a-48fc-b3e9-ee896db2cc37
begin
	cs_cofd = ["DKL_X", "LMS_X", "LMS_Y", "LMS_Z", "DKL_RG", "DKL_BY"]
	function loadcofd(dir;cs=cs_cofd)
		files = [c=>matchfile(Regex("cofd_map_$(c)\\w*.jld2");dir,join=true) for c in cs]
		ds = Dict(k=>isempty(v) ? nothing : load(v[1]) for (k,v) in files)
	end
end

# ╔═╡ 7eff433a-7a29-4ee5-90c6-65d2d5ef26bc
cofdmaps=Dict(r.subject=>loadcofd(joinpath(resultroot,r.subject,r.isisiteid)) for r in eachrow(alignment))

# ╔═╡ 689127b7-fcaa-4cc8-9f17-c948494d08eb
function cofdintensity(pdf,cofdmaps;cs=cs_cofd,s=50)
	t = select(pdf,[:subject,:siteid])
	cofdi = insertcols(t,(c=>
	[begin
		d = cofdmaps[r.subject][c]
		isnothing(d) ? missing : localaverage(d["cofdmagmap"],r.coord.-r[c];σ=s*r.ppum)
	end for r in eachrow(pdf)]
	for c in cs)...)

	color_cofd = insertcols(t,(c=>
	[begin
		d = cofdmaps[r.subject][c]
		isnothing(d) ? missing : d["minmaxcolor"]
	end for r in eachrow(pdf)]
	for c in cs)...)
	
	map_cofd = insertcols(t,(c=>
	[begin
		d = cofdmaps[r.subject][c]
		if isnothing(d)
			m=missing
		else
			m = alphamask(get(cgrad(d["minmaxcolor"]...),d["cofdmagmap"][r[c].roi...]))[1]
			m[map(c->c:c+1,floor.(Int,size(m)./2))...].=RGBA(0.5,0.5,0.5,1)
		end
		m
	end for r in eachrow(cofdi)]
	for c in cs)...)
	
	(;cofdi,color_cofd,map_cofd)
end

# ╔═╡ a8ed5c93-3f3f-40e4-8070-3098d2de7825
begin
	sigma_cofd = @bind s_cofd Slider(5:5:300,show_value=true,default=35)
	md"#### σ for COFDI: $sigma_cofd μm"
end

# ╔═╡ 2a257741-d8cd-437b-9e3a-9c761ef8af3a
cofdi,color_cofd,map_cofd = cofdintensity(pdf,cofdmaps,s=s_cofd);

# ╔═╡ cc17b590-6259-4610-b59e-dee40d7e04fd
md"## Achromatic"

# ╔═╡ b86938a5-c82e-4d6d-a38d-9c52c82dc4ea
plotlocalmap(cofdi.DKL_X,map_cofd.DKL_X,liname="COFD Intensity_A",dir=pmdir)

# ╔═╡ 208f63db-2d33-48e4-8651-c93fd56c1bd0
md"## L Cone"

# ╔═╡ 3ec51bfa-9f5f-40d2-b561-2b47fa4b98f6
plotlocalmap(cofdi.LMS_X,map_cofd.LMS_X,liname="COFD Intensity_L",dir=pmdir)

# ╔═╡ 5fd52e2c-2bef-45e2-b7fa-72f5d06d7428
md"## M Cone"

# ╔═╡ 92513f80-07b9-4218-837a-b4b70fde5665
plotlocalmap(cofdi.LMS_Y,map_cofd.LMS_Y,liname="COFD Intensity_M",dir=pmdir)

# ╔═╡ a4cab151-3c5e-4bfe-87c1-af3c415fde97
md"## S Cone"

# ╔═╡ c702f640-32d5-4641-bc7c-2f118e744739
plotlocalmap(cofdi.LMS_Z,map_cofd.LMS_Z,liname="COFD Intensity_S",dir=pmdir)

# ╔═╡ 5c1e0fcf-eb4b-43e1-bc9c-3e9210300ab5
md"## DKL Red/Green"

# ╔═╡ edc08795-e2bf-4710-9b61-4d82a58e8d83
plotlocalmap(cofdi.DKL_RG,map_cofd.DKL_RG,liname="COFD Intensity_RG",dir=pmdir)

# ╔═╡ a18461e2-a787-486c-9f20-0754ceb800c8
md"## DKL Blue/Yellow"

# ╔═╡ 4c6467b6-eecf-4353-a5ce-4a6a608aaae6
plotlocalmap(cofdi.DKL_BY,map_cofd.DKL_BY,liname="COFD Intensity_BY",dir=pmdir)

# ╔═╡ 09e38d52-b6a8-4dd6-9333-0c1683f81ffb
begin
	if saving
		foreach(c->
		foreach((id,m)->ismissing(m) || save(joinpath(pmdir,"$(id)_cofd_$(c).png"),m) ,pdf.siteid,map_cofd[!,c]),
		cs_cofd)
		
		scofdis = map(s->cofdintensity(pdf,cofdmaps;s).cofdi, srange)
		
		jldsave(joinpath(resultroot,"cofdi.jld2");srange,scofdis,color_cofd,siteid=pdf.siteid)
	end
end

# ╔═╡ Cell order:
# ╠═b6e9ff40-2e73-11ef-3453-d55ec196af6b
# ╟─c006ddce-951b-4395-9b0e-706c1ef2c753
# ╟─7729558e-2c19-4891-9527-4367146b4e0e
# ╟─fdc2399c-f601-49d0-af34-18a34e54ccb1
# ╠═84e0ab89-59db-4ca1-8a2f-5baa4b0f6528
# ╠═4a009372-619f-4986-b488-5df05e6aace6
# ╠═5b8b5618-cec1-47dd-bcff-f3c8ce845522
# ╟─c6c5f833-fe69-4d23-b4b0-78af456e780e
# ╟─f24ba5ea-a486-4065-b79e-cd6c1aba7233
# ╟─06f1add9-32d8-4c2b-815c-835e0f254318
# ╟─1b62a606-28a1-4c07-aa18-6382543c2b99
# ╟─5c9333e3-6543-4c16-b7a7-a8d5f79342f8
# ╟─d9395ae1-23f0-4805-a3ab-8002211954a5
# ╠═a4235ccb-8b31-4102-8b13-ce1c10177f10
# ╟─e692b13f-4bea-4d05-a043-dfdbc10eece1
# ╟─bc3f7c8a-22a1-41cb-a0d9-8e27ddf2d060
# ╠═53d7e60b-b541-46dd-8b02-0d4e649b4efb
# ╠═d8b4c153-dfc2-41a7-82d7-0b32eadd13b1
# ╟─48889a7f-12fe-4084-99da-37eb44a2f0b6
# ╟─43b9c22b-165f-4407-80d8-c9bef1b4bdbc
# ╟─fe4d198b-0ddc-4fd4-acb2-73d369592e2a
# ╟─2429d158-3dcf-45bc-b25a-9c2e4dddebab
# ╠═3a0cdbb0-d14a-48fc-b3e9-ee896db2cc37
# ╠═7eff433a-7a29-4ee5-90c6-65d2d5ef26bc
# ╟─689127b7-fcaa-4cc8-9f17-c948494d08eb
# ╠═2a257741-d8cd-437b-9e3a-9c761ef8af3a
# ╟─a8ed5c93-3f3f-40e4-8070-3098d2de7825
# ╟─cc17b590-6259-4610-b59e-dee40d7e04fd
# ╠═b86938a5-c82e-4d6d-a38d-9c52c82dc4ea
# ╟─208f63db-2d33-48e4-8651-c93fd56c1bd0
# ╠═3ec51bfa-9f5f-40d2-b561-2b47fa4b98f6
# ╟─5fd52e2c-2bef-45e2-b7fa-72f5d06d7428
# ╠═92513f80-07b9-4218-837a-b4b70fde5665
# ╟─a4cab151-3c5e-4bfe-87c1-af3c415fde97
# ╠═c702f640-32d5-4641-bc7c-2f118e744739
# ╟─5c1e0fcf-eb4b-43e1-bc9c-3e9210300ab5
# ╠═edc08795-e2bf-4710-9b61-4d82a58e8d83
# ╟─a18461e2-a787-486c-9f20-0754ceb800c8
# ╠═4c6467b6-eecf-4353-a5ce-4a6a608aaae6
# ╠═09e38d52-b6a8-4dd6-9333-0c1683f81ffb
