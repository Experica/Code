using FileIO
d=load("ucondidx.jld2")

cond = d["cond"]

cidx = d["cidx"]
uidx2= uidx[:,1]



cc=hcat(map(i->unique(i),cidx)...)'

unique(map(i->cond[i,:SpatialFreq],cc))
cond[12,:]
cond[1482,:]ys





ys,ns,ws,is=histrv(float.(p),0,1,nbins=20)

plot(ns)

findmax(ns)
mean(ys[3])
mean(ys[findmax(ns)[2]])



function hotellingt2test(X,mu,α = 0.05)
	n=size(X,1); p=size(X,2); m = mean(X,dims=1); S = cov(X);
	T2 = n*(m-mu)*inv(S)*(m-mu)'  # T-square statistics

	if n>= 50 # Chi-square approximation
		P = 1-chisqcdf(T2[1],p)
	else  # F approximation
		F = (n-p)/((n-1)*p)*T2[1]
		P = ccdf(FDist(p, n-p), F)
	end
	H = P<α
	return P
end

"""
OTFIT_CARANDINI Fits orientation curves like Carandini/Ferster 2000

[Rsp,Rp,Op,sigm,Rn,FITCURVE,ERR]=OTFIT_CARANDINI(ANGLES,...
     SPONTHINT, MAXRESPHINT, OTPREFHINT, WIDTHHINT,'DATA',DATA)

Finds the best fit to the function

R=Rsp+Rp*EXP(-(X-Op)^2)/(2*sig^2))+Rn*EXP(-(X-Op+180)^2/(2*sig^2))
where R is the response, Rsp is the spontaneous response, Rp is
the response at the preferred angle, Op is the preferred angle,
sigm is the tuning width, Rn is the firing rate at 180+Op.

FITCURVE is the fit function at 1 degree intervals (0:1:359).

Julia version of fitting function adapted from the Matlab code wrote by Van Hooser lab. Peichao
"""
# function otfit_carandini(angles,sponthint,maxresphint,otprefhint,widthhint,varargin)

# [Rsp,Rp,Ot,sigm,Rn,fitcurve,er,R2]=

spontfixed = NaN
assign(varargin{:})

Po = [sponthint maxresphint otprefhint widthhint maxresphint]

if ~isnan(spontfixed), Po = Po(2:end) end

Rsp_,Rp_,Ot_,sigm_,Rn_=otfit_carandini_conv("TOFITTING",Po,varargin{:});

Po = [ Rsp_ Rp_ Ot_ sigm_ Rn_]
if ~isnan(spontfixed), Po = Po(2:end) end

# starting point
options= optimset("Display","off","MaxFunEvals",10000,"TolX",1e-6);

assign(varargin{:});
searchArg = "@(x) otfit_carandini_err(x,angles,"
for i=1:2:length(varargin)
	searchArg = [searchArg "" varargin{i} "" "," varargin{i} ","]
end
needconvert = 1
searchArg = [searchArg '''needconvert'',needconvert)']

Pf = eval(['fminsearch(' searchArg ',Po,options);'])

[Rsp,Rp,Ot,sigm,Rn] = otfit_carandini_conv('TOREAL',Pf,varargin{:});

if ~isnan(spontfixed), er = otfit_carandini_err([Rp Ot sigm Rn],angles,varargin{:})
elseif er = otfit_carandini_err([Rsp Rp Ot sigm Rn],angles,varargin{:})
end

fitcurve = []
if nargout>5
	for i=1:length(varargin),
		if strcmp(varargin{i},'data'), break end
	end
	varargin = {varargin{setdiff(1:length(varargin),[i i+1])}}
	[d,fitcurve]=otfit_carandini_err(Pf,0:359,varargin{:},'needconvert',1)
end

R2 = 1 - sum(er)/(sum((data-mean(data)).^2))
# end
