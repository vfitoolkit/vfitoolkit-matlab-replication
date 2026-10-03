% HV2000 Model 4: Model 3 + iid temporary shock on log labor endowment.
% Uses experienceassetze (ebar' depends on z and e) and the M4 fn variants.

fprintf('\n========== HV2000 Model 4 (AR(1) + iid temporary shock) ==========\n');

%% Model-specific stochastic structure (Table 3)
Params.sigma2_y1=0.24;
Params.gamma=0.985;
Params.sigma2_eps1=0.02;
% Guarded: HV2000_Table6 varies this across {0.01,0.04,0.09} and restores 0.01
% afterwards. An unguarded assignment here silently clobbered its loop value.
if ~isfield(Params,'sigma2_eps2')
    Params.sigma2_eps2=0.01;
end
n_z=61;
n_e=3;

% Persistent AR(1) on log z. HV2000 footnote 14 says endpoints at +/- 6*sigma;
% interpreting as 6 STATIONARY stds (discretizeAR1_Tauchen scales by stationary std internally).
[zlog_grid,pi_z]=discretizeAR1_Tauchen(0,Params.gamma,sqrt(Params.sigma2_eps1),n_z,6);
% Note: Tanaka-Toda / Farmer-Toda would give a better discretization than Tauchen,
% but sticking with Tauchen to more closely follow the original paper.
z_grid=exp(zlog_grid);

% Iid log-e grid: 3 points on [-2*sigma_eps2, +2*sigma_eps2] (footnote 14).
[elog_grid,pi_e]=discretizeIID_Tauchen(0,sqrt(Params.sigma2_eps2),n_e,2);
% Note: Tanaka-Toda / Farmer-Toda would give a better discretization than Tauchen,
% but sticking with Tauchen to more closely follow the original paper.
e_grid=exp(elog_grid);

% Birth distribution: log(z_1) ~ N(0, sigma2_y1) onto the AR(1) zlog_grid
% (hand-discretized via midpoint normal CDF; Tauchen helpers build their own
% grids so we can't reuse them here).
mid=(zlog_grid(1:end-1)+zlog_grid(2:end))/2;
cdf_mid=normcdf(mid,0,sqrt(Params.sigma2_y1));
probs_z=zeros(n_z,1);
probs_z(1)=cdf_mid(1);
probs_z(end)=1-cdf_mid(end);
probs_z(2:end-1)=cdf_mid(2:end)-cdf_mid(1:end-1);

% Under the AR(1), Jensen's inequality on a growing log-variance makes the
% cohort mean of z drift with age. Build an age-varying z_grid_J via
% MarkovChainMoments_FHorz so that E[z_j]=1 at every j.
z_grid_J=repmat(z_grid,1,N_j);
pi_z_J=repmat(pi_z,1,1,N_j);
[zmean_J,~,~,~]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,probs_z);
for jj=1:N_j
    z_grid_J(:,jj)=z_grid_J(:,jj)/zmean_J(jj);
end

% Renormalize the iid e_grid so E[e]=1 under pi_e.
e_grid=e_grid/sum(pi_e.*e_grid);

jequaloneDist=zeros([n_a,n_z,n_e]);
% Loop variable is 'ee' not 'kk': this script runs in the caller's workspace, and
% HV2000_Table6 loops over variances with 'kk' (a clobber there silently wrote
% every variance into the same column of its results table).
for ee=1:n_e
    jequaloneDist(azeroindex,1,:,ee)=reshape(probs_z*pi_e(ee),[1,1,n_z,1]);
end

%% Build Model 4's vfoptions/simoptions from the shared bases.
% Model 4 swaps experienceassetz -> experienceassetze and adds the iid e channel.
vfoptions4=vfoptions;
vfoptions4.experienceassetz=0;
vfoptions4.experienceassetze=1;
vfoptions4.n_e=n_e;
vfoptions4.e_grid=e_grid;
vfoptions4.pi_e=pi_e;
vfoptions4.aprimeFn=@(d,ebar,z,e,agej,Jr,w,ybar_j,ybar_mean,g) HV2000_ebarprimeFn_M4(d,ebar,z,e,agej,Jr,w,ybar_j,ybar_mean,g);

simoptions4=simoptions;
simoptions4.experienceassetz=0;
simoptions4.experienceassetze=1;
simoptions4.n_e=n_e;
simoptions4.e_grid=e_grid;
simoptions4.pi_e=pi_e;
simoptions4.aprimeFn=vfoptions4.aprimeFn;
simoptions4.z_grid=z_grid_J;

%% Override ReturnFn and FnsToEvaluate to M4 variants (signatures include e)
ReturnFn_M4=@(d,aprime,a,ebar,z,e,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS)...
    HV2000_ReturnFn_M4(d,aprime,a,ebar,z,e,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS);

FnsToEvaluate_M4.K=@(d,aprime,a,ebar,z,e) a;
FnsToEvaluate_M4.L=@(d,aprime,a,ebar,z,e,agej,Jr,ybar_j) (agej<Jr)*z*e*ybar_j;
FnsToEvaluate_M4.Earnings=@(d,aprime,a,ebar,z,e,agej,Jr,ybar_j,w) (agej<Jr)*z*e*ybar_j*w;
FnsToEvaluate_M4.AccBeq=@(d,aprime,a,ebar,z,e,sj,r,tau) (1-sj)*aprime*(1+r*(1-tau));
FnsToEvaluate_M4.SSBenefits=@(d,aprime,a,ebar,z,e,agej,Jr,g,b_common,ybar_mean,w,haveSS)...
    HV2000_SSBenefitFn_M4(d,aprime,a,ebar,z,e,agej,Jr,g,b_common,ybar_mean,w,haveSS);
FnsToEvaluate_M4.Savings=@(d,aprime,a,ebar,z,e,g) (1+g)*aprime-a;
FnsToEvaluate_M4.Ebar=@(d,aprime,a,ebar,z,e) ebar;
FnsToEvaluate_M4.Consumption=@(d,aprime,a,ebar,z,e,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS)...
    HV2000_ConsumptionFn_M4(d,aprime,a,ebar,z,e,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS);
FnsToEvaluate_M4.Income=@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)...
    HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS);
FnsToEvaluate_M4.SavingRate=@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)...
    HV2000_SavingRateFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS);

%% This model is to big for my desktop, so turn on lowmemory=1
% vfoptions4.lowmemory=1;

%% Solve general equilibrium
% With per-age renormalization of z_grid_J (and iid renorm of e_grid),
% E[z_j*e]=1 at every j, so the precomputed Params.ybar_mean exactly equals
% mean earnings per worker divided by w. No ybar_mean GE price needed.
[p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d,n_a,n_z,N_j,[],pi_z_J,d_grid,a_grid,z_grid_J,ReturnFn_M4,FnsToEvaluate_M4,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions4,vfoptions4);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

%% Re-solve at GE prices and compute AllStats
tic;
vfoptions4.verbose=1;
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn_M4,Params,DiscountFactorParamNames,[],vfoptions4);
vftime=toc
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions4);
toc
tic;
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_M4,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions4);
toc

%% CheckTopOfGrid: mass near the top of the asset / ebar grids
% See Results.Model4_abar{0,w}.CheckTopOfGrid; large values => widen Params.amax / ebar_gridmax.
mass_by_a   = sum(StationaryDist,  2:ndims(StationaryDist));    % marginal over asset dim (dim 1)
mass_by_ebar= sum(StationaryDist, [1, 3:ndims(StationaryDist)]);% marginal over ebar dim (dim 2)
CheckTopOfGrid.MassAtTop10gridptsAssets = sum(mass_by_a(end-9:end));
CheckTopOfGrid.MassAtTop5gridptsAssets  = sum(mass_by_a(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsAssets  = sum(mass_by_a(end-1:end));
CheckTopOfGrid.MassAtTop5gridptsEbar    = sum(mass_by_ebar(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsEbar    = sum(mass_by_ebar(end-1:end));

% I ran out of memory during solving this model, so
clear V

%% Table 5: saving rates by income multiples (bin indicator folded into AggVars)
% Per HV2000 p. 380 (just below Table 5):
%   "multiples are calculated by taking a 10% band around each income multiple
%    and then dividing total saving of agents in the band by total income of
%    agents in the band. Income is defined as earnings after social security
%    taxes plus interest income and transfers."
% Bin = [0.9*m*Y, 1.1*m*Y] with Y=Params.IncomeMean (full-population mean);
% saving rate computed below as Savings.Mean / Income.Mean within the bin
% (toolkit's .Mean under conditionalrestriction is mass-weighted, so this
% ratio equals total saving in bin / total income in bin).
Params.IncomeMean=AllStats.Income.Mean;
HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];
% The bin saving rate is total saving in the band over total income in the band,
% which is a ratio of two AGGREGATES: the bin mass cancels, so no conditional
% restriction is needed. Folding the bin indicator into the functions and taking
% one AggVars pass avoids AllStats entirely -- AllStats sorts the whole grid to
% get medians/Lorenz/quantiles we never read, and that sort is what exhausts
% memory once the asset grid grows. An empty bin gives 0/0=NaN, which is exactly
% the value the old RestrictedSampleMass>0 guard produced.
% GPU arrayfun cannot compile closures that capture outer-workspace vars, so bake
% the bin edges in as numeric literals via str2func+sprintf.
binsig='@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ';
binincfn='HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
FnsToEvaluate_bins=struct();
for ii=1:length(HV2000_multiples)
    loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
    hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
    inbin=sprintf('((%s>=%.16g)*(%s<=%.16g))',binincfn,loVal,binincfn,hiVal);
    FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig '((1+g)*aprime-a)*' inbin]);
    FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
end
AggBins=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions4);
savingRates=zeros(length(HV2000_multiples),1);
for ii=1:length(HV2000_multiples)
    savingRates(ii)=gather(AggBins.(sprintf('Sav%d',ii)).Mean/AggBins.(sprintf('Inc%d',ii)).Mean);
end

%% Store
if Params.borrowFlag==0
    Results.Model4_abar0.Params=Params;
    Results.Model4_abar0.AllStats=AllStats;
    Results.Model4_abar0.GECondns=GECondns;
    Results.Model4_abar0.savingRates=savingRates;
    Results.Model4_abar0.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model4_abar0.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
elseif Params.borrowFlag==1
    Results.Model4_abarw.Params=Params;
    Results.Model4_abarw.AllStats=AllStats;
    Results.Model4_abarw.GECondns=GECondns;
    Results.Model4_abarw.savingRates=savingRates;
    Results.Model4_abarw.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model4_abarw.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
end

%% Within-age-group averaging saving rates (Tables 8/9; only when computeWithinAgeMean==1)
if isfield(Params,'computeWithinAgeMean') && Params.computeWithinAgeMean==1
    AgeStats_WAM=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_M4,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions4);
    % HV2000 Tables 8/9 impose an equal saving RATE within each age group, not an
    % equal saving LEVEL: agent i's counterfactual saving is sbar(j)*y_i, where
    % sbar(j) is age group j's own saving rate (its total saving over its total
    % income, the same aggregate-ratio definition used everywhere else in the paper).
    % The bin rate is then an income-weighted average of the sbar(j) of the agents in
    % the bin, so it stays inside the range of age-group saving rates and rises with
    % income because high-income agents sit in the high-saving-rate ages.
    % Using the age-group mean saving LEVEL as the numerator instead makes the ratio
    % decay like 1/(income multiple) -- the numerator is then capped near the largest
    % age-group mean saving while the denominator grows with the multiple -- which
    % inverts the profile the tables exist to show.
    Params.savingRateAgeMean=AgeStats_WAM.Savings.Mean./AgeStats_WAM.Income.Mean;
    % The bin saving rate is total saving in the band over total income in the band,
    % which is a ratio of two AGGREGATES: the bin mass cancels, so no conditional
    % restriction is needed. Folding the bin indicator into the functions and taking
    % one AggVars pass avoids AllStats entirely -- AllStats sorts the whole grid to
    % get medians/Lorenz/quantiles we never read, and that sort is what exhausts
    % memory once the asset grid grows. An empty bin gives 0/0=NaN, which is exactly
    % the value the old RestrictedSampleMass>0 guard produced.
    % GPU arrayfun cannot compile closures that capture outer-workspace vars, so bake
    % the bin edges in as numeric literals via str2func+sprintf.
    binsig='@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS,savingRateAgeMean) ';
    binincfn='HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
    FnsToEvaluate_bins=struct();
    for ii=1:length(HV2000_multiples)
        loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
        hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
        inbin=sprintf('((%s>=%.16g)*(%s<=%.16g))',binincfn,loVal,binincfn,hiVal);
        FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig 'savingRateAgeMean*' binincfn '*' inbin]);
        FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
    end
    AggBinsWAM=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions4);
    savingRatesWithinAge=zeros(length(HV2000_multiples),1);
    for ii=1:length(HV2000_multiples)
        savingRatesWithinAge(ii)=gather(AggBinsWAM.(sprintf('Sav%d',ii)).Mean/AggBinsWAM.(sprintf('Inc%d',ii)).Mean);
    end
    if Params.borrowFlag==0
        Results.Model4_abar0.savingRatesWithinAge=savingRatesWithinAge;
    elseif Params.borrowFlag==1
        Results.Model4_abarw.savingRatesWithinAge=savingRatesWithinAge;
    end
end

Y=Params.A*AllStats.K.Mean^Params.alpha*AllStats.L.Mean^(1-Params.alpha);
fprintf('Model 4 (borrowFlag=%d) GE: r=%.4f  K/Y=%.4f  S/Y=%.4f  L=%.4f  IncomeGini=%.3f\n',Params.borrowFlag,Params.r,AllStats.K.Mean/Y,AllStats.Savings.Mean/Y,AllStats.L.Mean,AllStats.Income.Gini);

%% Figures 4 and 7 (Model 4, a_underbar=0, baseline SS only). simoptions4 is still in scope.
if Params.borrowFlag==0 && Params.haveSS==1 && Params.doPart1figures==1
    HV2000_Figure4;
    HV2000_Figure7;
end
