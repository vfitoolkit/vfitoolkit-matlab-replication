% HV2000 Model 2: permanent shift in log labor endowment at birth.
% Each agent draws log endowment from N(0,sigma2_y1); fixed for life.
% Implemented as n_z=21 with identity pi_z and a lognormal jequaloneDist.

fprintf('\n========== HV2000 Model 2 (permanent shift at birth) ==========\n');

%% Model-specific stochastic structure (Table 3)
Params.sigma2_y1=0.45;
n_z=21;

% Log-z grid + birth probs: discretize N(0,sigma2_y1) on [-5*sigma_y1, +5*sigma_y1] (HV2000 footnote 14).
[zlog_grid,probs_z]=discretizeIID_Tauchen(0,sqrt(Params.sigma2_y1),n_z,5);
% Note: Tanaka-Toda / Farmer-Toda would give a better discretization than Tauchen,
% but sticking with Tauchen to more closely follow the original paper.
z_grid=exp(zlog_grid);                  % level multiplier exp(y_1-ybar_1)
pi_z=eye(n_z);                          % permanent shock: no transitions

% Renormalize z_grid so E[z]=1 under the birth distribution.
% pi_z is identity (permanent), so the cohort mean of z stays 1 at every age.
z_grid=z_grid/sum(probs_z.*z_grid);

%% Build Model 2's vfoptions/simoptions from the shared bases.
vfoptions2=vfoptions;
simoptions2=simoptions;
simoptions2.z_grid=z_grid;

jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(azeroindex,1,:)=reshape(probs_z,[1,1,n_z]); % (a=0, ebar=0, z drawn)

%% Solve general equilibrium
[p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d,n_a,n_z,N_j,[],pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions2,vfoptions2);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

%% Re-solve at GE prices and compute AllStats
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions2);

%% CheckTopOfGrid: mass near the top of the asset / ebar grids
% See Results.Model2_abar{0,w}.CheckTopOfGrid; large values => widen Params.amax / ebar_gridmax.
mass_by_a   = sum(StationaryDist,  2:ndims(StationaryDist));    % marginal over asset dim (dim 1)
mass_by_ebar= sum(StationaryDist, [1, 3:ndims(StationaryDist)]);% marginal over ebar dim (dim 2)
CheckTopOfGrid.MassAtTop10gridptsAssets = sum(mass_by_a(end-9:end));
CheckTopOfGrid.MassAtTop5gridptsAssets  = sum(mass_by_a(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsAssets  = sum(mass_by_a(end-1:end));
CheckTopOfGrid.MassAtTop5gridptsEbar    = sum(mass_by_ebar(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsEbar    = sum(mass_by_ebar(end-1:end));

AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions2);

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
binsig='@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ';
binincfn='HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
FnsToEvaluate_bins=struct();
for ii=1:length(HV2000_multiples)
    loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
    hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
    inbin=sprintf('((%s>=%.16g)*(%s<=%.16g))',binincfn,loVal,binincfn,hiVal);
    FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig '((1+g)*aprime-a)*' inbin]);
    FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
end
AggBins=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions2);
savingRates=zeros(length(HV2000_multiples),1);
for ii=1:length(HV2000_multiples)
    savingRates(ii)=gather(AggBins.(sprintf('Sav%d',ii)).Mean/AggBins.(sprintf('Inc%d',ii)).Mean);
end

%% Store
if Params.borrowFlag==0
    Results.Model2_abar0.Params=Params;
    Results.Model2_abar0.AllStats=AllStats;
    Results.Model2_abar0.GECondns=GECondns;
    Results.Model2_abar0.savingRates=savingRates;
    Results.Model2_abar0.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model2_abar0.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
elseif Params.borrowFlag==1
    Results.Model2_abarw.Params=Params;
    Results.Model2_abarw.AllStats=AllStats;
    Results.Model2_abarw.GECondns=GECondns;
    Results.Model2_abarw.savingRates=savingRates;
    Results.Model2_abarw.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model2_abarw.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
end

%% Within-age-group averaging saving rates (Tables 8/9; only when computeWithinAgeMean==1)
if isfield(Params,'computeWithinAgeMean') && Params.computeWithinAgeMean==1
    AgeStats_WAM=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions2);
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
    binsig='@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS,savingRateAgeMean) ';
    binincfn='HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
    FnsToEvaluate_bins=struct();
    for ii=1:length(HV2000_multiples)
        loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
        hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
        inbin=sprintf('((%s>=%.16g)*(%s<=%.16g))',binincfn,loVal,binincfn,hiVal);
        FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig 'savingRateAgeMean*' binincfn '*' inbin]);
        FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
    end
    AggBinsWAM=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions2);
    savingRatesWithinAge=zeros(length(HV2000_multiples),1);
    for ii=1:length(HV2000_multiples)
        savingRatesWithinAge(ii)=gather(AggBinsWAM.(sprintf('Sav%d',ii)).Mean/AggBinsWAM.(sprintf('Inc%d',ii)).Mean);
    end
    if Params.borrowFlag==0
        Results.Model2_abar0.savingRatesWithinAge=savingRatesWithinAge;
    elseif Params.borrowFlag==1
        Results.Model2_abarw.savingRatesWithinAge=savingRatesWithinAge;
    end
end

Y=Params.A*AllStats.K.Mean^Params.alpha*AllStats.L.Mean^(1-Params.alpha);
fprintf('Model 2 (borrowFlag=%d) GE: r=%.4f  K/Y=%.4f  S/Y=%.4f  L=%.4f  IncomeGini=%.3f\n',Params.borrowFlag,Params.r,AllStats.K.Mean/Y,AllStats.Savings.Mean/Y,AllStats.L.Mean,AllStats.Income.Gini);

%% Figure 3 (Model 2, a_underbar=0, baseline SS only)
if Params.borrowFlag==0 && Params.haveSS==1 && Params.doPart1figures==1
    HV2000_Figure3;
end

%% Figure 5 (Model 2, a_underbar=0, NO social security - same shape as Fig 3 but with haveSS=0)
if Params.borrowFlag==0 && Params.haveSS==0
    HV2000_Figure5;
end

% V, Policy and StationaryDist are no longer needed, and memory available is limited, so clear them
clear V Policy StationaryDist
