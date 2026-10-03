% HV2000 Model 3: persistent AR(1) on log labor endowment with regression
% to age-specific mean (gamma=0.985). Birth draw from N(0,sigma2_y1=0.24).

fprintf('\n========== HV2000 Model 3 (persistent AR(1)) ==========\n');

%% Model-specific stochastic structure (Table 3)
Params.sigma2_y1=0.24;
Params.gamma=0.985;
Params.sigma2_eps1=0.02;
n_z=61;

% AR(1) on log z (gamma persistence, innovation std sigma_eps1).
% HV2000 footnote 14 says grid endpoints at +/- 6*sigma; interpreting as 6 STATIONARY stds
% (the more common reading, and gives coverage ~99.7% of the invariant distribution).
% discretizeAR1_Tauchen scales by the stationary std internally, so Tauchen_q = 6.
[zlog_grid,pi_z]=discretizeAR1_Tauchen(0,Params.gamma,sqrt(Params.sigma2_eps1),n_z,6);
% Note: Tanaka-Toda / Farmer-Toda would give a better discretization than Tauchen,
% but sticking with Tauchen to more closely follow the original paper.
z_grid=exp(zlog_grid);

% Birth distribution: log(z_1) ~ N(0, sigma2_y1) discretized onto the
% AR(1)-built zlog_grid (Tauchen helpers each build their own grid, so we
% hand-discretize via midpoint normal CDF onto the AR(1) grid).
mid=(zlog_grid(1:end-1)+zlog_grid(2:end))/2;
cdf_mid=normcdf(mid,0,sqrt(Params.sigma2_y1));
probs_z=zeros(n_z,1);
probs_z(1)=cdf_mid(1);
probs_z(end)=1-cdf_mid(end);
probs_z(2:end-1)=cdf_mid(2:end)-cdf_mid(1:end-1);

% Under the AR(1), Jensen's inequality on a growing log-variance makes the
% cohort mean of z=exp(zlog) drift with age. Build an age-varying z_grid_J
% by computing the cohort mean at each age via MarkovChainMoments_FHorz and
% renormalizing column-by-column so E[z_j]=1 at every j.
z_grid_J=repmat(z_grid,1,N_j);
pi_z_J=repmat(pi_z,1,1,N_j);
[zmean_J,~,~,~]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,probs_z);
for jj=1:N_j
    z_grid_J(:,jj)=z_grid_J(:,jj)/zmean_J(jj);
end

%% Build Model 3's vfoptions/simoptions from the shared bases.
vfoptions3=vfoptions;
simoptions3=simoptions;
simoptions3.z_grid=z_grid_J;

jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(azeroindex,1,:)=reshape(probs_z,[1,1,n_z]);

%% Solve general equilibrium
% With per-age renormalization of z_grid_J, E[z_j]=1 at every j, so the
% precomputed Params.ybar_mean (per-worker mean of ybar_j) exactly equals
% mean earnings per worker divided by w. No ybar_mean GE price needed.
[p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d,n_a,n_z,N_j,[],pi_z_J,d_grid,a_grid,z_grid_J,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions3,vfoptions3);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

%% Re-solve at GE prices and compute AllStats
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions3);

%% CheckTopOfGrid: mass near the top of the asset / ebar grids
% See Results.Model3_abar{0,w}.CheckTopOfGrid; large values => widen Params.amax / ebar_gridmax.
mass_by_a   = sum(StationaryDist,  2:ndims(StationaryDist));    % marginal over asset dim (dim 1)
mass_by_ebar= sum(StationaryDist, [1, 3:ndims(StationaryDist)]);% marginal over ebar dim (dim 2)
CheckTopOfGrid.MassAtTop10gridptsAssets = sum(mass_by_a(end-9:end));
CheckTopOfGrid.MassAtTop5gridptsAssets  = sum(mass_by_a(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsAssets  = sum(mass_by_a(end-1:end));
CheckTopOfGrid.MassAtTop5gridptsEbar    = sum(mass_by_ebar(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsEbar    = sum(mass_by_ebar(end-1:end));

AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions3);

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
AggBins=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions3);
savingRates=zeros(length(HV2000_multiples),1);
for ii=1:length(HV2000_multiples)
    savingRates(ii)=gather(AggBins.(sprintf('Sav%d',ii)).Mean/AggBins.(sprintf('Inc%d',ii)).Mean);
end

%% Store
if Params.borrowFlag==0
    Results.Model3_abar0.Params=Params;
    Results.Model3_abar0.AllStats=AllStats;
    Results.Model3_abar0.GECondns=GECondns;
    Results.Model3_abar0.savingRates=savingRates;
    Results.Model3_abar0.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model3_abar0.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
elseif Params.borrowFlag==1
    Results.Model3_abarw.Params=Params;
    Results.Model3_abarw.AllStats=AllStats;
    Results.Model3_abarw.GECondns=GECondns;
    Results.Model3_abarw.savingRates=savingRates;
    Results.Model3_abarw.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model3_abarw.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
end

%% Within-age-group averaging saving rates (Tables 8/9; only when computeWithinAgeMean==1)
if isfield(Params,'computeWithinAgeMean') && Params.computeWithinAgeMean==1
    AgeStats_WAM=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions3);
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
    AggBinsWAM=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions3);
    savingRatesWithinAge=zeros(length(HV2000_multiples),1);
    for ii=1:length(HV2000_multiples)
        savingRatesWithinAge(ii)=gather(AggBinsWAM.(sprintf('Sav%d',ii)).Mean/AggBinsWAM.(sprintf('Inc%d',ii)).Mean);
    end
    if Params.borrowFlag==0
        Results.Model3_abar0.savingRatesWithinAge=savingRatesWithinAge;
    elseif Params.borrowFlag==1
        Results.Model3_abarw.savingRatesWithinAge=savingRatesWithinAge;
    end
end

Y=Params.A*AllStats.K.Mean^Params.alpha*AllStats.L.Mean^(1-Params.alpha);
fprintf('Model 3 (borrowFlag=%d) GE: r=%.4f  K/Y=%.4f  S/Y=%.4f  L=%.4f  IncomeGini=%.3f\n',Params.borrowFlag,Params.r,AllStats.K.Mean/Y,AllStats.Savings.Mean/Y,AllStats.L.Mean,AllStats.Income.Gini);

% V, Policy and StationaryDist are no longer needed, and memory available is limited, so clear them
clear V Policy StationaryDist
