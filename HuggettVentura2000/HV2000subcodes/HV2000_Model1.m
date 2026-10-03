% HV2000 Model 1: deterministic earnings (singleton productivity shock).
% Reads common setup from HuggettVentura2000.m workspace; writes results
% into Results.Model1.

fprintf('\n========== HV2000 Model 1 (deterministic earnings) ==========\n');

%% Model-specific stochastic structure
n_z=1;
z_grid=1;
pi_z=1;

%% Build Model 1's vfoptions/simoptions from the shared bases.
% Downstream scripts that need the Model 1 state can read simoptions1/vfoptions1.
vfoptions1=vfoptions;
simoptions1=simoptions;
simoptions1.z_grid=z_grid;

% Initial distribution at age 1: all newborns at (a=0, ebar=0, z=1)
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(azeroindex,1,1)=1;

%% Solve general equilibrium
[p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d,n_a,n_z,N_j,[],pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions1,vfoptions1);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

%% Re-solve at GE prices and compute AllStats
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);

%% CheckTopOfGrid: mass near the top of the asset / ebar grids
% See Results.Model1_abar{0,w}.CheckTopOfGrid; large values => widen Params.amax / ebar_gridmax.
mass_by_a   = sum(StationaryDist,  2:ndims(StationaryDist));    % marginal over asset dim (dim 1)
mass_by_ebar= sum(StationaryDist, [1, 3:ndims(StationaryDist)]);% marginal over ebar dim (dim 2)
CheckTopOfGrid.MassAtTop10gridptsAssets = sum(mass_by_a(end-9:end));
CheckTopOfGrid.MassAtTop5gridptsAssets  = sum(mass_by_a(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsAssets  = sum(mass_by_a(end-1:end));
CheckTopOfGrid.MassAtTop5gridptsEbar    = sum(mass_by_ebar(end-4:end));
CheckTopOfGrid.MassAtTop2gridptsEbar    = sum(mass_by_ebar(end-1:end));

AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

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
% one AggVars pass avoids AllStats entirely -- AllStats sorts the whole
% (a,ebar,z,j) grid to get medians/Lorenz/quantiles we never read, and that sort
% is what exhausts memory once the asset grid grows. An empty bin gives 0/0=NaN,
% which is exactly the value the old RestrictedSampleMass>0 guard produced.
% GPU arrayfun cannot compile closures that capture outer-workspace vars
% ("Use of functional workspace is not supported"); bake bin edges as numeric
% literals via str2func+sprintf so each closure has an empty workspace.
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
AggBins=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
savingRates=zeros(length(HV2000_multiples),1);
for ii=1:length(HV2000_multiples)
    savingRates(ii)=gather(AggBins.(sprintf('Sav%d',ii)).Mean/AggBins.(sprintf('Inc%d',ii)).Mean);
end

%% Store
if Params.borrowFlag==0
    Results.Model1_abar0.Params=Params;
    Results.Model1_abar0.AllStats=AllStats;
    Results.Model1_abar0.GECondns=GECondns;
    Results.Model1_abar0.savingRates=savingRates;
    Results.Model1_abar0.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model1_abar0.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
elseif Params.borrowFlag==1
    Results.Model1_abarw.Params=Params;
    Results.Model1_abarw.AllStats=AllStats;
    Results.Model1_abarw.GECondns=GECondns;
    Results.Model1_abarw.savingRates=savingRates;
    Results.Model1_abarw.CheckTopOfGrid=CheckTopOfGrid;
    Results.Model1_abarw.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
end

Y=Params.A*AllStats.K.Mean^Params.alpha*AllStats.L.Mean^(1-Params.alpha);
fprintf('Model 1 (borrowFlag=%d) GE: r=%.4f  K/Y=%.4f  S/Y=%.4f  L=%.4f  IncomeGini=%.3f\n',Params.borrowFlag,Params.r,AllStats.K.Mean/Y,AllStats.Savings.Mean/Y,AllStats.L.Mean,AllStats.Income.Gini);

% V, Policy and StationaryDist are no longer needed, and memory available is limited, so clear them
clear V Policy StationaryDist
