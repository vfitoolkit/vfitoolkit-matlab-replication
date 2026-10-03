% HV2000 Model 2 with permanent types (PType), used to produce HV2000 Figure 6.
%
% Background: HV2000 Section 5.3.2 proposes a counterfactual where accidental
% bequests are routed by permanent ability type rather than as a single
% population-wide lump sum. Newborns of type i inherit the average per-period
% bequests of type-i deceased. Fig 6 plots saving rates at multiples of mean
% income for Model 2 baseline ("Equal Transfers") vs this counterfactual
% ("Different Transfers"), in both borrowing-limit cases.
%
% Implementation: switch from the n_z=21 representation to N_i=21 permanent
% types. Each type i has its own ability level Params.alphai(i) and its own
% lump-sum transfer Params.T(i). The bequest-balance GE equation is now
% per-type, handled by heteroagentoptions.GEptype={'BequestBalance'}.

fprintf('\n========== HV2000 Model 2 PType (bequests by ability type, Fig 6) ==========\n');

% Save baseline state we will modify
T_orig=Params.T; % transfer
heteroagentoptions_orig=heteroagentoptions;

%% PType setup
N_i=21;
Params.sigma2_y1=0.45;
[alphai_log,probs_alphai]=discretizeIID_Tauchen(0,sqrt(Params.sigma2_y1),N_i,5);
alphai_levels=exp(alphai_log);
alphai_levels=alphai_levels/sum(probs_alphai.*alphai_levels); % renormalize so E[alphai]=1
Params.alphai=alphai_levels;
Params.alphai_dist=probs_alphai;
PTypeDistParamNames_pt={'alphai_dist'};

% Per-type lump-sum transfer (initial guess: same as the scalar baseline T for every type)
Params.T=T_orig*ones(N_i,1);

%% Within-type stochastic structure: trivial singleton z
n_z_pt=1;
z_grid_pt=1;
pi_z_pt=1;

%% Build Model 2 PType's vfoptions/simoptions from the shared bases (carried as
%% vfoptions5 / simoptions5 to match the per-model numbering of the other scripts).
% Substitute alphai (per-type ability) for z (the original z-grid level) in the underlying functions.
vfoptions5=vfoptions;
vfoptions5.aprimeFn=@(d,ebar,z,agej,Jr,w,ybar_j,ybar_mean,g,alphai) ...
    HuggettVentura2000_ebarprimeFn(d,ebar,alphai,agej,Jr,w,ybar_j,ybar_mean,g);
simoptions5=simoptions;
simoptions5.aprimeFn=vfoptions5.aprimeFn;
simoptions5.z_grid=z_grid_pt;
simoptions5.a_grid=a_grid;
simoptions5.d_grid=d_grid;

ReturnFn_PT=@(d,aprime,a,ebar,z,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS,alphai) ...
    HuggettVentura2000_ReturnFn(d,aprime,a,ebar,alphai,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS);

FnsToEvaluate_PT.K=@(d,aprime,a,ebar,z) a;
FnsToEvaluate_PT.L=@(d,aprime,a,ebar,z,agej,Jr,ybar_j,alphai) (agej<Jr)*alphai*ybar_j;
FnsToEvaluate_PT.Earnings=@(d,aprime,a,ebar,z,agej,Jr,ybar_j,w,alphai) (agej<Jr)*alphai*ybar_j*w;
FnsToEvaluate_PT.AccBeq=@(d,aprime,a,ebar,z,sj,r,tau) (1-sj)*aprime*(1+r*(1-tau));
FnsToEvaluate_PT.SSBenefits=@(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS) ...
    HV2000_SSBenefitFn(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS);
FnsToEvaluate_PT.Savings=@(d,aprime,a,ebar,z,g) (1+g)*aprime-a;
FnsToEvaluate_PT.Ebar=@(d,aprime,a,ebar,z) ebar;
FnsToEvaluate_PT.Consumption=@(d,aprime,a,ebar,z,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS,alphai) ...
    HuggettVentura2000_ConsumptionFn(d,aprime,a,ebar,alphai,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS);
FnsToEvaluate_PT.Income=@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS,alphai) ...
    HV2000_IncomeFn(d,aprime,a,ebar,alphai,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS);

%% Initial distribution at age 1: each type starts at (a=0, ebar=0)
jequaloneDist_PT=zeros([n_a,n_z_pt]);
jequaloneDist_PT(azeroindex,1,1)=1;

%% GE: BequestBalance is per-type
heteroagentoptions.GEptype={'BequestBalance'};
% Return the per-type T as an N_i-vector, not a struct keyed by ptype name. A
% 1x1 struct passes isscalar(), so HV2000_KotlikoffSummersStat would take its
% scalar branch and try to multiply a struct.
heteroagentoptions.GEptype_vectoroutput=1;

% fminsearch (fminalgo=1, the default) has 26 unknowns here (5 aggregate prices
% + 21 per-type T) and does not converge. All six GE eqns are written as
% "price - target", so the shooting map is p_new = p_old - factor*GEcondn:
% a factor of 1 would set the price straight to its implied value, so the 0.5
% factors below are half-steps toward it (damped, for stability).
heteroagentoptions.fminalgo=9; % Anderson acceleration of the shooting map
% Must be an nGEeqns-by-4 cell array (GEcondnName, price name, add, factor).
heteroagentoptions.fminalgo9.howtoupdate={ ...
    'CapitalMarket',  'r',        0, 0.3;           ... % r <- MPK-delta; this is the Aiyagari fixed point and oscillates undamped
    'LaborMarket',    'w',        0, 0.5;           ... % w <- MPL, definitional given K,L
    'TaxRate',        'tau',      0, 0.5;           ... % tau <- GoverY/(1-delta*K/Y), definitional
    'BequestBalance', 'T',        0, 1/(1+Params.n);... % gap is (1+n)*(T-AccBeq/(1+n)), so 1/(1+n) is the exact full step to T=AccBeq/(1+n)
    'SSBalance',      'theta',    0, 0.5;           ... % gap is w*L*(theta-SSBenefits/(w*L)); w*L~0.92 so 1 would be near a full step
    'bCommon',        'b_common', 0, 0.5};              % b_common <- 0.1242*Y, definitional
heteroagentoptions.anderson.memory=5;
heteroagentoptions.anderson.warmup=2;
heteroagentoptions.anderson.maxiter=200;
heteroagentoptions.verbose=1; % print the Anderson residual path, so the r damping can be judged from the diary

%% Loop over borrowFlag, solve GE, compute saving rates by income multiple
HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];

for borrowFlag=0:1
    Params.borrowFlag=borrowFlag;

    [p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z_pt,N_j,N_i,[],pi_z_pt,d_grid,a_grid,z_grid_pt,jequaloneDist_PT,ReturnFn_PT,FnsToEvaluate_PT,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightsParamNames,PTypeDistParamNames_pt,GEPriceParamNames,heteroagentoptions,simoptions5,vfoptions5);

    for pp=1:length(GEPriceParamNames)
        Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
    end
    
    %% Re-solve at GE prices and compute AllStats
    [V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_pt,N_j,N_i,d_grid,a_grid,z_grid_pt,pi_z_pt,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions5);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist_PT,AgeWeightsParamNames,PTypeDistParamNames_pt,Policy,n_d,n_a,n_z_pt,N_j,N_i,pi_z_pt,Params,simoptions5);
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_PT,Params,n_d,n_a,n_z_pt,N_j,N_i,d_grid,a_grid,z_grid_pt,simoptions5);

    %% CheckTopOfGrid (PType): weight per-type marginals by Params.alphai_dist
    % See Results.Model2_ptype_abar{0,w}.CheckTopOfGrid; large values => widen Params.amax / ebar_gridmax.
    mass_by_a=zeros(n_a(1),1); mass_by_ebar=zeros(n_a(2),1);
    for pp=1:N_i
        iistr=sprintf('ptype%03d',pp);
        SD_pp=StationaryDist.(iistr);
        mass_by_a    = mass_by_a    + Params.alphai_dist(pp) * squeeze(sum(SD_pp,  2:ndims(SD_pp)));
        mass_by_ebar = mass_by_ebar + Params.alphai_dist(pp) * squeeze(sum(SD_pp, [1, 3:ndims(SD_pp)]))';
    end
    CheckTopOfGrid.MassAtTop10gridptsAssets = sum(mass_by_a(end-9:end));
    CheckTopOfGrid.MassAtTop5gridptsAssets  = sum(mass_by_a(end-4:end));
    CheckTopOfGrid.MassAtTop2gridptsAssets  = sum(mass_by_a(end-1:end));
    CheckTopOfGrid.MassAtTop5gridptsEbar    = sum(mass_by_ebar(end-4:end));
    CheckTopOfGrid.MassAtTop2gridptsEbar    = sum(mass_by_ebar(end-1:end));

    % %% Check where the income dist mass is (as some of the savings rates by income multiple are coming out as having no mass)
    % ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate_PT,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

    %% Saving rates by income multiple (same AggVars bin pattern as Table 5)
    % Per HV2000 p. 380 (just below Table 5):
    %   "multiples are calculated by taking a 10% band around each income multiple
    %    and then dividing total saving of agents in the band by total income of
    %    agents in the band. Income is defined as earnings after social security
    %    taxes plus interest income and transfers."
    % Bin = [0.9*m*Y, 1.1*m*Y] with Y=Params.IncomeMean (full-population mean);
    % saving rate computed below as Savings.Mean / Income.Mean within the bin
    % (toolkit's .Mean under conditionalrestriction is mass-weighted, so this ratio equals total saving in bin / total income in bin).
    Params.IncomeMean=AllStats.Income.Mean;
    % The bin saving rate is total saving in the band over total income in the band,
    % which is a ratio of two AGGREGATES: the bin mass cancels, so no conditional
    % restriction is needed. That matters doubly here -- the PType grouped restricted
    % stats are the ones that were NaN/mis-weighted before the toolkit fix -- and it
    % avoids AllStats sorting the whole grid for stats we never read. An empty bin
    % gives 0/0=NaN, the value the old RestrictedSampleMass guard produced.
    % alphai stays in the signature so the toolkit auto-extracts the per-type ability
    % level from Params; GPU arrayfun cannot capture outer-workspace vars, so bake the
    % bin edges in as numeric literals via str2func+sprintf.
    binsig='@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS,alphai) ';
    binincfn='HV2000_IncomeFn(d,aprime,a,ebar,alphai,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
    FnsToEvaluate_bins=struct();
    for ii=1:length(HV2000_multiples)
        loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
        hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
        inbin=sprintf('((%s>=%.16g)*(%s<=%.16g))',binincfn,loVal,binincfn,hiVal);
        FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig '((1+g)*aprime-a)*' inbin]);
        FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
    end
    AggBins=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_bins,Params,n_d,n_a,n_z_pt,N_j,N_i,d_grid,a_grid,z_grid_pt,simoptions5);
    savingRates=zeros(length(HV2000_multiples),1);
    for ii=1:length(HV2000_multiples)
        savingRates(ii)=gather(AggBins.(sprintf('Sav%d',ii)).Mean/AggBins.(sprintf('Inc%d',ii)).Mean);
    end

    %% Store
    if borrowFlag==0
        Results.Model2_ptype_abar0.savingRates=savingRates;
        Results.Model2_ptype_abar0.AllStats=AllStats;
        Results.Model2_ptype_abar0.Params=Params;
        Results.Model2_ptype_abar0.GECondns=GECondns;
        Results.Model2_ptype_abar0.CheckTopOfGrid=CheckTopOfGrid;
        Results.Model2_ptype_abar0.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
    elseif borrowFlag==1
        Results.Model2_ptype_abarw.savingRates=savingRates;
        Results.Model2_ptype_abarw.AllStats=AllStats;
        Results.Model2_ptype_abarw.Params=Params;
        Results.Model2_ptype_abarw.GECondns=GECondns;
        Results.Model2_ptype_abarw.CheckTopOfGrid=CheckTopOfGrid;
        Results.Model2_ptype_abarw.TransferWealthPct=HV2000_KotlikoffSummersStat(Params,AllStats.K.Mean);
    end

    fprintf('Model 2 PType (borrowFlag=%d) GE solved\n',borrowFlag);
end

%% Restore baseline scalar T and original heteroagentoptions
Params.T=T_orig;
heteroagentoptions=heteroagentoptions_orig;
