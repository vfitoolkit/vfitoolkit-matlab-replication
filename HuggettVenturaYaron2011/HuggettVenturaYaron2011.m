% Huggett, Ventura & Yaron (2011) - Sources of Lifetime Inequality
% Two endogenous states, the first is assets, the second is human capital.
% Human capital is modelled as an 'experienceassetu' (u=uncertain due to an iid shock between periods)
%
% The authors' replication files are in ./data/ (from the AER website):
%   data/hvy_FIGS: PSID extracts and the code for Figures 1 and 2 (Master_Figures.m)
%   data/hvy_Table: code for Table 1 (estimation of the human capital shocks)
%   data/hvy-programfiles: the Fortran model code (hcapshocks-ge.f solves the model, distrib-estimation.f
%       estimates the initial distribution) and the calibration targets (mean-time-aer08, var-time-aer08, skew-time-aer08)
% We follow the authors' Fortran code wherever the paper leaves a detail open; these are noted where they occur.

%% Run infrastructure (diary, paths, output folders)
% Record the full run (incl. any error) to HVY2011_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./HVY2011_Diary.txt','file'), delete('./HVY2011_Diary.txt'); end
diary ./HVY2011_Diary.txt

% Subcodes (ReturnFn, ConsFn, NetTaxesFn, jequaloneDistFn, ability grid, data work) live in ./HVY2011subcodes/;
% MATLAB does not search subfolders of the working directory, so put it on
% the path. Note this does NOT change cwd, so the './SavedOutput/...' paths
% still resolve to here.
addpath('./HVY2011subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
if exist('./VFIToolkit-matlab','dir'); addpath(genpath('./VFIToolkit-matlab/')); end

% Output folders
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

Params.headlessFigures=1; % set to 1 to render figures offscreen (needed for headless MATLAB, e.g. the server); PNGs still saved
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end
CreateFigures=1;

doPart=[0,1,1,0,1,1,1,0,0,1,1,1,0];
% Default: reproduce all the results that need a GPU from the saved calibrations (SavedOutput/HVY2011calib.mat and
% HVY2011_Fig6.mat); doPart(7) only re-solves the saved alphas (AlphaToRun=0). The data work, Table 1 and the LaTeX tables
% (doPart 1, 9, 8) need no GPU and read the authors' files in ./data, and are run on a desktop.
% To run everything from scratch: (1) and (9) on a desktop; then (2), (3), (4) with CalibStartFromFig6=0 (from Table 2), (6),
% (7) with AlphaToRun=[1,1,1,1], then (4) again with CalibStartFromFig6=1 (from the best alpha of Fig 6), (5), (6), (7) with
% AlphaToRun=[0,0,0,0], (10), (11), (12); then (8) on a desktop.
% doPart(1): data work: PSID earnings age profiles (Fig 1), the calibration targets      -> SavedOutput/HVY2011_Data.mat
%            (no GPU needed; reads the PSID extracts in ./data/hvy_FIGS)
% doPart(2): solve the model at the HVY2011 Table 2 parameters (incl. beta=0.981): fit to the
%            targets, K/L, net taxes, grid checks                                       -> SavedOutput/HVY2011_Table2check.mat
% doPart(3): general eqm: beta such that K/L is that implied by r=0.042 (HVY2011 Appendix A),
%            at the Table 2 parameters                                                  -> SavedOutput/HVY2011_GE.mat
% doPart(4): calibrate the OLG (alpha_h and the initial dist over (h1,ability), with beta from
%            the general eqm), starting from the Table 2 values (CalibStartFromFig6=0) or the best
%            alpha of Fig 6 (CalibStartFromFig6=1)                                       -> SavedOutput/HVY2011calib.mat
% doPart(5): results at the calibrated parameters (Figs 2 and 3 of HVY2011)             -> SavedOutput/HVY2011full.mat
% doPart(6): the HVY2011 results (Figs 2-5, Tables 3, 5, 6) at BOTH parameter sets: our calibration
%            (HVY2011calib.mat) and HVY2011 Table 2 (with beta from the general eqm of doPart(3)) -> SavedOutput/HVY2011_Results.mat
% doPart(7): HVY2011 Figure 6: recalibrate the initial distribution at alpha=0.5,0.6,0.8,0.9 (about 2.5 hours each), and the
%            fraction of the variance of lifetime earnings due to initial conditions (needs doPart(4) and doPart(6) outputs) -> SavedOutput/HVY2011_Fig6.mat
% doPart(8): LaTeX tables for LatexDoc/App_HuggettVenturaYaron2011.tex, from the saved results (no GPU needed, but reads the authors' files in ./data) -> SavedOutput/LatexInputs/
% doPart(9): HVY2011 Table 1, estimation of human capital shocks from PSID wage data (no GPU needed; reads ./data/hvy_Table) -> SavedOutput/HVY2011_Table1.mat
% doPart(10): HVY2011 Table 4, RIP and HIP statistical earnings models estimated on 500 simulated panels of 200 agents, at both
%            parameter sets                                                             -> SavedOutput/HVY2011_Table4.mat
% doPart(11): the model with initial wealth differences (HVY2011 Section V.A): calibrate the trivariate initial distribution in
%            partial equilibrium, then Tables 3, 5 and 6, at both parameter sets                -> SavedOutput/HVY2011_Wealth.mat
% doPart(12): HVY2011 Section V.C sensitivity: cohort effects, low and high shocks, quantile ratio moments (recalibrate each,
%            then the Table 5 fractions and Table 6), relative to our calibration         -> SavedOutput/HVY2011_Sensitivity.mat
% doPart(13): identification of the initial wealth correlations: moment derivatives at the calibration of doPart(11), and the
%            objective at the authors' correlations (about 15 model solves)               -> SavedOutput/HVY2011_Identification.mat
% A part that is switched off loads what later parts need from its .mat (only named
% variables, never a whole workspace: .mat files store stale function handles).
% Setup (parameters, grids, functions) always runs, it is cheap.

figure_c=0;

%% Grid sizes
n_d=31; % fraction of time working, l (the rest, 1-l, is time learning)
n_a=[201,301]; % Assets, Human capital
n_z=0;
N_i=11; % number of ability levels (the authors' code uses an 11 point grid on ability)
n_u=5; % iid between period shocks to human capital accumulation (HVY2011 Appendix B says 5 points using Tauchen method from +-2 std dev)

vfoptions.experienceassetu=1; % Using an experience asset
simoptions.experienceassetu=1;

N_j=53; % total number of periods
Params.agejshifter=22; % starting age 23

%% Parameters (HVY2011 Table 2)

% Demographics
Params.J=N_j; % agents live from age 23 to 75
Params.agej=1:1:Params.J; % the model period
Params.Jr=39; % first period of retirement, age 61 (HVY2011 Table 2: J_R=39); agents work in periods 1 to 38, ages 23 to 60
Params.age=Params.agejshifter+Params.agej; % for graphs
Params.n=0.012; % population growth rate

% Preferences
Params.gamma=2; % CES utility parameter (HVY2011 call this rho)
Params.beta=0.981; % discount factor (gets determined in general eqm, to deliver r=0.042)

% Technology
Params.alpha_k=0.322; % capital share (HVY2011 call this gamma)
Params.r=0.042; % real interest rate (HVY2011 set prices, and choose beta so these are the eqm prices)
Params.KYratio=2.947; % capital-output ratio, 1959-2000 (HVY2011 footnote 17)
Params.delta_k=Params.alpha_k/Params.KYratio-Params.r; % depreciation, so that r=F_K-delta at K/Y=2.947 [HVY2011 Table 2 report delta=0.067; this is 0.0673]
Params.g=0.0019; % growth rate of labor-augmenting technology: the rental rate of human capital of a cohort grows at rate g with age
% Capital per unit of human capital labor (K/L) implied by r: r+delta=F_K gives K/(LA), and the rental rate
% R=A*F_L is normalized to one, which pins down A; together K/L=(alpha_k/(1-alpha_k))/(r+delta). (as in the authors' code)
Params.KLtarget=(Params.alpha_k/(1-Params.alpha_k))/(Params.r+Params.delta_k);

% Human capital
Params.alpha_h=0.7; % returns to scale of human capital production fn (gets determined by calibration)
% human capital shocks z~N(mu_z,sigma_z^2), HVY2011 Table 2 [sigma_z from their Table 1; mu_z set to match the fall in mean earnings at the end of working life]
Params.mu_z=-0.029;
Params.sigma_z=0.111;

% Initial distribution over (h1,ability), bivariate log-normal (HVY2011 Table 2; these get determined by calibration)
Params.mean_logh1=4.66;
Params.mean_logability=-1.12;
Params.stddev_logh1=sqrt(0.213);
Params.stddev_logability=sqrt(0.012);
Params.corr_logh1logability=0.041/(Params.stddev_logh1*Params.stddev_logability); % Table 2 gives the covariance, 0.041, so this is 0.811
% Our parametrization of the covariance matrix follows AH2021, corr=tanh(FTcorr) for a 2x2 matrix
Params.FTcorr_logh1logability=0.5*log((1+Params.corr_logh1logability)/(1-Params.corr_logh1logability));

% Tax system (HVY2011 Table 2 and Section III; the details are from the authors' code hcapshocks-ge.f)
Params.tau_ss=0.106; % social security payroll tax
% Income tax: average tax rate (in percent) of a*x-b*x^2, where x is income in tens of thousands of dollars
% (HVY2011 say it follows the CBO (2004) effective average federal income tax rates, implemented as in Huggett & Parra (2010))
Params.taxw1=2.1380; % before retirement (income is earnings plus interest)
Params.taxw2=0.04785;
Params.taxr1=1.8484; % in retirement (income is the transfer plus interest)
Params.taxr2=0.04273;
% The retirement transfer and the scale of the tax schedule depend on the data on mean earnings, set below

%% Grids
% Assets: borrowing is allowed (HVY2011 Appendix B impose ability-specific borrowing limits, and relax them until
% they do not bind; here the -Inf return to negative consumption, plus no debt at death, imposes the natural borrowing limit)
asset_grid=[-1500*(linspace(1,0,31)'.^2); 8000*(linspace(0,1,171)'.^3)];
asset_grid=unique(asset_grid); % the two pieces share the point zero, so this is 201 points
% Human capital: log-spaced, so the grid spacing is the same proportion (2.5 percent) everywhere. The experience asset puts
% next period human capital on the two neighbouring grid points, which adds noise to log human capital each period, so the
% spacing has to be fine relative to h (the authors' code instead interpolates the decision rules, on ability-specific grids)
h_grid=exp(linspace(log(2),log(4000),n_a(2))');
% Time working: the authors' code puts extra points near zero learning time
timelearning=[0; 0.004; 0.02; linspace(0.04,1,n_d-3)'];
l_grid=sort(1-timelearning);

% rename in toolkit notation
d_grid=l_grid;
a_grid=[asset_grid; h_grid];
z_grid=[];
pi_z=[];

%% Set up the permanent types
% Ability is correlated with initial human capital, but we just ignore this
% correlation while we set up the grid on ability
% (we will consider the correlation later when putting agents onto the grids)
% The grid on ability is parametrized by mean_logability and stddev_logability, which
% are calibrated. This idea of a PType parameter having a parametrized distribution
% across PTypes is common, and the way VFI Toolkit handles it is with PTypeParamFn.
Params.Number_i=N_i; % so we can pass it as an input to PTypeParamFn (and jequaloneDistFn)
PTypeParamFn=@(Params) HVY2011_AbilityGrid(Params);
% PTypeParamFn will always have Params as it's first and only entry
% Because any parameter that depends on PType will have to be kept in
% Params, the output of the function will also just be Params.
Params=PTypeParamFn(Params); % creates Params.ability (a parameter that depends on permanent type) and Params.abilitydist

PTypeDistParamNames={'abilitydist'};
% Note: the masses of the ptypes come from the initial distribution (jequaloneDist has ptype as a dimension, so its
% slice masses become the ptype masses), which accounts for the correlation between initial human capital and ability

%% Set up the experience asset, we need to aprimeFn which gives value of aprime given d2 and a2 (d2 is the decision variable relevant to experience asset, a2 is the experience asset)
% HVY2011: h'=exp(z')*H(h,s,a), H(h,s,a)=h+a*(h*s)^alpha, where s=1-l is time learning
vfoptions.aprimeFn=@(l,h,u,ability,alpha_h) exp(u)*(h+ability*(h*(1-l))^alpha_h);
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.a_grid=a_grid; % Needed when using experience asset
simoptions.d_grid=d_grid; % Needed when using experience asset

[u_grid,pi_u]=discretizeAR1_Tauchen(Params.mu_z,0,Params.sigma_z,n_u,2); % HVY2011 Appendix B: Tauchen, five equally spaced points on [mu-2sigma,mu+2sigma]
pi_u=pi_u(1,:)'; % iid

vfoptions.n_u=n_u;
vfoptions.u_grid=u_grid;
vfoptions.pi_u=pi_u;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;

%% Human capital production, to check the top of the human capital grid
if CreateFigures==1
    figure_c=figure_c+1;
    fig=figure(figure_c);
    abilityindexes=[1,ceil(N_i/2),N_i];
    abilitynames={'low ability','mid ability','high ability'};
    for pp=1:3
        subplot(3,1,pp);
        plot(h_grid,exp(u_grid(1))*(h_grid+Params.ability(abilityindexes(pp))*(h_grid*0.5).^Params.alpha_h))
        hold on
        plot(h_grid,exp(u_grid(ceil(n_u/2)))*(h_grid+Params.ability(abilityindexes(pp))*(h_grid*0.5).^Params.alpha_h))
        plot(h_grid,exp(u_grid(end))*(h_grid+Params.ability(abilityindexes(pp))*(h_grid*0.5).^Params.alpha_h))
        plot(h_grid,h_grid)
        hold off
        legend('low u shock','mid u shock','high u shock','45 degree','Location','northwest')
        title(['human capital production, half time learning: ',abilitynames{pp}])
    end
    saveas(fig,'./SavedOutput/Graphs/HVY2011_HumanCapitalProduction.png')
end

%% Initial agent distribution as a function, as we need that so that the parameters of this can be part of our calibration
% jequaloneDist as a function has first four inputs (a_grid,z_grid,n_a,n_z,...) followed by any parameters
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability,mu_z,sigma_z,ability,Number_i) HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability,mu_z,sigma_z,ability,Number_i);
% To use jequaloneDist as a function you must include a_grid and z_grid in simoptions
simoptions.a_grid=a_grid;
simoptions.z_grid=[]; % Even though we don't use z_grid, it is still needed in simoptions to be able to call jequaloneDist as a function
% Params.ability and Params.Number_i are passed as inputs to jequaloneDistFn, We already set these above

%% Age weights
% No mortality before age 75 (HVY2011), so the mass of each age just reflects population growth
Params.mewj=(1/(1+Params.n)).^(0:1:Params.J-1);
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

%% doPart(1): Data work, the earnings age profiles of HVY2011 Figure 1 (these are the calibration targets)
if doPart(1)==1
    HVY2011_DataWork
    save ./SavedOutput/HVY2011_Data.mat agedata Fig1_Mean_time Fig1_VarLog_time Fig1_Gini_time Fig1_Skew_time Fig1_Mean_cohort Fig1_VarLog_cohort Fig1_Gini_cohort Fig1_Skew_cohort CellN Data_P9010_time Data_P7525_time Data_P9010_cohort Data_P7525_cohort
else
    fprintf('Working directory: %s \n',pwd)
    if ~exist('./SavedOutput/HVY2011_Data.mat','file')
        error('./SavedOutput/HVY2011_Data.mat not found (see working directory above); run doPart(1) or copy it there')
    end
    load('./SavedOutput/HVY2011_Data.mat','agedata','Fig1_Mean_time','Fig1_VarLog_time','Fig1_Gini_time','Fig1_Skew_time','Fig1_Mean_cohort','Fig1_VarLog_cohort','Fig1_Gini_cohort','Fig1_Skew_cohort','CellN','Data_P9010_time','Data_P7525_time','Data_P9010_cohort','Data_P7525_cohort')
end

% Parameters that depend on the data on mean earnings (the authors' code)
% Retirement transfer: 40 percent of mean earnings at the end of the working life (HVY2011 Section III). The time-effects
% mean earnings profile is normalized to 100 at age 60, and the rental rate is one at age 23, so this is 0.4*100 (the authors' code).
Params.ssbenefit=0.4*Fig1_Mean_time(end);
% Income tax schedule: model units are converted to dollars by setting mean earnings (weighted by the working-age
% population) equal to $32,922 (the authors' code)
popweights_working=Params.mewj(1:Params.Jr-1)/sum(Params.mewj(1:Params.Jr-1));
Params.taxfactor=32922/sum(popweights_working'.*Fig1_Mean_time);

%% Value function and return function
DiscountFactorParamNames={'beta'};

ReturnFn=@(l,aprime,a,h,r,gamma,agej,Jr,J,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)...
    HuggettVenturaYaron2011_ReturnFn(l,aprime,a,h,r,gamma,agej,Jr,J,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2);

% Solution method: divide-and-conquer plus grid interpolation layer
vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=50;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

vfoptions.verbose=0;
simoptions.verbose=0;

%% Functions to evaluate
% Earnings are Rj*h*l, where Rj=(1+g)^(j-1) is the rental rate of human capital for the cohort at age j
FnsToEvaluate.earnings=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr);
FnsToEvaluate.logearnings=@(l,aprime,a,h,agej,Jr,g) log(max(((1+g)^(agej-1))*h*l,10^(-3)))*(agej<Jr); % (zero in retirement)
FnsToEvaluate.humancapital=@(l,aprime,a,h) h; % h is human capital
FnsToEvaluate.assets=@(l,aprime,a,h) a; % assets
FnsToEvaluate.timestudying=@(l,aprime,a,h,agej,Jr) (1-l)*(agej<Jr); % 1-l is fraction of time invested in human capital
FnsToEvaluate.zeroearnings=@(l,aprime,a,h,agej,Jr) (l==0)*(agej<Jr); % workers who spend all their time learning (their log earnings is log(10^(-3)))
FnsToEvaluate.consumption=@(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2);
FnsToEvaluate.logconsumption=@(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) log(max(HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2),10^(-10))); % HVY2011 Fig 5 (the max only matters at grid points nobody reaches, where consumption at the policy can be negative)

%% General eqm (HVY2011 Appendix A)
% Given prices (r=0.042, and the rental rate), find the discount factor beta so that the capital-labor ratio of the
% households equals that implied by the prices. Aggregates are for the cross-section at a point in time: assets
% and taxes of age j are divided by (1+g)^(j-1); labor is human capital times time working.
GEPriceParamNames={'beta'};

FnsToEvaluate_GE.K=@(l,aprime,a,h,agej,g) a/((1+g)^(agej-1));
FnsToEvaluate_GE.L=@(l,aprime,a,h,agej,Jr) h*l*(agej<Jr);

GeneralEqmEqns.capitallabor=@(K,L,KLtarget) K/L-KLtarget;

heteroagentoptions.fminalgo=8; % Matlab lsqnonlin()
heteroagentoptions.verbose=1;

% HVY2011 also require that net taxes are positive (the government consumes the difference)
FnsToEvaluate_Check=FnsToEvaluate_GE;
FnsToEvaluate_Check.nettaxes=@(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) HuggettVenturaYaron2011_NetTaxesFn(l,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2);

%% Calibration targets (HVY2011 Appendix D): mean earnings, variance of log earnings, and mean/median earnings,
% ages 23-60 (model periods 1-38), time-effects view; the objective is the sum of squared log ratios
TargetMoments.AgeConditionalStats.earnings.Mean=[Fig1_Mean_time', nan(1,N_j-length(Fig1_Mean_time))];
TargetMoments.AgeConditionalStats.logearnings.Variance=[Fig1_VarLog_time', nan(1,N_j-length(Fig1_VarLog_time))];
TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian=[Fig1_Skew_time', nan(1,N_j-length(Fig1_Skew_time))];
Nworkages=Params.Jr-1; % 38 ages with targets

% Keep the HVY2011 Table 2 parameters (later parts overwrite Params with the general eqm beta and the calibrated values)
Params_Table2=Params;


%% doPart(2): Solve the model at the HVY2011 Table 2 parameters
if doPart(2)==1
    tic;
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    vftime=toc;
    tic;
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistFn,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
    disttime=toc;
    tic;
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_Check,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    statstime=toc;
    fprintf('Run times: value fn %8.1fs, agent dist %8.1fs, stats %8.1fs \n',vftime,disttime,statstime)

    % Fit to the calibration targets
    Table2check.age=Params.age(1:Nworkages);
    Table2check.earningsMean=gather(AgeConditionalStats.earnings.Mean(1:Nworkages));
    Table2check.logearningsVariance=gather(AgeConditionalStats.logearnings.Variance(1:Nworkages));
    Table2check.earningsRatioMeanToMedian=gather(AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages));
    Table2check.earningsGini=gather(AgeConditionalStats.earnings.Gini(1:Nworkages));
    Table2check.timestudying=gather(AgeConditionalStats.timestudying.Mean(1:Nworkages));
    Table2check.humancapital=gather(AgeConditionalStats.humancapital.Mean(1:Nworkages));
    Table2check.zeroearnings=gather(AgeConditionalStats.zeroearnings.Mean(1:Nworkages));
    Table2check.minlogearnings=gather(AgeConditionalStats.logearnings.Minimum(1:Nworkages));
    Table2check.objective1=sum(log(Table2check.earningsMean./Fig1_Mean_time').^2);
    Table2check.objective2=sum(log(Table2check.logearningsVariance./Fig1_VarLog_time').^2);
    Table2check.objective3=sum(log(Table2check.earningsRatioMeanToMedian./Fig1_Skew_time').^2);
    Table2check.objective=Table2check.objective1+Table2check.objective2+Table2check.objective3;
    Table2check.KLratio=gather(AggVars.K.Mean/AggVars.L.Mean);
    Table2check.nettaxes=gather(AggVars.nettaxes.Mean);

    fprintf('\nModel at the HVY2011 Table 2 parameters (beta=%6.4f, alpha_h=%4.2f) \n',Params.beta,Params.alpha_h)
    fprintf('Calibration objective (sum of squared log ratios, ages 23-60): %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',Table2check.objective,Table2check.objective1,Table2check.objective2,Table2check.objective3)
    fprintf('K/L: model %8.4f, implied by r=0.042: %8.4f (if the model K/L is above the target, beta is too high) \n',Table2check.KLratio,Params.KLtarget)
    fprintf('Net taxes (per person, cross-section): %8.4f (HVY2011 require this to be positive) \n',Table2check.nettaxes)
    fprintf('    age   mean earnings (model, data)   var log earnings (model, data)   mean/median (model, data)   time learning   human capital   fraction zero earnings   min log earnings \n')
    for jj=1:Nworkages
        fprintf('    %3i   %8.2f %8.2f   %8.4f %8.4f   %8.4f %8.4f   %8.4f   %8.2f   %8.2e   %8.3f \n',Table2check.age(jj),Table2check.earningsMean(jj),Fig1_Mean_time(jj),Table2check.logearningsVariance(jj),Fig1_VarLog_time(jj),Table2check.earningsRatioMeanToMedian(jj),Fig1_Skew_time(jj),Table2check.timestudying(jj),Table2check.humancapital(jj),Table2check.zeroearnings(jj),Table2check.minlogearnings(jj))
    end
    % The authors' model output (data/hvy_FIGS/ModeloutputR2.xls, sheet time-alpha0.70), for comparison
    fprintf('HVY2011 model (authors'' output, ModeloutputR2.xls), ages 23, 30, 40, 50, 60: mean earnings 86.66, 95.14, 109.62, 119.23, 106.13; var log earnings 0.2637, 0.3434, 0.4038, 0.4252, 0.4867 \n')
    fprintf('This model, same ages:                                              mean earnings %6.2f, %6.2f, %6.2f, %6.2f, %6.2f; var log earnings %6.4f, %6.4f, %6.4f, %6.4f, %6.4f \n',Table2check.earningsMean([1,8,18,28,38]),Table2check.logearningsVariance([1,8,18,28,38]))

    % The authors' code (hcapshocks-ge.f) hard-codes an initial distribution that differs slightly from Table 2:
    % muh=4.67, sigmah=0.45, sigmaa=0.12, corr=0.78 (Table 2: 4.66, sqrt(0.213)=0.4615, sqrt(0.012)=0.1095, 0.811).
    % Their model output is presumably at these values, so solve the model there too, for a like-for-like comparison.
    Params_authorscode=Params;
    Params_authorscode.mean_logh1=4.67;
    Params_authorscode.stddev_logh1=0.45;
    Params_authorscode.stddev_logability=0.12;
    Params_authorscode.FTcorr_logh1logability=0.5*log((1+0.78)/(1-0.78));
    Params_authorscode=PTypeParamFn(Params_authorscode); % ability grid for stddev_logability=0.12
    [~, Policy_authorscode]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params_authorscode, DiscountFactorParamNames, vfoptions);
    StationaryDist_authorscode=StationaryDist_Case1_FHorz_PType(jequaloneDistFn,AgeWeightsParamNames,PTypeDistParamNames,Policy_authorscode,n_d,n_a,n_z,N_j,N_i,[],Params_authorscode,simoptions);
    AgeConditionalStats_authorscode=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_authorscode,Policy_authorscode,FnsToEvaluate,Params_authorscode,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    Table2check.authorscode_earningsMean=gather(AgeConditionalStats_authorscode.earnings.Mean(1:Nworkages));
    Table2check.authorscode_logearningsVariance=gather(AgeConditionalStats_authorscode.logearnings.Variance(1:Nworkages));
    Table2check.authorscode_earningsRatioMeanToMedian=gather(AgeConditionalStats_authorscode.earnings.RatioMeanToMedian(1:Nworkages));
    fprintf('This model at the initial distribution in the authors'' code:            mean earnings %6.2f, %6.2f, %6.2f, %6.2f, %6.2f; var log earnings %6.4f, %6.4f, %6.4f, %6.4f, %6.4f \n',Table2check.authorscode_earningsMean([1,8,18,28,38]),Table2check.authorscode_logearningsVariance([1,8,18,28,38]))
    clear Policy_authorscode StationaryDist_authorscode AgeConditionalStats_authorscode

    % Check the grids: mass at the bottom of the asset grid (borrowing limit), near the top of the asset grid,
    % near the top of the human capital grid, and workers who spend all their time learning (zero earnings)
    massbottoma=zeros(1,N_j); masstopa=zeros(1,N_j); masstoph=zeros(1,N_j);
    for ii=1:N_i
        Dist_ii=gather(StationaryDist.(sprintf('ptype%03i',ii))); % [n_a(1),n_a(2),N_j]
        Dist_ii=Dist_ii./max(sum(sum(Dist_ii,1),2),eps); % distribution at each age, within ptype
        massbottoma=massbottoma+StationaryDist.ptweights(ii)*reshape(sum(sum(Dist_ii(1:2,:,:),1),2),[1,N_j]);
        masstopa=masstopa+StationaryDist.ptweights(ii)*reshape(sum(sum(Dist_ii(end-4:end,:,:),1),2),[1,N_j]);
        masstoph=masstoph+StationaryDist.ptweights(ii)*reshape(sum(sum(Dist_ii(:,end-4:end,:),1),2),[1,N_j]);
    end
    Table2check.maxmassbottoma=max(massbottoma);
    Table2check.maxmasstopa=max(masstopa);
    Table2check.maxmasstoph=max(masstoph);
    Table2check.minearnings=gather(min(AgeConditionalStats.earnings.Minimum(1:Nworkages)));
    fprintf('Grid checks (max over ages): mass at the bottom two asset points %8.2e, top five asset points %8.2e, top five human capital points %8.2e; min earnings of workers %8.4f \n',Table2check.maxmassbottoma,Table2check.maxmasstopa,Table2check.maxmasstoph,Table2check.minearnings)
    fprintf('Ptype masses (from the initial distribution): min %8.2e, max %8.2e \n',min(StationaryDist.ptweights),max(StationaryDist.ptweights))

    if CreateFigures==1
        figure_c=figure_c+1;
        fig=figure(figure_c);
        fig.Position=[100 100 900 700];
        subplot(2,2,1); plot(Table2check.age,Table2check.earningsMean,'-k+',Table2check.age,Fig1_Mean_time,'-r')
        xlim([23,60]); ylim([0,140]); xlabel('Age'); title('Mean earnings'); legend({'Model','Data-time'},'Location','southeast')
        subplot(2,2,2); plot(Table2check.age,Table2check.logearningsVariance,'-k+',Table2check.age,Fig1_VarLog_time,'-r')
        xlim([23,60]); ylim([0,0.9]); xlabel('Age'); title('Variance of log earnings')
        subplot(2,2,3); plot(Table2check.age,Table2check.earningsGini,'-k+',Table2check.age,Fig1_Gini_time,'-r')
        xlim([23,60]); ylim([0,0.6]); xlabel('Age'); title('Earnings Gini')
        subplot(2,2,4); plot(Table2check.age,Table2check.earningsRatioMeanToMedian,'-k+',Table2check.age,Fig1_Skew_time,'-r')
        xlim([23,60]); ylim([0.8,1.6]); xlabel('Age'); title('Earnings skewness (mean/median)')
        saveas(fig,'./SavedOutput/Graphs/HVY2011_Table2check.png')

        figure_c=figure_c+1;
        fig=figure(figure_c);
        subplot(3,2,1); plot(Params.age, AgeConditionalStats.humancapital.Mean)
        title('human capital'); xlim([23,75])
        subplot(3,2,2); plot(Params.age, AgeConditionalStats.assets.Mean)
        title('assets'); xlim([23,75])
        subplot(3,2,3); plot(Params.age, AgeConditionalStats.timestudying.Mean)
        title('time studying'); xlim([23,75])
        subplot(3,2,4); plot(Params.age, AgeConditionalStats.earnings.Mean)
        title('earnings'); xlim([23,75])
        subplot(3,2,5); plot(Params.age, AgeConditionalStats.consumption.Mean)
        title('consumption'); xlim([23,75])
        subplot(3,2,6); plot(Params.age, AgeConditionalStats.assets.Minimum)
        title('minimum assets'); xlim([23,75])
        saveas(fig,'./SavedOutput/Graphs/HVY2011_Table2check_LifeCycle.png')
    end

    save ./SavedOutput/HVY2011_Table2check.mat Table2check
end


%% doPart(3): General eqm (beta) at the Table 2 parameters
if doPart(3)==1
    tic;
    [p_eqm,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDistFn, ReturnFn, FnsToEvaluate_GE, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    getime=toc;
    fprintf('General eqm at the Table 2 parameters: beta=%8.5f (HVY2011 Table 2: 0.981), K/L-target=%8.2e, run time %8.1fs \n',p_eqm.beta,GeneralEqmConditions.capitallabor,getime)
    save ./SavedOutput/HVY2011_GE.mat p_eqm GeneralEqmConditions
elseif any(doPart(4:5)==1)
    load('./SavedOutput/HVY2011_GE.mat','p_eqm','GeneralEqmConditions')
end
if any(doPart(3:5)==1)
    Params.beta=p_eqm.beta; % the initial guess for beta in the calibration
end


%% Calibrate the OLG

% To do this using the toolkit there are two things we need to setup
% First, just say all the model parameters we want to calibrate
% For exercise of HVY2011 these are alpha_h and all the parameters that define the initial distribution over
% human capital and ability (beta is determined by the general eqm, as part of the calibration).
CalibParamNames={'alpha_h','mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTcorr_logh1logability'};
% CalibParamNames gives the names of all the model parameters that will be
% calibrated, these parameters must all appear in Params, and the values in
% Params will be used as the initial values.
% All other parameters in Params will remain fixed.

caliboptions.constrain0to1={'alpha_h'};
% The standard deviations must be positive
caliboptions.constrainpositive={'stddev_logh1','stddev_logability'};
% The correlation is parametrized as corr=tanh(FTcorr) (AH2021 for a 2x2 matrix), so any
% real FTcorr gives a valid correlation. We just keep it from running off: FTcorr in [-3,3] means corr in (-0.995,0.995)
caliboptions.constrainAtoB={'FTcorr_logh1logability'};
caliboptions.constrainAtoBlimits.FTcorr_logh1logability=[-3,3];
% Note: constraints are given by parameter name (a binary vector is silently ignored by the current toolkit)

% Second, we need to say which model statistics we want to target (TargetMoments, set above)
% Note: When setting up TargetMoments there are some rules you must follow
% There are two options TargetMomements.AgeConditionalStats and
% TargetMoments.AllStats (you can use both). Within these you must
% follow the structure that you get when you run the commands
% AgeConditionalStats=LifeCycleProfiles_FHorz_Case1()
% and
% AllStats=EvalFnOnAgentDist_AggVars_FHorz_Case1()
% [If you are using PType, then these will be the PType equivalents]

% So we want a FnsToEvaluate which is just what we need for calibration (including general eqm eqns)
FnsToEvaluate_Calib=FnsToEvaluate_GE; % things we need for general eqm
FnsToEvaluate_Calib.earnings=FnsToEvaluate.earnings; % the calibration targets are about earnings
FnsToEvaluate_Calib.logearnings=FnsToEvaluate.logearnings;

% The objective (HVY2011 Appendix D): sum of squared log ratios of model to data moments
caliboptions.metric='sum_logratiosquared';

%% Set calibrate OLG options
caliboptions.relativeGEweight=50;
caliboptions.jointoptimization=1; % joint-optimization or nested-optimization, default is joint

% How CalibrateOLGModel_PType works:
%      caliboptions.jointoptimization=1 means do a joint-optimization (as opposed to nested), so it will
%  Minimize Obj=Obj1+caliboptions.relativeGEweight*Obj2
%    where Obj1=sum of log-ratio squared (as we set caliboptions.metric, default would be sum of squares) with
%               weights based on caliboptions.weights, which is ones by default)
%          Obj2=sum of squares of the general eqm eqns (weighted using heteroagentoptions.multiGEweights, which
%               is ones by default)

%% doPart(4): Done. Now just run the calibrate OLG model command
% The first calibration, starting from the Table 2 values, stopped at alpha_h=0.71 (objective 0.152), but recalibrating the
% initial distribution at fixed alpha (Fig 6, doPart(7)) found a better fit at alpha_h=0.8 (objective 0.094): the objective is
% U-shaped in alpha with its bottom near 0.8. So the calibration starts from the best alpha of Fig 6 (needs HVY2011_Fig6.mat).
CalibStartFromFig6=1; % =0 starts from the Table 2 values
if doPart(4)==1
    if CalibStartFromFig6==1
        tempFig6=load('./SavedOutput/HVY2011_Fig6.mat','Fig6');
        [~,bestalphaindex]=min(tempFig6.Fig6.objective_targets);
        startfields=fieldnames(tempFig6.Fig6.CalibParams{bestalphaindex}); % includes beta
        for pp=1:length(startfields)
            Params.(startfields{pp})=tempFig6.Fig6.CalibParams{bestalphaindex}.(startfields{pp});
        end
        Params.alpha_h=tempFig6.Fig6.alpha(bestalphaindex);
        Params=PTypeParamFn(Params);
        fprintf('Calibration starts from the Fig 6 calibration at alpha_h=%4.2f (objective %8.4f) \n',Params.alpha_h,tempFig6.Fig6.objective_targets(bestalphaindex))
    end
    tic;
    [CalibParams,calibsummary]=CalibrateOLGModel_PType(CalibParamNames,TargetMoments,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, jequaloneDistFn, AgeWeightsParamNames, GEPriceParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate_Calib, GeneralEqmEqns, heteroagentoptions, caliboptions, vfoptions,simoptions);
    calibtime=toc;
    fprintf('Calibration run time %8.1fs \n',calibtime)
    save ./SavedOutput/HVY2011calib.mat CalibParams calibsummary
elseif doPart(5)==1
    load('./SavedOutput/HVY2011calib.mat','CalibParams','calibsummary')
end

if any(doPart(4:5)==1)
    % From now on, use our calibrated parameter values
    for pp=1:length(CalibParamNames)
        Params.(CalibParamNames{pp})=CalibParams.(CalibParamNames{pp});
    end
    if isfield(CalibParams,'beta')
        Params.beta=CalibParams.beta; % with joint optimization the general eqm price is also returned
    end
    Params=PTypeParamFn(Params); % ability grid for the calibrated mean_logability and stddev_logability
    fprintf('\nCalibrated parameters (HVY2011 Table 2 in brackets): \n')
    fprintf('beta=%8.4f (0.981), alpha_h=%8.4f (0.70) \n',Params.beta,Params.alpha_h)
    fprintf('mean_logh1=%8.4f (4.66), mean_logability=%8.4f (-1.12) \n',Params.mean_logh1,Params.mean_logability)
    fprintf('var_logh1=%8.4f (0.213), var_logability=%8.4f (0.012), cov=%8.4f (0.041) [corr of logs %6.3f (0.811)] \n',Params.stddev_logh1^2,Params.stddev_logability^2,tanh(Params.FTcorr_logh1logability)*Params.stddev_logh1*Params.stddev_logability,tanh(Params.FTcorr_logh1logability))
end


%% doPart(5): Produce some model outputs
if doPart(5)==1
[p_eqm,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDistFn, ReturnFn, FnsToEvaluate_GE, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
fprintf('General eqm at the calibrated parameters: beta=%8.5f, K/L-target=%8.2e \n',Params.beta,GeneralEqmConditions.capitallabor)

[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistFn,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_Check,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
fprintf('Net taxes (per person, cross-section): %8.4f (HVY2011 require this to be positive) \n',AggVars.nettaxes.Mean)

% Fit to the calibration targets at the calibrated parameters (gathered, so the .mat can be read without a GPU)
CalibFit.age=Params.age(1:Nworkages);
CalibFit.earningsMean=gather(AgeConditionalStats.earnings.Mean(1:Nworkages));
CalibFit.logearningsVariance=gather(AgeConditionalStats.logearnings.Variance(1:Nworkages));
CalibFit.earningsRatioMeanToMedian=gather(AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages));
CalibFit.earningsGini=gather(AgeConditionalStats.earnings.Gini(1:Nworkages));
CalibFit.timestudying=gather(AgeConditionalStats.timestudying.Mean(1:Nworkages));
CalibFit.humancapital=gather(AgeConditionalStats.humancapital.Mean(1:Nworkages));
CalibFit.objective1=sum(log(CalibFit.earningsMean./Fig1_Mean_time').^2);
CalibFit.objective2=sum(log(CalibFit.logearningsVariance./Fig1_VarLog_time').^2);
CalibFit.objective3=sum(log(CalibFit.earningsRatioMeanToMedian./Fig1_Skew_time').^2);
CalibFit.objective=CalibFit.objective1+CalibFit.objective2+CalibFit.objective3;
fprintf('Calibration objective at the calibrated parameters (sum of squared log ratios, ages 23-60): %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',CalibFit.objective,CalibFit.objective1,CalibFit.objective2,CalibFit.objective3)
fprintf('    age   mean earnings (model, data)   var log earnings (model, data)   mean/median (model, data)   time learning   human capital \n')
for jj=1:Nworkages
    fprintf('    %3i   %8.2f %8.2f   %8.4f %8.4f   %8.4f %8.4f   %8.4f   %8.2f \n',CalibFit.age(jj),CalibFit.earningsMean(jj),Fig1_Mean_time(jj),CalibFit.logearningsVariance(jj),Fig1_VarLog_time(jj),CalibFit.earningsRatioMeanToMedian(jj),Fig1_Skew_time(jj),CalibFit.timestudying(jj),CalibFit.humancapital(jj))
end

%% Take a look specifically at the target and calibrated model moments
% Figure 2 of HVY2011
figure_c=figure_c+1;
fig=figure(figure_c);
fig.Position=[100 100 900 700];
agemodel=Params.age(1:Nworkages);
subplot(2,2,1); plot(agemodel,AgeConditionalStats.earnings.Mean(1:Nworkages),'-k+',agedata,Fig1_Mean_time,'-r')
xlim([23,60]); ylim([0,140]); xlabel('Age'); title('Panel A: Mean earnings'); legend({'Model','Data-time'},'Location','southeast')
subplot(2,2,2); plot(agemodel,AgeConditionalStats.logearnings.Variance(1:Nworkages),'-k+',agedata,Fig1_VarLog_time,'-r')
xlim([23,60]); ylim([0,0.9]); xlabel('Age'); title('Panel B: Variance of log earnings')
subplot(2,2,3); plot(agemodel,AgeConditionalStats.earnings.Gini(1:Nworkages),'-k+',agedata,Fig1_Gini_time,'-r')
xlim([23,60]); ylim([0,0.6]); xlabel('Age'); title('Panel C: Earnings Gini')
subplot(2,2,4); plot(agemodel,AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages),'-k+',agedata,Fig1_Skew_time,'-r')
xlim([23,60]); ylim([0.8,1.6]); xlabel('Age'); title('Panel D: Earnings skewness (mean/median)')
saveas(fig,'./SavedOutput/Graphs/HVY2011_Fig2.png')

% Figure 3 of HVY2011
figure_c=figure_c+1;
fig=figure(figure_c);
subplot(1,2,1); plot(agemodel,AgeConditionalStats.timestudying.Mean(1:Nworkages))
xlim([23,60]); xlabel('Age'); title('Mean time in human capital')
subplot(1,2,2); plot(agemodel,AgeConditionalStats.humancapital.Mean(1:Nworkages))
xlim([23,60]); xlabel('Age'); title('Mean human capital profile')
saveas(fig,'./SavedOutput/Graphs/HVY2011_Fig3.png')

save ./SavedOutput/HVY2011full.mat Params CalibParams calibsummary p_eqm GeneralEqmConditions AgeConditionalStats CalibFit
end


%% doPart(6): HVY2011 results (Figs 2-5, Tables 3, 5, 6) at both parameter sets: our calibration, and HVY2011 Table 2
% Everything here holds prices fixed at those of the parameter set (r=0.042 and the rental rate are fixed throughout,
% and beta is the general eqm beta of the parameter set), as HVY2011 do for their counterfactuals.
if doPart(6)==1
tempGE=load('./SavedOutput/HVY2011_GE.mat','p_eqm'); % general eqm beta at the Table 2 parameters (doPart(3))
tempCalib=load('./SavedOutput/HVY2011calib.mat','CalibParams'); % our calibration (doPart(4)), includes its general eqm beta
ParamSetNames={'calibrated','table2'};
ParamSetLabels={'our calibration','HVY2011 Table 2'};
[~,zeroassetindex]=min(abs(asset_grid));
u_gridc=gather(u_grid);
pi_uc=gather(pi_u);
EarningsReturnFn=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr); % for the expected present value of lifetime earnings
Results=struct();
for ss=1:2
    setname=ParamSetNames{ss};
    ParamsSet=Params_Table2;
    if ss==1
        calibfields=fieldnames(tempCalib.CalibParams);
        for pp=1:length(calibfields)
            ParamsSet.(calibfields{pp})=tempCalib.CalibParams.(calibfields{pp});
        end
    else
        ParamsSet.beta=tempGE.p_eqm.beta;
    end
    ParamsSet=PTypeParamFn(ParamsSet); % ability grid
    ParamsSet.discountearnings=1/(1+ParamsSet.r); % lifetime earnings are discounted at the interest rate (HVY2011 footnote 24)
    % The initial distribution, and the ptype masses implied by it (the toolkit uses these masses anyway, as jequaloneDist has ptype as a dimension; set abilitydist to match)
    jequaloneDistSet=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,ParamsSet.mean_logh1,ParamsSet.mean_logability,ParamsSet.stddev_logh1,ParamsSet.stddev_logability,ParamsSet.FTcorr_logh1logability,ParamsSet.mu_z,ParamsSet.sigma_z,ParamsSet.ability,N_i);
    ParamsSet.abilitydist=gather(reshape(sum(sum(jequaloneDistSet,1),2),[N_i,1]));
    var_h=ParamsSet.stddev_logh1^2;
    var_a=ParamsSet.stddev_logability^2;
    corr_ha=tanh(ParamsSet.FTcorr_logh1logability);
    cov_ha=corr_ha*ParamsSet.stddev_logh1*ParamsSet.stddev_logability;
    fprintf('\n=== HVY2011 results at the parameters of %s: beta=%8.5f, alpha_h=%6.4f, mean_logh1=%6.4f, mean_logability=%7.4f, var_logh1=%6.4f, var_logability=%6.4f, cov=%7.4f (corr %6.3f) \n',ParamSetLabels{ss},ParamsSet.beta,ParamsSet.alpha_h,ParamsSet.mean_logh1,ParamsSet.mean_logability,var_h,var_a,cov_ha,corr_ha)
    Results.(setname).Params=ParamsSet;

    %% Benchmark model: value fn, agent dist, life-cycle profiles
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsSet, DiscountFactorParamNames, vfoptions);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistSet,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],ParamsSet,simoptions);
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,ParamsSet,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    Results.(setname).age=Params.age(1:Nworkages);
    Results.(setname).earningsMean=gather(AgeConditionalStats.earnings.Mean(1:Nworkages));
    Results.(setname).logearningsVariance=gather(AgeConditionalStats.logearnings.Variance(1:Nworkages));
    Results.(setname).earningsGini=gather(AgeConditionalStats.earnings.Gini(1:Nworkages));
    Results.(setname).earningsRatioMeanToMedian=gather(AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages));
    Results.(setname).timestudying=gather(AgeConditionalStats.timestudying.Mean(1:Nworkages));
    Results.(setname).humancapital=gather(AgeConditionalStats.humancapital.Mean(1:Nworkages));
    Results.(setname).logconsumptionVariance=gather(AgeConditionalStats.logconsumption.Variance(1:Nworkages));
    Results.(setname).objective1=sum(log(Results.(setname).earningsMean./Fig1_Mean_time').^2);
    Results.(setname).objective2=sum(log(Results.(setname).logearningsVariance./Fig1_VarLog_time').^2);
    Results.(setname).objective3=sum(log(Results.(setname).earningsRatioMeanToMedian./Fig1_Skew_time').^2);
    Results.(setname).objective=Results.(setname).objective1+Results.(setname).objective2+Results.(setname).objective3;
    fprintf('Calibration objective: %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',Results.(setname).objective,Results.(setname).objective1,Results.(setname).objective2,Results.(setname).objective3)
    fprintf('Mean time learning at age 23: %6.3f (HVY2011: approximately 0.25), and after age 55 at most %6.3f (HVY2011: less than 0.05) \n',Results.(setname).timestudying(1),max(Results.(setname).timestudying(34:end)))

    %% Table 3: properties of the initial distribution (bivariate log-normal; h1 is initial human capital before the age-23 shock)
    Table3.meanability=exp(ParamsSet.mean_logability+var_a/2);
    Table3.cvability=sqrt(exp(var_a)-1);
    Table3.meanh1=exp(ParamsSet.mean_logh1+var_h/2);
    Table3.cvh1=sqrt(exp(var_h)-1);
    Table3.corr=(exp(cov_ha)-1)/sqrt((exp(var_h)-1)*(exp(var_a)-1));
    Results.(setname).Table3=Table3;
    fprintf('Table 3 (HVY2011 benchmark in brackets): mean ability %6.3f (0.329), CV ability %6.3f (0.112), mean h1 %6.1f (116.9), CV h1 %6.3f (0.487), corr(a,h1) %6.3f (0.746) \n',Table3.meanability,Table3.cvability,Table3.meanh1,Table3.cvh1,Table3.corr)

    %% Table 5: fraction of the variance of lifetime utility and lifetime earnings due to initial conditions
    % Var(X)=Var(E[X|x1])+E[Var(X|x1)], where the initial condition x1=(h1,a) (zero initial assets). Everything is exact, no simulation.
    % E[X|x1]: for lifetime utility it is the value fn, and for lifetime earnings (discounted at r) it is the value fn of the policy
    % with earnings as the return fn. As in the authors' code, h1 is hit by the age-23 shock z1, so E[X|h1,a] averages X at age 23
    % over the five values of z1.
    % Var(X): X=sum_j d^(j-1) f_j, and W_j(s)=E[sum_{k>=j} d^(k-j) f_k | s_j=s]. The martingale M_j=E[X|age j info] has increments
    % M_j-M_{j-1}=d^(j-1)*(W_j(s_j)-E[W_j|s_{j-1}]), with E[W_j|s_{j-1}]=(W_{j-1}(s_{j-1})-f_{j-1}(s_{j-1}))/d, so
    %   Var(X)=Var_1(W_1) + sum_{j=2}^J d^(2(j-1)) * ( E_j[W_j^2] - E_{j-1}[((W_{j-1}-f_{j-1})/d)^2] )
    % where E_j is over the agent distribution at age j, and f_j (period earnings, period utility) is evaluated on the grid.
    % (E_j[W_j^2] is a moment of the value fn, which is not a function of (d,aprime,a), so it is summed over the agent distribution directly)
    [Wearn,~]=ValueFnFromPolicy_FHorz_PType(Policy,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,ParamsSet,{'discountearnings'},vfoptions);
    FnsToEvaluate_f=struct();
    FnsToEvaluate_f.earnings=FnsToEvaluate.earnings;
    FnsToEvaluate_f.utility=@(l,aprime,a,h,r,gamma,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) (HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)^(1-gamma))/(1-gamma);
    PeriodValues=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsToEvaluate_f,ParamsSet,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
    P1=gather(MVNormal_ProbabilitiesOnGrid([log(h_grid);log(ParamsSet.ability)],[ParamsSet.mean_logh1,ParamsSet.mean_logability],[var_h,cov_ha;cov_ha,var_a],[n_a(2),N_i])); % distribution of (h1,a) on the grids
    P1=P1/sum(P1(:));
    EU_x1=zeros(n_a(2),N_i);
    Eearn_x1=zeros(n_a(2),N_i);
    for ii=1:N_i
        Vii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        V1=reshape(Vii(zeroassetindex,:,1),[n_a(2),1]);
        W1=reshape(Wii(zeroassetindex,:,1),[n_a(2),1]);
        for zz=1:n_u
            hz=min(max(exp(u_gridc(zz))*h_grid,h_grid(1)),h_grid(end));
            EU_x1(:,ii)=EU_x1(:,ii)+pi_uc(zz)*interp1(h_grid,V1,hz);
            Eearn_x1(:,ii)=Eearn_x1(:,ii)+pi_uc(zz)*interp1(h_grid,W1,hz);
        end
    end
    meanU_exact=sum(P1(:).*EU_x1(:));
    meanEarn_exact=sum(P1(:).*Eearn_x1(:));
    varU_initial=sum(P1(:).*(EU_x1(:)-meanU_exact).^2);
    varEarn_initial=sum(P1(:).*(Eearn_x1(:)-meanEarn_exact).^2);
    % Moments by age over the agent distribution: E_j[W_j], E_j[W_j^2], E_j[(W_j-f_j)/d], E_j[((W_j-f_j)/d)^2]
    % (sums are over states with positive mass only, as W is -Inf at states nobody reaches)
    momEarn=zeros(4,N_j);
    momU=zeros(4,N_j);
    for ii=1:N_i
        Dii=reshape(gather(StationaryDist.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Dii=StationaryDist.ptweights(ii)*Dii./max(sum(Dii,1),eps); % distribution at each age, times the ptype mass
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Fii=reshape(gather(PeriodValues.earnings.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Uii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Gii=reshape(gather(PeriodValues.utility.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Wii(Dii==0)=0; Fii(Dii==0)=0; Uii(Dii==0)=0; Gii(Dii==0)=0;
        Cii=(Wii-Fii)/ParamsSet.discountearnings; % E[W_{j+1}|s_j]
        CUii=(Uii-Gii)/ParamsSet.beta;
        momEarn=momEarn+[sum(Dii.*Wii,1); sum(Dii.*Wii.^2,1); sum(Dii.*Cii,1); sum(Dii.*Cii.^2,1)];
        momU=momU+[sum(Dii.*Uii,1); sum(Dii.*Uii.^2,1); sum(Dii.*CUii,1); sum(Dii.*CUii.^2,1)];
    end
    discE=ParamsSet.discountearnings.^(2*(1:1:N_j-1)); % d^(2(j-1)) for j=2,...,J
    discU=ParamsSet.beta.^(2*(1:1:N_j-1));
    varEarn_total=(momEarn(2,1)-momEarn(1,1)^2)+sum(discE.*(momEarn(2,2:end)-momEarn(4,1:end-1)));
    varU_total=(momU(2,1)-momU(1,1)^2)+sum(discU.*(momU(2,2:end)-momU(4,1:end-1)));
    Table5.utility=varU_initial/varU_total;
    Table5.earnings=varEarn_initial/varEarn_total;
    Table5.wealth=Table5.earnings; % lifetime wealth is lifetime earnings plus initial wealth, which is zero in the benchmark model
    Table5.utility_inclz1=(momU(2,1)-momU(1,1)^2)/varU_total; % if the age-23 shock z1 is counted as an initial condition
    Table5.earnings_inclz1=(momEarn(2,1)-momEarn(1,1)^2)/varEarn_total;
    Results.(setname).Table5=Table5;
    fprintf('Table 5, fraction of variance due to initial conditions (HVY2011 benchmark in brackets): lifetime utility %6.3f (0.640), lifetime earnings %6.3f (0.615), lifetime wealth %6.3f (0.615) \n',Table5.utility,Table5.earnings,Table5.wealth)
    fprintf('   counting the age-23 shock as an initial condition: lifetime utility %6.3f, lifetime earnings %6.3f \n',Table5.utility_inclz1,Table5.earnings_inclz1)
    fprintf('   checks: E_j[W_j]-E_{j-1}[E[W_j|s_{j-1}]] should be zero, max abs over ages: earnings %8.2e (relative to mean %8.2e), utility %8.2e (relative to mean %8.2e) \n',max(abs(momEarn(1,2:end)-momEarn(3,1:end-1))),max(abs(momEarn(1,2:end)-momEarn(3,1:end-1))./abs(momEarn(1,2:end))),max(abs(momU(1,2:end)-momU(3,1:end-1))),max(abs(momU(1,2:end)-momU(3,1:end-1))./abs(momU(1,2:end))))
    fprintf('          each age''s variance increment should be non-negative, min: earnings %8.2e, utility %8.2e \n',min(momEarn(2,2:end)-momEarn(4,1:end-1)),min(momU(2,2:end)-momU(4,1:end-1)))
    fprintf('          mean of lifetime earnings and utility: over (h1,a) %8.2f and %10.5f, over the age-23 distribution %8.2f and %10.5f (should be close) \n',meanEarn_exact,meanU_exact,momEarn(1,1),momU(1,1))
    clear PeriodValues

    %% Table 6: equivalent variation and change in expected lifetime wealth from a one standard deviation change in an initial condition
    % Baseline: median h1 and median ability (exp of the means of the logs), zero assets. Ability is a ptype, so solve the
    % model for three ability levels: median, and plus/minus one standard deviation of log ability.
    Params6=ParamsSet;
    Params6.ability=exp(ParamsSet.mean_logability+[-1;0;1]*ParamsSet.stddev_logability);
    Params6.abilitydist=ones(3,1)/3;
    [V6, Policy6]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,3, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params6, DiscountFactorParamNames, vfoptions);
    [W6,~]=ValueFnFromPolicy_FHorz_PType(Policy6,n_d,n_a,n_z,N_j,3,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,Params6,{'discountearnings'},vfoptions);
    h1values=exp(ParamsSet.mean_logh1+[-1,0,1]*ParamsSet.stddev_logh1);
    EV6=zeros(3,3); % (h1 value, ability value): expected lifetime utility
    EW6=zeros(3,3); % (h1 value, ability value): expected lifetime earnings
    for kk=1:3
        Vkk=reshape(gather(V6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        Wkk=reshape(gather(W6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        V1=reshape(Vkk(zeroassetindex,:,1),[n_a(2),1]);
        W1=reshape(Wkk(zeroassetindex,:,1),[n_a(2),1]);
        for zz=1:n_u
            EV6(:,kk)=EV6(:,kk)+pi_uc(zz)*interp1(h_grid,V1,exp(u_gridc(zz))*h1values');
            EW6(:,kk)=EW6(:,kk)+pi_uc(zz)*interp1(h_grid,W1,exp(u_gridc(zz))*h1values');
        end
    end
    % Equivalent variation: the percent change in consumption at all ages that gives the baseline agent the same expected
    % utility (with CES utility, scaling consumption by lambda scales expected utility by lambda^(1-gamma))
    Table6.EV=100*([EV6(3,2);EV6(1,2);EV6(2,3);EV6(2,1)]/EV6(2,2)).^(1/(1-ParamsSet.gamma))-100; % h1 +1sd, h1 -1sd, ability +1sd, ability -1sd
    Table6.lifetimewealth=100*[EW6(3,2);EW6(1,2);EW6(2,3);EW6(2,1)]/EW6(2,2)-100;
    Results.(setname).Table6=Table6;
    fprintf('Table 6 (HVY2011 benchmark in brackets): equivalent variation (percent): h1 +1sd %6.1f (39.3), h1 -1sd %6.1f (-28.3), ability +1sd %6.1f (5.7), ability -1sd %6.1f (-2.6) \n',Table6.EV)
    fprintf('                                         change in lifetime wealth (percent): h1 +1sd %6.1f (47.5), h1 -1sd %6.1f (-31.7), ability +1sd %6.1f (8.1), ability -1sd %6.1f (-3.9) \n',Table6.lifetimewealth)
    clear V6 Policy6 W6 Wearn

    %% Figure 4: earnings dispersion without ability differences, and without shocks (prices held fixed)
    % (A) every agent has the median ability, initial human capital unchanged
    ParamsA=ParamsSet;
    ParamsA.ability=exp(ParamsSet.mean_logability)*ones(N_i,1);
    [~, PolicyA]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsA, DiscountFactorParamNames, vfoptions);
    StationaryDistA=StationaryDist_Case1_FHorz_PType(jequaloneDistSet,AgeWeightsParamNames,PTypeDistParamNames,PolicyA,n_d,n_a,n_z,N_j,N_i,[],ParamsA,simoptions);
    AgeConditionalStatsA=LifeCycleProfiles_FHorz_Case1_PType(StationaryDistA,PolicyA,FnsToEvaluate,ParamsA,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    Results.(setname).noability_logearningsVariance=gather(AgeConditionalStatsA.logearnings.Variance(1:Nworkages));
    clear PolicyA StationaryDistA AgeConditionalStatsA
    % (B) no shocks: sigma_z=0, with mu_z adjusted to keep the mean shock exp(mu_z+sigma_z^2/2) unchanged; initial (h1,a) unchanged
    ParamsB=ParamsSet;
    ParamsB.mu_z=ParamsSet.mu_z+ParamsSet.sigma_z^2/2;
    ParamsB.sigma_z=0;
    vfoptionsB=vfoptions;
    vfoptionsB.u_grid=ParamsB.mu_z*ones(n_u,1); % every shock is the mean shock
    simoptionsB=simoptions;
    simoptionsB.u_grid=vfoptionsB.u_grid;
    jequaloneDistB=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,ParamsB.mean_logh1,ParamsB.mean_logability,ParamsB.stddev_logh1,ParamsB.stddev_logability,ParamsB.FTcorr_logh1logability,ParamsB.mu_z,ParamsB.sigma_z,ParamsB.ability,N_i);
    ParamsB.abilitydist=gather(reshape(sum(sum(jequaloneDistB,1),2),[N_i,1]));
    [~, PolicyB]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsB, DiscountFactorParamNames, vfoptionsB);
    StationaryDistB=StationaryDist_Case1_FHorz_PType(jequaloneDistB,AgeWeightsParamNames,PTypeDistParamNames,PolicyB,n_d,n_a,n_z,N_j,N_i,[],ParamsB,simoptionsB);
    AgeConditionalStatsB=LifeCycleProfiles_FHorz_Case1_PType(StationaryDistB,PolicyB,FnsToEvaluate,ParamsB,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptionsB);
    Results.(setname).noshocks_logearningsVariance=gather(AgeConditionalStatsB.logearnings.Variance(1:Nworkages));
    Results.(setname).noshocks_earningsMean=gather(AgeConditionalStatsB.earnings.Mean(1:Nworkages));
    clear PolicyB StationaryDistB AgeConditionalStatsB

    %% Figures
    if CreateFigures==1
        agemodel=Results.(setname).age;
        % Figure 2: model and data (time effects)
        figure_c=figure_c+1;
        fig=figure(figure_c);
        fig.Position=[100 100 900 700];
        subplot(2,2,1); plot(agemodel,Results.(setname).earningsMean,'-k+',agedata,Fig1_Mean_time,'-r')
        xlim([23,60]); ylim([0,140]); xlabel('Age'); title('Panel A: Mean earnings'); legend({'Model','Data-time'},'Location','southeast')
        subplot(2,2,2); plot(agemodel,Results.(setname).logearningsVariance,'-k+',agedata,Fig1_VarLog_time,'-r')
        xlim([23,60]); ylim([0,0.9]); xlabel('Age'); title('Panel B: Variance of log earnings')
        subplot(2,2,3); plot(agemodel,Results.(setname).earningsGini,'-k+',agedata,Fig1_Gini_time,'-r')
        xlim([23,60]); ylim([0,0.6]); xlabel('Age'); title('Panel C: Earnings Gini')
        subplot(2,2,4); plot(agemodel,Results.(setname).earningsRatioMeanToMedian,'-k+',agedata,Fig1_Skew_time,'-r')
        xlim([23,60]); ylim([0.8,1.6]); xlabel('Age'); title('Panel D: Earnings skewness (mean/median)')
        saveas(fig,['./SavedOutput/Graphs/HVY2011_Fig2_',setname,'.png'])
        % Figure 3: human capital by age
        figure_c=figure_c+1;
        fig=figure(figure_c);
        fig.Position=[100 100 900 380];
        subplot(1,2,1); plot(agemodel,Results.(setname).timestudying,'-r','LineWidth',1)
        xlim([23,60]); ylim([0,0.5]); xlabel('Age'); title('Panel A: Mean time in human capital')
        subplot(1,2,2); plot(agemodel,Results.(setname).humancapital,'-r','LineWidth',1)
        xlim([23,60]); ylim([0,140]); xlabel('Age'); title('Panel B: Mean human capital profile')
        saveas(fig,['./SavedOutput/Graphs/HVY2011_Fig3_',setname,'.png'])
        % Figure 4: dispersion of earnings by age
        figure_c=figure_c+1;
        fig=figure(figure_c);
        fig.Position=[100 100 900 380];
        subplot(1,2,1); plot(agemodel,Results.(setname).noability_logearningsVariance,'-k',agemodel,Results.(setname).logearningsVariance,'--b','LineWidth',1)
        xlim([23,60]); ylim([0,0.9]); xlabel('Age'); ylabel('Var(log earnings)'); title('Panel A: Eliminating ability differences')
        legend({'No ability differences','Benchmark: time'},'Location','southeast')
        subplot(1,2,2); plot(agemodel,Results.(setname).noshocks_logearningsVariance,'-k',agemodel,Results.(setname).logearningsVariance,'--b','LineWidth',1)
        xlim([23,60]); ylim([0,0.9]); xlabel('Age'); ylabel('Var(log earnings)'); title('Panel B: Eliminating shocks')
        legend({'No shocks','Benchmark: time'},'Location','northeast')
        saveas(fig,['./SavedOutput/Graphs/HVY2011_Fig4_',setname,'.png'])
        % Figure 5: rise in consumption dispersion from age 25 (the data series in HVY2011 Fig 5 are from other papers and are not reproduced)
        figure_c=figure_c+1;
        fig=figure(figure_c);
        fig.Position=[100 100 640 480];
        age25index=find(agemodel==25);
        plot(agemodel(age25index:end),Results.(setname).logconsumptionVariance(age25index:end)-Results.(setname).logconsumptionVariance(age25index),'-k')
        xlim([25,60]); ylim([-0.05,0.3]); xlabel('Age'); title('Rise in consumption dispersion: model (time effects)')
        legend({'Model-time effects'},'Location','northwest')
        saveas(fig,['./SavedOutput/Graphs/HVY2011_Fig5_',setname,'.png'])
    end
    fprintf('Figure 5: rise in the variance of log consumption from age 25 to 60: %6.3f (HVY2011: between Primiceri and van Rens (2009) and Aguiar and Hurst (2008), less than Deaton and Paxson (1994)) \n',Results.(setname).logconsumptionVariance(end)-Results.(setname).logconsumptionVariance(find(Results.(setname).age==25)))

    save ./SavedOutput/HVY2011_Results.mat Results % save after each parameter set
    clear V Policy StationaryDist AgeConditionalStats
end
end


%% doPart(7): HVY2011 Figure 6, lifetime inequality and the elasticity parameter alpha
% For each alpha, set the distribution of initial conditions to best match the earnings facts (beta from the general eqm, as in
% doPart(4)), then compute the fraction of the variance of lifetime earnings due to initial conditions (as Table 5 in doPart(6)).
% Our calibration (doPart(4), where alpha is also chosen) gives the highlighted point, from doPart(6).
if doPart(7)==1
AlphaList=[0.5,0.6,0.8,0.9];
tempCalib=load('./SavedOutput/HVY2011calib.mat','CalibParams'); % starting point for each alpha: our calibration
tempResults=load('./SavedOutput/HVY2011_Results.mat','Results'); % the highlighted point (doPart(6))
[~,zeroassetindex]=min(abs(asset_grid));
u_gridc=gather(u_grid);
pi_uc=gather(pi_u);
EarningsReturnFn=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr);
CalibParamNames7={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTcorr_logh1logability'}; % alpha_h is fixed
caliboptions7=caliboptions;
caliboptions7.constrain0to1={};
AlphaToRun=[0,0,0,0]; % which alphas to calibrate in this run; the others use their calibration saved in HVY2011_Fig6.mat by an earlier run
if all(AlphaToRun==1)
    Fig6=struct();
    Fig6.alpha=AlphaList;
    Fig6.earnings=nan(1,length(AlphaList));
    Fig6.utility=nan(1,length(AlphaList));
    Fig6.objective=nan(1,length(AlphaList));
    Fig6.CalibParams=cell(1,length(AlphaList));
else
    load('./SavedOutput/HVY2011_Fig6.mat','Fig6')
end
Fig6.objective_targets=nan(1,length(AlphaList)); % the calibration objective itself (sum of squared log ratios; the toolkit's objvalue also scales it and adds the general eqm term)
for aa=1:length(AlphaList)
    Params7=Params_Table2;
    if AlphaToRun(aa)==1
        % Initial guess: the calibration at the nearest alpha already done (our calibration, or an earlier alpha), with the mean of
        % log ability shifted so that ability*(h*s)^alpha is unchanged at a typical h*s of 36 (otherwise, at a high alpha the initial
        % guess has so much human capital production that median earnings at some age is zero, and the objective is infinite)
        startalpha=tempCalib.CalibParams.alpha_h;
        startparams=tempCalib.CalibParams;
        for bb=1:length(AlphaList)
            if ~isempty(Fig6.CalibParams{bb}) && abs(AlphaList(bb)-AlphaList(aa))<abs(startalpha-AlphaList(aa))
                startalpha=AlphaList(bb);
                startparams=Fig6.CalibParams{bb};
            end
        end
        calibfields=fieldnames(startparams);
        for pp=1:length(calibfields)
            Params7.(calibfields{pp})=startparams.(calibfields{pp});
        end
        Params7.alpha_h=AlphaList(aa);
        Params7.mean_logability=startparams.mean_logability-(AlphaList(aa)-startalpha)*log(36);
        Params7=PTypeParamFn(Params7);
        fprintf('\n=== Figure 6: calibrating the initial distribution with alpha_h=%4.2f (initial guess from alpha_h=%4.2f, mean_logability %7.4f) \n',Params7.alpha_h,startalpha,Params7.mean_logability)
        tic;
        [CalibParams7,calibsummary7]=CalibrateOLGModel_PType(CalibParamNames7,TargetMoments,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params7, DiscountFactorParamNames, jequaloneDistFn, AgeWeightsParamNames, GEPriceParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate_Calib, GeneralEqmEqns, heteroagentoptions, caliboptions7, vfoptions,simoptions);
        calibtime7=toc;
        Fig6.CalibParams{aa}=CalibParams7;
        Fig6.objective(aa)=calibsummary7.objvalue;
        fprintf('alpha_h=%4.2f: calibration run time %8.1fs \n',AlphaList(aa),calibtime7)
    else
        CalibParams7=Fig6.CalibParams{aa};
        fprintf('\n=== Figure 6: alpha_h=%4.2f, using the calibration from an earlier run \n',AlphaList(aa))
    end
    calibfields=fieldnames(CalibParams7); % includes beta (joint optimization)
    for pp=1:length(calibfields)
        Params7.(calibfields{pp})=CalibParams7.(calibfields{pp});
    end
    Params7.alpha_h=AlphaList(aa);
    Params7=PTypeParamFn(Params7);
    Params7.discountearnings=1/(1+Params7.r);
    var_h=Params7.stddev_logh1^2;
    var_a=Params7.stddev_logability^2;
    corr_ha=tanh(Params7.FTcorr_logh1logability);
    cov_ha=corr_ha*Params7.stddev_logh1*Params7.stddev_logability;
    fprintf('alpha_h=%4.2f: beta=%8.5f, mean_logh1=%6.4f, mean_logability=%7.4f, var_logh1=%6.4f, var_logability=%6.4f, corr=%6.3f \n',Params7.alpha_h,Params7.beta,Params7.mean_logh1,Params7.mean_logability,var_h,var_a,corr_ha)

    % Fraction of the variance of lifetime earnings (and utility) due to initial conditions, exactly as for Table 5 in doPart(6)
    jequaloneDist7=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,Params7.mean_logh1,Params7.mean_logability,Params7.stddev_logh1,Params7.stddev_logability,Params7.FTcorr_logh1logability,Params7.mu_z,Params7.sigma_z,Params7.ability,N_i);
    Params7.abilitydist=gather(reshape(sum(sum(jequaloneDist7,1),2),[N_i,1]));
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params7, DiscountFactorParamNames, vfoptions);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist7,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params7,simoptions);
    % The calibration objective (as in doPart(4)), so the fit can be compared across alpha
    AgeConditionalStats7=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_Calib,Params7,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    objective1=sum(log(gather(AgeConditionalStats7.earnings.Mean(1:Nworkages))./Fig1_Mean_time').^2);
    objective2=sum(log(gather(AgeConditionalStats7.logearnings.Variance(1:Nworkages))./Fig1_VarLog_time').^2);
    objective3=sum(log(gather(AgeConditionalStats7.earnings.RatioMeanToMedian(1:Nworkages))./Fig1_Skew_time').^2);
    Fig6.objective_targets(aa)=objective1+objective2+objective3;
    fprintf('alpha_h=%4.2f: calibration objective %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',AlphaList(aa),Fig6.objective_targets(aa),objective1,objective2,objective3)
    clear AgeConditionalStats7
    [Wearn,~]=ValueFnFromPolicy_FHorz_PType(Policy,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,Params7,{'discountearnings'},vfoptions);
    FnsToEvaluate_f=struct();
    FnsToEvaluate_f.earnings=FnsToEvaluate.earnings;
    FnsToEvaluate_f.utility=@(l,aprime,a,h,r,gamma,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) (HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)^(1-gamma))/(1-gamma);
    PeriodValues=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsToEvaluate_f,Params7,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
    % Var(E[X|x1]), x1=(h1,a) before the age-23 shock
    P1=gather(MVNormal_ProbabilitiesOnGrid([log(h_grid);log(Params7.ability)],[Params7.mean_logh1,Params7.mean_logability],[var_h,cov_ha;cov_ha,var_a],[n_a(2),N_i]));
    P1=P1/sum(P1(:));
    EU_x1=zeros(n_a(2),N_i);
    Eearn_x1=zeros(n_a(2),N_i);
    for ii=1:N_i
        Vii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        V1=reshape(Vii(zeroassetindex,:,1),[n_a(2),1]);
        W1=reshape(Wii(zeroassetindex,:,1),[n_a(2),1]);
        for zz=1:n_u
            hz=min(max(exp(u_gridc(zz))*h_grid,h_grid(1)),h_grid(end));
            EU_x1(:,ii)=EU_x1(:,ii)+pi_uc(zz)*interp1(h_grid,V1,hz);
            Eearn_x1(:,ii)=Eearn_x1(:,ii)+pi_uc(zz)*interp1(h_grid,W1,hz);
        end
    end
    varU_initial=sum(P1(:).*(EU_x1(:)-sum(P1(:).*EU_x1(:))).^2);
    varEarn_initial=sum(P1(:).*(Eearn_x1(:)-sum(P1(:).*Eearn_x1(:))).^2);
    % Var(X), by the martingale decomposition over ages
    momEarn=zeros(4,N_j);
    momU=zeros(4,N_j);
    for ii=1:N_i
        Dii=reshape(gather(StationaryDist.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Dii=StationaryDist.ptweights(ii)*Dii./max(sum(Dii,1),eps);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Fii=reshape(gather(PeriodValues.earnings.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Uii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Gii=reshape(gather(PeriodValues.utility.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Wii(Dii==0)=0; Fii(Dii==0)=0; Uii(Dii==0)=0; Gii(Dii==0)=0;
        Cii=(Wii-Fii)/Params7.discountearnings;
        CUii=(Uii-Gii)/Params7.beta;
        momEarn=momEarn+[sum(Dii.*Wii,1); sum(Dii.*Wii.^2,1); sum(Dii.*Cii,1); sum(Dii.*Cii.^2,1)];
        momU=momU+[sum(Dii.*Uii,1); sum(Dii.*Uii.^2,1); sum(Dii.*CUii,1); sum(Dii.*CUii.^2,1)];
    end
    varEarn_total=(momEarn(2,1)-momEarn(1,1)^2)+sum(Params7.discountearnings.^(2*(1:1:N_j-1)).*(momEarn(2,2:end)-momEarn(4,1:end-1)));
    varU_total=(momU(2,1)-momU(1,1)^2)+sum(Params7.beta.^(2*(1:1:N_j-1)).*(momU(2,2:end)-momU(4,1:end-1)));
    Fig6.earnings(aa)=varEarn_initial/varEarn_total;
    Fig6.utility(aa)=varU_initial/varU_total;
    fprintf('alpha_h=%4.2f: fraction of variance due to initial conditions: lifetime earnings %6.3f, lifetime utility %6.3f \n',Params7.alpha_h,Fig6.earnings(aa),Fig6.utility(aa))
    fprintf('   check: E_j[W_j]-E_{j-1}[E[W_j|s_{j-1}]], max relative over ages: earnings %8.2e, utility %8.2e \n',max(abs(momEarn(1,2:end)-momEarn(3,1:end-1))./abs(momEarn(1,2:end))),max(abs(momU(1,2:end)-momU(3,1:end-1))./abs(momU(1,2:end))))
    save ./SavedOutput/HVY2011_Fig6.mat Fig6 % save after each alpha
    clear V Policy StationaryDist Wearn PeriodValues
end
% The highlighted point: our calibration
Fig6.alpha_calibrated=tempResults.Results.calibrated.Params.alpha_h;
Fig6.earnings_calibrated=tempResults.Results.calibrated.Table5.earnings;
Fig6.utility_calibrated=tempResults.Results.calibrated.Table5.utility;
save ./SavedOutput/HVY2011_Fig6.mat Fig6
Fig6.objective_calibrated=tempResults.Results.calibrated.objective;
save ./SavedOutput/HVY2011_Fig6.mat Fig6
fprintf('\nFigure 6: fraction of the variance of lifetime earnings due to initial conditions, and the calibration objective, by alpha: \n')
alphaall=[Fig6.alpha,Fig6.alpha_calibrated];
earningsall=[Fig6.earnings,Fig6.earnings_calibrated];
objectiveall=[Fig6.objective_targets,Fig6.objective_calibrated];
[alphaall,sortindex]=sort(alphaall);
earningsall=earningsall(sortindex);
objectiveall=objectiveall(sortindex);
for aa=1:length(alphaall)
    fprintf('   alpha=%5.3f: lifetime earnings %6.3f, calibration objective %8.4f \n',alphaall(aa),earningsall(aa),objectiveall(aa))
end
if CreateFigures==1
    figure_c=figure_c+1;
    fig=figure(figure_c);
    fig.Position=[100 100 640 480];
    plot(alphaall,earningsall,':b','LineWidth',1.5); hold on;
    plot(Fig6.alpha_calibrated,Fig6.earnings_calibrated,'ob','MarkerSize',8); hold off;
    xlim([0.5,0.9]); ylim([0,1]); xlabel('Elasticity parameter, \alpha');
    legend({'Lifetime earnings'},'Location','southeast')
    title('Lifetime inequality and the elasticity parameter \alpha')
    saveas(fig,'./SavedOutput/Graphs/HVY2011_Fig6.png')
end
end


%% doPart(9): HVY2011 Table 1, estimation of the human capital shocks (data work, no GPU needed)
% (comes before doPart(8) in the code, as doPart(8) writes its LaTeX table)
if doPart(9)==1
    HVY2011_Table1
end


%% doPart(10): HVY2011 Table 4, statistical models of earnings (RIP and HIP) estimated on data simulated from the model
% HVY2011 Section IV.D and footnote 21, following the authors' code (data/hvy-programfiles/gmm-estimation.f, gmm-estimation-hip.f).
% For each parameter set (our calibration, and HVY2011 Table 2), simulate 500 panels of 200 agents (SimPanelValues), each followed
% over ages 23-60 [the authors instead draw 200 agents at random for each age from one simulation of 40,000 agents]. For each panel:
% measurement error N(0,0.15^2) is added to log earnings (each earnings value is multiplied by exp(N(0,0.15^2)) and floored at
% 0.0001, as in the authors' code); the age means are removed; and the variances (ages 26-60) and autocovariances (lags 1 to LAG)
% of log earnings are computed, dividing by N. The RIP and HIP models are estimated by minimizing the equally weighted squared
% distance to these moments, with the authors' parametrization, bounds and initial guesses (fminsearch from their first simplex
% vertex; they use AMOEBA). LAG=1 is 'cov=1' and LAG=5 is 'cov=5'. Table 4 is the mean and standard deviation across the 500 panels.
% Also the autocorrelations of earnings growth rates (HVY2011 Section IV.D; the authors' INDIC=1, from the LAG=5 moments).
if doPart(10)==1
tempGE=load('./SavedOutput/HVY2011_GE.mat','p_eqm');
tempCalib=load('./SavedOutput/HVY2011calib.mat','CalibParams');
ParamSetNames={'calibrated','table2'};
nMonte=500;
nSize=200;
RETX=Nworkages; % 38 working ages
p0RIP=[-3.7;-1.531525;-0.852252;-4.044717];
p0HIP=[-3.87;-1.339970;0.314;-4.655511;-7.0074;0.0231];
optsT4=optimset('MaxIter',20000,'MaxFunEvals',20000,'TolX',1e-8,'TolFun',1e-12);
LagList=[1,5];
FnsToEvaluate_Sim=struct();
FnsToEvaluate_Sim.earnings=FnsToEvaluate.earnings;
simoptionsSim=simoptions;
simoptionsSim.numbersims=nSize;
Table4=struct();
for ss=1:2
    setname=ParamSetNames{ss};
    ParamsSet=Params_Table2;
    if ss==1
        calibfields=fieldnames(tempCalib.CalibParams);
        for pp=1:length(calibfields)
            ParamsSet.(calibfields{pp})=tempCalib.CalibParams.(calibfields{pp});
        end
    else
        ParamsSet.beta=tempGE.p_eqm.beta;
    end
    ParamsSet=PTypeParamFn(ParamsSet);
    jequaloneDistSet=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,ParamsSet.mean_logh1,ParamsSet.mean_logability,ParamsSet.stddev_logh1,ParamsSet.stddev_logability,ParamsSet.FTcorr_logh1logability,ParamsSet.mu_z,ParamsSet.sigma_z,ParamsSet.ability,N_i);
    ParamsSet.abilitydist=gather(reshape(sum(sum(jequaloneDistSet,1),2),[N_i,1])); % SimPanelValues uses abilitydist for the number of agents of each type
    [~, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsSet, DiscountFactorParamNames, vfoptions);

    % estimates: rows are panels; columns rho, sigma2_alpha, sigma2_eta, sigma2_epsilon (RIP), plus sigma2_beta, sigma_alphabeta (HIP)
    estRIP=zeros(nMonte,4,length(LagList));
    estHIP=zeros(nMonte,6,length(LagList));
    corrgrowth=zeros(nMonte,4);
    tic;
    for mm=1:nMonte
        SimPanel=SimPanelValues_FHorz_Case1_PType(jequaloneDistSet,PTypeDistParamNames,Policy,FnsToEvaluate_Sim,ParamsSet,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,simoptionsSim);
        e=gather(SimPanel.earnings(1:RETX,:))'; % [200 agents, 38 ages]
        le=log(max(e.*exp(0.15*randn(nSize,RETX)),0.0001));
        dle=le-mean(le,1);
        VARLOGE=mean(dle.^2,1)';
        for ll=1:length(LagList)
            LAG=LagList(ll);
            COVARLOGE=zeros(RETX,LAG);
            for H=4+LAG:RETX
                for IK=1:LAG
                    COVARLOGE(H,IK)=mean(dle(:,H).*dle(:,H-IK));
                end
            end
            g=fminsearch(@(guess) HVY2011_Table4_Objective(guess,0,VARLOGE,COVARLOGE,LAG),p0RIP,optsT4);
            estRIP(mm,:,ll)=[max(0.91,1-0.1*exp(g(3))), exp(g(2)), exp(g(1)), exp(g(4))];
            g=fminsearch(@(guess) HVY2011_Table4_Objective(guess,1,VARLOGE,COVARLOGE,LAG),p0HIP,optsT4);
            estHIP(mm,:,ll)=[max(0.40,1-0.1*exp(g(3))), exp(g(2)), exp(g(1)), exp(g(4)), exp(g(5)), exp(g(6))-1];
        end
        % Autocorrelation of earnings growth rates (with LAG=5): GR(:,:,IP) is the growth from age H-IP to H+1-IP, H=4+LAG,...,38
        LAG=5;
        GR=zeros(nSize,RETX-(4+LAG)+1,LAG);
        for IP=1:LAG
            GR(:,:,IP)=dle(:,(4+LAG:RETX)+1-IP)-dle(:,(4+LAG:RETX)-IP);
        end
        STDEVGR=reshape(sqrt(mean(mean(GR.^2,1),2)),[LAG,1]);
        for IP=1:LAG-1
            corrgrowth(mm,IP)=mean(mean(GR(:,:,1).*GR(:,:,IP+1)))/(STDEVGR(IP)*STDEVGR(IP+1));
        end
    end
    fprintf('Table 4 at the parameters of %s: %i panels simulated and estimated in %8.1fs \n',setname,nMonte,toc)
    Table4.(setname).RIPmean=reshape(mean(estRIP,1),[4,length(LagList)])';
    Table4.(setname).RIPstd=reshape(sqrt(mean(estRIP.^2,1)-mean(estRIP,1).^2),[4,length(LagList)])';
    Table4.(setname).HIPmean=reshape(mean(estHIP,1),[6,length(LagList)])';
    Table4.(setname).HIPstd=reshape(sqrt(mean(estHIP.^2,1)-mean(estHIP,1).^2),[6,length(LagList)])';
    Table4.(setname).corrgrowth=mean(corrgrowth,1);
    Table4.(setname).estRIP=estRIP;
    Table4.(setname).estHIP=estHIP;
    fprintf('\nTable 4 at the parameters of %s (mean across %i panels of %i agents, std dev in brackets): \n',setname,nMonte,nSize)
    fprintf('            rho, sigma2_alpha, sigma2_eta, sigma2_epsilon [, sigma2_beta, sigma_alphabeta] \n')
    for ll=1:length(LagList)
        fprintf('RIP, cov=%i: %s \n',LagList(ll),sprintf('%7.4f (%6.4f)  ',[Table4.(setname).RIPmean(ll,:);Table4.(setname).RIPstd(ll,:)]))
    end
    for ll=1:length(LagList)
        fprintf('HIP, cov=%i: %s \n',LagList(ll),sprintf('%9.6f (%8.6f)  ',[Table4.(setname).HIPmean(ll,:);Table4.(setname).HIPstd(ll,:)]))
    end
    fprintf('Autocorrelation of earnings growth rates, lags 1 to 4: %s (HVY2011: first order -0.386, higher orders nearly zero) \n',sprintf('%7.3f ',Table4.(setname).corrgrowth))
    save ./SavedOutput/HVY2011_Table4.mat Table4 % save after each parameter set
    clear Policy
end
fprintf('HVY2011 Table 4 (benchmark): RIP cov=1: 0.964, 0.283, 0.013, 0.025; cov=5: 0.963, 0.271, 0.014, 0.026; HIP cov=1: 0.860, 0.264, 0.032, 0.006, 0.000060, 0.0003; cov=5: 0.862, 0.268, 0.029, 0.010, 0.00006, 0.0005 \n')
end


%% doPart(11): the model with initial wealth differences (HVY2011 Section V.A; second columns of Tables 3 and 5, Table 6)
% Partial equilibrium (HVY2011 footnote 25): all parameters, and the prices and beta, are those of the parameter set; only the
% initial distribution changes. Initial assets are k1=(exp(x)-0.381)*(mean earnings at age 23), x~N(-0.277,0.921^2) (PSID, as
% in the authors' code), and (log h1, log a, x) is trivariate normal; the means and variances of (log h1, log a) and the three
% correlations are calibrated to the same earnings targets (as in the authors' distrib-estimation-initw.f), starting from the
% benchmark parameter set and the authors' correlations of wealth with ability (0.466) and with human capital (0.375).
% Then Table 3 (properties of the initial distribution), Table 5 (fractions of the variance due to initial conditions, where
% lifetime wealth is lifetime earnings plus initial wealth) and Table 6 (including +-1 st. dev. of log shifted initial wealth).
if doPart(11)==1
tempGE=load('./SavedOutput/HVY2011_GE.mat','p_eqm');
tempCalib=load('./SavedOutput/HVY2011calib.mat','CalibParams');
ParamSetNames={'calibrated','table2'};
[~,zeroassetindex]=min(abs(asset_grid));
u_gridc=gather(u_grid);
pi_uc=gather(pi_u);
EarningsReturnFn=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr);
jequaloneDistFnW=@(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTw_ha,FTw_hk,FTw_ak,mean_logwealth,stddev_logwealth,wealthshift,earn23,mu_z,sigma_z,ability,Number_i) HuggettVenturaYaron2011_jequaloneDistFn_wealth(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTw_ha,FTw_hk,FTw_ak,mean_logwealth,stddev_logwealth,wealthshift,earn23,mu_z,sigma_z,ability,Number_i);
CalibParamNamesW={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTw_ha','FTw_hk','FTw_ak'};
caliboptionsW=struct();
caliboptionsW.metric='sum_logratiosquared';
caliboptionsW.constrainpositive={'stddev_logh1','stddev_logability'};
caliboptionsW.constrainAtoB={'FTw_ha','FTw_hk','FTw_ak'};
caliboptionsW.constrainAtoBlimits.FTw_ha=[-3,3];
caliboptionsW.constrainAtoBlimits.FTw_hk=[-3,3];
caliboptionsW.constrainAtoBlimits.FTw_ak=[-3,3];
FnsToEvaluate_CalibW=struct();
FnsToEvaluate_CalibW.earnings=FnsToEvaluate.earnings;
FnsToEvaluate_CalibW.logearnings=FnsToEvaluate.logearnings;
WealthModel=struct();
for ss=1:2
    setname=ParamSetNames{ss};
    ParamsW=Params_Table2;
    if ss==1
        calibfields=fieldnames(tempCalib.CalibParams);
        for pp=1:length(calibfields)
            ParamsW.(calibfields{pp})=tempCalib.CalibParams.(calibfields{pp});
        end
    else
        ParamsW.beta=tempGE.p_eqm.beta;
    end
    ParamsW.discountearnings=1/(1+ParamsW.r);
    % Initial wealth (HVY2011 Section V.A and the authors' code)
    ParamsW.mean_logwealth=-0.277;
    ParamsW.stddev_logwealth=0.921412;
    ParamsW.wealthshift=0.381;
    ParamsW.earn23=Fig1_Mean_time(1); % data mean earnings at age 23 (the authors' TARGETEARN(1))
    % Initial guess for the correlations: the benchmark corr(log h1, log a), and the authors' corr(log h1,x)=0.37452, corr(log a,x)=0.466
    CorrW0=[1, tanh(ParamsW.FTcorr_logh1logability), 0.37452; tanh(ParamsW.FTcorr_logh1logability), 1, 0.466; 0.37452, 0.466, 1];
    LogCorrW0=logm(CorrW0);
    ParamsW.FTw_ha=LogCorrW0(2,1);
    ParamsW.FTw_hk=LogCorrW0(3,1);
    ParamsW.FTw_ak=LogCorrW0(3,2);
    ParamsW=PTypeParamFn(ParamsW);
    fprintf('\n=== Model with initial wealth differences, starting from the parameters of %s (alpha_h=%6.4f, beta=%8.5f) \n',setname,ParamsW.alpha_h,ParamsW.beta)

    %% Calibrate the initial distribution (partial equilibrium)
    tic;
    [CalibParamsW,calibsummaryW]=CalibrateLifeCycleModel_PType(CalibParamNamesW,TargetMoments,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsW, DiscountFactorParamNames, jequaloneDistFnW, AgeWeightsParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate_CalibW, caliboptionsW, vfoptions,simoptions);
    fprintf('Calibration of the initial distribution with initial wealth: run time %8.1fs \n',toc)
    for pp=1:length(CalibParamNamesW)
        ParamsW.(CalibParamNamesW{pp})=CalibParamsW.(CalibParamNamesW{pp});
    end
    ParamsW=PTypeParamFn(ParamsW);
    [CorrW,~]=GFT_inverse_mapping([ParamsW.FTw_ha;ParamsW.FTw_hk;ParamsW.FTw_ak],10^(-9));
    stddevW=[ParamsW.stddev_logh1; ParamsW.stddev_logability; ParamsW.stddev_logwealth];
    CoVarW=CorrW.*(stddevW*stddevW');
    fprintf('mean_logh1=%6.4f, mean_logability=%7.4f, var_logh1=%6.4f, var_logability=%6.4f; correlations of logs: (h1,a) %6.3f, (h1,x) %6.3f, (a,x) %6.3f \n',ParamsW.mean_logh1,ParamsW.mean_logability,CoVarW(1,1),CoVarW(2,2),CorrW(1,2),CorrW(1,3),CorrW(2,3))
    WealthModel.(setname).Params=ParamsW;
    WealthModel.(setname).CorrMatrix=CorrW;

    %% Solve the model with this initial distribution
    jequaloneDistW=HuggettVenturaYaron2011_jequaloneDistFn_wealth(a_grid,z_grid,n_a,n_z,ParamsW.mean_logh1,ParamsW.mean_logability,ParamsW.stddev_logh1,ParamsW.stddev_logability,ParamsW.FTw_ha,ParamsW.FTw_hk,ParamsW.FTw_ak,ParamsW.mean_logwealth,ParamsW.stddev_logwealth,ParamsW.wealthshift,ParamsW.earn23,ParamsW.mu_z,ParamsW.sigma_z,ParamsW.ability,N_i);
    ParamsW.abilitydist=gather(reshape(sum(sum(jequaloneDistW,1),2),[N_i,1]));
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsW, DiscountFactorParamNames, vfoptions);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistW,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],ParamsW,simoptions);
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_CalibW,ParamsW,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    objective1=sum(log(gather(AgeConditionalStats.earnings.Mean(1:Nworkages))./Fig1_Mean_time').^2);
    objective2=sum(log(gather(AgeConditionalStats.logearnings.Variance(1:Nworkages))./Fig1_VarLog_time').^2);
    objective3=sum(log(gather(AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages))./Fig1_Skew_time').^2);
    WealthModel.(setname).objective=objective1+objective2+objective3;
    fprintf('Calibration objective: %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',WealthModel.(setname).objective,objective1,objective2,objective3)

    %% Table 3: properties of the initial distribution of (a,h1) (log-normal)
    Table3.meanability=exp(ParamsW.mean_logability+CoVarW(2,2)/2);
    Table3.cvability=sqrt(exp(CoVarW(2,2))-1);
    Table3.meanh1=exp(ParamsW.mean_logh1+CoVarW(1,1)/2);
    Table3.cvh1=sqrt(exp(CoVarW(1,1))-1);
    Table3.corr=(exp(CoVarW(1,2))-1)/sqrt((exp(CoVarW(1,1))-1)*(exp(CoVarW(2,2))-1));
    WealthModel.(setname).Table3=Table3;
    fprintf('Table 3, with initial wealth (HVY2011 in brackets): mean ability %6.3f (0.328), CV ability %6.3f (0.124), mean h1 %6.1f (117.5), CV h1 %6.3f (0.476), corr(a,h1) %6.3f (0.655) \n',Table3.meanability,Table3.cvability,Table3.meanh1,Table3.cvh1,Table3.corr)

    %% Table 5: fractions of the variance due to initial conditions x1=(h1,a,k1), computed exactly as in doPart(6)
    % E[X|x1]: the value fn (utility) or the value fn of the policy with earnings as return fn, at age 23 at (k1, exp(z1)*h1),
    % averaged over the five values of z1 (linear interpolation in assets and human capital). Lifetime wealth is lifetime earnings
    % plus k1, and k1 is part of the initial condition, so Var(wealth)=Var(E[wealth|x1]) + (Var(earnings)-Var(E[earnings|x1])).
    [Wearn,~]=ValueFnFromPolicy_FHorz_PType(Policy,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,ParamsW,{'discountearnings'},vfoptions);
    FnsToEvaluate_f=struct();
    FnsToEvaluate_f.earnings=FnsToEvaluate.earnings;
    FnsToEvaluate_f.utility=@(l,aprime,a,h,r,gamma,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) (HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)^(1-gamma))/(1-gamma);
    PeriodValues=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsToEvaluate_f,ParamsW,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
    % The distribution of x1=(h1,a,x): (h1,a) on the grids, and x on nodes given (log h1, log a), as in the jequaloneDistFn
    % (but with h1 before the age-23 shock)
    P1=gather(MVNormal_ProbabilitiesOnGrid([log(h_grid);log(ParamsW.ability)],[ParamsW.mean_logh1,ParamsW.mean_logability],CoVarW(1:2,1:2),[n_a(2),N_i]));
    P1=P1/sum(P1(:));
    Bx=CoVarW(3,1:2)/CoVarW(1:2,1:2);
    condsd=sqrt(max(CoVarW(3,3)-Bx*CoVarW(3,1:2)',10^(-12)));
    [LogH,LogA]=ndgrid(log(h_grid),log(ParamsW.ability));
    condmean=ParamsW.mean_logwealth+Bx(1)*(LogH-ParamsW.mean_logh1)+Bx(2)*(LogA-ParamsW.mean_logability);
    nx=21;
    xnodes=ParamsW.mean_logwealth+ParamsW.stddev_logwealth*linspace(-3,3,nx)';
    xedges=[-Inf; (xnodes(1:end-1)+xnodes(2:end))/2; Inf];
    k1nodes=(exp(xnodes)-ParamsW.wealthshift)*ParamsW.earn23;
    EU_x1=zeros(n_a(2),N_i,nx);
    Eearn_x1=zeros(n_a(2),N_i,nx);
    Px1=zeros(n_a(2),N_i,nx);
    for kk=1:nx
        Px1(:,:,kk)=P1.*(0.5*erfc(-(xedges(kk+1)-condmean)/(condsd*sqrt(2)))-0.5*erfc(-(xedges(kk)-condmean)/(condsd*sqrt(2))));
    end
    Px1=Px1/sum(Px1(:));
    for ii=1:N_i
        Vii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        V1=Vii(:,:,1);
        W1=Wii(:,:,1);
        for zz=1:n_u
            hz=min(max(exp(u_gridc(zz))*h_grid,h_grid(1)),h_grid(end));
            for kk=1:nx
                EU_x1(:,ii,kk)=EU_x1(:,ii,kk)+pi_uc(zz)*interp2(h_grid',asset_grid,V1,hz,k1nodes(kk)*ones(n_a(2),1));
                Eearn_x1(:,ii,kk)=Eearn_x1(:,ii,kk)+pi_uc(zz)*interp2(h_grid',asset_grid,W1,hz,k1nodes(kk)*ones(n_a(2),1));
            end
        end
    end
    Ewealth_x1=Eearn_x1+reshape(k1nodes,[1,1,nx]);
    % (states with zero probability can have a value fn of -Inf; leave them out)
    keep=(Px1(:)>0);
    varU_initial=sum(Px1(keep).*(EU_x1(keep)-sum(Px1(keep).*EU_x1(keep))).^2);
    varEarn_initial=sum(Px1(keep).*(Eearn_x1(keep)-sum(Px1(keep).*Eearn_x1(keep))).^2);
    varWealth_initial=sum(Px1(keep).*(Ewealth_x1(keep)-sum(Px1(keep).*Ewealth_x1(keep))).^2);
    % Var(X), by the martingale decomposition over ages (as in doPart(6))
    momEarn=zeros(4,N_j);
    momU=zeros(4,N_j);
    for ii=1:N_i
        Dii=reshape(gather(StationaryDist.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Dii=StationaryDist.ptweights(ii)*Dii./max(sum(Dii,1),eps);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Fii=reshape(gather(PeriodValues.earnings.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Uii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Gii=reshape(gather(PeriodValues.utility.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Wii(Dii==0)=0; Fii(Dii==0)=0; Uii(Dii==0)=0; Gii(Dii==0)=0;
        Cii=(Wii-Fii)/ParamsW.discountearnings;
        CUii=(Uii-Gii)/ParamsW.beta;
        momEarn=momEarn+[sum(Dii.*Wii,1); sum(Dii.*Wii.^2,1); sum(Dii.*Cii,1); sum(Dii.*Cii.^2,1)];
        momU=momU+[sum(Dii.*Uii,1); sum(Dii.*Uii.^2,1); sum(Dii.*CUii,1); sum(Dii.*CUii.^2,1)];
    end
    varEarn_total=(momEarn(2,1)-momEarn(1,1)^2)+sum(ParamsW.discountearnings.^(2*(1:1:N_j-1)).*(momEarn(2,2:end)-momEarn(4,1:end-1)));
    varU_total=(momU(2,1)-momU(1,1)^2)+sum(ParamsW.beta.^(2*(1:1:N_j-1)).*(momU(2,2:end)-momU(4,1:end-1)));
    Table5.utility=varU_initial/varU_total;
    Table5.earnings=varEarn_initial/varEarn_total;
    Table5.wealth=varWealth_initial/(varWealth_initial+varEarn_total-varEarn_initial);
    WealthModel.(setname).Table5=Table5;
    fprintf('Table 5, with initial wealth (HVY2011 in brackets): lifetime utility %6.3f (0.661), lifetime earnings %6.3f (0.613), lifetime wealth %6.3f (0.626) \n',Table5.utility,Table5.earnings,Table5.wealth)
    fprintf('   checks: E_j[W_j]-E_{j-1}[E[W_j|s_{j-1}]], max relative over ages: earnings %8.2e, utility %8.2e; mean lifetime earnings over x1 %8.2f and over the age-23 distribution %8.2f (should be close) \n',max(abs(momEarn(1,2:end)-momEarn(3,1:end-1))./abs(momEarn(1,2:end))),max(abs(momU(1,2:end)-momU(3,1:end-1))./abs(momU(1,2:end))),sum(Px1(keep).*Eearn_x1(keep)),momEarn(1,1))
    clear Wearn PeriodValues StationaryDist AgeConditionalStats

    %% Table 6: one standard deviation changes in log h1, log a and x (log shifted initial wealth), from their means
    Params6=ParamsW;
    Params6.ability=exp(ParamsW.mean_logability+[-1;0;1]*ParamsW.stddev_logability);
    Params6.abilitydist=ones(3,1)/3;
    [V6, Policy6]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,3, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params6, DiscountFactorParamNames, vfoptions);
    [W6,~]=ValueFnFromPolicy_FHorz_PType(Policy6,n_d,n_a,n_z,N_j,3,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,Params6,{'discountearnings'},vfoptions);
    h1values=exp(ParamsW.mean_logh1+[-1,0,1]*ParamsW.stddev_logh1);
    k1values=(exp(ParamsW.mean_logwealth+[-1,0,1]*ParamsW.stddev_logwealth)-ParamsW.wealthshift)*ParamsW.earn23;
    EV6=zeros(3,3,3); % (h1 value, ability value, k1 value): expected lifetime utility
    EW6=zeros(3,3,3); % expected lifetime wealth (lifetime earnings plus k1)
    for kk=1:3
        Vkk=reshape(gather(V6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        Wkk=reshape(gather(W6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        for hh=1:3
            for ww=1:3
                for zz=1:n_u
                    EV6(hh,kk,ww)=EV6(hh,kk,ww)+pi_uc(zz)*interp2(h_grid',asset_grid,Vkk(:,:,1),exp(u_gridc(zz))*h1values(hh),k1values(ww));
                    EW6(hh,kk,ww)=EW6(hh,kk,ww)+pi_uc(zz)*interp2(h_grid',asset_grid,Wkk(:,:,1),exp(u_gridc(zz))*h1values(hh),k1values(ww));
                end
                EW6(hh,kk,ww)=EW6(hh,kk,ww)+k1values(ww);
            end
        end
    end
    % rows: h1 +1sd, h1 -1sd, ability +1sd, ability -1sd, wealth +1sd, wealth -1sd
    EVnew=[EV6(3,2,2);EV6(1,2,2);EV6(2,3,2);EV6(2,1,2);EV6(2,2,3);EV6(2,2,1)];
    EWnew=[EW6(3,2,2);EW6(1,2,2);EW6(2,3,2);EW6(2,1,2);EW6(2,2,3);EW6(2,2,1)];
    Table6.EV=100*(EVnew/EV6(2,2,2)).^(1/(1-ParamsW.gamma))-100;
    Table6.lifetimewealth=100*EWnew/EW6(2,2,2)-100;
    WealthModel.(setname).Table6=Table6;
    fprintf('Table 6, with initial wealth (HVY2011 in brackets): equivalent variation (percent): h1 +1sd %6.1f (39.3), h1 -1sd %6.1f (-28.3), ability +1sd %6.1f (5.7), ability -1sd %6.1f (-2.6), wealth +1sd %6.1f (7.1), wealth -1sd %6.1f (-1.6) \n',Table6.EV)
    fprintf('                                               change in lifetime wealth (percent): h1 +1sd %6.1f (47.5), h1 -1sd %6.1f (-31.7), ability +1sd %6.1f (8.1), ability -1sd %6.1f (-3.9), wealth +1sd %6.1f (5.0), wealth -1sd %6.1f (-1.3) \n',Table6.lifetimewealth)
    clear V Policy V6 Policy6 W6

    save ./SavedOutput/HVY2011_Wealth.mat WealthModel % save after each parameter set
end
end


%% doPart(12): HVY2011 Section V.C sensitivity (other than Fig 6), relative to our calibration
% Four scenarios, each recalibrated (with beta from the general eqm, as in doPart(4)), starting from our calibration:
%  'cohort':    the cohort effects view of the earnings facts; alpha_h and the initial distribution are recalibrated. The
%               retirement transfer and the scale of the income tax use the cohort-effects mean earnings profile [the paper does
%               not say; mu_z is kept]
%  'lowshock':  sigma_z=0.104 (the estimate minus its standard error); mu_z adjusted so the mean shock exp(mu_z+sigma_z^2/2)
%               is unchanged (so mean earnings fall at the same rate late in life); alpha_h fixed (HVY2011 footnote 28)
%  'highshock': sigma_z=0.118, likewise
%  'quantiles': the targets are the mean earnings profile and the profiles of the 90-10 and 75-25 quantile ratios of earnings
%               (time effects; computed in the data work), the ratios via CustomModelStats; alpha_h fixed
% For each, the fractions of the variance of lifetime utility and earnings due to initial conditions (as Table 5, doPart(6)),
% and the equivalent variations and changes in lifetime wealth of Table 6.
if doPart(12)==1
tempCalib=load('./SavedOutput/HVY2011calib.mat','CalibParams');
ScenarioNames={'cohort','lowshock','highshock','quantiles'};
[~,zeroassetindex]=min(abs(asset_grid));
EarningsReturnFn=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr);
Sensitivity=struct();
for sc=1:length(ScenarioNames)
    scname=ScenarioNames{sc};
    Params12=Params_Table2;
    calibfields=fieldnames(tempCalib.CalibParams); % includes beta
    for pp=1:length(calibfields)
        Params12.(calibfields{pp})=tempCalib.CalibParams.(calibfields{pp});
    end
    vfoptions12=vfoptions;
    simoptions12=simoptions;
    TargetMoments12=TargetMoments;
    caliboptions12=caliboptions;
    caliboptions12.constrain0to1={};
    FnsToEvaluate_Calib12=FnsToEvaluate_Calib;
    CalibParamNames12={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTcorr_logh1logability'};
    if strcmp(scname,'cohort')
        TargetMoments12.AgeConditionalStats.earnings.Mean=[Fig1_Mean_cohort', nan(1,N_j-length(Fig1_Mean_cohort))];
        TargetMoments12.AgeConditionalStats.logearnings.Variance=[Fig1_VarLog_cohort', nan(1,N_j-length(Fig1_VarLog_cohort))];
        TargetMoments12.AgeConditionalStats.earnings.RatioMeanToMedian=[Fig1_Skew_cohort', nan(1,N_j-length(Fig1_Skew_cohort))];
        Params12.ssbenefit=0.4*Fig1_Mean_cohort(end);
        Params12.taxfactor=32922/sum(popweights_working'.*Fig1_Mean_cohort);
        CalibParamNames12=[{'alpha_h'},CalibParamNames12];
        caliboptions12.constrain0to1={'alpha_h'};
    elseif strcmp(scname,'lowshock') || strcmp(scname,'highshock')
        if strcmp(scname,'lowshock')
            sigma12=0.104;
        else
            sigma12=0.118;
        end
        Params12.mu_z=Params_Table2.mu_z+(Params_Table2.sigma_z^2-sigma12^2)/2;
        Params12.sigma_z=sigma12;
        [u_grid12,pi_u12]=discretizeAR1_Tauchen(Params12.mu_z,0,Params12.sigma_z,n_u,2);
        pi_u12=pi_u12(1,:)';
        vfoptions12.u_grid=u_grid12;
        vfoptions12.pi_u=pi_u12;
        simoptions12.u_grid=u_grid12;
        simoptions12.pi_u=pi_u12;
    elseif strcmp(scname,'quantiles')
        TargetMoments12=struct();
        TargetMoments12.AgeConditionalStats.earnings.Mean=TargetMoments.AgeConditionalStats.earnings.Mean;
        TargetMoments12.CustomModelStats.p9010=Data_P9010_time';
        TargetMoments12.CustomModelStats.p7525=Data_P7525_time';
        caliboptions12.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,caliboptions,vfoptions,simoptions) HuggettVenturaYaron2011_QuantileRatios(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,caliboptions,vfoptions,simoptions);
    end
    Params12=PTypeParamFn(Params12);
    fprintf('\n=== Sensitivity, %s: calibrating %s (start alpha_h=%6.4f) \n',scname,strjoin(CalibParamNames12,', '),Params12.alpha_h)

    %% Calibrate
    tic;
    [CalibParams12,calibsummary12]=CalibrateOLGModel_PType(CalibParamNames12,TargetMoments12,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params12, DiscountFactorParamNames, jequaloneDistFn, AgeWeightsParamNames, GEPriceParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate_Calib12, GeneralEqmEqns, heteroagentoptions, caliboptions12, vfoptions12,simoptions12);
    fprintf('Sensitivity, %s: calibration run time %8.1fs \n',scname,toc)
    calibfields=fieldnames(CalibParams12); % includes beta (joint optimization)
    for pp=1:length(calibfields)
        Params12.(calibfields{pp})=CalibParams12.(calibfields{pp});
    end
    Params12=PTypeParamFn(Params12);
    Params12.discountearnings=1/(1+Params12.r);
    var_h=Params12.stddev_logh1^2;
    var_a=Params12.stddev_logability^2;
    corr_ha=tanh(Params12.FTcorr_logh1logability);
    cov_ha=corr_ha*Params12.stddev_logh1*Params12.stddev_logability;
    fprintf('Sensitivity, %s: beta=%8.5f, alpha_h=%6.4f, mean_logh1=%6.4f, mean_logability=%7.4f, var_logh1=%6.4f, var_logability=%6.4f, corr=%6.3f \n',scname,Params12.beta,Params12.alpha_h,Params12.mean_logh1,Params12.mean_logability,var_h,var_a,corr_ha)
    Sensitivity.(scname).Params=Params12;
    Sensitivity.(scname).calibsummary=calibsummary12;

    %% Solve the model at the recalibrated parameters
    jequaloneDist12=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,Params12.mean_logh1,Params12.mean_logability,Params12.stddev_logh1,Params12.stddev_logability,Params12.FTcorr_logh1logability,Params12.mu_z,Params12.sigma_z,Params12.ability,N_i);
    Params12.abilitydist=gather(reshape(sum(sum(jequaloneDist12,1),2),[N_i,1]));
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params12, DiscountFactorParamNames, vfoptions12);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist12,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params12,simoptions12);
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_Calib12,Params12,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions12);
    if strcmp(scname,'quantiles')
        CustomStats12=HuggettVenturaYaron2011_QuantileRatios(V,Policy,StationaryDist,Params12,FnsToEvaluate_Calib12,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],[],caliboptions12,vfoptions12,simoptions12);
        objective1=sum(log(gather(AgeConditionalStats.earnings.Mean(1:Nworkages))./Fig1_Mean_time').^2);
        objective2=sum(log(CustomStats12.p9010./Data_P9010_time').^2);
        objective3=sum(log(CustomStats12.p7525./Data_P7525_time').^2);
        fprintf('Sensitivity, %s: calibration objective %8.4f = mean %8.4f + 90-10 %8.4f + 75-25 %8.4f \n',scname,objective1+objective2+objective3,objective1,objective2,objective3)
    else
        objective1=sum(log(gather(AgeConditionalStats.earnings.Mean(1:Nworkages))./TargetMoments12.AgeConditionalStats.earnings.Mean(1:Nworkages)).^2);
        objective2=sum(log(gather(AgeConditionalStats.logearnings.Variance(1:Nworkages))./TargetMoments12.AgeConditionalStats.logearnings.Variance(1:Nworkages)).^2);
        objective3=sum(log(gather(AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages))./TargetMoments12.AgeConditionalStats.earnings.RatioMeanToMedian(1:Nworkages)).^2);
        fprintf('Sensitivity, %s: calibration objective %8.4f = mean %8.4f + var(log) %8.4f + mean/median %8.4f \n',scname,objective1+objective2+objective3,objective1,objective2,objective3)
    end
    Sensitivity.(scname).objective=objective1+objective2+objective3;
    clear AgeConditionalStats

    %% Table 5: fractions of the variance due to initial conditions (exactly as in doPart(6), with this scenario's shocks)
    u_gridc12=gather(vfoptions12.u_grid);
    pi_uc12=gather(vfoptions12.pi_u);
    [Wearn,~]=ValueFnFromPolicy_FHorz_PType(Policy,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,Params12,{'discountearnings'},vfoptions12);
    FnsToEvaluate_f=struct();
    FnsToEvaluate_f.earnings=FnsToEvaluate.earnings;
    FnsToEvaluate_f.utility=@(l,aprime,a,h,r,gamma,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2) (HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)^(1-gamma))/(1-gamma);
    PeriodValues=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsToEvaluate_f,Params12,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions12);
    P1=gather(MVNormal_ProbabilitiesOnGrid([log(h_grid);log(Params12.ability)],[Params12.mean_logh1,Params12.mean_logability],[var_h,cov_ha;cov_ha,var_a],[n_a(2),N_i]));
    P1=P1/sum(P1(:));
    EU_x1=zeros(n_a(2),N_i);
    Eearn_x1=zeros(n_a(2),N_i);
    for ii=1:N_i
        Vii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1),n_a(2),N_j]);
        V1=reshape(Vii(zeroassetindex,:,1),[n_a(2),1]);
        W1=reshape(Wii(zeroassetindex,:,1),[n_a(2),1]);
        for zz=1:n_u
            hz=min(max(exp(u_gridc12(zz))*h_grid,h_grid(1)),h_grid(end));
            EU_x1(:,ii)=EU_x1(:,ii)+pi_uc12(zz)*interp1(h_grid,V1,hz);
            Eearn_x1(:,ii)=Eearn_x1(:,ii)+pi_uc12(zz)*interp1(h_grid,W1,hz);
        end
    end
    varU_initial=sum(P1(:).*(EU_x1(:)-sum(P1(:).*EU_x1(:))).^2);
    varEarn_initial=sum(P1(:).*(Eearn_x1(:)-sum(P1(:).*Eearn_x1(:))).^2);
    momEarn=zeros(4,N_j);
    momU=zeros(4,N_j);
    for ii=1:N_i
        Dii=reshape(gather(StationaryDist.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Dii=StationaryDist.ptweights(ii)*Dii./max(sum(Dii,1),eps);
        Wii=reshape(gather(Wearn.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Fii=reshape(gather(PeriodValues.earnings.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Uii=reshape(gather(V.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Gii=reshape(gather(PeriodValues.utility.(sprintf('ptype%03i',ii))),[n_a(1)*n_a(2),N_j]);
        Wii(Dii==0)=0; Fii(Dii==0)=0; Uii(Dii==0)=0; Gii(Dii==0)=0;
        Cii=(Wii-Fii)/Params12.discountearnings;
        CUii=(Uii-Gii)/Params12.beta;
        momEarn=momEarn+[sum(Dii.*Wii,1); sum(Dii.*Wii.^2,1); sum(Dii.*Cii,1); sum(Dii.*Cii.^2,1)];
        momU=momU+[sum(Dii.*Uii,1); sum(Dii.*Uii.^2,1); sum(Dii.*CUii,1); sum(Dii.*CUii.^2,1)];
    end
    varEarn_total=(momEarn(2,1)-momEarn(1,1)^2)+sum(Params12.discountearnings.^(2*(1:1:N_j-1)).*(momEarn(2,2:end)-momEarn(4,1:end-1)));
    varU_total=(momU(2,1)-momU(1,1)^2)+sum(Params12.beta.^(2*(1:1:N_j-1)).*(momU(2,2:end)-momU(4,1:end-1)));
    Table5.utility=varU_initial/varU_total;
    Table5.earnings=varEarn_initial/varEarn_total;
    Sensitivity.(scname).Table5=Table5;
    fprintf('Sensitivity, %s: fraction of variance due to initial conditions: lifetime utility %6.3f, lifetime earnings %6.3f \n',scname,Table5.utility,Table5.earnings)
    fprintf('   check: E_j[W_j]-E_{j-1}[E[W_j|s_{j-1}]], max relative over ages: earnings %8.2e, utility %8.2e \n',max(abs(momEarn(1,2:end)-momEarn(3,1:end-1))./abs(momEarn(1,2:end))),max(abs(momU(1,2:end)-momU(3,1:end-1))./abs(momU(1,2:end))))
    clear Wearn PeriodValues StationaryDist

    %% Table 6: one standard deviation changes in log h1 and log a (as in doPart(6))
    Params6=Params12;
    Params6.ability=exp(Params12.mean_logability+[-1;0;1]*Params12.stddev_logability);
    Params6.abilitydist=ones(3,1)/3;
    [V6, Policy6]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,3, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params6, DiscountFactorParamNames, vfoptions12);
    [W6,~]=ValueFnFromPolicy_FHorz_PType(Policy6,n_d,n_a,n_z,N_j,3,d_grid,a_grid,z_grid,pi_z,EarningsReturnFn,Params6,{'discountearnings'},vfoptions12);
    h1values=exp(Params12.mean_logh1+[-1,0,1]*Params12.stddev_logh1);
    EV6=zeros(3,3);
    EW6=zeros(3,3);
    for kk=1:3
        Vkk=reshape(gather(V6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        Wkk=reshape(gather(W6.(sprintf('ptype%03i',kk))),[n_a(1),n_a(2),N_j]);
        V1=reshape(Vkk(zeroassetindex,:,1),[n_a(2),1]);
        W1=reshape(Wkk(zeroassetindex,:,1),[n_a(2),1]);
        for zz=1:n_u
            EV6(:,kk)=EV6(:,kk)+pi_uc12(zz)*interp1(h_grid,V1,exp(u_gridc12(zz))*h1values');
            EW6(:,kk)=EW6(:,kk)+pi_uc12(zz)*interp1(h_grid,W1,exp(u_gridc12(zz))*h1values');
        end
    end
    Table6.EV=100*([EV6(3,2);EV6(1,2);EV6(2,3);EV6(2,1)]/EV6(2,2)).^(1/(1-Params12.gamma))-100; % h1 +1sd, h1 -1sd, ability +1sd, ability -1sd
    Table6.lifetimewealth=100*[EW6(3,2);EW6(1,2);EW6(2,3);EW6(2,1)]/EW6(2,2)-100;
    Sensitivity.(scname).Table6=Table6;
    fprintf('Sensitivity, %s: equivalent variation (percent): h1 +1sd %6.1f, h1 -1sd %6.1f, ability +1sd %6.1f, ability -1sd %6.1f \n',scname,Table6.EV)
    clear V Policy V6 Policy6 W6

    save ./SavedOutput/HVY2011_Sensitivity.mat Sensitivity % save after each scenario
end
fprintf('\nHVY2011 Section V.C: cohort effects: initial conditions 60.3%% (earnings) and 64.9%% (utility), equivalent variations ability +1sd 8.0%%, h1 +1sd 34.1%%; low shocks 0.656 and 0.695; high shocks 0.570 and 0.620; quantile ratio moments 53.8%% and 57.3%% \n')
end


%% doPart(13): identification of the initial wealth correlations (model with initial wealth differences, our calibration)
% Are the correlations of initial wealth with initial human capital and with learning ability pinned down by the earnings
% targets? At the calibrated parameters of doPart(11) (our calibration):
%  (i) moment derivatives: central differences of the 114 targeted log moments (log mean earnings, log variance of log
%      earnings, log mean/median, ages 23-60) with respect to the seven calibrated quantities, with the three correlations
%      perturbed directly (and mapped to the Archakov-Hansen parameters exactly, gamma=vecl(logm(C)));
%  (ii) for each correlation, the share of its derivative column not explained by the other six columns (1-R^2; small means its
%      effect on the moments can be undone by the other parameters);
%  (iii) the Gauss-Newton approximation to the rise in the objective when the two wealth correlations move to the authors'
%      values (0.375 with h1, 0.466 with a) and the other parameters re-optimize locally;
%  (iv) the objective at the authors' wealth correlations with all other parameters at our values (an upper bound on the rise
%      in the objective; re-optimizing the other parameters could only lower it).
% Only the ability parameters change the value fn (through the ability grid); the others only change the initial distribution.
if doPart(13)==1
tempWealth=load('./SavedOutput/HVY2011_Wealth.mat','WealthModel');
ParamsI=tempWealth.WealthModel.calibrated.Params;
CorrI=tempWealth.WealthModel.calibrated.CorrMatrix;
% The seven calibrated quantities, in natural units: means and std devs of (log h1, log a), and the correlations (h1,a), (h1,x), (a,x)
IdentNames={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','corr_ha','corr_hx','corr_ax'};
theta0=[ParamsI.mean_logh1; ParamsI.mean_logability; ParamsI.stddev_logh1; ParamsI.stddev_logability; CorrI(2,1); CorrI(3,1); CorrI(3,2)];
steps=[0.02; 0.01; 0.01; 0.005; 0.05; 0.05; 0.05];
thetaAuthors=theta0;
thetaAuthors(6)=0.37452; % the authors' corr(log h1, x)
thetaAuthors(7)=0.466;   % the authors' corr(log a, x)
% Parameter vectors to evaluate: the calibration, +-step in each of the seven, and the authors' wealth correlations
nIdent=length(theta0);
ThetaList=[theta0, theta0*ones(1,2*nIdent), thetaAuthors];
for kk=1:nIdent
    ThetaList(kk,1+2*kk-1)=theta0(kk)+steps(kk);
    ThetaList(kk,1+2*kk)=theta0(kk)-steps(kk);
end
targetvecI=log([Fig1_Mean_time; Fig1_VarLog_time; Fig1_Skew_time]); % the 114 targeted log moments
logmomentsI=zeros(3*Nworkages,size(ThetaList,2));
FnsToEvaluate_I=struct();
FnsToEvaluate_I.earnings=FnsToEvaluate.earnings;
FnsToEvaluate_I.logearnings=FnsToEvaluate.logearnings;
[~, Policy0]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, PTypeParamFn(ParamsI), DiscountFactorParamNames, vfoptions);
for tt=1:size(ThetaList,2)
    ParamsT=ParamsI;
    ParamsT.mean_logh1=ThetaList(1,tt);
    ParamsT.mean_logability=ThetaList(2,tt);
    ParamsT.stddev_logh1=ThetaList(3,tt);
    ParamsT.stddev_logability=ThetaList(4,tt);
    CorrT=[1, ThetaList(5,tt), ThetaList(6,tt); ThetaList(5,tt), 1, ThetaList(7,tt); ThetaList(6,tt), ThetaList(7,tt), 1];
    if min(eig(CorrT))<=0
        error('doPart(13): the correlation matrix at evaluation %i is not positive definite',tt)
    end
    LogCorrT=logm(CorrT);
    ParamsT.FTw_ha=LogCorrT(2,1);
    ParamsT.FTw_hk=LogCorrT(3,1);
    ParamsT.FTw_ak=LogCorrT(3,2);
    ParamsT=PTypeParamFn(ParamsT); % ability grid
    jequaloneDistT=HuggettVenturaYaron2011_jequaloneDistFn_wealth(a_grid,z_grid,n_a,n_z,ParamsT.mean_logh1,ParamsT.mean_logability,ParamsT.stddev_logh1,ParamsT.stddev_logability,ParamsT.FTw_ha,ParamsT.FTw_hk,ParamsT.FTw_ak,ParamsT.mean_logwealth,ParamsT.stddev_logwealth,ParamsT.wealthshift,ParamsT.earn23,ParamsT.mu_z,ParamsT.sigma_z,ParamsT.ability,N_i);
    ParamsT.abilitydist=gather(reshape(sum(sum(jequaloneDistT,1),2),[N_i,1]));
    if ThetaList(2,tt)~=theta0(2) || ThetaList(4,tt)~=theta0(4) % the ability grid changed, so the value fn changes
        [~, PolicyT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, ParamsT, DiscountFactorParamNames, vfoptions);
    else
        PolicyT=Policy0;
    end
    StationaryDistT=StationaryDist_Case1_FHorz_PType(jequaloneDistT,AgeWeightsParamNames,PTypeDistParamNames,PolicyT,n_d,n_a,n_z,N_j,N_i,[],ParamsT,simoptions);
    AgeConditionalStatsT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDistT,PolicyT,FnsToEvaluate_I,ParamsT,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
    logmomentsI(:,tt)=log([gather(AgeConditionalStatsT.earnings.Mean(1:Nworkages))'; gather(AgeConditionalStatsT.logearnings.Variance(1:Nworkages))'; gather(AgeConditionalStatsT.earnings.RatioMeanToMedian(1:Nworkages))']);
    clear PolicyT StationaryDistT AgeConditionalStatsT
end
clear Policy0
objectiveI=sum((logmomentsI-targetvecI).^2,1); % the calibration objective at each evaluation
% (i) moment derivatives, [114 moments, 7 parameters]
JacI=zeros(3*Nworkages,nIdent);
for kk=1:nIdent
    JacI(:,kk)=(logmomentsI(:,1+2*kk-1)-logmomentsI(:,1+2*kk))/(2*steps(kk));
end
% (ii) share of each column not explained by the other six
unexplained=zeros(nIdent,1);
for kk=1:nIdent
    others=JacI(:,setdiff(1:nIdent,kk));
    resid=JacI(:,kk)-others*(others\JacI(:,kk));
    unexplained(kk)=sum(resid.^2)/sum(JacI(:,kk).^2);
end
% (iii) Gauss-Newton: objective(theta0+d) ~ objective(theta0)+2*r'*J*d+d'*J'*J*d; moving the two wealth correlations by D and
% re-optimizing the other five locally raises the objective by about D'*S*D, S the Schur complement of the other five in J'*J
JJ=JacI'*JacI;
wealthidx=[6,7];
otheridx=1:5;
SchurW=JJ(wealthidx,wealthidx)-JJ(wealthidx,otheridx)*(JJ(otheridx,otheridx)\JJ(otheridx,wealthidx));
DeltaW=thetaAuthors(wealthidx)-theta0(wealthidx);
riseGN=DeltaW'*SchurW*DeltaW;
riseGNfixed=DeltaW'*JJ(wealthidx,wealthidx)*DeltaW; % the same, with the other five held fixed
Identification=struct();
Identification.names=IdentNames;
Identification.theta0=theta0;
Identification.steps=steps;
Identification.Jacobian=JacI;
Identification.columnnorm=sqrt(sum(JacI.^2,1))';
Identification.unexplained=unexplained;
Identification.objective0=objectiveI(1);
Identification.objectiveAuthors=objectiveI(end);
Identification.riseGN=riseGN;
Identification.riseGNfixed=riseGNfixed;
Identification.singularvalues=svd(JacI);
fprintf('\nIdentification of the initial wealth correlations (model with initial wealth differences, our calibration) \n')
fprintf('   parameter            value     norm of d(log moments)/d(param)   share not explained by the other parameters \n')
for kk=1:nIdent
    fprintf('   %-18s %8.4f   %10.4f   %10.4f \n',IdentNames{kk},theta0(kk),Identification.columnnorm(kk),unexplained(kk))
end
fprintf('   singular values of the moment derivatives: %s \n',sprintf('%8.4f ',Identification.singularvalues))
tempResultsI=load('./SavedOutput/HVY2011_Results.mat','Results');
fprintf('   objective at the calibration %8.4f (doPart(11) found %8.4f, should agree); at the authors'' wealth correlations (0.375, 0.466), others fixed: %8.4f; benchmark model without initial wealth: %8.4f \n',objectiveI(1),tempWealth.WealthModel.calibrated.objective,objectiveI(end),tempResultsI.Results.calibrated.objective)
fprintf('   Gauss-Newton rise in the objective from moving the wealth correlations to the authors'' values: %8.4f with the other parameters fixed, %8.4f with them re-optimized locally \n',riseGNfixed,riseGN)
save ./SavedOutput/HVY2011_Identification.mat Identification
end


%% doPart(8): write the LaTeX tables (for LatexDoc/App_HuggettVenturaYaron2011.tex), from the saved results
% Uses only the .mat files saved by doPart(1), (6) and (7), plus the authors' model output, so it needs no GPU.
% The document shows each HVY2011 original just above the replication, so the replication tables contain only our numbers.
if doPart(8)==1
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
LTX=load('./SavedOutput/HVY2011_Results.mat','Results');
ParamSetNames={'calibrated','table2'};
for ss=1:2
    setname=ParamSetNames{ss};
    R=LTX.Results.(setname);
    P=R.Params;
    var_h=P.stddev_logh1^2;
    var_a=P.stddev_logability^2;
    cov_ha=tanh(P.FTcorr_logh1logability)*P.stddev_logh1*P.stddev_logability;

    % Table 2: parameter values
    FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table2_',setname,'.tex'],'w');
    fprintf(FID,'\\begin{tabular}{p{3.4cm}p{5.2cm}p{6.8cm}} \n \\hline \\hline \n');
    fprintf(FID,'Category & Symbol & Parameter value \\\\ \\hline \n');
    fprintf(FID,'Demographics & $(J,J_R,n)$ & $(J,J_R,n)=(%i,%i,%1.3f)$ \\\\ \n',P.J,P.Jr,P.n);
    fprintf(FID,'Preferences & $\\beta$, $u(c)=c^{1-\\rho}/(1-\\rho)$ & $(\\beta,\\rho)=(%1.4f,%i)$ \\\\ \n',P.beta,P.gamma);
    fprintf(FID,'Technology & $(\\gamma,\\delta,g)$ & $(\\gamma,\\delta,g)=(%1.3f,%1.4f,%1.4f)$ \\\\ \n',P.alpha_k,P.delta_k,P.g);
    fprintf(FID,'Tax system & $T_j=T_j^{ss}+T_j^{inc}$ & $T_j^{ss}(e_j)=%1.3fe_j$ for $j<J_R$, $T_j^{ss}=-%1.0f$ otherwise; $T_j^{inc}$ see text \\\\ \n',P.tau_ss,P.ssbenefit);
    fprintf(FID,'Human capital shocks & $z\\sim N(\\mu,\\sigma^2)$ & $(\\mu,\\sigma)=(%1.3f,%1.3f)$ \\\\ \n',P.mu_z,P.sigma_z);
    fprintf(FID,'Human capital technology & $h''=\\exp(z'')H(h,s,a)$, $H=h+a(hs)^{\\alpha}$ & $\\alpha=%1.3f$ \\\\ \n',P.alpha_h);
    fprintf(FID,'Initial conditions & $\\psi=LN(\\mu_x,\\Sigma)$ & $\\mu_x=(\\mu_h,\\mu_a)=(%1.3f,%1.3f)$ \\\\ \n',P.mean_logh1,P.mean_logability);
    fprintf(FID,' & & $(\\sigma_h^2,\\sigma_a^2,\\sigma_{ha})=(%1.3f,%1.4f,%1.4f)$ \\\\ \n',var_h,var_a,cov_ha);
    fprintf(FID,'\\hline \n \\end{tabular} \n');
    fclose(FID);

    % Table 3: properties of the initial distribution
    FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table3_',setname,'.tex'],'w');
    fprintf(FID,'\\begin{tabular*}{0.70\\textwidth}{@{\\extracolsep{\\fill}}lr} \n \\hline \\hline \n');
    fprintf(FID,'Statistic & Benchmark model \\\\ \\hline \n');
    fprintf(FID,'Mean learning ability ($a$) & %1.3f \\\\ \n',R.Table3.meanability);
    fprintf(FID,'Coefficient of variation ($a$) & %1.3f \\\\ \n',R.Table3.cvability);
    fprintf(FID,'Mean initial human capital ($h_1$) & %3.1f \\\\ \n',R.Table3.meanh1);
    fprintf(FID,'Coefficient of variation ($h_1$) & %1.3f \\\\ \n',R.Table3.cvh1);
    fprintf(FID,'Correlation ($a,h_1$) & %1.3f \\\\ \n',R.Table3.corr);
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);

    % Table 5: sources of lifetime inequality
    FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table5_',setname,'.tex'],'w');
    fprintf(FID,'\\begin{tabular*}{0.90\\textwidth}{@{\\extracolsep{\\fill}}lr} \n \\hline \\hline \n');
    fprintf(FID,'Statistic & Benchmark model \\\\ \\hline \n');
    fprintf(FID,'Fraction of variance in lifetime utility due to initial conditions & %1.3f \\\\ \n',R.Table5.utility);
    fprintf(FID,'Fraction of variance in lifetime earnings due to initial conditions & %1.3f \\\\ \n',R.Table5.earnings);
    fprintf(FID,'Fraction of variance in lifetime wealth due to initial conditions & %1.3f \\\\ \n',R.Table5.wealth);
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);

    % Table 6: changes in initial conditions
    FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table6_',setname,'.tex'],'w');
    fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llrr} \n \\hline \\hline \n');
    fprintf(FID,'Variable & Change in variable & Equivalent variation (\\%%) & Change in lifetime wealth (\\%%) \\\\ \\hline \n');
    fprintf(FID,'Human capital & $+1$ st. deviation & %3.1f & %3.1f \\\\ \n',R.Table6.EV(1),R.Table6.lifetimewealth(1));
    fprintf(FID,' & $-1$ st. deviation & %3.1f & %3.1f \\\\ \n',R.Table6.EV(2),R.Table6.lifetimewealth(2));
    fprintf(FID,'Learning ability & $+1$ st. deviation & %3.1f & %3.1f \\\\ \n',R.Table6.EV(3),R.Table6.lifetimewealth(3));
    fprintf(FID,' & $-1$ st. deviation & %3.1f & %3.1f \\\\ \n',R.Table6.EV(4),R.Table6.lifetimewealth(4));
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);
end

% The calibration objective: our calibration, the HVY2011 Table 2 parameters in our model, and the authors' own model output
% (data/hvy_FIGS/ModeloutputR2.xls, sheet time-alpha0.70: their model at the Table 2 parameters)
AuthorsModel_Mean=readmatrix('./data/hvy_FIGS/ModeloutputR2.xls','Sheet','time-alpha0.70','Range','D6:D43');
AuthorsModel_VarLog=readmatrix('./data/hvy_FIGS/ModeloutputR2.xls','Sheet','time-alpha0.70','Range','H6:H43');
AuthorsModel_Skew=readmatrix('./data/hvy_FIGS/ModeloutputR2.xls','Sheet','time-alpha0.70','Range','L6:L43');
AuthorsObj1=sum(log(AuthorsModel_Mean./Fig1_Mean_time).^2);
AuthorsObj2=sum(log(AuthorsModel_VarLog./Fig1_VarLog_time).^2);
AuthorsObj3=sum(log(AuthorsModel_Skew./Fig1_Skew_time).^2);
FID=fopen('./SavedOutput/LatexInputs/HVY2011_Objective.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lrrrr} \n \\hline \\hline \n');
fprintf(FID,' & Mean earnings & Var. log earnings & Mean/median & Total \\\\ \\hline \n');
fprintf(FID,'Our calibration & %1.4f & %1.4f & %1.4f & %1.4f \\\\ \n',LTX.Results.calibrated.objective1,LTX.Results.calibrated.objective2,LTX.Results.calibrated.objective3,LTX.Results.calibrated.objective);
fprintf(FID,'Table 2 parameters, our model & %1.4f & %1.4f & %1.4f & %1.4f \\\\ \n',LTX.Results.table2.objective1,LTX.Results.table2.objective2,LTX.Results.table2.objective3,LTX.Results.table2.objective);
fprintf(FID,'Table 2 parameters, authors'' model & %1.4f & %1.4f & %1.4f & %1.4f \\\\ \n',AuthorsObj1,AuthorsObj2,AuthorsObj3,AuthorsObj1+AuthorsObj2+AuthorsObj3);
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

% Table 5 alternative: counting the age-23 shock as an initial condition
FID=fopen('./SavedOutput/LatexInputs/HVY2011_Table5_inclz1.tex','w');
fprintf(FID,'\\begin{tabular*}{0.90\\textwidth}{@{\\extracolsep{\\fill}}lrr} \n \\hline \\hline \n');
fprintf(FID,'Statistic & Our calibration & Table 2 parameters \\\\ \\hline \n');
fprintf(FID,'Lifetime utility & %1.3f & %1.3f \\\\ \n',LTX.Results.calibrated.Table5.utility_inclz1,LTX.Results.table2.Table5.utility_inclz1);
fprintf(FID,'Lifetime earnings & %1.3f & %1.3f \\\\ \n',LTX.Results.calibrated.Table5.earnings_inclz1,LTX.Results.table2.Table5.earnings_inclz1);
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

% Figure 6 values (once doPart(7) has been run)
if exist('./SavedOutput/HVY2011_Fig6.mat','file')
    LTX7=load('./SavedOutput/HVY2011_Fig6.mat','Fig6');
    alphaall=[LTX7.Fig6.alpha,LTX7.Fig6.alpha_calibrated];
    earningsall=[LTX7.Fig6.earnings,LTX7.Fig6.earnings_calibrated];
    [alphaall,sortindex]=sort(alphaall);
    earningsall=earningsall(sortindex);
    FID=fopen('./SavedOutput/LatexInputs/HVY2011_Fig6values.tex','w');
    fprintf(FID,'\\begin{tabular}{l%s} \n \\hline \\hline \n',repmat('r',1,length(alphaall)));
    fprintf(FID,'$\\alpha$%s \\\\ \\hline \n',sprintf(' & %1.3f',alphaall));
    fprintf(FID,'Lifetime earnings%s \\\\ \n',sprintf(' & %1.3f',earningsall));
    fprintf(FID,'\\hline \n \\end{tabular} \n');
    fclose(FID);
end
% Section V.C sensitivity (once doPart(12) has been run). The paper reports these numbers only in its text (there is no
% original table), so the HVY2011 numbers are shown alongside.
if exist('./SavedOutput/HVY2011_Sensitivity.mat','file')
    LTX12=load('./SavedOutput/HVY2011_Sensitivity.mat','Sensitivity');
    FID=fopen('./SavedOutput/LatexInputs/HVY2011_Sensitivity.tex','w');
    fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrr} \n \\hline \\hline \n');
    fprintf(FID,' & & \\multicolumn{2}{c}{Replication} & \\multicolumn{2}{c}{HVY2011} & \\\\ \n');
    fprintf(FID,'Case & $\\alpha$ & Earnings & Utility & Earnings & Utility & Objective \\\\ \\hline \n');
    fprintf(FID,'Benchmark (our calibration) & %1.3f & %1.3f & %1.3f & 0.615 & 0.640 & %1.4f \\\\ \n',LTX.Results.calibrated.Params.alpha_h,LTX.Results.calibrated.Table5.earnings,LTX.Results.calibrated.Table5.utility,LTX.Results.calibrated.objective);
    sclabels={'Cohort effects view','Low shocks ($\sigma=0.104$)','High shocks ($\sigma=0.118$)','Quantile ratio moments'};
    scnames={'cohort','lowshock','highshock','quantiles'};
    paperE=[0.603,0.656,0.570,0.538];
    paperU=[0.649,0.695,0.620,0.573];
    for sc=1:4
        if isfield(LTX12.Sensitivity,scnames{sc})
            S12=LTX12.Sensitivity.(scnames{sc});
            fprintf(FID,'%s & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.4f \\\\ \n',sclabels{sc},S12.Params.alpha_h,S12.Table5.earnings,S12.Table5.utility,paperE(sc),paperU(sc),S12.objective);
        end
    end
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);
end

% Tables 3, 5 and 6 with the model with initial wealth differences (once doPart(11) has been run)
if exist('./SavedOutput/HVY2011_Wealth.mat','file')
    LTX11=load('./SavedOutput/HVY2011_Wealth.mat','WealthModel');
    Wsetnames=fieldnames(LTX11.WealthModel);
    for ss=1:length(Wsetnames)
        RB=LTX.Results.(Wsetnames{ss}); % benchmark
        RW=LTX11.WealthModel.(Wsetnames{ss}); % with initial wealth differences
        % Table 3: two columns
        FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table3W_',Wsetnames{ss},'.tex'],'w');
        fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lrr} \n \\hline \\hline \n');
        fprintf(FID,'Statistic & Benchmark model & \\parbox[b]{4cm}{\\raggedleft Benchmark model with initial wealth differences} \\\\ \\hline \n');
        fprintf(FID,'Mean learning ability ($a$) & %1.3f & %1.3f \\\\ \n',RB.Table3.meanability,RW.Table3.meanability);
        fprintf(FID,'Coefficient of variation ($a$) & %1.3f & %1.3f \\\\ \n',RB.Table3.cvability,RW.Table3.cvability);
        fprintf(FID,'Mean initial human capital ($h_1$) & %3.1f & %3.1f \\\\ \n',RB.Table3.meanh1,RW.Table3.meanh1);
        fprintf(FID,'Coefficient of variation ($h_1$) & %1.3f & %1.3f \\\\ \n',RB.Table3.cvh1,RW.Table3.cvh1);
        fprintf(FID,'Correlation ($a,h_1$) & %1.3f & %1.3f \\\\ \n',RB.Table3.corr,RW.Table3.corr);
        fprintf(FID,'\\hline \n \\end{tabular*} \n');
        fclose(FID);
        % Table 5: two columns
        FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table5W_',Wsetnames{ss},'.tex'],'w');
        fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}p{8.5cm}rr} \n \\hline \\hline \n');
        fprintf(FID,'Statistic & \\parbox[b]{2cm}{\\raggedleft Benchmark model} & \\parbox[b]{3.5cm}{\\raggedleft Benchmark model with initial wealth differences} \\\\ \\hline \n');
        fprintf(FID,'Fraction of variance in lifetime utility due to initial conditions & %1.3f & %1.3f \\\\ \n',RB.Table5.utility,RW.Table5.utility);
        fprintf(FID,'Fraction of variance in lifetime earnings due to initial conditions & %1.3f & %1.3f \\\\ \n',RB.Table5.earnings,RW.Table5.earnings);
        fprintf(FID,'Fraction of variance in lifetime wealth due to initial conditions & %1.3f & %1.3f \\\\ \n',RB.Table5.wealth,RW.Table5.wealth);
        fprintf(FID,'\\hline \n \\end{tabular*} \n');
        fclose(FID);
        % Table 6: model with initial wealth differences (all six rows)
        FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table6W_',Wsetnames{ss},'.tex'],'w');
        fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llrr} \n \\hline \\hline \n');
        fprintf(FID,'Variable & Change in variable & Equivalent variation (\\%%) & Change in lifetime wealth (\\%%) \\\\ \\hline \n');
        rownames={'Human capital','','Learning ability','','Initial wealth',''};
        changenames={'$+1$ st. deviation','$-1$ st. deviation'};
        for rr=1:6
            fprintf(FID,'%s & %s & %3.1f & %3.1f \\\\ \n',rownames{rr},changenames{2-mod(rr,2)},RW.Table6.EV(rr),RW.Table6.lifetimewealth(rr));
        end
        fprintf(FID,'\\hline \n \\end{tabular*} \n');
        fclose(FID);
    end
end

% Table 4: statistical models of earnings (once doPart(10) has been run)
if exist('./SavedOutput/HVY2011_Table4.mat','file')
    LTX10=load('./SavedOutput/HVY2011_Table4.mat','Table4');
    T4setnames=fieldnames(LTX10.Table4);
    for ss=1:length(T4setnames)
        T4=LTX10.Table4.(T4setnames{ss});
        FID=fopen(['./SavedOutput/LatexInputs/HVY2011_Table4_',T4setnames{ss},'.tex'],'w');
        fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \n \\hline \\hline \n');
        fprintf(FID,'Statistical model & $\\rho$ & $\\sigma^2_{\\alpha}$ & $\\sigma^2_{\\eta}$ & $\\sigma^2_{\\epsilon}$ & $\\sigma^2_{\\beta}$ & $\\sigma_{\\alpha\\beta}$ \\\\ \\hline \n');
        fprintf(FID,'RIP model & & & & & & \\\\ \n');
        lagnames={'cov = 1','cov = 5'};
        for ll=1:2
            fprintf(FID,'\\quad %s & %1.3f & %1.3f & %1.3f & %1.3f & --- & --- \\\\ \n',lagnames{ll},T4.RIPmean(ll,:));
            fprintf(FID,' & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & --- & --- \\\\ \n',T4.RIPstd(ll,:));
        end
        fprintf(FID,'HIP model & & & & & & \\\\ \n');
        for ll=1:2
            fprintf(FID,'\\quad %s & %1.3f & %1.3f & %1.3f & %1.3f & %1.6f & %1.4f \\\\ \n',lagnames{ll},T4.HIPmean(ll,:));
            fprintf(FID,' & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.6f) & (%1.4f) \\\\ \n',T4.HIPstd(ll,:));
        end
        fprintf(FID,'\\hline \n \\end{tabular*} \n');
        fclose(FID);
    end
end

% Table 1: estimation of human capital shocks (once doPart(9) has been run)
if exist('./SavedOutput/HVY2011_Table1.mat','file')
    LTX9=load('./SavedOutput/HVY2011_Table1.mat','Table1');
    T1=LTX9.Table1;
    FID=fopen('./SavedOutput/LatexInputs/HVY2011_Table1.tex','w');
    fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrrr} \n \\hline \\hline \n');
    fprintf(FID,'Age range & Period & $e_{min}$ & $e_{max}$ & $N$ & $\\sigma$ & SE($\\sigma$) & $\\sigma_{\\epsilon}$ & SE($\\sigma_{\\epsilon}$) \\\\ \\hline \n');
    for ss=1:size(T1.specs,1)
        if T1.specs(ss,3)==69 && T1.specs(ss,4)==100
            periodstr='1969--2004';
        elseif T1.specs(ss,3)==69
            periodstr='1969--1981';
        else
            periodstr='1982--2004';
        end
        fprintf(FID,'%i--%i & %s & %i & 1.8 M & %5.1f & %1.3f & %1.3f & %1.3f & %1.3f \\\\',T1.specs(ss,1),T1.specs(ss,2),periodstr,T1.specs(ss,5),T1.N(ss),T1.sigma(ss),T1.se_sigma(ss),T1.sigma_u(ss),T1.se_sigma_u(ss));
        if ss==2 || ss==6 || ss==8
            fprintf(FID,'[1ex] \n');
        else
            fprintf(FID,' \n');
        end
    end
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);
end
fprintf('doPart(8): LaTeX tables written to ./SavedOutput/LatexInputs/ \n')
end


%% End of run
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','on'); % do not leave an interactive session silently offscreen
end
fprintf('HuggettVenturaYaron2011: run finished (doPart=%s)\n',mat2str(doPart))
diary off
