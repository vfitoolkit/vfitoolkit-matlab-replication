% Huggett & Ventura (2000) - Understanding why high income households save more than low income households
% https://doi.org/10.1016/S0304-3932(99)00058-6
%
% Four subscripts: HV2000_Model1, HV2000_Model2, HV2000_Model3, HV2000_Model4.
% Solve the four models.
%
% This script does setup, calls the four model codes, and then create
% results.
%
% This does not attempt a full replication, it is instead a demo of how
% easy AI Coding an OLG model is using Clause + VFI Toolkit.
%
% We just produce half of Table 4, specifically the parts for alowerbar=0
%
% ebar (average past indexed earnings) is implemented as an experience asset:
% Models 1-3 use experienceassetz with aprime(d,ebar,z);
% Model 4 uses experienceassetze with aprime(d,ebar,z,e).
% HV2000 has no labor choice, so the decision dimension d is a singleton dummy.

% Record the full run (incl. any error) to HV2000_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./HV2000_Diary.txt','file'), delete('./HV2000_Diary.txt'); end
diary ./HV2000_Diary.txt

% Subcodes (model scripts, tables, figures, and the ReturnFn/ConsumptionFn
% family) live in ./HV2000subcodes/; MATLAB does not search subfolders of the
% working directory, so put it on the path. Note this does NOT change cwd, so
% the './SavedOutput/...' paths inside the subcodes still resolve to here.
addpath('./HV2000subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
addpath(genpath('./VFIToolkit-matlab/'));

doPart=[0,0,0,0,1];
% doPart(1): main results
% doPart(2): sensitivity - shocks
% doPart(3): sensitivity - borrowing limits
% doPart(4): rep. agent: average within-age
% doPart(5): alternative treatment of bequests

close all % close any figures, make sure they are all cleanly built from scratch
% Make sure the subfolders to save output exist
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

%% Sizes (asset / ebar; z and possibly e are set per-model)
n_d=1;
n_a=[251,75]; % [assets, ebar]
% Asset points raised 201->251 alongside amax 50->150 below: the cubic spacing puts
% half the points below amax/8, so tripling amax on the old 201 would have coarsened
% the grid exactly where most of the mass sits.
N_j=81;       % ages 20..100



%% Parameters - HV2000 Table 2
Params.agejshifter=19;
Params.J=N_j;
Params.agej=(1:1:N_j)';
Params.Jr=46;

% Preferences
Params.beta=1.011;
Params.sigma=1.12;

% Production
Params.A=0.89594;
Params.alpha=0.36;
Params.delta=0.06;
Params.g=0.018;

% Demographics
Params.n=0.012;

% Government
Params.GoverY=0.195;

%% US 1990 data: survival probabilities sj and age-earnings profile ybar_j
HuggettVentura2000_USdata;
Params.sj=HV2000_sj;
Params.ybar_j=HV2000_ybar;
clear HV2000_sj HV2000_ybar

%% Age weights mewj
AgeWeightsParamNames={'mewj'};

Params.mewj=ones(N_j,1);
for jj=2:N_j
    Params.mewj(jj)=Params.mewj(jj-1)*Params.sj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj/sum(Params.mewj);
Params.mewj=Params.mewj';

%% Initial guesses for GE prices (carried over between models)
Params.r=0.06;

KoverL=(Params.alpha*Params.A/(Params.r+Params.delta))^(1/(1-Params.alpha));

Params.w=(1-Params.alpha)*Params.A*KoverL^Params.alpha;

KoverY=KoverL^(1-Params.alpha)/Params.A;

Params.tau=Params.GoverY/(1-Params.delta*KoverY);
Params.theta=0.10;
Params.T=0.05;

fractionworking=(Params.agej<Params.Jr);

% Population-weighted mean of ybar_j over working ages. Constant (no GE
% dependence), so we precompute. Bend points and the AIME cap inside
% ReturnFn/ebarprimeFn/etc are written as 0.20*w*ybar_mean, 1.24*w*ybar_mean,
% 2.47*w*ybar_mean so they scale with the GE-determined wage w.
Params.ybar_mean=sum(Params.mewj(fractionworking).*Params.ybar_j(fractionworking)')/sum(Params.mewj(fractionworking));
Params.b_common=0.1242*Params.A*KoverL^Params.alpha;


%% Discount factor
Params.beta_g=Params.beta*(1+Params.g)^(1-Params.sigma);
DiscountFactorParamNames={'beta_g','sj'};

%% Grids (asset, ebar, dummy d)
d_grid=0;

Params.amax=150; % was 50, which bound: 0.37-0.58% of mass sat in the top ten asset
% gridpoints across Tables 6-9 (tolerance 1e-5), i.e. a_max was cutting off households
% that wanted to hold more. With the cubic spacing the old top ten points covered
% a in [42.5,50], so the constrained mass was real accumulation rather than rounding.
% NB: changing this changes every number in the replication. The cached
% SavedOutput/HV2000_AllModels{0,1}.mat and SavedOutput/Intermediates/*.mat must be
% deleted before re-running, or those stages will be skipped and silently mix grids.
Params.amin=-1.5; % covers a̲=-w with margin: w lies in [0.94, 0.98] across the four ā=0 GE solutions
Params.borrowFlag=1; % binary: 0 => a̲=0 (no borrowing); 1 => a̲=-w (HV2000 Table 2 column ā)
Params.haveSS=1; % binary: 1 => social security on; 0 => no SS (used by HV2000_Table7)
Params.headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end
Params.doPart1figures=0; % 1 during the doPart(1) main loop, 0 elsewhere; guards Fig 3, 4, 7 creation
asset_grid=[linspace(Params.amin,0,31)'; Params.amax*linspace(1/(n_a(1)-31),1,n_a(1)-31)'.^3];
azeroindex=find(asset_grid==0,1); % newborn asset-index: jequaloneDist should start mass at a=0, not at a=amin

Params.ebar_gridmax=5;
ebar_grid=linspace(0,Params.ebar_gridmax,n_a(2))';

a_grid=[asset_grid;ebar_grid];

%% experienceassetz setup -- BASE simoptions/vfoptions
% Each model script copies these into its own per-model variant:
%   HV2000_Model1     -> simoptions1 / vfoptions1
%   HV2000_Model2     -> simoptions2 / vfoptions2
%   HV2000_Model3     -> simoptions3 / vfoptions3
%   HV2000_Model4     -> simoptions4 / vfoptions4 (flips to experienceassetze, adds e-channel)
%   HV2000_Model2_ptype -> simoptions5 / vfoptions5 (PType variant of Model 2)
% Downstream Tables/Figures pick the variant that matches the model they extend
% (e.g. HV2000_Table6 and HV2000_Figure4/Figure7 use simoptions4).
vfoptions.experienceassetz=1;
simoptions.experienceassetz=1;
vfoptions.aprimeFn=@(d,ebar,z,agej,Jr,w,ybar_j,ybar_mean,g) HuggettVentura2000_ebarprimeFn(d,ebar,z,agej,Jr,w,ybar_j,ybar_mean,g);
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

%% VFI solver options
vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=50;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Return function (Models 1-3; Model 4 overrides to HV2000_ReturnFn_M4)
ReturnFn=@(d,aprime,a,ebar,z,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS)...
    HuggettVentura2000_ReturnFn(d,aprime,a,ebar,z,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS);

%% FnsToEvaluate (Models 1-3; Model 4 builds its own with e in signatures)
FnsToEvaluate.K=@(d,aprime,a,ebar,z) a;
FnsToEvaluate.L=@(d,aprime,a,ebar,z,agej,Jr,ybar_j) (agej<Jr)*z*ybar_j;
FnsToEvaluate.Earnings=@(d,aprime,a,ebar,z,agej,Jr,ybar_j,w) (agej<Jr)*z*ybar_j*w;
FnsToEvaluate.AccBeq=@(d,aprime,a,ebar,z,sj,r,tau) (1-sj)*aprime*(1+r*(1-tau));
FnsToEvaluate.SSBenefits=@(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS)...
    HV2000_SSBenefitFn(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS);
FnsToEvaluate.Savings=@(d,aprime,a,ebar,z,g) (1+g)*aprime-a;
FnsToEvaluate.Ebar=@(d,aprime,a,ebar,z) ebar;
FnsToEvaluate.Consumption=@(d,aprime,a,ebar,z,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS)...
    HuggettVentura2000_ConsumptionFn(d,aprime,a,ebar,z,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS);
FnsToEvaluate.Income=@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)...
    HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS);
FnsToEvaluate.SavingRate=@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)...
    HV2000_SavingRateFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS);

%% General equilibrium setup (HV2000 conditions 1-7, Section 3.3)
GEPriceParamNames={'r','w','tau','T','theta','b_common'};

heteroagentoptions.verbose=1;

heteroagentoptions.intermediateEqns.Y=@(A,alpha,K,L) A*K^alpha*L^(1-alpha); % Cobb-Douglas output (HV2000 eqn 3, transformed-variables form). Used by TaxRate and bCommon.

GeneralEqmEqns.CapitalMarket=@(r,A,alpha,K,L,delta) r-(alpha*A*K^(alpha-1)*L^(1-alpha)-delta); % Eqm condition 2 (Section 3.3): r=F_1(K,L)-delta, the firm's FOC on capital.
GeneralEqmEqns.LaborMarket=@(w,A,alpha,K,L) w-(1-alpha)*A*K^alpha*L^(-alpha);                  % Eqm condition 2 (Section 3.3): w=F_2(K,L), the firm's FOC on labor.
GeneralEqmEqns.TaxRate=@(tau,GoverY,delta,K,Y) tau-GoverY/(1-delta*K/Y);                       % Section 4.1: tau=0.195/(1-delta*K/Y), from eqm cond 5 (gov budget) given G/Y=0.195 average 1959-1993.
GeneralEqmEqns.BequestBalance=@(T,n,AccBeq) T*(1+n)-AccBeq;                                    % Eqm condition 7 (eqn for T-hat, Section 3.3): lump-sum transfer T equals per-capita aggregate after-tax bequests, scaled by 1+n for newborn cohort.
GeneralEqmEqns.SSBalance=@(theta,w,L,SSBenefits) theta*w*L-SSBenefits;                         % Eqm condition 6 (Section 3.3): SS payroll-tax revenue theta*w*L equals aggregate benefits paid to retirees (pay-as-you-go).
GeneralEqmEqns.bCommon=@(b_common,Y) b_common-0.1242*Y;                                        % Section 4.3 (p. 376): common medical/hospital SS benefit b = 0.1242*Y (= 7.72%+4.70% of GDP per person over 20, 1990-1994 average from SS Bulletin 1996).
% Since this model has exogenous labor and exogenous labor productivity, I
% expect you could replace SSBalance with w*SSBalance (by redefining
% SSBalance) and thus eliminate the GeneralEqmEqns.SSBalance. But model is
% easy enough to solve, so not bothering.

%% Run each model (warm-starts GE prices from previous model's solution)
Results=struct();
if doPart(1)==1

    Params.doPart1figures=doPart(1); % 1 during the doPart(1) main loop
    % Solve every model twice: once with a̲=0 (no borrowing) and once with a̲=-w
    % (borrow up to one period's wage), matching the two columns in HV2000 Table 2.
    % Each borrowing limit saves to its own .mat and is skipped when that .mat
    % already exists. This lets the two run as separate jobs and lets a
    % crashed or timed-out run resume without redoing the finished limit.
    % NB: delete the .mat files to force a fresh solve after changing parameters.
    for borrowFlag=0:1
        AllModelsFile=sprintf('./SavedOutput/HV2000_AllModels%d.mat',borrowFlag);
        if ~exist(AllModelsFile,'file')
            Params.borrowFlag=borrowFlag;

            HV2000_Model1;

            HV2000_Model2;

            HV2000_Model3;

            HV2000_Model4;

            % V, Policy and StationaryDist are no longer needed, and memory available is limited, so clear them
            % (this clear lives here rather than at the end of HV2000_Model4, because
            % HV2000_Table6 calls that script and then re-uses Policy and StationaryDist)
            clear V Policy StationaryDist

            save(AllModelsFile,'Results','Params','n_d','n_a','N_j','d_grid','a_grid');
        end
    end
    Params.doPart1figures=0;

    % Combine the two borrowing limits, then write Tables 2-5 (Table 4 has both an
    % a̲=0 and an a̲=-w column, so both files are needed).
    if exist('./SavedOutput/HV2000_AllModels0.mat','file') && exist('./SavedOutput/HV2000_AllModels1.mat','file')
        R0=load('./SavedOutput/HV2000_AllModels0.mat','Results');
        R1=load('./SavedOutput/HV2000_AllModels1.mat','Results');
        Results=R0.Results;
        ResultsFields=fieldnames(R1.Results);
        for ii=1:length(ResultsFields)
            Results.(ResultsFields{ii})=R1.Results.(ResultsFields{ii});
        end
        clear R0 R1 ResultsFields

        HV2000_WriteTables;

        save('./SavedOutput/HV2000_Results_Part1.mat','Results')
    else
        fprintf('doPart(1): need both HV2000_AllModels0.mat and HV2000_AllModels1.mat for HV2000_WriteTables; re-run to finish the missing borrowing limit\n');
    end
end

%% Tables 6 and 7 (sensitivity exercises) - extra GE solves.
if doPart(2)==1
    HV2000_Table6; % Model 4 x 3 temp-shock variances (3 extra GE solves)
end
if doPart(3)==1
    HV2000_Table7; % Models 2-4 x 2 borrowing limits, no SS (6 extra GE solves)
end

%% Tables 8 and 9 (within-age-group averaging counterfactual)
if doPart(4)==1
    HV2000_Tables89; % within-age-group averaging counterfactual (12 extra GE solves)

    save('./SavedOutput/HV2000_Results_Part3.mat','Results')

    save('./SavedOutput/HV2000_doPart3.mat','-v7.3')
end


%% Figure 6 (alternative bequest treatments in Model 2) - extra GE solves.
if doPart(5)==1
    % Switches Model 2 to N_i=21 permanent types and routes accidental bequests
    % by ability type (per-type T) via heteroagentoptions.GEptype={'BequestBalance'}.
    % Adds 2 extra GE solves (one per borrowing limit).
    HV2000_Model2_ptype; % PType variant of Model 2: bequests by ability type
    HV2000_Figure6;      % plot Figure 6 ("Equal Transfers" vs "Different Transfers")

    save('./SavedOutput/HV2000_Results_Part4.mat','Results')

    save('./SavedOutput/HV2000_doPart4.mat','-v7.3')
end


%% CheckTopOfGrid report: mass near the top of the asset / ebar grids across all solved models
% Threshold: 1e-5 (0.001% of population mass); any value above this triggers a warning
% suggesting Params.amax (asset grid ceiling) or Params.ebar_gridmax (ebar grid ceiling) be widened.
checkTol=1e-5;
CheckTopOfGrid_models={'Model1_abar0','Model1_abarw','Model2_abar0','Model2_abarw', ...
                       'Model3_abar0','Model3_abarw','Model4_abar0','Model4_abarw', ...
                       'Model2_ptype_abar0','Model2_ptype_abarw'};
fprintf('\n========== CheckTopOfGrid across all solved models (threshold = %.0e) ==========\n',checkTol);
fprintf('%-22s %14s %14s %14s %14s %14s\n','Model', ...
    'Top10Assets','Top5Assets','Top2Assets','Top5Ebar','Top2Ebar');
for mm=1:length(CheckTopOfGrid_models)
    nm=CheckTopOfGrid_models{mm};
    if isfield(Results,nm) && isfield(Results.(nm),'CheckTopOfGrid')
        C=Results.(nm).CheckTopOfGrid;
        fprintf('%-22s %14.4e %14.4e %14.4e %14.4e %14.4e\n',nm, ...
            C.MassAtTop10gridptsAssets,C.MassAtTop5gridptsAssets,C.MassAtTop2gridptsAssets, ...
            C.MassAtTop5gridptsEbar,C.MassAtTop2gridptsEbar);
        if C.MassAtTop10gridptsAssets>checkTol
            warning('%s: %.4e mass in top 10 asset gridpts (>%.0e); widen Params.amax',nm,C.MassAtTop10gridptsAssets,checkTol);
        end
        if C.MassAtTop5gridptsEbar>checkTol
            warning('%s: %.4e mass in top 5 ebar gridpts (>%.0e); widen Params.ebar_gridmax',nm,C.MassAtTop5gridptsEbar,checkTol);
        end
    end
end

%% End of run
diary off
