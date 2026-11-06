% De Nardi, French & Jones (2010) - Why do the Elderly Save? The Role of Medical Expenses

% Because VFI Toolkit does not presently have a way to allow survival
% probabilities to depend on the state of the model, we instead have to
% implement this via adding a 'dead' value to the health state, so it takes 
% three values instead of two. To be able to use this probability to
% determine the warm-glow of bequests, we have to further differentiate
% 'death' from 'dead', so the health state now takes four values, and
% 'death' indicates the first period of being dead in which a warm-glow is
% received.

% A few notational changes: I use j for age (DFJ2010 use t). I use q_c to
% index the income quintiles and q1,q2,...,q5 as names for them (DFJ2010 use I).

% Asset profiles look wrong.
% When I disable the warm-glow of bequests (theta=0), things are okay. 
% [Female medical expenses shocks near end of life are still large enough 
% to create hump-shape for assets, but this does not happen for males. Not 
% clear if the female are incorrect, as DFJ2010 paper never plots results 
% for just females.]

% I believe that the code here is fine and that the problem with bequests
% is as follows: DFJ2010 solve the model by the FOCs/Euler eqns, hence
% their code does not use the warm-glow of bequests function, only it's
% derivative. Hence their parameters are all fine, but the warm-glow is
% missing a constant term and so here where I use the warm-glow, rather
% than it's derivative, this missing constant means the warm-glow is way to
% big and this is what is making things look wrong.
% [UPDATE: the C code (wealth20_Singles.c line 1084) does use the warm-glow
% level, phi0*(b+K0)^(1-nu)/(1-nu), so there is no missing constant. The more
% likely culprit was k=215: Table 3 reports k in thousands, so k=215,000.
% Now fixed below; the comments above are kept as a record until it is re-run.]

% psi=(zeta+xi) [in the medical expenses] has unit variance by construction in DFJ2010 (Table 1's
% variances 0.050/(1-0.922^2)+0.665=0.999); the 'normalizepsi' parameter in the codes makes this explicit (=1.0).

%% Three parameter values not in the paper (tau_e, estateexemption, r) are set
% in the DFJ2010 GAUSS master program wealth_singles/singlsav15_wave8.gau
% (lines 220-222): tauBeq=0.15, exBeq=6e5, mu_r=1/0.98-1. The C solver reads
% them from asstvec.fmt, which gmmcrit25.g builds as cfloor|tauBeq|exBeq|mu_r|sigma_r.
% [Earlier versions of this code used placeholders tau_e=0.05, estateexemption=1e6, r=0.05.]

%% Run infrastructure (diary, paths, output folders)
% Record the full run (incl. any error) to DFJ2010_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./DFJ2010_Diary.txt','file'), delete('./DFJ2010_Diary.txt'); end
diary ./DFJ2010_Diary.txt

close all % close any figures, make sure they are all cleanly built from scratch

% Subcodes (ReturnFn and the AfterTaxIncome/Consumption/GovTransfers/MedicalExpenses
% FnsToEvaluate) live in ./DFJ2010subcodes/; MATLAB does not search subfolders of
% the working directory, so put it on the path. Note this does NOT change cwd, so
% the './SavedOutput/...' paths still resolve to here.
addpath('./DFJ2010subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
if exist('./VFIToolkit-matlab','dir'); addpath(genpath('./VFIToolkit-matlab/')); end

% Make sure the subfolders to save output exist
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

Params.headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end

doPart=[0,0,0,0,0,0,0,0,0,0,0,0,0,1];
% doPart(1): data work (DFJ2010_DataWork: RAND HRS extract, GMM targets, initial distribution; writes SavedOutput/DFJ2010_DataWork.mat)
% doPart(2): solve the model and the descriptive figures (value fn, policy fn, the youngest cohort's life-cycle profiles)
% doPart(3): model vs data at the current parameters (each cohort from its 1996 distribution, the GMM objective, Figure 5 style plot, runtime estimate)
% doPart(4): GMM estimation
% doPart(5): counterfactuals at the current parameters (DFJ2010 Figures 9-11: no medical expenses, no medical expense risk, 80 percent consumption floor)
% doPart(6): the model with ENDOGENOUS medical expenses (DFJ2010 Sec. IX) at its published parameters, compared with the data (DFJ2010_EndogMedex1 and 2 in the subcodes; Figures 12 and A4); overwrites the exogenous model objects, so it runs last
% doPart(7): GMM estimation of the endogenous medical-expense model (nu, omega, beta, eMin) on its targets (asset medians, mean and 90th pct of 2-year medex, AC1 and AC2 of log 2-year medex), the autocorrelations in the TOOLKIT definition (AutoCorr targets; DFJ2010_EndogMedex3 with ACdefinition=toolkit; needs part 6 in the same run)
% doPart(8): the same estimation with the autocorrelations in DFJ2010's definition (CustomModelStats targets, DFJ2010_EM_CustomStats; DFJ2010_EndogMedex3 with ACdefinition=DFJ; needs part 6; prints the two side by side when part 7 ran too)
% doPart(9): the endogenous model with medical goods as a decision variable (n_d=101) vs the bisection split: runtime, GPU memory, solution (DFJ2010_EndogMedex4; needs part 6)
% doPart(10): the endogenous model's figures at the estimates (DFJ2010 Figures 12, 13, 14 and A4; DFJ2010_EndogMedex5): needs SavedOutput/DFJ2010_EndogMedex_Estimation_AC<DFJ|toolkit>.mat from parts 8/7, runs without parts 6-9 (Figure 14 panel A is the exogenous model)
% doPart(11): DFJ2010 Tables 1 and 2, Figure 8 panels A and B, appendix Figures A1 and A2 (DFJ2010_AppendixTablesFigs; writes SavedOutput/LatexInputs/DFJ2010_Table1.tex and Table2.tex): runs without parts 6-10; needs DFJ2010_DataWork.mat with Targets.MedexObs (data work of 2026-11-01)
% doPart(12): DFJ2010 Table A1 column 2: the benchmark re-estimated with the consumption floor fixed at $5,000 (upsilon and beta estimated; a few hours): runs without parts 6-11
% doPart(13): the benchmark (Table 3 column 1: upsilon, beta, cfloor) re-estimated from doPart(12)'s estimates with cfloor 5000, as doPart(12)'s objective is below doPart(4)'s: needs SavedOutput/DFJ2010_Estimation_TableA1column2.mat, runs without parts 6-11
% doPart(14): map of the benchmark objective over upsilon (3 to 7.5) and cfloor (1500 to 10000) at beta .965, one evaluation per point (about an hour): runs without parts 6-11
% Parts 2-4 need SavedOutput/DFJ2010_DataWork.mat (from part 1, or an earlier run of DFJ2010subcodes/DFJ2010_DataWork.m on its own); parts 3-4 solve the value fn themselves if part 2 is off.

%% Part 1: data work
if doPart(1)==1
    calledfromDFJmain=1; % tells DFJ2010_DataWork not to open its own diary (it can also be run on its own)
    DFJ2010_DataWork
    clear calledfromDFJmain
end

%% Setup the model action and state space
n_d=0;
DFJ2010_AssetGrid; % sets n_a and a_grid (subcode shared with DFJ2010_DataWork, which places the initial distribution on the grid)
n_z=[4,15]; % h (health status, 0=good, 1=bad, 2=death, 3=dead) and zeta (markov for medical expense shocks)
n_e=7; % xi (iid for medical expense shocks)
N_j=33; % 31 is 70 to 100
% Permanent types: 2 genders, 5 income groups (quintiles of permanent income)
Names_i={'femaleq1','femaleq2','femaleq3','femaleq4','femaleq5',...
    'maleq1','maleq2','maleq3','maleq4','maleq5'}; % q1=first (bottom) quintile, q5=fifth (top) quintile

% Note: DFJ2010 used 8 points for zeta and 8 points for xi (lines 68,69 of their wealth20_Singles.c)

% Near the beginning of DFJ2010 wealth20_Singles.c code, there is a comment
%    "- We need to check if the difference in indexes exists also for singles (as it did for couples)
%    Note the different indexes for males and females: for a male age t=1, 
%    means age=70,  for a female t=1 means age=67. Both 
%    live up to age 100 (t=31 for males and 34 for females)."
% Reading through the code it seems to use same index for age in various
% matrices regardless of gender, so pretty sure this was resolved and model
% is just ages 70 to 100.

figure_c=0; % counter for figures

%% Parameters
% Six are just initial guesses as they will be GMM estimated later:
% upsilon, delta, theta, k, cfloor, 

% Demographics
Params.agejshifter=69;
Params.J=N_j;
Params.agej=1:1:N_j;

% Preference parameters: which specification of DFJ2010 Table 3 to start from (and, in part 4, to estimate)
Table3column=[1,2,3,4]; % 1 benchmark (no bequest motive, delta=0); 2 health-dependent preferences; 3 bequests; 4 all
% A vector runs part 4 (GMM) for each column in turn, in one run; parts 2 and 3 use the first entry.
% The values are the GAUSS estimates in singlsav15_wave8.gau lines 330-337 (Table 3 rounds them), in the order there:
% GAUSS columns '+delta', '+delta+beq', '+beq', 'benchmark' are Table 3 columns 2, 4, 3, 1. cfloor is stored as its log.
% upsilon=CES utility param (nu in DFJ2010), beta=discount factor, delta=utility cost of bad health,
% theta=intensity of bequest motive, k=curvature of the warm-glow (bequests are a luxury), cfloor=consumption floor
Table3.upsilon=[3.81476, 3.74871,  3.83941,  3.6606];
Table3.beta=   [0.97106, 0.971026, 0.96536,  0.971094];
Table3.delta=  [0,       -0.206486, 0,       -0.360916];
Table3.phi0=   [1,       1,        0.119341, 0.108545]; % GAUSS keeps phi0, the 'MPC out of terminal wealth'; phi0=1 means no bequest motive
Table3.k=      [215000,  215000,   272975,   214583];
Table3.logcfloor=[7.88704, 7.8836, 7.8879,   7.88358];
Table3.EstimParamNames={{'upsilon','beta','cfloor'},{'upsilon','beta','delta','cfloor'},{'upsilon','beta','theta','k','cfloor'},{'upsilon','beta','delta','theta','k','cfloor'}};
Params.upsilon=Table3.upsilon(Table3column(1));
Params.beta=Table3.beta(Table3column(1));
Params.delta=Table3.delta(Table3column(1));
Params.k=Table3.k(Table3column(1));
Params.cfloor=exp(Table3.logcfloor(Table3column(1)));
phi0=Table3.phi0(Table3column(1));
EstimParamNames=Table3.EstimParamNames{Table3column(1)};
% theta is not stored directly: GAUSS keeps phi0 and gmmcrit25.g line 102 maps it as
% theta=((R*(1-phi0)/phi0)^upsilon)/(R*beta) with R=1/0.98 (this gives 2419 for column 4, as in Table 3)
if phi0==1
    Params.theta=0;
else
    Params.theta=(((1/0.98)*(1-phi0)/phi0)^Params.upsilon)/((1/0.98)*Params.beta);
end
% Note: Table 3 of DFJ2010 reports k "in 000s" (=215), and gmmcrit25.g line 100 does K0=parmvec[6]*1000, so k in dollars is 215,000.
% [Earlier versions of this code used k=215, which makes the warm-glow enormous at low estates.]
% Note: theta (phi0 in the codes) is not estimated directly; singlsav15_wave8.gau line 336 stores the 'MPC out of terminal
% wealth' (0.108545) and gmmcrit25.g line 102 maps it as phi0=((R*(1-phi0)/phi0)^nu)/(R*beta), which gives 2419 as in Table 3.
% The C code (wealth20_Singles.c line 1084) uses phi0*(netbequest+K0)^(1-nu)/(1-nu), with no '-1' in the numerator.

%% First-stage parameters (shared subcode): interest rate, taxes, health, survival, income, medical-expense process
DFJ2010_FirstStageParams

%% Grids and shock transitions (shared subcode)
DFJ2010_Grids

%% ReturnFn
DiscountFactorParamNames={'beta'}; 
% Note: conditional survival probabilities are handled by pi_h_J
% Specifically the transtions from h=0 & h=1 (good and bad health) to h=2 (death).

ReturnFn=@(aprime,a,h,zeta,xi,r,upsilon,delta,theta,k,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, tau_e, estateexemption, beta, agej, J, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)...
    DeNardiFrenchJones2010_ReturnFn(aprime,a,h,zeta,xi,r,upsilon,delta,theta,k,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, tau_e, estateexemption, beta, agej, J, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6);


%% Solve value fn and policy fn (parts 2-5 and 10-14)
if any(doPart([2:5,10:14])==1) % (parts 10-14 need the exogenous model: Figure 14 panel A, Figures 8 and A2 panel A, Table A1 column 2, the benchmark from the column-2 estimates, the objective map)
vfoptions.verbose=1;
vfoptions.verboseparams=1; % show params being sent to each ptype
vfoptions.divideandconquer=1;
% vfoptions.level1n=9; % you could reduce runtimes through playing with this a bit
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=10;
tic;
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
vftime=toc
% Note: the level of accuracy being used on assets is very high (501 points, plus 10
% interpolation points between each two consecutive grid points). 
% DFJ2010 appear to have used 110-to-235 points on cash-on-hand (very similar to using 
% them on assets), although their interpolation would not be limited to
% grid.

% When using grid interpolation layer, you have to tell simoptions
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
end


%% Part 2: descriptive figures of the value fn and policy fn
if doPart(2)==1
%% Quick look a value fn in final period
figure_c=figure_c+1;
figure(figure_c);
subplot(2,1,1); surf(squeeze(V.maleq1(:,:,1,1,end)))
% The warm-glow of bequests is so immense it is the only thing you can see
subplot(2,1,2); surf(squeeze(V.maleq1(:,1:2,1,1,end)))
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])
% Can see that in the first two 'health' states, utility is positive, and then becomes zero when dead

% NOTE: According to eqns in DFJ2010 paper, the warm-glow function phi(e) does not have a '-1' in numerator, while 
% utility fn does. Seems odd given their choice to use same upsilon for both. This left the warm-glow as being very 
% close to zero, so I have put the -1 in here so that there is actually a warm-glow worth keeping assets for.
% [This should be in their code somewhere, but I have not yet tried to find it.]

%% Plot the policy fn, take a look at it
PolicyVals=PolicyInd2Val_FHorz_PType(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,simoptions);
% The only policy is aprime
% First, for the middle (third) quintile, look at male and female at different ages
% row: male/female
% column: good health (1st index), bad health (2nd index)
figure_c=figure_c+1;
figure(figure_c);
subplot(2,2,1); plot(a_grid,PolicyVals.maleq3(1,:,1,8,4,1),a_grid,PolicyVals.maleq3(1,:,1,8,4,10),a_grid,PolicyVals.maleq3(1,:,1,8,4,20),a_grid,PolicyVals.maleq3(1,:,1,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: male, good health; mid values for med shocks')
subplot(2,2,2); plot(a_grid,PolicyVals.maleq3(1,:,2,8,4,1),a_grid,PolicyVals.maleq3(1,:,2,8,4,10),a_grid,PolicyVals.maleq3(1,:,2,8,4,20),a_grid,PolicyVals.maleq3(1,:,2,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: male, bad health; mid values for med shocks')
subplot(2,2,3); plot(a_grid,PolicyVals.femaleq3(1,:,1,8,4,1),a_grid,PolicyVals.femaleq3(1,:,1,8,4,10),a_grid,PolicyVals.femaleq3(1,:,1,8,4,20),a_grid,PolicyVals.femaleq3(1,:,1,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: female, good health; mid values for med shocks')
subplot(2,2,4); plot(a_grid,PolicyVals.femaleq3(1,:,2,8,4,1),a_grid,PolicyVals.femaleq3(1,:,2,8,4,10),a_grid,PolicyVals.femaleq3(1,:,2,8,4,20),a_grid,PolicyVals.femaleq3(1,:,2,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: female, bad health; mid values for med shocks')
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])

% DFJ2010 don't appear to report any policy functions, so not sure if these
% are reasonable or not.
% The fact there is a 'cliff' at low assets looks odd, but maybe this is
% the impact of the consumption floor? (I've not solved a consumption floor
% model before; you have to have aprime=0 to qualify for the consumption floor 
% transfers.) I doubt it, cliff is around $1million, and consumption floor is
% just 2653

end % of doPart(2) figures

%% Initial distribution, cohorts and quintile restrictions (shared subcode)
DFJ2010_InitialDistCohorts

%% Part 2: agent distribution, model statistics and descriptive figures
if doPart(2)==1
%% Agent distribution
% For the descriptive figures below, follow the youngest cohort (aged 72-76 in 1996, entering at 74) from its 1996
% distribution. The GMM section does the same for every cohort. Ages before entry have zero mass, so their statistics are NaN.
c_desc=1;
jequaloneDist_desc=struct();
cohortmass_desc=zeros(1,length(Names_i));
for ii=1:length(Names_i)
    temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c_desc)-Params.agejshifter);
    cohortmass_desc(ii)=gather(sum(temp(:)));
    if cohortmass_desc(ii)>0
        jequaloneDist_desc.(Names_i{ii})=temp/cohortmass_desc(ii);
    else
        jequaloneDist_desc.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
    end
end
Params_desc=Params;
Params_desc.ptypemasses=Params.ptypemasses(:)'.*cohortmass_desc;
Params_desc.ptypemasses=Params_desc.ptypemasses/sum(Params_desc.ptypemasses); % ptype weights within this cohort
Params_desc.mewj=Params.mewj(cohortentryage(c_desc)-Params.agejshifter)*ones(1,N_j); % the cohort keeps its entry-age weight
simoptions_desc=simoptions;
simoptions_desc.jequaloneDistAge=cohortentryage(c_desc)-Params.agejshifter;
simoptions_desc.agemass_withCohort=1; % non-empty: the age weights are a cohort mass, not expected to sum to one
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist_desc,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_desc,simoptions_desc);
disttime=toc

%% Calculate the model statistics that are relevant to the GMM estimation targets
FnsToEvaluate.assets=@(aprime,a,h,zeta,xi) a;
simoptions.conditionalrestrictions.alive=@(aprime,a,h,zeta,xi) (h==0) || (h==1);

% Note: only actually use the median in the estimation so internally the GMM command will set simoptions.whichstats=[0,1,0,0,0,0,0]
tic;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
statstime=toc

%% Calculate the various other model stats so we can create some graphs and better understand what the model is doing
FnsToEvaluate.assets=@(aprime,a,h,zeta,xi) a;
FnsToEvaluate.fractionliving=@(aprime,a,h,zeta,xi) (h==0) || (h==1);
FnsToEvaluate.goodhealth=@(aprime,a,h,zeta,xi) (h==0);
FnsToEvaluate.earnings=@(aprime,a,h,zeta,xi,earnings) earnings;
FnsToEvaluate.aftertaxincome=@(aprime,a,h,zeta,xi,r,earnings, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)...
    DeNardiFrenchJones2010_AfterTaxIncome(aprime,a,h,zeta,xi,r,earnings, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6);
FnsToEvaluate.medicalexpenses=@(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi)...
    DeNardiFrenchJones2010_MedicalExpenses(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi);
FnsToEvaluate.govtransfers=@(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)...
    DeNardiFrenchJones2010_GovTransfers(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6);
FnsToEvaluate.receivegovtransfers=@(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)...
   (DeNardiFrenchJones2010_GovTransfers(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)>0);
FnsToEvaluate.consumption=@(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)...
    DeNardiFrenchJones2010_Consumption(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6);
FnsToEvaluate.avgbequestsize=@(aprime,a,h,zeta,xi) a*(h==2); % note, this is in the period after it was left
FnsToEvaluate.avgbequestsize2=@(aprime,a,h,zeta,xi,agej,J) (agej==J)*((h==0) || (h==1))*aprime; % just the final period

simoptions.conditionalrestrictions.alive=@(aprime,a,h,zeta,xi) (h==0) || (h==1);
simoptions.conditionalrestrictions.death=@(aprime,a,h,zeta,xi) (h==2);
simoptions.conditionalrestrictions.ingoodhealth=@(aprime,a,h,zeta,xi) (h==0);
simoptions.conditionalrestrictions.inbadhealth=@(aprime,a,h,zeta,xi) (h==1);

tic;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
statstime2=toc

% DFJ2010 Figure 7: median consumption by permanent-income quintile, simulations (the cohort aged 74 in 1996, followed to 100)
figure_c=figure_c+1;
figure(figure_c);
hold on
for q=1:5
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.(['aliveq',num2str(q)]).consumption.Median/1000,'-','Color',[0 0 0.4+0.12*q],'LineWidth',1)
end
hold off
xlim([74,100]); xlabel('Age'); ylabel('Median consumption (thousands of 1998 dollars)')
title('Median consumption by PI quintile: model (cohort entering at 74)')
legend('q1','q2','q3','q4','q5','Location','northeast')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig7_consumption.png')

%% 
% Plot of assets by earnings quintile (DFJ2010 do a lot of versions of this)
figure_c=figure_c+1;
figure(figure_c);
subplot(1,2,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.maleq5.Mean)
hold off
title('Mean assets, male')
legend('income q1','income q2','income q3','income q4','income q5')
subplot(1,2,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.femaleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.assets.femaleq5.Mean)
hold off
title('Mean assets, female')
legend('income q1','income q2','income q3','income q4','income q5')
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])

% See how much of the assets at end of life is just because the people with
% bad health all died and mostly just good health left (as seen in previous figure).
figure_c=figure_c+1;
figure(figure_c);
subplot(1,2,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.maleq5.Mean)
hold off
title('Mean assets, male in good health')
legend('income q1','income q2','income q3','income q4','income q5')
subplot(1,2,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.femaleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.ingoodhealth.assets.femaleq5.Mean)
hold off
title('Mean assets, female in good health')
legend('income q1','income q2','income q3','income q4','income q5')
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])

% Look at some things relating to health and death
figure_c=figure_c+1;
figure(figure_c);
subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.fractionliving.Mean)
title('fraction of the youngest cohort (enters at 74) still alive')
subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.goodhealth.Mean)
title('fraction in good health')
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])

% Look at earnings, after tax income, medical expenses, gov transfers, and consumption
% Note: DFJ2010 show the medical expenses data in their Figure 3 (data, not model, so should be similar but not same)
% All use same legend, so only put it into the 4th
figure_c=figure_c+1;
figure(figure_c);
subplot(3,2,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.maleq5.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.femaleq1.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.earnings.femaleq5.Mean)
hold off
title('Earnings (Mean)')
% legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')
subplot(3,2,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.maleq5.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.femaleq1.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.aftertaxincome.femaleq5.Mean)
title('After-Tax Income (Mean)')
% legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')
subplot(3,2,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.maleq5.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.femaleq1.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.medicalexpenses.femaleq5.Mean)
title('Medical expenses (Mean)')
% legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')
subplot(3,2,4); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.maleq1.Mean./AgeConditionalStats.alive.receivegovtransfers.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.maleq2.Mean./AgeConditionalStats.alive.receivegovtransfers.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.maleq3.Mean./AgeConditionalStats.alive.receivegovtransfers.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.maleq4.Mean./AgeConditionalStats.alive.receivegovtransfers.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.maleq5.Mean./AgeConditionalStats.alive.receivegovtransfers.maleq5.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.femaleq1.Mean./AgeConditionalStats.alive.receivegovtransfers.femaleq1.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.femaleq2.Mean./AgeConditionalStats.alive.receivegovtransfers.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.femaleq3.Mean./AgeConditionalStats.alive.receivegovtransfers.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.femaleq4.Mean./AgeConditionalStats.alive.receivegovtransfers.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.govtransfers.femaleq5.Mean./AgeConditionalStats.alive.receivegovtransfers.femaleq5.Mean)
title('Gov Transfers (Mean, conditional on receiving them)')
legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')
subplot(3,2,5); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.maleq1.Mean)
hold on
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.maleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.maleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.maleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.maleq5.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.femaleq1.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.femaleq2.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.femaleq3.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.femaleq4.Mean)
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.alive.consumption.femaleq5.Mean)
title('Consumption (Mean)')
% legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')
subplot(3,2,6); plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.maleq1.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.maleq1.Mean(end)])
hold on % note, shift avgbequestsize one period, to when it was left, then include those left at end of final period
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.maleq2.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.maleq2.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.maleq3.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.maleq3.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.maleq4.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.maleq4.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.maleq5.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.maleq5.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.femaleq1.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.femaleq1.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.femaleq2.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.femaleq2.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.femaleq3.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.femaleq3.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.femaleq4.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.femaleq4.Mean(end)])
plot(Params.agejshifter+(1:1:N_j),[AgeConditionalStats.death.avgbequestsize.femaleq5.Mean(2:end),AgeConditionalStats.alive.avgbequestsize2.femaleq5.Mean(end)])
title('Bequest (left at end of period; mean among those who die)')
saveas(figure(figure_c),['./SavedOutput/Graphs/DFJ2010_Figure',num2str(figure_c),'.png'])
% legend('femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5','Location','northwest')




%% Runtimes summary
fprintf(' \n')
fprintf('When doing GMM, you have solve the model a lot of times, as a rule of thumb this has to take at most a few minutes each time \n')
fprintf('So lets look at runtimes: \n')
fprintf('  Value fn took:    %4.1f seconds \n',vftime)
fprintf('  Agent dist took:  %4.1f seconds \n',disttime)
fprintf('  Model stats took: %4.1f seconds \n',statstime)
fprintf('So total runtime was: %4.1f seconds \n', vftime+disttime+statstime)
fprintf('  (The full set of model stats used for the graphs took: %4.1f seconds) \n',statstime2)
fprintf(' \n')
% Even just on laptop it is easy fast enough you can GMM estimate it.

% Note: currently calculating way more model stats that are needed in the
% estimation. If you strip it back to just those for estimation the runtime
% for those stats falls substantially.

%% Save the descriptive results
save ./SavedOutput/DFJ2010_Results.mat Params AgeConditionalStats vftime disttime statstime statstime2
end % of doPart(2)

%% GMM targets: median assets by cohort, quintile and year (shared subcode)
DFJ2010_AssetTargets

%% Part 3: model vs data at the current parameter values: each cohort's distribution from its entry age (as the estimator does it),
if doPart(3)==1
% model vs data medians, the GMM objective, and a figure in the style of DFJ2010 Figure 5
modelmedian=nan(Jdata,1);
Params_c=Params;
simoptions_c=simoptions_gmm;
simoptions_c.agemass_withCohort=1;
% For DFJ2010 Figures 6 and 8: median assets by cohort among everyone alive and among the survivors to 2006 (all quintiles
% pooled), and the distribution of assets in the period of death (deaths between 1996 and 2006, all cohorts pooled)
Fig6_median=nan(ncohorts,N_j); Fig6_median_survivors=nan(ncohorts,N_j);
deathdist=zeros(n_a,1);
cohortshare=zeros(1,ncohorts);
tic; % time the cohort loop: together with the value fn solve this is what one objective evaluation costs
for c=1:ncohorts
    jequaloneDist_c=struct();
    cohortmass_c=zeros(1,length(Names_i));
    for ii=1:length(Names_i)
        temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c)-Params.agejshifter);
        cohortmass_c(ii)=gather(sum(temp(:)));
        if cohortmass_c(ii)>0
            jequaloneDist_c.(Names_i{ii})=temp/cohortmass_c(ii);
        else
            jequaloneDist_c.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
        end
    end
    Params_c.ptypemasses=Params.ptypemasses(:)'.*cohortmass_c;
    Params_c.ptypemasses=Params_c.ptypemasses/sum(Params_c.ptypemasses);
    Params_c.mewj=Params.mewj(cohortentryage(c)-Params.agejshifter)*ones(1,N_j);
    simoptions_c.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
    StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_c);
    ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy,FnsToEvaluate_gmm,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_c);
    % (assets, health, age) marginal of this cohort's distribution, ptypes pooled with their within-cohort weights; and the
    % same weighted by the probability of still being alive in 2006 (from the health-and-survival chain pi_h_J, which does
    % not depend on assets), which gives the survivors' distribution at every earlier age
    jentry=cohortentryage(c)-Params.agejshifter;
    jend=min(jentry+10,N_j); % 2006
    cohortshare(c)=sum(Params.ptypemasses(:)'.*cohortmass_c);
    Dh_all=zeros(n_a,n_z(1),N_j); Dh_surv=zeros(n_a,n_z(1),N_j);
    for ii=1:length(Names_i)
        Dh=gather(reshape(sum(reshape(StationaryDist_c.(Names_i{ii}),[n_a,n_z(1),n_z(2)*n_e,N_j]),3),[n_a,n_z(1),N_j]));
        S=zeros(n_z(1),N_j); S(1:2,jend)=1; % prob of being alive at jend, given health h at age j
        for jj=jend-1:-1:jentry
            S(:,jj)=pi_h_J.(Names_i{ii})(:,:,jj)*S(:,jj+1);
        end
        Dh_all=Dh_all+Params_c.ptypemasses(ii)*Dh;
        Dh_surv=Dh_surv+Params_c.ptypemasses(ii)*Dh.*reshape(S,[1,n_z(1),N_j]);
    end
    for jj=jentry:jend
        w=sum(Dh_all(:,1:2,jj),2); cw=cumsum(w); Fig6_median(c,jj)=a_grid(find(cw>=0.5*cw(end),1));
        w=sum(Dh_surv(:,1:2,jj),2); cw=cumsum(w); Fig6_median_survivors(c,jj)=a_grid(find(cw>=0.5*cw(end),1));
    end
    deathdist=deathdist+cohortshare(c)*sum(Dh_all(:,3,jentry+1:jend),3); % h=2 (index 3) is the period of death; assets then are those held one period before
    rows=find(MomentList.cohort==c);
    for rr=rows'
        modelmedian(rr)=ACS_c.(['aliveq',num2str(MomentList.PIq(rr))]).assets.Median(cohortentryage(c)+MomentList.year(rr)-1996-Params.agejshifter);
    end
    % Independent check of the toolkit's restricted grouped medians: the weighted median of assets among the alive of
    % quintile q (female and male ptypes pooled with their within-cohort weights) computed by hand from the (a,h) marginals
    % [added 2026-09-29 after toolkit commits 617cf5a1/ed156e23 changed the median definition and the PType restricted stats]
    Dh_q=zeros(n_a,n_z(1),N_j,5);
    for ii=1:length(Names_i)
        Dh=gather(reshape(sum(reshape(StationaryDist_c.(Names_i{ii}),[n_a,n_z(1),n_z(2)*n_e,N_j]),3),[n_a,n_z(1),N_j]));
        q=Params.quintile.(Names_i{ii});
        Dh_q(:,:,:,q)=Dh_q(:,:,:,q)+Params_c.ptypemasses(ii)*Dh;
    end
    rows=find(MomentList.cohort==c);
    maxdiff_c=0;
    for rr=rows'
        jj=cohortentryage(c)+MomentList.year(rr)-1996-Params.agejshifter;
        w=sum(Dh_q(:,1:2,jj,MomentList.PIq(rr)),2); cw=cumsum(w);
        byhand=a_grid(find(cw>=0.5*cw(end),1));
        maxdiff_c=max(maxdiff_c,abs(byhand-modelmedian(rr)));
    end
    fprintf('   cohort %d: max |by-hand weighted median - toolkit aliveq median| over its %d cells = %.1f dollars (a_grid step there is 140-750)\n',c,numel(rows),maxdiff_c)
end
cohortlooptime=toc;
fprintf('\nRuntime of one GMM objective evaluation: value fn %.1f s + five cohort distributions and their medians %.1f s = %.1f s \n',vftime,cohortlooptime,vftime+cohortlooptime)
nEstim=length(EstimParamNames);
fprintf('   (lsqnonlin with %d parameters needs about %d evaluations per iteration, so roughly %.1f minutes per iteration; the std errors afterwards need a further %d evaluations, about %.1f minutes) \n',nEstim,nEstim+1,(nEstim+1)*(vftime+cohortlooptime)/60,1+4*2*nEstim,(1+4*2*nEstim)*(vftime+cohortlooptime)/60)
fprintf('\nModel vs data median assets by cohort, PI quintile and year (thousands of 1998 dollars), current parameters (Table 3 column %d: upsilon %.3f, beta %.4f, delta %.3f, theta %.0f, k %.0f, cfloor %.0f):\n',Table3column(1),Params.upsilon,Params.beta,Params.delta,Params.theta,Params.k,Params.cfloor)
fprintf('   cohort  q  year  age     data    model   (data se)     n\n')
for rr=1:Jdata
    fprintf('     %d     %d  %d  %3d   %6.0f   %6.0f    (%5.0f)  %4d\n',MomentList.cohort(rr),MomentList.PIq(rr),MomentList.year(rr),cohortentryage(MomentList.cohort(rr))+MomentList.year(rr)-1996,MomentList.median(rr)/1000,modelmedian(rr)/1000,MomentList.se(rr)/1000,MomentList.n(rr))
end
gap=modelmedian(tkorder)-MomentList.median(tkorder);
if exist('chi2inv','file')
    fprintf('GMM objective at current parameters (diagonal W): %.1f  [chi2(%d) 5%% critical value %.0f]\n',gap'*WeightingMatrix*gap,Jdata-length(EstimParamNames),chi2inv(0.95,Jdata-length(EstimParamNames)))
else
    fprintf('GMM objective at current parameters (diagonal W): %.1f  [%d degrees of freedom]\n',gap'*WeightingMatrix*gap,Jdata-length(EstimParamNames))
end

% DFJ2010 Figure 5: median assets by PI quintile for the cohorts aged 74 and 84 in 1996 (cohorts 1 and 3), data (solid) and
% model (dashed), 1996-2006. In 1996 the model's distribution is the data's (the initial distribution), so both lines start
% from the 1996 data median; 1998-2006 are the targets.
figure_c=figure_c+1;
figure(figure_c);
hold on
for c=[1,3]
    for q=1:5
        rows=find(MomentList.cohort==c & MomentList.PIq==q);
        [~,yy]=sort(MomentList.year(rows)); rows=rows(yy);
        ages=[cohortentryage(c); cohortentryage(c)+MomentList.year(rows(:))-1996];
        plot(ages,[DataWork.Targets.MedianAssets(c,q,1); MomentList.median(rows(:))]/1000,'k-',ages,[DataWork.Targets.MedianAssets(c,q,1); modelmedian(rows(:))]/1000,'k--')
    end
end
hold off
xlim([70,102]); xlabel('Age'); ylabel('Assets (000s of 1998 dollars)')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig5_fit.png')

% DFJ2010 Figure 6: median assets by birth cohort, everyone alive (solid) vs survivors to 2006 (dashed), model
figure_c=figure_c+1;
figure(figure_c);
hold on
for c=1:ncohorts
    jj=find(~isnan(Fig6_median(c,:)));
    plot(Params.agejshifter+jj,Fig6_median(c,jj)/1000,'-','Color',[0 0 0.8],'LineWidth',1)
    plot(Params.agejshifter+jj,Fig6_median_survivors(c,jj)/1000,'--','Color',[0.8 0 0],'LineWidth',1)
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Median assets by birth cohort: model, everyone alive (solid) vs survivors to 2006 (dashed)')
legend({'everyone alive','survivors to 2006'},'Location','northeast')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig6_survivors.png')
fprintf('\nDFJ2010 Figure 6: cohort entering at 74, median assets at entry: everyone %.0fk, survivors to 2006 %.0fk; in 2006: %.0fk (DFJ2010 data: 60k, 84k, 44k)\n',Fig6_median(1,cohortentryage(1)-Params.agejshifter)/1000,Fig6_median_survivors(1,cohortentryage(1)-Params.agejshifter)/1000,Fig6_median_survivors(1,cohortentryage(1)-Params.agejshifter+10)/1000)

% DFJ2010 Figure 8: CDF of assets held one period before death, data (solid) and model (dashed)
AssetsBeforeDeath=DataWork.Targets.AssetsBeforeDeath;
[abd_sorted]=sort(AssetsBeforeDeath);
figure_c=figure_c+1;
figure(figure_c);
plot(abd_sorted/1000,(1:numel(abd_sorted))/numel(abd_sorted),'k-','LineWidth',1)
hold on
plot(a_grid/1000,cumsum(deathdist)/sum(deathdist),'b--','LineWidth',1)
hold off
xlim([0,500]); ylim([0,1]); xlabel('Assets (thousands of 1998 dollars)'); ylabel('CDF')
title('Assets held one period before death: data (solid) and model (dashed)')
legend({'data (HRS)','model'},'Location','southeast')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig8_assetsbeforedeath.png')
cdfmodel=cumsum(deathdist)/sum(deathdist);
fprintf('Assets one period before death: data quartiles %s thousand, model quartiles %s thousand; share with a<=0 (data) or at the zero grid point (model): %.3f vs %.3f\n',mat2str(round(quantile(AssetsBeforeDeath,[0.25,0.5,0.75])/1000)),mat2str(round(a_grid(arrayfun(@(p) find(cdfmodel>=p,1),[0.25,0.5,0.75]))'/1000)),mean(AssetsBeforeDeath<=0),cdfmodel(1))

save ./SavedOutput/DFJ2010_Validation.mat MomentList modelmedian cohortentryage tkorder cohortlooptime Fig6_median Fig6_median_survivors deathdist cohortshare
end % of doPart(3)

%% Part 4: GMM estimation
% EstimParamNames follows Table3column (set with the preference parameters above): the benchmark estimates
% upsilon, beta and cfloor with delta=0 and theta=0 fixed; the other columns add delta and/or theta and k.
estimoptions.verbose=1;
% Standard errors: the standard GMM sandwich formula, with the diagonal W and the full bootstrap covariance of the medians (CoVarMatrixDataMoments) (estimoptions.efficientW=0 default)
estimoptions.cohortagejshifter=1; % jequaloneDist has an N_j dimension: the entry ages are the slices with mass (the five cohorts)
% Finite-difference step for the Jacobian used in the standard errors: the moments are medians of assets, which sit on the
% 501-point asset grid, so for the default step (eedefault=3, about 1e-4 of the parameter) 94-97 percent of the moments have
% exactly zero derivative and the rest jump by a whole grid step, giving meaningless (tiny) standard errors. The largest step
% (eedefault=4, about 1e-2) gives std errors of 0.36, 0.034 and 341 for the benchmark, in line with DFJ2010's 0.50, 0.04, 346.
% [Found from estsummary.doublechecks of the 2026-09-19 run; the point estimates are unaffected.]
estimoptions.eedefault=4;

% Each objective evaluation is one solve of the model plus five cohort distributions and their medians (see the runtime estimate printed by part 3)
% One column per entry of Table3column, each started from DFJ's published estimates for that column and saved to its own file.
if doPart(4)==1
    vfoptions_gmm=vfoptions;
    vfoptions_gmm.verbose=0; % the value fn is solved at every objective evaluation; keep the diary readable
    vfoptions_gmm.verboseparams=0;
    Table3results=struct(); % summary across the columns run
    for cc=1:length(Table3column)
        col=Table3column(cc);
        Params.upsilon=Table3.upsilon(col);
        Params.beta=Table3.beta(col);
        Params.delta=Table3.delta(col);
        Params.k=Table3.k(col);
        Params.cfloor=exp(Table3.logcfloor(col));
        phi0=Table3.phi0(col);
        if phi0==1
            Params.theta=0;
        else
            Params.theta=(((1/0.98)*(1-phi0)/phi0)^Params.upsilon)/((1/0.98)*Params.beta);
        end
        EstimParamNames=Table3.EstimParamNames{col};
        estimoptions.constrainpositive=intersect({'upsilon','theta','k','cfloor'},EstimParamNames,'stable');
        estimoptions.constrain0to1=intersect({'beta'},EstimParamNames,'stable');
        fprintf('\nGMM estimation: Table 3 column %d, estimating %s \n',col,strjoin(EstimParamNames,', '))
        fprintf('   starting from DFJ2010 estimates: upsilon %.4f, beta %.4f, delta %.4f, theta %.1f, k %.0f, cfloor %.1f \n',Params.upsilon,Params.beta,Params.delta,Params.theta,Params.k,Params.cfloor)
        tic;
        [EstimParams,EstimParamsConfInts,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TargetMoments,WeightingMatrix,CoVarMatrixDataMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_gmm,estimoptions,vfoptions_gmm,simoptions_gmm);
        estimtime=toc;
        fprintf('GMM estimation took %.1f minutes \n',estimtime/60)
        fprintf('Estimated parameters (Table 3 column %d), standard errors in brackets:\n',col)
        for pp=1:length(EstimParamNames)
            fprintf('   %-8s %12.4f  (%.4f) \n',EstimParamNames{pp},EstimParams.(EstimParamNames{pp}),estsummary.EstimParamsStdDev.(EstimParamNames{pp}))
        end
        fprintf('Objective at the estimates: %.2f  [%d moments, %d parameters]\n',estsummary.objectivefnval,Jdata,length(EstimParamNames))
        Table3column_this=col;
        save(['./SavedOutput/DFJ2010_Estimation_column',num2str(col),'.mat'],'EstimParams','EstimParamsConfInts','estsummary','Table3column_this','EstimParamNames','estimtime')
        Table3results(cc).column=col;
        Table3results(cc).EstimParams=EstimParams;
        Table3results(cc).StdDev=estsummary.EstimParamsStdDev;
        Table3results(cc).objective=estsummary.objectivefnval;
        Table3results(cc).minutes=estimtime/60;
    end
    % Summary of the columns run (Table 3 layout: rows are parameters, columns are specifications)
    fprintf('\nSummary of the GMM estimates (standard errors in brackets), Table 3 columns %s \n',mat2str(Table3column))
    allnames={'upsilon','beta','delta','theta','k','cfloor'};
    fprintf('   %-8s',''); for cc=1:length(Table3column), fprintf('  %26s',['column ',num2str(Table3column(cc))]); end; fprintf('\n')
    for pp=1:length(allnames)
        fprintf('   %-8s',allnames{pp})
        for cc=1:length(Table3column)
            if isfield(Table3results(cc).EstimParams,allnames{pp})
                fprintf('  %14.4f (%9.4f)',Table3results(cc).EstimParams.(allnames{pp}),Table3results(cc).StdDev.(allnames{pp}))
            else
                fprintf('  %26s','fixed')
            end
        end
        fprintf('\n')
    end
    fprintf('   %-8s',''); for cc=1:length(Table3column), fprintf('  %26s',sprintf('obj %.1f, %.0f min',Table3results(cc).objective,Table3results(cc).minutes)); end; fprintf('\n')
    save ./SavedOutput/DFJ2010_Estimation.mat Table3results Table3column
    % Put the preference parameters back to the FIRST column's published values, so that parts 5 and 6 (which use
    % Params) run at the same specification as parts 2 and 3 when several parts run together
    Params.upsilon=Table3.upsilon(Table3column(1));
    Params.beta=Table3.beta(Table3column(1));
    Params.delta=Table3.delta(Table3column(1));
    Params.k=Table3.k(Table3column(1));
    Params.cfloor=exp(Table3.logcfloor(Table3column(1)));
    phi0=Table3.phi0(Table3column(1));
    if phi0==1
        Params.theta=0;
    else
        Params.theta=(((1/0.98)*(1-phi0)/phi0)^Params.upsilon)/((1/0.98)*Params.beta);
    end
    EstimParamNames=Table3.EstimParamNames{Table3column(1)};
end

%% Part 5: counterfactuals (DFJ2010 Figures 9-11), cohort aged 74 in 1996 followed to age 100 with no attrition
% DFJ2010 Sec. VIII: fix the parameters, change one feature of the model, and compare median assets by PI quintile for the
% cohort with average age 74 in 1996. 'We rule out attrition and assume that every individual lives to age 100. People
% will face the same risks and have the same expectations as before, but in their realized lives they do not die.' So the
% policy is solved with mortality risk, but the distribution is iterated with the health transitions conditional on survival.
if doPart(5)==1
pi_z_J_noattrition=pi_z_J;
for ii=1:length(Names_i)
    pih=pi_h_J.(Names_i{ii});
    for jj=1:N_j
        for hh=1:2
            pih(hh,:,jj)=[pih(hh,1:2,jj)/sum(pih(hh,1:2,jj)),0,0];
        end
        pi_z_J_noattrition.(Names_i{ii})(:,:,jj)=kron(pi_zeta,pih(:,:,jj));
    end
end
% the cohort entering at 74 (as in the descriptive part 2)
c_cf=1;
jequaloneDist_cf=struct();
cohortmass_cf=zeros(1,length(Names_i));
for ii=1:length(Names_i)
    temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c_cf)-Params.agejshifter);
    cohortmass_cf(ii)=gather(sum(temp(:)));
    if cohortmass_cf(ii)>0
        jequaloneDist_cf.(Names_i{ii})=temp/cohortmass_cf(ii);
    else
        jequaloneDist_cf.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
    end
end
simoptions_cf=simoptions_gmm; % aliveq1..5 restrictions, median only
simoptions_cf.jequaloneDistAge=cohortentryage(c_cf)-Params.agejshifter;
simoptions_cf.agemass_withCohort=1;
vfoptions_cf=vfoptions;
vfoptions_cf.verbose=0;
vfoptions_cf.verboseparams=0;
cfnames={'baseline','nomedex','nomedexrisk','cfloor80'};
CF=struct();
tic;
for cc=1:length(cfnames)
    Params_cf=Params;
    Params_cf.ptypemasses=Params.ptypemasses(:)'.*cohortmass_cf;
    Params_cf.ptypemasses=Params_cf.ptypemasses/sum(Params_cf.ptypemasses);
    Params_cf.mewj=Params.mewj(cohortentryage(c_cf)-Params.agejshifter)*ones(1,N_j);
    switch cfnames{cc}
        case 'nomedex' % zero out all out-of-pocket medical expenses (m=exp(m_coeff+sigma_coeff*psi), so m_coeff=-100 gives m=0)
            for ii=1:length(Names_i)
                Params_cf.m_coeff_healthgood.(Names_i{ii})=-100*ones(size(Params.m_coeff_healthgood.(Names_i{ii})));
                Params_cf.m_coeff_healthbad.(Names_i{ii})=-100*ones(size(Params.m_coeff_healthbad.(Names_i{ii})));
                Params_cf.sigma_coeff_healthgood.(Names_i{ii})=zeros(size(Params.sigma_coeff_healthgood.(Names_i{ii})));
                Params_cf.sigma_coeff_healthbad.(Names_i{ii})=zeros(size(Params.sigma_coeff_healthbad.(Names_i{ii})));
            end
        case 'nomedexrisk' % no shocks, but the same mean medical expense conditional on age, ptype and health: E[exp(sigma*psi)]=exp(sigma^2/2)
            for ii=1:length(Names_i)
                Params_cf.m_coeff_healthgood.(Names_i{ii})=Params.m_coeff_healthgood.(Names_i{ii})+0.5*(Params.sigma_coeff_healthgood.(Names_i{ii})/Params.normalizepsi).^2;
                Params_cf.m_coeff_healthbad.(Names_i{ii})=Params.m_coeff_healthbad.(Names_i{ii})+0.5*(Params.sigma_coeff_healthbad.(Names_i{ii})/Params.normalizepsi).^2;
                Params_cf.sigma_coeff_healthgood.(Names_i{ii})=zeros(size(Params.sigma_coeff_healthgood.(Names_i{ii})));
                Params_cf.sigma_coeff_healthbad.(Names_i{ii})=zeros(size(Params.sigma_coeff_healthbad.(Names_i{ii})));
            end
        case 'cfloor80' % consumption floor at 80 percent of its value
            Params_cf.cfloor=0.8*Params.cfloor;
    end
    [~,Policy_cf]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_cf,DiscountFactorParamNames,vfoptions_cf);
    StationaryDist_cf=StationaryDist_Case1_FHorz_PType(jequaloneDist_cf,AgeWeightsParamNames,PTypeDistParamNames,Policy_cf,n_d,n_a,n_z,N_j,Names_i,pi_z_J_noattrition,Params_cf,simoptions_cf);
    ACS_cf=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_cf,Policy_cf,FnsToEvaluate_gmm,Params_cf,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cf);
    CF.(cfnames{cc})=zeros(5,N_j);
    for q=1:5
        CF.(cfnames{cc})(q,:)=ACS_cf.(['aliveq',num2str(q)]).assets.Median;
    end
    if strcmp(cfnames{cc},'baseline')
        % DFJ2010 Figures 3 and 4: average medical expenses and average (non-asset) income by PI quintile, simulated for this
        % cohort with no attrition from the first-stage profiles (paper Sec. VI.B); they depend on the first-stage parameters,
        % the 1996 distribution and the health dynamics only
        FnsToEvaluate_34.medicalexpenses=@(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi)...
            DeNardiFrenchJones2010_MedicalExpenses(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi);
        FnsToEvaluate_34.income=@(aprime,a,h,zeta,xi,earnings) earnings;
        simoptions_34=simoptions_cf;
        simoptions_34.whichstats=[1,0,0,0,0,0,0]; % mean only
        ACS_34=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_cf,Policy_cf,FnsToEvaluate_34,Params_cf,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_34);
        Fig34.medicalexpenses=zeros(5,N_j); Fig34.income=zeros(5,N_j);
        for q=1:5
            Fig34.medicalexpenses(q,:)=ACS_34.(['aliveq',num2str(q)]).medicalexpenses.Mean;
            Fig34.income(q,:)=ACS_34.(['aliveq',num2str(q)]).income.Mean;
        end
        fig34names={'medicalexpenses','income'}; fig34titles={'Average medical expenses by PI quintile (cohort aged 74 in 1996, no attrition)','Average income by PI quintile (cohort aged 74 in 1996, no attrition)'};
        fig34files={'DFJ2010_Fig3_medex.png','DFJ2010_Fig4_income.png'}; fig34ylabels={'Mean out-of-pocket medical expenses (thousands of 1998 dollars)','Mean non-asset income (thousands of 1998 dollars)'};
        for ff=1:2
            figure_c=figure_c+1;
            figure(figure_c);
            hold on
            for q=1:5
                plot(74:100,Fig34.(fig34names{ff})(q,(74:100)-Params.agejshifter)/1000,'-','Color',[0 0 0.4+0.12*q],'LineWidth',1)
            end
            hold off
            xlabel('Age'); ylabel(fig34ylabels{ff}); title(fig34titles{ff})
            legend('q1','q2','q3','q4','q5','Location','northwest')
            saveas(figure(figure_c),['./SavedOutput/Graphs/',fig34files{ff}])
        end
        fprintf('DFJ2010 Figures 3 and 4: mean medical expenses at 74/90/100 by quintile (thousands): %s; mean income: %s (DFJ2010: top quintile income about 20k)\n',mat2str(round(Fig34.medicalexpenses(:,[74,90,100]-Params.agejshifter)/1000,1)),mat2str(round(Fig34.income(:,[74,90,100]-Params.agejshifter)/1000,1)))
    end
    fprintf('Counterfactual %-12s solved (%.0f s so far): median assets at ages 74/90/100 by quintile (thousands): %s\n',cfnames{cc},toc,mat2str(round(CF.(cfnames{cc})(:,[74,90,100]-Params.agejshifter)/1000)))
end
cftime=toc;
% Figures 9-11: baseline (dashed) vs counterfactual (solid), five quintiles, ages 74-100
ages_cf=74:100;
jj_cf=ages_cf-Params.agejshifter;
cftitles={'','Figure 9: no medical expenses','Figure 10: no medical expense risk','Figure 11: 80 percent of the consumption floor'};
cffiles={'','DFJ2010_Fig9_nomedex.png','DFJ2010_Fig10_nomedexrisk.png','DFJ2010_Fig11_cfloor80.png'};
for cc=2:4
    figure_c=figure_c+1;
    figure(figure_c);
    hold on
    for q=1:5
        plot(ages_cf,CF.baseline(q,jj_cf)/1000,'--','Color',[0 0 0.4+0.12*q],'LineWidth',1)
        plot(ages_cf,CF.(cfnames{cc})(q,jj_cf)/1000,'-','Color',[0 0 0.4+0.12*q],'LineWidth',1)
    end
    hold off
    xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
    title(['Median assets by PI quintile, cohort aged 74 in 1996: baseline (dashed) vs ',cftitles{cc}])
    legend({'q1 baseline','q1','q2 baseline','q2','q3 baseline','q3','q4 baseline','q4','q5 baseline','q5'},'Location','northeast','NumColumns',1)
    saveas(figure(figure_c),['./SavedOutput/Graphs/',cffiles{cc}])
end
fprintf('\nCounterfactuals took %.0f seconds. DFJ2010: with no medical expenses the top quintile depletes its net worth by 94 (baseline keeps ~40k at 100);\n',cftime)
fprintf('   no medical-expense risk has only a small effect; the 80 percent floor raises 90-year-olds in q5 from 100k to 120k and q4 from 40k to 50k\n')
fprintf('   here: top quintile at 94: baseline %.0fk, no medex %.0fk, no risk %.0fk; at 100: %.0fk / %.0fk / %.0fk\n',CF.baseline(5,94-Params.agejshifter)/1000,CF.nomedex(5,94-Params.agejshifter)/1000,CF.nomedexrisk(5,94-Params.agejshifter)/1000,CF.baseline(5,100-Params.agejshifter)/1000,CF.nomedex(5,100-Params.agejshifter)/1000,CF.nomedexrisk(5,100-Params.agejshifter)/1000)
fprintf('   q5 and q4 at 90: baseline %.0fk / %.0fk, floor at 80 percent %.0fk / %.0fk\n',CF.baseline(5,90-Params.agejshifter)/1000,CF.baseline(4,90-Params.agejshifter)/1000,CF.cfloor80(5,90-Params.agejshifter)/1000,CF.cfloor80(4,90-Params.agejshifter)/1000)
save ./SavedOutput/DFJ2010_Counterfactuals.mat CF cfnames cftime cohortentryage Fig34
end % of doPart(5)

%% Part 6: the model with endogenous medical expenses (DFJ2010 Sec. IX), stage 1: published parameters vs data
if doPart(6)==1
    DFJ2010_EndogMedex1 % parameters, grids, return function, value function (overwrites the exogenous model's)
    DFJ2010_EndogMedex2 % initial distribution and cohorts, targets (DFJ2010_EMTargets), model vs data (Figures 12 and A4), autocorrelations, the full objective
end

%% Part 7: GMM estimation of the endogenous medical-expense model (stage 4 of the proposal), autocorrelation targets in the toolkit definition
if doPart(7)==1
    if doPart(6)==0
        error('doPart(7) (estimation of the endogenous model) needs doPart(6) in the same run: it uses the model, cohorts and targets that part 6 sets up')
    end
    ACdefinition='toolkit';
    DFJ2010_EndogMedex3
end

%% Part 8: the same estimation with the autocorrelation targets in DFJ2010's definition (CustomModelStats)
if doPart(8)==1
    if doPart(6)==0
        error('doPart(8) (estimation of the endogenous model, DFJ autocorrelation definition) needs doPart(6) in the same run: it uses the model, cohorts and targets that part 6 sets up')
    end
    ACdefinition='DFJ';
    DFJ2010_EndogMedex3
end

%% Part 9: the endogenous model with medical goods as a decision variable on a grid, vs the bisection split
if doPart(9)==1
    if doPart(6)==0
        error('doPart(9) needs doPart(6) in the same run: it compares against the bisection solve of the endogenous model')
    end
    DFJ2010_EndogMedex4
end

%% Part 10: the endogenous medical-expense model's figures at the estimates (DFJ2010 Figures 12, 13, 14 and A4)
if doPart(10)==1
    if any(doPart(6:9)==1)
        error('doPart(10) needs the exogenous model in the workspace (Figure 14 panel A), which parts 6-9 overwrite: run it without them')
    end
    EMfigdefinition='DFJ'; % the estimates to use: 'DFJ' (doPart(8), autocorrelations in DFJ2010's definition) or 'toolkit' (doPart(7))
    DFJ2010_EndogMedex5
end

%% Part 11: DFJ2010 Tables 1 and 2, Figure 8 (both panels), and appendix Figures A1 and A2
if doPart(11)==1
    if any(doPart(6:10)==1)
        error('doPart(11) needs the exogenous model in the workspace (Figures 8 and A2 panel A), which parts 6-10 overwrite: run it without them')
    end
    EMfigdefinition='DFJ'; % Figure A2 panel B at these estimates: 'DFJ' (doPart(8)) or 'toolkit' (doPart(7))
    DFJ2010_AppendixTablesFigs
end

%% Part 12: DFJ2010 Table A1 column 2: the benchmark with the consumption floor fixed at $5,000
% DFJ2010 online appendix C: re-estimating the benchmark (Table 3 column 1: no bequest motive, delta=0) with the floor set
% exogenously to $5,000 tests whether risk aversion and the floor are separately identified. Their estimates: upsilon 6.04
% (1.01), beta .96 (.06), objective 110.2 on 99 df against 82.3 on 98 for the benchmark; the difference of 17.9 exceeds the
% chi-squared(1) critical value of 3.8, so the restriction is rejected and the floor is identified. Same estimator, targets,
% weighting and options as doPart(4); started from DFJ2010's column-2 estimates (rounded, as the appendix reports them).
if doPart(12)==1
    if any(doPart(6:11)==1)
        error('doPart(12) needs the exogenous model in the workspace, which parts 6-11 overwrite: run it without them')
    end
    vfoptions_gmm=vfoptions;
    vfoptions_gmm.verbose=0;
    vfoptions_gmm.verboseparams=0;
    Params_A1c2=Params; % the benchmark specification (Table 3 column 1, set above): delta=0, no bequest motive
    Params_A1c2.upsilon=6.04; % DFJ2010 Table A1 column 2
    Params_A1c2.beta=0.96;
    Params_A1c2.delta=Table3.delta(1);
    Params_A1c2.k=Table3.k(1);
    Params_A1c2.theta=0;
    Params_A1c2.cfloor=5000; % fixed, not estimated
    EstimParamNames_A1c2={'upsilon','beta'};
    estimoptions_A1c2=estimoptions; % as doPart(4): cohorts, eedefault=4, verbose
    estimoptions_A1c2.constrainpositive={'upsilon'};
    estimoptions_A1c2.constrain0to1={'beta'};
    % Diagnostics (2026-11-03): the first run of this part gave model medians of exactly zero in every cell at every parameter
    % value lsqnonlin tried (objective 1391.3 throughout). One evaluation each (skipestimation=1, skipjacobian=1) separates the
    % floor from risk aversion, and the benchmark at DFJ2010's published values checks this route against the write-up (98.1).
    diagcases={'benchmark at DFJ2010 Table 3 column 1',Table3.upsilon(1),Table3.beta(1),exp(Table3.logcfloor(1)); ...
               'column-2 start: upsilon 6.04, beta .96, cfloor 5000',6.04,0.96,5000; ...
               'benchmark upsilon and beta, cfloor 5000',Table3.upsilon(1),Table3.beta(1),5000; ...
               'upsilon 6.04 and beta .96, benchmark cfloor',6.04,0.96,exp(Table3.logcfloor(1))};
    fprintf('\nTable A1 column 2 diagnostics: the objective and the model medians at four points (one evaluation each)\n')
    for dd=1:size(diagcases,1)
        Params_d=Params_A1c2;
        Params_d.upsilon=diagcases{dd,2}; Params_d.beta=diagcases{dd,3}; Params_d.cfloor=diagcases{dd,4};
        estimoptions_d=estimoptions_A1c2;
        estimoptions_d.skipestimation=1;
        estimoptions_d.skipjacobian=1;
        estimoptions_d.verbose=0;
        [~,~,estsummary_d]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames_A1c2,TargetMoments,WeightingMatrix,CoVarMatrixDataMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_d,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_gmm,estimoptions_d,vfoptions_gmm,simoptions_gmm);
        mv=estsummary_d.currentmomentvec(:); tv=estsummary_d.targetmomentvec(:);
        fprintf('   %-52s objective %8.1f; model medians: %3d of %d exactly zero, mean %7.0f (data mean %7.0f), top-quintile cells mean %7.0f\n',diagcases{dd,1},estsummary_d.objectivefnval,sum(mv==0),numel(mv),mean(mv),mean(tv),mean(mv(MomentList.PIq(tkorder)==5)))
    end
    fprintf('\nGMM estimation: DFJ2010 Table A1 column 2, the benchmark with cfloor fixed at %.0f, estimating %s\n',Params_A1c2.cfloor,strjoin(EstimParamNames_A1c2,', '))
    fprintf('   starting from DFJ2010 Table A1 column 2: upsilon %.2f, beta %.2f\n',Params_A1c2.upsilon,Params_A1c2.beta)
    tic;
    [EstimParams_A1c2,EstimParamsConfInts_A1c2,estsummary_A1c2]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames_A1c2,TargetMoments,WeightingMatrix,CoVarMatrixDataMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_A1c2,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_gmm,estimoptions_A1c2,vfoptions_gmm,simoptions_gmm);
    estimtime_A1c2=toc;
    fprintf('GMM estimation (Table A1 column 2) took %.1f minutes\n',estimtime_A1c2/60)
    fprintf('Estimated parameters (DFJ2010 Table A1 column 2 in brackets: upsilon 6.04 (1.01), beta .96 (.06)), standard errors in brackets:\n')
    for pp=1:length(EstimParamNames_A1c2)
        fprintf('   %-8s %12.4f  (%.4f)\n',EstimParamNames_A1c2{pp},EstimParams_A1c2.(EstimParamNames_A1c2{pp}),estsummary_A1c2.EstimParamsStdDev.(EstimParamNames_A1c2{pp}))
    end
    Est_c1=load('./SavedOutput/DFJ2010_Estimation_column1.mat','estsummary'); % the benchmark of doPart(4)
    fprintf('Objective at the estimates: %.2f [%d moments, %d parameters, %d df]; the benchmark (doPart(4), cfloor estimated): %.2f [%d df]; difference %.2f against the chi-squared(1) 5 percent critical value 3.84 (DFJ2010: 110.2 vs 82.3, difference 17.9)\n',estsummary_A1c2.objectivefnval,Jdata,length(EstimParamNames_A1c2),Jdata-length(EstimParamNames_A1c2),Est_c1.estsummary.objectivefnval,Jdata-3,estsummary_A1c2.objectivefnval-Est_c1.estsummary.objectivefnval)
    save ./SavedOutput/DFJ2010_Estimation_TableA1column2.mat EstimParams_A1c2 EstimParamsConfInts_A1c2 estsummary_A1c2 EstimParamNames_A1c2 estimtime_A1c2 Params_A1c2
end

%% Part 13: the benchmark re-estimated from the Table A1 column 2 estimates
% doPart(12) (2026-11-04) gave objective 86.35 with cfloor fixed at 5000 (upsilon 5.86, beta .963), below the 96.55 of the
% benchmark of doPart(4), which estimates cfloor too and started from DFJ2010's Table 3 column 1. A restricted model cannot fit
% better than the unrestricted one at its minimum, so doPart(4) stopped at a local minimum. This re-estimates the benchmark
% (upsilon, beta, cfloor; same estimator, targets, weighting and options as doPart(4)) starting from doPart(12)'s estimates.
if doPart(13)==1
    if any(doPart(6:11)==1)
        error('doPart(13) needs the exogenous model in the workspace, which parts 6-11 overwrite: run it without them')
    end
    vfoptions_gmm=vfoptions;
    vfoptions_gmm.verbose=0;
    vfoptions_gmm.verboseparams=0;
    Est_A1c2=load('./SavedOutput/DFJ2010_Estimation_TableA1column2.mat','EstimParams_A1c2','estsummary_A1c2','Params_A1c2','EstimParamNames_A1c2'); % from doPart(12)
    Params_c1b=Params; % the benchmark specification (Table 3 column 1, set above): delta=0, no bequest motive
    Params_c1b.upsilon=Est_A1c2.EstimParams_A1c2.upsilon;
    Params_c1b.beta=Est_A1c2.EstimParams_A1c2.beta;
    Params_c1b.delta=Table3.delta(1);
    Params_c1b.k=Table3.k(1);
    Params_c1b.theta=0;
    Params_c1b.cfloor=Est_A1c2.Params_A1c2.cfloor; % 5000, now estimated
    EstimParamNames_c1b=Table3.EstimParamNames{1}; % upsilon, beta, cfloor
    estimoptions_c1b=estimoptions; % as doPart(4): cohorts, eedefault=4, verbose
    estimoptions_c1b.constrainpositive=intersect({'upsilon','theta','k','cfloor'},EstimParamNames_c1b,'stable');
    estimoptions_c1b.constrain0to1=intersect({'beta'},EstimParamNames_c1b,'stable');
    fprintf('\nGMM estimation: the benchmark (Table 3 column 1), estimating %s\n',strjoin(EstimParamNames_c1b,', '))
    fprintf('   starting from the Table A1 column 2 estimates of doPart(12): upsilon %.4f, beta %.4f, cfloor %.1f\n',Params_c1b.upsilon,Params_c1b.beta,Params_c1b.cfloor)
    tic;
    [EstimParams_c1b,EstimParamsConfInts_c1b,estsummary_c1b]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames_c1b,TargetMoments,WeightingMatrix,CoVarMatrixDataMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_c1b,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_gmm,estimoptions_c1b,vfoptions_gmm,simoptions_gmm);
    estimtime_c1b=toc;
    fprintf('GMM estimation (benchmark from the column-2 estimates) took %.1f minutes\n',estimtime_c1b/60)
    Est_c1=load('./SavedOutput/DFJ2010_Estimation_column1.mat','EstimParams','estsummary'); % the benchmark of doPart(4)
    fprintf('Estimated parameters, standard errors in brackets; doPart(4) benchmark (started from DFJ2010) alongside:\n')
    for pp=1:length(EstimParamNames_c1b)
        fprintf('   %-8s %12.4f  (%.4f)    doPart(4): %12.4f  (%.4f)\n',EstimParamNames_c1b{pp},EstimParams_c1b.(EstimParamNames_c1b{pp}),estsummary_c1b.EstimParamsStdDev.(EstimParamNames_c1b{pp}),Est_c1.EstimParams.(EstimParamNames_c1b{pp}),Est_c1.estsummary.EstimParamsStdDev.(EstimParamNames_c1b{pp}))
    end
    fprintf('Objective at the estimates: %.2f [%d df]; doPart(4) benchmark: %.2f; Table A1 column 2 (doPart(12), cfloor fixed at 5000): %.2f [%d df]; column 2 minus this benchmark %.2f against the chi-squared(1) 5 percent critical value 3.84 (DFJ2010: 110.2 vs 82.3, difference 17.9)\n',estsummary_c1b.objectivefnval,Jdata-length(EstimParamNames_c1b),Est_c1.estsummary.objectivefnval,Est_A1c2.estsummary_A1c2.objectivefnval,Jdata-length(Est_A1c2.EstimParamNames_A1c2),Est_A1c2.estsummary_A1c2.objectivefnval-estsummary_c1b.objectivefnval)
    save ./SavedOutput/DFJ2010_Estimation_column1_fromA1c2start.mat EstimParams_c1b EstimParamsConfInts_c1b estsummary_c1b EstimParamNames_c1b estimtime_c1b Params_c1b
end

%% Part 14: map of the benchmark objective over upsilon and cfloor
% doPart(13) (2026-11-05) found the benchmark objective has at least two local minima: 96.55 at upsilon 3.72, cfloor 2673
% (doPart(4), started from DFJ2010) and 86.33 at upsilon 5.86, cfloor 5000 (started from doPart(12)); in both beta is about
% .965, and in doPart(13) cfloor did not move at all (the median moments make the objective piecewise flat). One evaluation
% at each point of a grid (skipestimation=1, skipjacobian=1, as the doPart(12) diagnostics), beta fixed at .965, maps where
% the lowest objective lies and how flat the objective is in cfloor. About a minute per evaluation, 49 evaluations.
if doPart(14)==1
    if any(doPart(6:11)==1)
        error('doPart(14) needs the exogenous model in the workspace, which parts 6-11 overwrite: run it without them')
    end
    vfoptions_gmm=vfoptions;
    vfoptions_gmm.verbose=0;
    vfoptions_gmm.verboseparams=0;
    Params_map=Params; % the benchmark specification (Table 3 column 1, set above): delta=0, no bequest motive
    Params_map.beta=0.965;
    Params_map.delta=Table3.delta(1);
    Params_map.k=Table3.k(1);
    Params_map.theta=0;
    EstimParamNames_map=Table3.EstimParamNames{1}; % upsilon, beta, cfloor (only used to set up the objective; nothing is estimated)
    estimoptions_map=estimoptions;
    estimoptions_map.constrainpositive=intersect({'upsilon','theta','k','cfloor'},EstimParamNames_map,'stable');
    estimoptions_map.constrain0to1=intersect({'beta'},EstimParamNames_map,'stable');
    estimoptions_map.skipestimation=1;
    estimoptions_map.skipjacobian=1;
    estimoptions_map.verbose=0;
    upsilon_map=3:0.75:7.5;
    cfloor_map=[1500,2500,3500,5000,6500,8000,10000];
    objective_map=nan(length(upsilon_map),length(cfloor_map));
    nzero_map=nan(length(upsilon_map),length(cfloor_map)); % number of model medians that are exactly zero
    fprintf('\nMap of the benchmark objective (beta fixed at %.3f, delta=0, no bequest motive): rows upsilon, columns cfloor\n',Params_map.beta)
    tic;
    for uu=1:length(upsilon_map)
        for cc=1:length(cfloor_map)
            Params_map.upsilon=upsilon_map(uu);
            Params_map.cfloor=cfloor_map(cc);
            [~,~,estsummary_map]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames_map,TargetMoments,WeightingMatrix,CoVarMatrixDataMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_map,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_gmm,estimoptions_map,vfoptions_gmm,simoptions_gmm);
            objective_map(uu,cc)=estsummary_map.objectivefnval;
            nzero_map(uu,cc)=sum(estsummary_map.currentmomentvec(:)==0);
        end
    end
    maptime=toc;
    fprintf('   %-8s',''); fprintf('  %9.0f',cfloor_map); fprintf('   (cfloor)\n')
    for uu=1:length(upsilon_map)
        fprintf('   %-8.2f',upsilon_map(uu)); fprintf('  %9.1f',objective_map(uu,:)); fprintf('\n')
    end
    fprintf('Number of the %d model medians that are exactly zero, same layout:\n',Jdata)
    for uu=1:length(upsilon_map)
        fprintf('   %-8.2f',upsilon_map(uu)); fprintf('  %9d',nzero_map(uu,:)); fprintf('\n')
    end
    [minobj_map,minind_map]=min(objective_map(:));
    [uumin,ccmin]=ind2sub(size(objective_map),minind_map);
    fprintf('Lowest objective on the map: %.2f at upsilon %.2f, cfloor %.0f (doPart(4): 96.55 at 3.72, 2673; doPart(13): 86.33 at 5.86, 5000; both with beta about .965). The map took %.1f minutes\n',minobj_map,upsilon_map(uumin),cfloor_map(ccmin),maptime/60)
    figure_c=figure_c+1;
    figure(figure_c);
    contourf(cfloor_map,upsilon_map,objective_map,20)
    colorbar
    hold on
    plot(2672.5,3.7225,'wo','MarkerFaceColor','w')
    plot(5000,5.8575,'ws','MarkerFaceColor','w')
    hold off
    xlabel('consumption floor (dollars)')
    ylabel('risk aversion \upsilon')
    title(sprintf('Benchmark objective, beta=%.3f (circle: doPart(4); square: doPart(13))',Params_map.beta))
    saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_ObjectiveMap.png')
    save ./SavedOutput/DFJ2010_ObjectiveMap.mat upsilon_map cfloor_map objective_map nzero_map Params_map maptime
end

%% End of run
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','on'); % do not leave an interactive session silently offscreen
end
diary off
