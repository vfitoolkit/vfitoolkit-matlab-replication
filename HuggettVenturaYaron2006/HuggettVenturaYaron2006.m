% Huggett, Ventura & Yaron (2006) - Human Capital and Earnings Distribution Dynamics
% This model is simple, just a deterministic life-cycle model.
% The main complication is that there are a lot of permanent types.
% Note that the permanent types are part of the calibration, and moreover the initial distribution 
% is part of the calibration and is a joint-distribution over the ability/type and the human capital levels.

% Human capital, h, is modelled as an experience asset: hprime(l,h) (see vfoptions.aprimeFn below).

%% Run infrastructure (diary, paths, output folders)
% Record the full run (incl. any error) to HVY2006_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./HVY2006_Diary.txt','file'), delete('./HVY2006_Diary.txt'); end
diary ./HVY2006_Diary.txt

% Subcodes (ReturnFn, jequaloneDistFn, ability grid, data work) live in
% ./HVY2006subcodes/; MATLAB does not search subfolders of the working
% directory, so put it on the path. Note this does NOT change cwd, so the
% './SavedOutput/...' paths (and the PSID files) still resolve to here.
addpath('./HVY2006subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
if exist('./VFIToolkit-matlab','dir'); addpath(genpath('./VFIToolkit-matlab/')); end

% Output folders
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

Params.headlessFigures=1; % set to 1 to render figures offscreen (needed for headless MATLAB, e.g. the server); PNGs still saved
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end

doPart=[0,0,0,1,1,1,1];
% doPart(1): pre-calibration check: solve model at the initial-guess parameters (VFI, agent dist, life-cycle profiles)
% doPart(2): data work: import PSID, earnings age profiles (Figs 1-3)       -> SavedOutput/HVY2006_Data.mat
% doPart(3): calibrate the initial dist over (h1,ability), for alpha=0.5,...,0.9 x start age 10,20 -> SavedOutput/HVY2006calib.mat
% doPart(4): results for all cases: Tables 1,2,3,5,6, Fig 5, cross-section, Prop 1 check -> SavedOutput/HVY2006_Results.mat, HVY2006full.mat
% doPart(5): the non-parametric case: Tables 2,3 (non-parametric rows), Table 4, Fig 4 -> SavedOutput/HVY2006_NonParametric.mat
% doPart(6): fixed initial human capital (common h1), non-parametric ability: Fig 6 -> SavedOutput/HVY2006_FixedH1.mat
% doPart(7): LaTeX tables for LatexDoc/App_HuggettVenturaYaron2006.tex, from the saved results (no GPU needed) -> SavedOutput/LatexInputs/
% A part that is switched off loads what later parts need from its .mat (only named
% variables, never a whole workspace: .mat files store stale function handles).
% Setup (parameters, grids, functions) always runs, it is cheap.

%% Grid sizes
n_d=101; % share of time to invest in new human capital
n_h=501; % Human capital [I call this n_h, as if I call it n_a I got confused about it being human capital not ability]
n_z=0;
N_i=21; % number of different ability levels
% Note: HVY2006, pg 277 says initial distribution was based on a 
% rectangular grid with 20 points on each of human capital and learning
% ability. Here I use 21 points for ability as I just prefer to use odd
% numbers when discretizing. [I use way more points for human capital.]

vfoptions.experienceasset=1; % Using an experience asset
simoptions.experienceasset=1;

N_j=39; % 39 or 49, depending on exercise [starting age 20, starting age 10, respectively]
if N_j==39
    Params.agejshifter=19; % starting age 20, to correspond to 39 periods
elseif N_j==49
    Params.agejshifter=9; % starting age 10, to correspond to 49 periods
end

figure_c=0;
CreateFigures=1 % Set to 0 to skip creating figures (so can be run headless)

% If you are running this with doPart(2)=1, you need to download the PSID
% data set that it uses:
% Download the data from PSID, then import
% Go to PSID and look for previous "carts"
% https://simba.isr.umich.edu/DC/c.aspx
% Search for and download the cart called 
% "325247 - HVY2006" associated with the email address "robertdkirkby@gmail.com"
% This will download a zip file, unzip it and put the files in the same
% folder as this current m-file.
% Rename the file which is .sps to _sps.txt
% Rename the file which is _format.sps to _formate_sps.txt
% (if you don't rename these ImportPSIDdata will just give a warning telling you to do so)


%% Parameters
Params.J=N_j;
Params.agej=1:1:Params.J; % the model period

% Human capital
Params.delta=0.0114; % depreciation rate of human capital
Params.alpha=0.7; % returns to scale of human capital production fn
% Note: HVY2006 try five different values of alpha: 0.5,0.6,0.7,0.8,0.9

% Interest rate
Params.r=0.04; % real interest rate (used to determine the discount factor)
Params.beta=1/(1+Params.r); % the discount factor (HVY2006 just call it 1/(1+r) )

% Wage growth
Params.g=0.0014;
Params.w=(1+Params.g).^(Params.agej-1); % wage rate (per unit of human-capital-hour)

%% Parameters for ability and for initial dist (on ability and human capital)

Params.mean_h1=92.3;
Params.coeffofvariation_h1=0.481;
Params.mean_ability=0.209;
Params.coeffofvariation_ability=0.347;
Params.corr_h1ability=0.781;
% These parameters are taken from Table 5; they are for parametric case with initial age 20 and alpha=0.7

% From these we can get the standard deviations
Params.stddev_h1=Params.coeffofvariation_h1*Params.mean_h1; % coeff. of variation=std. dev./mean, is the definition of the coeff. of variation
Params.stddev_ability=Params.coeffofvariation_ability*Params.mean_ability;

% But everything is really done around log(ability) and log(human capital)
% So we actually need the mean, std dev and correlations for the logs.
% There is no way to do this for std dev and correlation, so just have to
% stick in some guesses instead (we are anyway going to calibrate all these
% parameters later so no big deal).

% We want to discretize the log-normal distribution, currently we have mean
% and std. dev. and corr for the variables, we need to figure out the same
% for the logs (as that is what we will discretize as normal).
% Cannot just do this from formulas (maybe you can, you could with a single
% log-normal, the bivariate complicates things, I am not aware of formula,
% maybe there is one but I am fairly confident there is not).
Params.mean_logh1=log(Params.mean_h1);
Params.stddev_logh1=0.1; % PLACEHOLDER
Params.mean_logability=log(Params.mean_ability);
Params.stddev_logability=0.1; % PLACEHOLDER
Params.corr_logh1logability=0.5; % PLACEHOLDER

%%
% maxh is going to be coming from delta*h=ability*(h*l)^alpha
% Rearrange gives h=(ability/delta)^(1/(1-alpha))
% max of ability is roughly 0.4 (take a look at ability_grid below)
hmax=(0.4/Params.delta)^(1/(1-(Params.alpha)));
% Doing this gives hmax of around 10000, but solving model a few times is
% clear noone ever goes anywhere near this so I just use hmax=1000
h_grid=linspace(0.1,1000,n_h)'; % grid on human capital

l_grid=linspace(0,1,n_d)'; % grid on human capital investment time

% rename in toolkit notation
n_a=n_h;
d_grid=l_grid;
a_grid=h_grid;
z_grid=[];
pi_z=[];

%%  Set up the permanent types
% Ability is correlated with initial human capital, but we just ignore this correlation while we set up the grid on ability
% (we will consider the correlation later when putting agents onto the grids)

[ability_grid,pi_ability]=discretizeAR1_FarmerToda(Params.mean_logability,0,Params.stddev_logability,N_i);
ability_grid=exp(ability_grid);
pi_ability=pi_ability(1,:)'; % iid
% Comment: originally I just created a grid on ability, rather than log of ability, but this gave some points a negative 
% value for ability, so switched to current setup of creating the grid on log of ability, then take exponent.
% Note that the initial distribution is anyway log-normal, so working with a grid on log-ability makes more sense for that anyway.

% Ability is a parameter that depends on permanent type
Params.ability=ability_grid;

PTypeDistParamNames={'abilitydist'};
Params.abilitydist=pi_ability; % Note: This is not really used, as is replaced with the masses from the initial distribution which accounts for the correlation between initial human capital and initial assets


%% Set up the experience asset, we need to aprimeFn which gives value of aprime given d2 and a2 (d2 is the decision variable relevant to experience asset, a2 is the experience asset)
% f_hl=ability*(h*l)^alpha; % new production of human capital
% hprime=f_hl+h*(1-delta);

% aprimeFn: hprime(l,h,parameters)
vfoptions.aprimeFn=@(l,h,ability,alpha,delta) ability*(h*l)^alpha+h*(1-delta);
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.a_grid=a_grid; % Needed when using experience asset
simoptions.d_grid=d_grid; % Needed when using experience asset

%% Figure 1: next period human capital h' against h, when all time is invested in human capital (l=1)
% Shows how far h can grow: h stops growing where h' meets the 45 degree line, which for high
% ability is far above the h grid, so in practice it is the finite N_j that bounds h.
if CreateFigures==1
    figure_c=figure_c+1;
    fig1=figure(figure_c);
    hprime_lowability=Params.ability(1)*h_grid.^Params.alpha+h_grid*(1-Params.delta); % aprimeFn with l=1
    hprime_highability=Params.ability(end)*h_grid.^Params.alpha+h_grid*(1-Params.delta);
    plot(h_grid,hprime_lowability,'-b','LineWidth',1); hold on;
    plot(h_grid,hprime_highability,'-r','LineWidth',1);
    plot(h_grid,h_grid,'--k');
    hold off;
    xlabel('h (human capital)'); ylabel('h'' (next period human capital)');
    legend({'lowest ability','highest ability','45 degree line'},'Location','northwest');
    title('Human capital production, all time invested (l=1)');
    saveas(fig1,'./SavedOutput/Graphs/HVY2006_HumanCapitalProduction.png');
end

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(l,h,w)...
    HuggettVenturaYaron2006_ReturnFn(l,h,w);

%% Initial agent distribution as a function, as we need that so that the parameters of this can be part of our calibration
% jequaloneDist as a function has first four inputs (a_grid,z_grid,n_a,n_z,...) followed by any parameters
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability, ability, Number_i) HuggettVenturaYaron2006_jequaloneDistFn(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability, ability, Number_i);
% To use jequaloneDist as a function you must include a_grid and z_grid in simoptions
simoptions.a_grid=h_grid;
simoptions.z_grid=[]; % Even though we don't use z_grid, it is still needed in simoptions to be able to call jequaloneDist as a function
% Params.ability=ability_grid; % so we can pass it as an input to jequaloneDistFn, We already set this above
Params.Number_i=N_i; % so we can pass it as an input to jequaloneDistFn

% Our parametrization of the covariance matrix follows AH2021:
% So we need to set up initial guess for
Params.FTcorr_logh1logability=0.5*log((1+Params.corr_logh1logability)/(1-Params.corr_logh1logability)); % This formula is exact when using AH2021 for 2x2 Covar matrix

%% Age weights
% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age


%% doPart(1): Pre-calibration check
if doPart(1)==1
%% Now solve the value function iteration problem, just to check that things are working before we go to calibration
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
vftime=toc


%% HVY2006 consider both a parametric and a non-parametric initial distribution over initial human  capital, h1 and ability.
% Here we will just do the parametric, which uses a bi-variate log-normal distribution.

Params.logh1logability_Mean=[Params.mean_logh1,Params.mean_logability];

% Create covariance matrix for std devs and correlation
CorrMatrix=[1, Params.corr_logh1logability; Params.corr_logh1logability, 1];
Params.logh1logability_CoVarMatrix = corr2cov_homemade([Params.stddev_logh1,Params.stddev_logability],CorrMatrix);

% Now we can find the probabilities of a bivariate log-normal distribution over the existing grids
tic;
P=MVNormal_ProbabilitiesOnGrid([log(h_grid);log(ability_grid)],Params.logh1logability_Mean, Params.logh1logability_CoVarMatrix, [n_h,N_i]);
% Note: P is on defined on (h,ability).
mvntime=toc

tic;
P=gpuArray(MVNormal_ProbabilitiesOnGrid(gather([log(h_grid);log(ability_grid)]),gather(Params.logh1logability_Mean), gather(Params.logh1logability_CoVarMatrix), [n_h,N_i]));
% Note: P is on defined on (h,ability).
mvntime2=toc % This is actually way way faster for some reason

% Make sure all the agent ptypes have positive (non-zero) probabilities [as otherwise this would cause problems]
sum(P,1) % with N_i=21 the highest and lowest are roughly 10 to minus eight
% sum(P,1)==0 shows they are all non-zero (only just!)

% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
% jequaloneDist=zeros(n_a,N_i,'gpuArray'); % Put no households anywhere on grid (I just leave this here to show what size it is)
jequaloneDist=P; % joint log-normal distribution onto our existing grids
% (the agent distribution below uses jequaloneDistFn, the same distribution but parametrized with FTcorr)

%% We now compute the 'stationary distribution' of households
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistFn,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
statdisttime=toc

%% Calculate the life-cycle profiles
FnsToEvaluate.humancapital=@(l,h) h; % h is human capital
FnsToEvaluate.timel=@(l,h) l; % l is fraction of time invested in human capital

tic;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
acstime=toc

if CreateFigures==1
    figure_c=figure_c+1;
    fig=figure(figure_c);
    subplot(2,1,1); plot(Params.agejshifter+Params.agej, AgeConditionalStats.humancapital.Mean)
    title('Human capital: mean life-cycle profile')
    subplot(2,1,2); plot(Params.agejshifter+Params.agej, AgeConditionalStats.timel.Mean)
    title('Fraction of time learning: mean life-cycle profile')
    saveas(fig,'./SavedOutput/Graphs/HVY2006_PreCalibCheck.png')
end
% Note, we can check our distribution of the initial distribuiton by
% looking at the means and correlation for h and ability
end

%% doPart(2): Let's do some data work to prepare for the model calibration
% Following code imports PSID data, then estimates the life-cycle profiles
% for earnings that we use for model calibration
if doPart(2)==1
    HVY2006_DataWork
    save ./SavedOutput/HVY2006_Data.mat Fig1_Mean Fig1_Median Fig1_Gini B ptiles ptilemass % B is the 23 earnings percentiles by age (Fig 2), used for the cross-section
else
    if ~exist('./SavedOutput/HVY2006_Data.mat','file')
        % Report where MATLAB is actually looking, then stop
        fprintf('Working directory (pwd) is: %s\n',pwd)
        fprintf('Contents of ./SavedOutput/ :\n')
        dir('./SavedOutput/')
        error('./SavedOutput/HVY2006_Data.mat not found (see working directory above); run doPart(2) or copy it there')
    end
    load('./SavedOutput/HVY2006_Data.mat','Fig1_Mean','Fig1_Median','Fig1_Gini','B','ptiles','ptilemass')
end

%% Calibration exercise of HVY2006
% Choose the parameters that parametrize the initial distribution over
% human capital and ability to minimize the sum of squares of the log
% ratios for the life-cycle profile of earnings in mean, Gini, and
% ratio-of-mean-to-median (HVY2006 footnote 18; they use median/mean, the
% 'inverse skewness', which gives exactly the same squared log ratio).

%% One more thing, we need to deal with the 'ability'
% ability is a parameter that depends on PType
% But it is also determined/parametrized by mean_logability and stddev_logability which
% are parameters to be calibrated.
% This idea of a PType parameter having a parametrized distribution across PTypes
% is common, and the way VFI Toolkit handles it is with PTypeParamFn

PTypeParamFn=@(Params) HVY2006_AbilityGrid(Params);
% PTypeParamFn will always have Params as it's first and only entry
% Because any parameter that depends on PType will have to be kept in
% Params, the output of the function will also just be Params.

%% Now, we calibrate the model parameters to target the life cycle profile of earnings in mean, mean/median, and Gini

% To do this using the toolkit there are two things we need to setup
% First, just say all the model parameters we want to calibrate
% For exercise of HVY2006 these are all the parameters that define the
% initial distribution over human capital and ability.
CalibParamNames={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTcorr_logh1logability'};
% CalibParamNames gives the names of all the model parameters that will be
% calibrated, these parameters must all appear in Params, and the values in
% Params will be used as the initial values.
% All other parameters in Params will remain fixed.

% Because the third and fourth parameters are standard deviations we want to constrain them to be positive
caliboptions.constrainpositive={'stddev_logh1','stddev_logability'};
% The fifth parameter parametrizes the correlation, corr=tanh(FTcorr) (AH2021 for a 2x2 matrix), so any
% real FTcorr gives a valid correlation. We just keep it from running off: FTcorr in [-3,3] means corr in (-0.995,0.995)
% [Note: constraining FTcorr to 0 to 1 would cap corr at tanh(1)=0.76, below the ~0.80 implied by HVY2006 Table 5]
caliboptions.constrainAtoB={'FTcorr_logh1logability'};
caliboptions.constrainAtoBlimits.FTcorr_logh1logability=[-3,3];
% Note: constraints are given by parameter name (a binary vector is silently ignored by the current toolkit)

% Second, we need to say which model statistics we want to target (set for each case below, as
% TargetMoments.AgeConditionalStats.earnings.Mean/RatioMeanToMedian/Gini)
% Note: When setting up TargetMoments there are some rules you must follow
% There are two options TargetMomements.AgeConditionalStats and
% TargetMoments.AllStats (you can use both). Within these you must
% follow the structure that you get when you run the commands
% AgeConditionalStats=LifeCycleProfiles_FHorz_Case1()
% and
% AllStats=EvalFnOnAgentDist_AggVars_FHorz_Case1()
% [If you are using PType, then these will be the PType equivalents]

% So we want a FnsToEvaluate which is just earnings (this will be faster as it only includes what we actually need)
clear FnsToEvaluate
FnsToEvaluate.earnings=@(l,h,w) w*h*(1-l);

% We want to set a few options, one relates to the objective function
caliboptions.metric='sum_logratiosquared';


%% The cases: HVY2006 calibrate for alpha=0.5,0.6,0.7,0.8,0.9, with human capital accumulation starting at age 10 or at age 20
alpha_vec=[0.5,0.6,0.7,0.8,0.9];
startage_vec=[10,20];
% The baseline case (alpha=0.7, accumulation from age 20) goes last, so its results are what
% remains in the workspace for the figures at the end of doPart(4)
CaseList=zeros(0,2); % rows are [alpha, startage]
for ss=1:length(startage_vec)
    for aa=1:length(alpha_vec)
        if ~(alpha_vec(aa)==0.7 && startage_vec(ss)==20)
            CaseList=[CaseList; alpha_vec(aa), startage_vec(ss)];
        end
    end
end
CaseList=[CaseList; 0.7, 20];
nCases=size(CaseList,1);

% HVY2006 Table 5 (parametric case), used as initial guesses for our calibration (and to compare our results with)
% rows: mean(a), CV(a), skewness(a), mean(h1), CV(h1), skewness(h1), corr(a,h1); columns: alpha=0.5,0.6,0.7,0.8,0.9
Table5paper.startage10=[0.499 0.322 0.207 0.139 0.089; 0.514 0.436 0.353 0.235 0.198; 1.125 1.092 1.061 1.027 1.010; 64.0 69.2 74.7 75.1 78.6; 0.454 0.453 0.434 0.403 0.184; 1.100 1.100 1.090 1.077 1.071; 0.070 0.145 0.171 0.333 0.351];
Table5paper.startage20=[0.467 0.321 0.209 0.136 0.088; 0.613 0.474 0.347 0.257 0.158; 1.191 1.109 1.058 1.033 1.012; 86.7 89.5 92.3 96.6 100.1; 0.427 0.439 0.481 0.468 0.459; 1.088 1.092 1.109 1.105 1.100; 0.600 0.621 0.781 0.792 0.796];
% HVY2006 Table 2 and Table 3, parametric rows (to compare our results with)
Table2paper.startage10=[7.5 6.4 5.9 5.2 6.2];
Table2paper.startage20=[6.8 7.0 5.2 5.0 6.4];
Table3paper.startage10=[77.8 72.9 68.0 77.1 79.2];
Table3paper.startage20=[73.0 65.4 72.8 70.8 78.3];


%% doPart(3): calibrate each case
if doPart(3)==1
    save ./SavedOutput/HVY2006precalib.mat

    CalibParams=struct();
    calibsummary=struct();
    for cc=1:nCases
        %% Set up this case
        alpha_c=CaseList(cc,1);
        startage_c=CaseList(cc,2);
        casename=sprintf('alpha%02i_start%02i',round(10*alpha_c),startage_c);
        N_j=58-startage_c+1; % 49 or 39 model periods, ending at age 58
        Params.J=N_j;
        Params.agej=1:1:N_j;
        Params.agejshifter=startage_c-1;
        Params.w=(1+Params.g).^(Params.agej-1); % wage rate (per unit of human-capital-hour)
        Params.mewj=ones(1,N_j)/N_j; % equal mass at each age
        Params.alpha=alpha_c;
        % The data targets are for ages 20-58; with accumulation from age 10 the first ten model ages have no data (NaN)
        nodata=nan(1,N_j-length(Fig1_Mean));
        TargetMoments.AgeConditionalStats.earnings.Mean=[nodata, Fig1_Mean'];
        TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian=[nodata, Fig1_Mean'./Fig1_Median'];
        TargetMoments.AgeConditionalStats.earnings.Gini=[nodata, Fig1_Gini'];

        %% Initial guess: the joint log-normal implied by the HVY2006 Table 5 values for this case
        % For a log-normal X with mean m and coeff of variation cv: var(log X)=log(1+cv^2), E[log X]=log(m)-var(log X)/2,
        % and corr(X1,X2)=(exp(rho*s1*s2)-1)/sqrt((exp(s1^2)-1)*(exp(s2^2)-1)) where rho is the correlation of the logs
        paperT5=Table5paper.(sprintf('startage%i',startage_c));
        colc=find(alpha_vec==alpha_c);
        s2_a=log(1+paperT5(2,colc)^2);
        s2_h=log(1+paperT5(5,colc)^2);
        Params.mean_logability=log(paperT5(1,colc))-s2_a/2;
        Params.stddev_logability=sqrt(s2_a);
        Params.mean_logh1=log(paperT5(4,colc))-s2_h/2;
        Params.stddev_logh1=sqrt(s2_h);
        rho_log=log(1+paperT5(7,colc)*sqrt((exp(s2_a)-1)*(exp(s2_h)-1)))/sqrt(s2_a*s2_h);
        Params.FTcorr_logh1logability=atanh(min(max(rho_log,-0.99),0.99));

        fprintf('\n========== Calibrating case %i of %i: alpha=%1.1f, accumulation starts at age %i ==========\n',cc,nCases,alpha_c,startage_c)
        [CalibParams_c,calibsummary_c]=CalibrateLifeCycleModel_PType(CalibParamNames,TargetMoments,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, jequaloneDistFn,AgeWeightParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate, caliboptions, vfoptions,simoptions);
        % CalibParams is the calibrated parameter values
        % calibsummary is a structure containing various info on how the calibration went
        CalibParams.(casename)=CalibParams_c;
        calibsummary.(casename)=calibsummary_c;

        save ./SavedOutput/HVY2006calib.mat CalibParams calibsummary CaseList % save after each case, so a run that stops early keeps the cases it finished
    end
else
    load('./SavedOutput/HVY2006calib.mat','CalibParams','calibsummary','CaseList')
    nCases=size(CaseList,1);
end


%% doPart(4): Report the things in HVY2006
if doPart(4)==1
% We already replicated Figure 1 (in HVY2006_DataWork)
% We already replicated Figure 2 (in HVY2006_DataWork)
% We already replicated Figure 3 (in HVY2006_DataWork: cohort dummies, time dummies, restricted time dummies)

%% Table 1 is just the parameters
fprintf('\nTable 1 of HVY2006: parameter values \n')
fprintf('Model periods J=49 (accumulation from age 10) and J=39 (from age 20); interest rate r=%1.2f; rental growth g=%1.4f; depreciation delta=%1.4f; alpha in %s \n', Params.r, Params.g, Params.delta, mat2str(alpha_vec))
fprintf('End of outputs relating to Table 1 \n')

% Functions to evaluate: earnings for Fig 5 and Table 2, ability and human capital for Table 5, time learning to look at
FnsToEvaluate.humancapital=@(l,h) h; % h is human capital
FnsToEvaluate.timel=@(l,h) l; % l is fraction of time invested in human capital
FnsToEvaluate.ability=@(l,h,ability) ability;

Results=struct();
for cc=1:nCases
    %% Set up this case (same as in doPart(3))
    alpha_c=CaseList(cc,1);
    startage_c=CaseList(cc,2);
    casename=sprintf('alpha%02i_start%02i',round(10*alpha_c),startage_c);
    N_j=58-startage_c+1; % 49 or 39 model periods, ending at age 58
    Params.J=N_j;
    Params.agej=1:1:N_j;
    Params.agejshifter=startage_c-1;
    Params.w=(1+Params.g).^(Params.agej-1); % wage rate (per unit of human-capital-hour)
    Params.mewj=ones(1,N_j)/N_j; % equal mass at each age
    Params.alpha=alpha_c;
    % The data targets are for ages 20-58; with accumulation from age 10 the first ten model ages have no data (NaN)
    nodata=nan(1,N_j-length(Fig1_Mean));
    TargetMoments.AgeConditionalStats.earnings.Mean=[nodata, Fig1_Mean'];
    TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian=[nodata, Fig1_Mean'./Fig1_Median'];
    TargetMoments.AgeConditionalStats.earnings.Gini=[nodata, Fig1_Gini'];

    % Use the calibrated parameter values for this case
    for pp=1:length(CalibParamNames)
        Params.(CalibParamNames{pp})=CalibParams.(casename).(CalibParamNames{pp});
    end
    Params=PTypeParamFn(Params); % ability grid for the calibrated mean_logability and stddev_logability

    %% Solve the model at the calibrated parameters
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, [], [], ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    jequaloneDist=jequaloneDistFn(a_grid,z_grid,n_a,n_z,Params.mean_logh1,Params.mean_logability,Params.stddev_logh1,Params.stddev_logability,Params.FTcorr_logh1logability, Params.ability, Params.Number_i);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);

    Results.(casename).alpha=alpha_c;
    Results.(casename).startage=startage_c;
    Results.(casename).CalibParams=CalibParams.(casename);
    Results.(casename).age=Params.agejshifter+(1:1:N_j);
    Results.(casename).earningsMean=gather(AgeConditionalStats.earnings.Mean);
    Results.(casename).earningsMedian=gather(AgeConditionalStats.earnings.Median);
    Results.(casename).earningsGini=gather(AgeConditionalStats.earnings.Gini);
    Results.(casename).earningsRatioMeanToMedian=gather(AgeConditionalStats.earnings.RatioMeanToMedian);

    %% Table 2: goodness of fit, HVY2006 footnote 20: sum over the J ages with data of |log ratio| of mean, dispersion and skewness, divided by 3J
    dataages=~isnan(TargetMoments.AgeConditionalStats.earnings.Mean);
    AbsLogMean=abs(log(TargetMoments.AgeConditionalStats.earnings.Mean(dataages)./Results.(casename).earningsMean(dataages)));
    AbsLogRatioMeanToMedian=abs(log(TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian(dataages)./Results.(casename).earningsRatioMeanToMedian(dataages)));
    AbsLogGini=abs(log(TargetMoments.AgeConditionalStats.earnings.Gini(dataages)./Results.(casename).earningsGini(dataages)));
    Results.(casename).Table2=100*(sum(AbsLogMean)+sum(AbsLogRatioMeanToMedian)+sum(AbsLogGini))/(3*sum(dataages));
    % Not in Table 2, but the same thing in the metric of the calibration (squared rather than absolute log ratios)
    Results.(casename).SumSqLogDeviation=100*(sum(AbsLogMean.^2)+sum(AbsLogRatioMeanToMedian.^2)+sum(AbsLogGini.^2))/(3*sum(dataages));

    %% Table 3: percentage of the variance of the present value of earnings due to learning ability differences
    % HVY2006 footnote 22: 100*var(E[PV|a])/var(PV), with the present value of earnings 'as of age 20' (HVY2006 p281),
    % for both starting ages. The value function at age 20 is exactly the present value of earnings from age 20
    % [the ReturnFn is earnings, discount 1/(1+r)], so weight V at age 20 by the agent distribution at age 20.
    % [With accumulation from age 20 this is the initial condition (a,h1); from age 10 it is the (a,h) people reach by age 20.]
    P1=gather(jequaloneDist); P1=P1/sum(P1(:)); % initial distribution on (h1,ability) (also used for Tables 5 and 6)
    jage20=20-Params.agejshifter; % model period of age 20
    ptweights=gather(StationaryDist.ptweights);
    PV=zeros(n_h,N_i);
    P20=zeros(n_h,N_i); % distribution at age 20 on (h,ability)
    for ii=1:N_i
        V_ii=gather(V.(sprintf('ptype%03i',ii))); % [n_h,N_j]
        PV(:,ii)=V_ii(:,jage20);
        Dist_ii=gather(StationaryDist.(sprintf('ptype%03i',ii))); % [n_h,N_j]
        P20(:,ii)=ptweights(ii)*Dist_ii(:,jage20)/max(sum(Dist_ii(:,jage20)),eps);
    end
    P20=P20/sum(P20(:));
    Mean_PV=sum(P20(:).*PV(:));
    Var_PV=sum(P20(:).*(PV(:)-Mean_PV).^2);
    mass_i=sum(P20,1); % mass of each ability level
    CondMean_PV=sum(P20.*PV,1)./max(mass_i,eps);
    Var_CondMean_PV=sum(mass_i.*(CondMean_PV-Mean_PV).^2);
    Results.(casename).Table3=100*Var_CondMean_PV/Var_PV;

    %% Table 5: ability and human capital at birth (the first model age)
    Results.(casename).Table5=zeros(7,1); % rows as in HVY2006 Table 5
    Results.(casename).Table5(1)=gather(AgeConditionalStats.ability.Mean(1));
    Results.(casename).Table5(2)=gather(AgeConditionalStats.ability.StdDeviation(1)/AgeConditionalStats.ability.Mean(1)); % Coeff of variation=std dev/mean
    Results.(casename).Table5(3)=gather(AgeConditionalStats.ability.RatioMeanToMedian(1)); % HVY2006 measure skewness as mean/median
    Results.(casename).Table5(4)=gather(AgeConditionalStats.humancapital.Mean(1));
    Results.(casename).Table5(5)=gather(AgeConditionalStats.humancapital.StdDeviation(1)/AgeConditionalStats.humancapital.Mean(1));
    Results.(casename).Table5(6)=gather(AgeConditionalStats.humancapital.RatioMeanToMedian(1));
    % Correlation of (ability, h1): computed directly from the initial distribution, which is on (h,ability)
    h1vals=gather(h_grid); % column, so varies along the h dimension of P1
    abilityvals=gather(Params.ability); abilityvals=abilityvals(:)'; % row, so varies along the ability dimension of P1
    Mean_h1=sum(sum(P1.*h1vals)); Mean_ability=sum(sum(P1.*abilityvals));
    Cov_h1ability=sum(sum(P1.*(h1vals-Mean_h1).*(abilityvals-Mean_ability)));
    StdDev_h1=sqrt(sum(sum(P1.*(h1vals-Mean_h1).^2))); StdDev_ability=sqrt(sum(sum(P1.*(abilityvals-Mean_ability).^2)));
    Results.(casename).Table5(7)=Cov_h1ability/(StdDev_h1*StdDev_ability);

    %% Check: how much mass is near the top of the human capital grid (if it is not negligible, h_grid needs a higher max)
    massnearmax=zeros(1,N_j);
    for ii=1:N_i
        Dist_ii=gather(StationaryDist.(sprintf('ptype%03i',ii))); % [n_h,N_j]
        Dist_ii=Dist_ii./max(sum(Dist_ii,1),eps); % distribution over h at each age, within ptype
        massnearmax=massnearmax+ptweights(ii)*sum(Dist_ii(end-4:end,:),1);
    end
    Results.(casename).maxmassnearhmax=max(massnearmax);

    %% Table 6: persistence of individual earnings (HVY2006 report it for alpha=0.7)
    % The model is deterministic given (a,h1), so follow each initial condition (h1 on the grid, each
    % ability level) forward. Earnings growth rates are very sensitive to the discreteness of the l grid
    % (l moves in steps of 0.01 as h crosses grid points, which is noise of about 1% in earnings growth),
    % so use the model's closed-form solution: the value function is linear in h, with slope
    % lambda_j=w_j+beta*(1-delta)*lambda_{j+1} (lambda_{J+1}=0), so the optimal human capital put into
    % producing new human capital, x=h*l, is x_j=min(h,(beta*lambda_{j+1}*alpha*a/w_j)^(1/(1-alpha))).
    % [This reproduces the toolkit's mean earnings profile to within 0.3% at every age.]
    if alpha_c==0.7
        lambda=zeros(1,N_j+1);
        for jj=N_j:-1:1
            lambda(jj)=Params.w(jj)+Params.beta*(1-Params.delta)*lambda(jj+1);
        end
        hpath=repmat(gather(h_grid),1,N_i); % initial human capital, one column per ability level
        Epath=zeros(n_h,N_i,N_j);
        for jj=1:N_j
            xstar=(Params.beta*lambda(jj+1)*Params.alpha*abilityvals/Params.w(jj)).^(1/(1-Params.alpha)); % 1-by-N_i
            xnow=min(hpath,xstar); % l=1 (all time learning) when xstar>=h
            Epath(:,:,jj)=Params.w(jj)*(hpath-xnow); % earnings w*h*(1-l)
            hpath=abilityvals.*xnow.^Params.alpha+hpath*(1-Params.delta);
        end
        Table6ages=[45,40];
        Table6lags=[1,5,10];
        Results.(casename).Table6=zeros(6,2); % rows: levels lag 1,5,10, then growth rates lag 1,5,10; columns: age 45, age 40
        for tt=1:2
            jj=Table6ages(tt)-Params.agejshifter; % model period of this age
            for ll=1:3
                jlag=jj-Table6lags(ll);
                % levels: corr(E_j, E_j-lag), weighted by the initial distribution
                X1=reshape(Epath(:,:,jj),[],1); X2=reshape(Epath(:,:,jlag),[],1); wts=P1(:);
                X1m=sum(wts.*X1); X2m=sum(wts.*X2);
                Results.(casename).Table6(ll,tt)=sum(wts.*(X1-X1m).*(X2-X2m))/sqrt(sum(wts.*(X1-X1m).^2)*sum(wts.*(X2-X2m).^2));
                % growth rates: z_j=log(E_j/E_j-1); corr(z_j, z_j-lag), using the initial conditions with positive earnings at all four ages
                E_j=reshape(Epath(:,:,jj),[],1); E_jm1=reshape(Epath(:,:,jj-1),[],1);
                E_l=reshape(Epath(:,:,jlag),[],1); E_lm1=reshape(Epath(:,:,jlag-1),[],1);
                positive=(E_j>0 & E_jm1>0 & E_l>0 & E_lm1>0);
                wts=P1(:).*positive; wts=wts/sum(wts);
                Z1=zeros(size(E_j)); Z2=zeros(size(E_j));
                Z1(positive)=log(E_j(positive)./E_jm1(positive)); Z2(positive)=log(E_l(positive)./E_lm1(positive));
                Z1m=sum(wts.*Z1); Z2m=sum(wts.*Z2);
                Results.(casename).Table6(3+ll,tt)=sum(wts.*(Z1-Z1m).*(Z2-Z2m))/sqrt(sum(wts.*(Z1-Z1m).^2)*sum(wts.*(Z2-Z2m).^2));
            end
        end

        %% Cross-section earnings distribution (HVY2006 Section 5.1, end): pool ages 20-58, with population
        % weights for 1% population growth, and each age's earnings distribution (the exact paths above,
        % weighted by the initial distribution). The data side is done the same way, from the percentiles.
        ages_cs=20:1:58;
        popweight=1.01.^(-(ages_cs-20)); popweight=popweight/sum(popweight);
        csvalues=zeros(n_h*N_i,length(ages_cs)); csweights=zeros(n_h*N_i,length(ages_cs));
        for aa=1:length(ages_cs)
            jj=ages_cs(aa)-Params.agejshifter; % model period of this age
            csvalues(:,aa)=reshape(Epath(:,:,jj),[],1);
            csweights(:,aa)=popweight(aa)*P1(:);
        end
        csvalues=csvalues(:); csweights=csweights(:);
        keep=(csweights>0); csvalues=csvalues(keep); csweights=csweights(keep)/sum(csweights(keep));
        [csvalues,sortindex]=sort(csvalues); csweights=csweights(sortindex);
        cumpop=[0; cumsum(csweights)]; cuminc=[0; cumsum(csweights.*csvalues)]; cuminc=cuminc/cuminc(end);
        csGini=1-sum((cuminc(2:end)+cuminc(1:end-1)).*csweights); % exact Gini of the discrete distribution
        csMedian=csvalues(find(cumpop(2:end)>=0.5,1));
        [cumpop_u,iu]=unique(cumpop,'last'); % (tiny weights can leave the cumulative population share unchanged in floating point)
        csTopShares=1-interp1(cumpop_u,cuminc(iu),1-[0.2,0.1,0.05,0.01]); % shares of earnings of the top 20%, 10%, 5% and 1%
        Results.(casename).CrossSection=[csGini, sum(csweights.*csvalues)/csMedian, 100*csTopShares];
    end
end
% The baseline case (alpha=0.7, accumulation from age 20) was last, so the workspace now holds its
% Params, N_j, TargetMoments, V, Policy, jequaloneDist, StationaryDist and AgeConditionalStats.

%% Cross-section earnings distribution in the data, constructed in the same way from the 23 percentiles by age (B, from HVY2006_DataWork)
ages_cs=20:1:58;
popweight=1.01.^(-(ages_cs-20)); popweight=popweight/sum(popweight);
csvalues=B(:); % B is ages-by-percentiles
csweights=reshape(popweight'*ptilemass,[],1); % each percentile stands for the population mass ptilemass, at each age
[csvalues,sortindex]=sort(csvalues); csweights=csweights(sortindex);
cumpop=[0; cumsum(csweights)]; cuminc=[0; cumsum(csweights.*csvalues)]; cuminc=cuminc/cuminc(end);
csGini=1-sum((cuminc(2:end)+cuminc(1:end-1)).*csweights);
csMedian=csvalues(find(cumpop(2:end)>=0.5,1));
[cumpop_u,iu]=unique(cumpop,'last'); % (tiny weights can leave the cumulative population share unchanged in floating point)
csTopShares=1-interp1(cumpop_u,cuminc(iu),1-[0.2,0.1,0.05,0.01]);
CrossSectionData=[csGini, sum(csweights.*csvalues)/csMedian, 100*csTopShares];

%% Check of HVY2006 Proposition 1(iii): with a common ability, the cohort Lorenz curves for human capital
% and earnings become more equal with age (so differences in ability are essential for rising dispersion).
% Baseline case (alpha=0.7, from age 20), with every agent given the mean ability and keeping the initial
% human capital distribution. Uses the closed-form solution (as for Table 6).
Prop1ability=Mean_ability; % mean ability of the baseline initial distribution
Prop1h1mass=sum(P1,2); % marginal distribution of initial human capital
lambda=zeros(1,N_j+1);
for jj=N_j:-1:1
    lambda(jj)=Params.w(jj)+Params.beta*(1-Params.delta)*lambda(jj+1);
end
Prop1Astar=(Params.beta*lambda(2:N_j+1)*Params.alpha*Prop1ability./Params.w).^(1/(1-Params.alpha)); % the cutoff A_j(a): full-time learning below it
% Proposition 1(ii)-(iii) condition: a*A_j^alpha+(1-delta)*A_j >= A_j+1
Prop1condition=all(Prop1ability*Prop1Astar(1:end-1).^Params.alpha+(1-Params.delta)*Prop1Astar(1:end-1)>=Prop1Astar(2:end));
hnow=gather(h_grid);
Prop1h=zeros(n_h,N_j); Prop1E=zeros(n_h,N_j);
for jj=1:N_j
    xnow=min(hnow,Prop1Astar(jj));
    Prop1h(:,jj)=hnow;
    Prop1E(:,jj)=Params.w(jj)*(hnow-xnow);
    hnow=Prop1ability*xnow.^Params.alpha+hnow*(1-Params.delta);
end
% Lorenz curves at population shares 0.01,...,0.99 for each age, and the Gini by age
popshares=0.01:0.01:0.99;
Prop1LorenzH=zeros(length(popshares),N_j); Prop1LorenzE=zeros(length(popshares),N_j);
Prop1GiniH=zeros(1,N_j); Prop1GiniE=zeros(1,N_j);
keep=(Prop1h1mass>0); wts=Prop1h1mass(keep)/sum(Prop1h1mass(keep));
for jj=1:N_j
    for vv=1:2
        if vv==1
            vals=Prop1h(keep,jj);
        else
            vals=Prop1E(keep,jj);
        end
        [vals,sortindex]=sort(vals); w_s=wts(sortindex);
        cumpop=[0; cumsum(w_s)]; cuminc=[0; cumsum(w_s.*vals)]; cuminc=cuminc/cuminc(end);
        [cumpop_u,iu]=unique(cumpop,'last');
        if vv==1
            Prop1LorenzH(:,jj)=interp1(cumpop_u,cuminc(iu),popshares');
            Prop1GiniH(jj)=1-sum((cuminc(2:end)+cuminc(1:end-1)).*w_s);
        else
            Prop1LorenzE(:,jj)=interp1(cumpop_u,cuminc(iu),popshares');
            Prop1GiniE(jj)=1-sum((cuminc(2:end)+cuminc(1:end-1)).*w_s);
        end
    end
end
Prop1allpositive=all(Prop1E(keep,:)>0,1); % Proposition 1(iii) assumes everyone has strictly positive earnings
% Lorenz curve 'more equal' from age j to j+1: L_j+1(p)>=L_j(p) at every population share p
Prop1violationH=max(max(Prop1LorenzH(:,1:end-1)-Prop1LorenzH(:,2:end))); % positive means a violation
jjE=find(Prop1allpositive,1):N_j; % ages from which all earnings are positive
Prop1violationE=max(max(Prop1LorenzE(:,jjE(1:end-1))-Prop1LorenzE(:,jjE(2:end))));

%% Print Tables 2, 3, 5 and 6 (ours, with the HVY2006 values alongside)
alphastr=sprintf('%8.1f',alpha_vec);
for ss=1:length(startage_vec)
    fprintf('\nTable 2 of HVY2006: mean absolute deviation (%%), parametric, accumulation starts at age %i \n',startage_vec(ss))
    fprintf('alpha               %s \n',alphastr)
    ours=zeros(1,length(alpha_vec)); oursSq=zeros(1,length(alpha_vec));
    for aa=1:length(alpha_vec)
        casename=sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss));
        ours(aa)=Results.(casename).Table2;
        oursSq(aa)=Results.(casename).SumSqLogDeviation;
    end
    fprintf('Ours                %s \n',sprintf('%8.1f',ours))
    fprintf('HVY2006             %s \n',sprintf('%8.1f',Table2paper.(sprintf('startage%i',startage_vec(ss)))))
    fprintf('(ours, mean squared log deviation, the calibration metric: %s) \n',sprintf('%8.2f',oursSq))
end
fprintf('End of outputs relating to Table 2 \n')

for ss=1:length(startage_vec)
    fprintf('\nTable 3 of HVY2006: percentage of variance in PV of earnings due to learning ability differences, parametric, accumulation starts at age %i \n',startage_vec(ss))
    fprintf('alpha               %s \n',alphastr)
    ours=zeros(1,length(alpha_vec));
    for aa=1:length(alpha_vec)
        casename=sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss));
        ours(aa)=Results.(casename).Table3;
    end
    fprintf('Ours                %s \n',sprintf('%8.1f',ours))
    fprintf('HVY2006             %s \n',sprintf('%8.1f',Table3paper.(sprintf('startage%i',startage_vec(ss)))))
end
fprintf('(PV of earnings as of age 20, for both starting ages, as in HVY2006) \n')
fprintf('End of outputs relating to Table 3 \n')

Table5rownames={'Mean (a)','Coef. of variation (a)','Skewness (a)','Mean (h1)','Coef. of variation (h1)','Skewness (h1)','Correlation (a,h1)'};
Table5panelnames={'A','B'};
for ss=1:length(startage_vec)
    fprintf('\nTable 5 of HVY2006: ability and human capital at birth, parametric, Panel %s: accumulation starts at age %i (ours, then HVY2006 in brackets) \n',Table5panelnames{ss},startage_vec(ss))
    fprintf('alpha                      %s \n',sprintf('%18.1f',alpha_vec))
    paperT5=Table5paper.(sprintf('startage%i',startage_vec(ss)));
    for rr=1:7
        rowstr='';
        for aa=1:length(alpha_vec)
            casename=sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss));
            rowstr=[rowstr,sprintf('%9.3f (%6.3f)',Results.(casename).Table5(rr),paperT5(rr,aa))];
        end
        fprintf('%-26s %s \n',Table5rownames{rr},rowstr)
    end
end
fprintf('End of outputs relating to Table 5 \n')

fprintf('\nTable 6 of HVY2006: persistence in individual earnings, alpha=0.7, parametric (HVY2006 report non-parametric, accumulation from age 10) \n')
Table6paper=[0.9999 0.9997; 0.9966 0.9854; 0.9679 0.8671; 0.9995 0.9994; 0.9960 0.9652; 0.9750 0.5229]; % rows as below; columns age 45, age 40
Table6rownames={'Corr(E_j,E_j-1)','Corr(E_j,E_j-5)','Corr(E_j,E_j-10)','Corr(z_j,z_j-1)','Corr(z_j,z_j-5)','Corr(z_j,z_j-10)'};
fprintf('%-20s %12s %12s %12s %12s %14s %14s \n','','start10 j=45','start10 j=40','start20 j=45','start20 j=40','HVY2006 j=45','HVY2006 j=40')
for rr=1:6
    fprintf('%-20s %12.4f %12.4f %12.4f %12.4f %14.4f %14.4f \n',Table6rownames{rr},Results.alpha07_start10.Table6(rr,1),Results.alpha07_start10.Table6(rr,2),Results.alpha07_start20.Table6(rr,1),Results.alpha07_start20.Table6(rr,2),Table6paper(rr,1),Table6paper(rr,2))
end
fprintf('(E_j is earnings and z_j=log(E_j/E_j-1) is earnings growth) \n')
fprintf('End of outputs relating to Table 6 \n')

CrossSectionpaper.data=[0.33, 1.16, 40.2, 25.1, 15.5, 4.7]; % HVY2006 Section 5.1
CrossSectionpaper.model=[0.327, 1.18, 40.9, 27.0, 17.5, 6.1]; % HVY2006 Section 5.1, non-parametric, alpha=0.7
fprintf('\nCross-section earnings distribution (HVY2006 Section 5.1): ages 20-58 pooled with 1%% population growth \n')
fprintf('%-46s %8s %9s %8s %8s %8s %8s \n','','Gini','Mean/Med','Top20%','Top10%','Top5%','Top1%')
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Data, ours',CrossSectionData)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Data, HVY2006',CrossSectionpaper.data)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, ours, parametric, alpha=0.7, from age 10',Results.alpha07_start10.CrossSection)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, ours, parametric, alpha=0.7, from age 20',Results.alpha07_start20.CrossSection)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, HVY2006, non-parametric, alpha=0.7',CrossSectionpaper.model)
fprintf('(Top shares are the percentage of all earnings received by the top 20%%, 10%%, 5%% and 1%% of earners) \n')
fprintf('End of outputs relating to the cross-section \n')

fprintf('\nCheck of HVY2006 Proposition 1(iii): baseline (alpha=0.7, from age 20) but everyone has the mean ability, a=%6.4f \n',Prop1ability)
fprintf('Condition a*A_j^alpha+(1-delta)*A_j>=A_j+1 at all ages (needed for the proposition): %i \n',Prop1condition)
fprintf('Gini of human capital at ages 20, 30, 40, 50, 58: %s \n',sprintf('%7.4f',Prop1GiniH([1,11,21,31,39])))
fprintf('Gini of earnings at ages 20, 30, 40, 50, 58:      %s \n',sprintf('%7.4f',Prop1GiniE([1,11,21,31,39])))
fprintf('Everyone has strictly positive earnings from age %i \n',Params.agejshifter+jjE(1))
fprintf('Largest amount by which a Lorenz curve becomes LESS equal from one age to the next (<=0, up to rounding, means the proposition holds): \n')
fprintf('   human capital, all ages: %9.2e;  earnings, from age %i: %9.2e \n',Prop1violationH,Params.agejshifter+jjE(1),Prop1violationE)
fprintf('End of outputs relating to Proposition 1 \n')

fprintf('\nCalibrated parameters, and check that the top of the h grid is not binding (max over ages of the mass within 5 grid points of h max): \n')
for cc=1:nCases
    casename=sprintf('alpha%02i_start%02i',round(10*CaseList(cc,1)),CaseList(cc,2));
    fprintf('%s: mean_logh1=%8.4f, mean_logability=%8.4f, stddev_logh1=%8.4f, stddev_logability=%8.4f, FTcorr=%8.4f (corr of logs %6.3f); mass near h max %8.2e \n',casename,Results.(casename).CalibParams.mean_logh1,Results.(casename).CalibParams.mean_logability,Results.(casename).CalibParams.stddev_logh1,Results.(casename).CalibParams.stddev_logability,Results.(casename).CalibParams.FTcorr_logh1logability,tanh(Results.(casename).CalibParams.FTcorr_logh1logability),Results.(casename).maxmassnearhmax)
end

save ./SavedOutput/HVY2006_Results.mat Results CaseList alpha_vec startage_vec Table2paper Table3paper Table5paper Table6paper CrossSectionData CrossSectionpaper Prop1ability Prop1condition Prop1GiniH Prop1GiniE Prop1LorenzH Prop1LorenzE Prop1violationH Prop1violationE


if CreateFigures==1
    % HVY2006 Figure 5: data, and the model (alpha=0.7) with accumulation starting at age 10 and at age 20,
    % age profiles of mean, Gini and mean/median earnings; axes and shape as HVY2006 Fig 5
    figure_c=figure_c+1;
    fig5=figure(figure_c);
    fig5.Position=[100 100 420 1190]; % same shape as HVY2006 Fig 5 (three stacked panels, each about 1.27 times as wide as tall)
    agedata=20:1:58;
    age10=Results.alpha07_start10.age; age20=Results.alpha07_start20.age;
    keep10=(age10>=20); % plot ages 20-58, where there is data
    subplot(3,1,1);
    plot(agedata,Fig1_Mean,'ok','MarkerSize',3); hold on;
    plot(age10(keep10),Results.alpha07_start10.earningsMean(keep10),'-b','LineWidth',1);
    plot(age20,Results.alpha07_start20.earningsMean,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,120]);
    xlabel('Age');
    title('Panel A: Mean Earnings');
    legend({'Data','Model, start age 10','Model, start age 20'},'Location','southeast');
    subplot(3,1,2);
    plot(agedata,Fig1_Gini,'ok','MarkerSize',3); hold on;
    plot(age10(keep10),Results.alpha07_start10.earningsGini(keep10),'-b','LineWidth',1);
    plot(age20,Results.alpha07_start20.earningsGini,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,0.45]);
    xlabel('Age');
    title('Panel B: Dispersion (Gini Coefficient)');
    subplot(3,1,3);
    plot(agedata,Fig1_Mean./Fig1_Median,'ok','MarkerSize',3); hold on;
    plot(age10(keep10),Results.alpha07_start10.earningsRatioMeanToMedian(keep10),'-b','LineWidth',1);
    plot(age20,Results.alpha07_start20.earningsRatioMeanToMedian,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,1.4]);
    xlabel('Age');
    title('Panel C: Skewness (Mean/Median)');
    saveas(fig5,'./SavedOutput/Graphs/HVY2006_Fig5_parametric.png');

    % Check of HVY2006 Proposition 1(iii): with a common ability, earnings dispersion falls with age
    figure_c=figure_c+1;
    figP1=figure(figure_c);
    figP1.Position=[100 100 640 480]; % panel shape as the HVY2006 figures
    plot(agedata,Fig1_Gini,'ok','MarkerSize',3); hold on;
    plot(age20,Results.alpha07_start20.earningsGini,'-b','LineWidth',1);
    plot(Params.agejshifter+(1:1:N_j),Prop1GiniE,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,0.45]); % Gini axis as HVY2006 Figs 1, 3, 4, 5
    xlabel('Age'); ylabel('Gini of earnings');
    title('Proposition 1(iii): dispersion with a common ability');
    legend({'Data','Model, baseline','Model, everyone has the mean ability'},'Location','southeast');
    saveas(figP1,'./SavedOutput/Graphs/HVY2006_Prop1.png');

    % The figures below are not in HVY2006; they are for the baseline case (alpha=0.7, accumulation from age 20)

    % Do a plot of the initial dist to be able to see it
    figure_c=figure_c+1;
    fig=figure(figure_c);
    subplot(2,1,1); plot(Params.ability,cumsum(Params.abilitydist),Params.ability,cumsum(sum(jequaloneDist,1)))
    title('ability dist')
    legend('from Params', 'from jequaloneDist')
    subplot(2,1,2); plot(a_grid,cumsum(sum(jequaloneDist,2)))
    title('initial human capital dist')

    saveas(fig,'./SavedOutput/Graphs/HVY2006_InitialDist.png')
    % Plot mean earnings conditional on ptype
    figure_c=figure_c+1;
    fig=figure(figure_c);
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype021.Mean)
    hold off
    title('Mean earnings conditional on ptype')
    % Plot median earnings conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype021.Median)
    hold off
    title('Median earnings conditional on ptype')


    saveas(fig,'./SavedOutput/Graphs/HVY2006_EarningsByAbility.png')
    % Plot mean timel conditional on ptype
    figure_c=figure_c+1;
    fig=figure(figure_c);
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype021.Mean)
    hold off
    title('Mean timel (time in human capital creation) conditional on ptype')
    % Plot median timel conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype021.Median)
    hold off
    title('Median timel (time in human capital creation) conditional on ptype')


    saveas(fig,'./SavedOutput/Graphs/HVY2006_TimeLearningByAbility.png')
    % Plot mean humancapital conditional on ptype
    figure_c=figure_c+1;
    fig=figure(figure_c);
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype021.Mean)
    hold off
    title('Mean human capital conditional on ptype')
    % Plot median humancapital conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype021.Median)
    hold off
    title('Median human capital conditional on ptype')

    saveas(fig,'./SavedOutput/Graphs/HVY2006_HumanCapitalByAbility.png')
    % So issue is coming from the human capital distribution
    figure_c=figure_c+1;
    fig=figure(figure_c);
    plot(h_grid,cumsum(StationaryDist.ptype001(:,1)))
    hold on
    for jj=2:N_j
        plot(h_grid,cumsum(StationaryDist.ptype001(:,jj)))
    end
    hold off
    saveas(fig,'./SavedOutput/Graphs/HVY2006_HumanCapitalDistByAge.png')
end


% Double-check: following is zero
cumsum(StationaryDist.ptweights)-cumsum(sum(jequaloneDist,1))'

save ./SavedOutput/HVY2006full.mat
end





%% doPart(5): the non-parametric case (HVY2006 Section 5.1, Fig 4, Table 4, and the non-parametric rows of Tables 2 and 3)
if doPart(5)==1
% The initial distribution can be any histogram on a rectangular grid of 20 points for initial human capital
% times 20 points for ability (400 bins). The model solution does not depend on the initial distribution, and
% with the closed-form solution (HVY2006 Proposition 1(ii), which reproduces the toolkit solution, see Table 6)
% each bin (h1,a) gives one deterministic earnings path. The earnings distribution at each age is then the 400
% paths weighted by the histogram, so we calibrate the 400 weights directly: they are parametrized as
% weights=exp(theta)/sum(exp(theta)), and theta minimizes the same objective as the parametric case (HVY2006
% footnote 18), using lsqnonlin (Levenberg-Marquardt, as there are more parameters than targets).
% The grids span +-2.5 standard deviations of log(h1) and log(a) from the parametric calibration of the same
% case, and the parametric calibration (evaluated on the grid) is the initial guess.
nh1_np=20; na_np=20; nbins=nh1_np*na_np;
optsNP=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','off','MaxFunctionEvaluations',1e6,'MaxIterations',5000,'FunctionTolerance',1e-10,'StepTolerance',1e-12); % FunctionTolerance 1e-10: with 1e-12 one case crawled to the iteration limit with no change in the fit

% HVY2006 Table 4 (non-parametric case); rows as in Table 5; columns: alpha=0.5,0.6,0.7,0.8,0.9
Table4paper.startage10=[0.466 0.319 0.209 0.139 0.087; 0.601 0.463 0.358 0.243 0.212; 1.303 1.190 1.183 1.168 1.103; 69.6 71.4 74.9 76.0 83.5; 0.456 0.453 0.422 0.397 0.261; 1.152 1.146 1.151 1.155 1.142; 0.10 0.205 0.305 0.397 0.418];
Table4paper.startage20=[0.453 0.320 0.210 0.134 0.089; 0.669 0.504 0.365 0.324 0.168; 1.251 1.188 1.147 1.131 1.111; 86.8 88.1 93.4 94.5 99.6; 0.475 0.486 0.510 0.457 0.501; 1.148 1.163 1.167 1.135 1.124; 0.621 0.689 0.781 0.792 0.741];
% HVY2006 Table 2 and Table 3, non-parametric rows
Table2paperNP.startage10=[3.5 3.2 2.6 2.5 2.8];
Table2paperNP.startage20=[3.1 3.5 2.8 3.9 3.8];
Table3paperNP.startage10=[78.2 73.8 70.5 68.4 73.1];
Table3paperNP.startage20=[62.7 63.1 63.5 64.4 65.1];

% Data targets, ages 20-58
Tmean=Fig1_Mean'; Tgini=Fig1_Gini'; Tmeanmedian=Fig1_Mean'./Fig1_Median';
nd=length(Tmean);

NonParam=struct();
for cc=1:nCases
    %% Set up this case
    alpha_c=CaseList(cc,1);
    startage_c=CaseList(cc,2);
    casename=sprintf('alpha%02i_start%02i',round(10*alpha_c),startage_c);
    N_j=58-startage_c+1; % 49 or 39 model periods, ending at age 58
    Params.J=N_j;
    Params.agej=1:1:N_j;
    Params.agejshifter=startage_c-1;
    Params.w=(1+Params.g).^(Params.agej-1); % wage rate (per unit of human-capital-hour)
    Params.alpha=alpha_c;
    jdata=(20:1:58)-Params.agejshifter; % model periods of the data ages 20-58
    jage20=20-Params.agejshifter;

    %% Grids for (h1,a), from the parametric calibration of this case
    cp=CalibParams.(casename);
    logh1grid_np=linspace(cp.mean_logh1-2.5*cp.stddev_logh1,cp.mean_logh1+2.5*cp.stddev_logh1,nh1_np)'; % column
    logagrid_np=linspace(cp.mean_logability-2.5*cp.stddev_logability,cp.mean_logability+2.5*cp.stddev_logability,na_np); % row
    h1grid_np=exp(logh1grid_np);
    agrid_np=exp(logagrid_np);

    %% Earnings path of each bin, closed-form solution
    lambda=zeros(1,N_j+1);
    for jj=N_j:-1:1
        lambda(jj)=Params.w(jj)+Params.beta*(1-Params.delta)*lambda(jj+1);
    end
    hnow=repmat(h1grid_np,1,na_np);
    Enp=zeros(nh1_np,na_np,N_j);
    for jj=1:N_j
        xstar=(Params.beta*lambda(jj+1)*Params.alpha*agrid_np/Params.w(jj)).^(1/(1-Params.alpha)); % 1-by-na_np
        xnow=min(hnow,xstar); % l=1 (all time learning) when xstar>=h
        Enp(:,:,jj)=Params.w(jj)*(hnow-xnow);
        hnow=agrid_np.*xnow.^Params.alpha+hnow*(1-Params.delta);
    end
    PV20=zeros(nh1_np,na_np); % present value of earnings as of age 20 (for Table 3)
    for jj=jage20:N_j
        PV20=PV20+Params.beta^(jj-jage20)*Enp(:,:,jj);
    end

    %% Earnings at each data age, sorted once (the sort order does not depend on the weights)
    Esorted=zeros(nbins,nd); sortindex=zeros(nbins,nd);
    for dd=1:nd
        [Esorted(:,dd),sortindex(:,dd)]=sort(reshape(Enp(:,:,jdata(dd)),[],1));
    end
    SortMat=sparse(1:nbins*nd,sortindex(:),1,nbins*nd,nbins); % SortMat*weights gives the weights in the sorted order of each age (stacked)

    %% Model statistics as functions of theta
    % Weights, sorted by earnings at each age: Ws (nbins-by-nd). Mean and Gini of the discrete distribution, and the
    % median by linear interpolation within the bin where the cumulative weight crosses 0.5 (so it is continuous in theta).
    WeightsFn=@(theta) exp(theta-max(theta))/sum(exp(theta-max(theta)));
    WsFn=@(theta) reshape(SortMat*WeightsFn(theta),nbins,nd);
    MeanFromWs=@(Ws) sum(Ws.*Esorted,1);
    GiniFromWs=@(Ws) 1-sum(Ws.*(2*cumsum(Ws.*Esorted,1)-Ws.*Esorted),1)./sum(Ws.*Esorted,1);
    MedianFromWs=@(Ws,cumWs,ix) Esorted(max(ix-1,1)+(0:nd-1)*nbins)+(ix>1).*(0.5-cumWs(max(ix-1,1)+(0:nd-1)*nbins))./Ws(ix+(0:nd-1)*nbins).*(Esorted(ix+(0:nd-1)*nbins)-Esorted(max(ix-1,1)+(0:nd-1)*nbins));
    MedianFn=@(theta) MedianFromWs(WsFn(theta),cumsum(WsFn(theta),1),sum(cumsum(WsFn(theta),1)<0.5,1)+1);
    % Residuals: log ratios of data to model for mean, dispersion and skewness at each age (HVY2006 footnote 18)
    ResidualsFn=@(theta) [log(Tmean./MeanFromWs(WsFn(theta))), log(Tgini./GiniFromWs(WsFn(theta))), log(Tmeanmedian./(MeanFromWs(WsFn(theta))./MedianFn(theta)))]';

    %% Initial guess: the parametric calibration (joint log-normal) evaluated on the grid
    rho=tanh(cp.FTcorr_logh1logability);
    zh=(logh1grid_np-cp.mean_logh1)/cp.stddev_logh1; za=(logagrid_np-cp.mean_logability)/cp.stddev_logability;
    density=exp(-(zh.^2-2*rho*zh.*za+za.^2)/(2*(1-rho^2)));
    weights0=density(:)/sum(density(:));
    theta0=log(max(weights0,1e-12));

    fprintf('\n========== Non-parametric calibration, case %i of %i: alpha=%1.1f, accumulation starts at age %i ==========\n',cc,nCases,alpha_c,startage_c)
    tic;
    [thetahat,resnorm,~,exitflag]=lsqnonlin(ResidualsFn,theta0,[],[],optsNP);
    fprintf('lsqnonlin exitflag %i, objective (sum of squared log ratios) %8.5f, initial guess %8.5f, time %5.1fs \n',exitflag,resnorm,sum(ResidualsFn(theta0).^2),toc)

    %% Results for this case
    weightshat=reshape(WeightsFn(thetahat),nh1_np,na_np); % on (h1,a)
    NonParam.(casename).h1grid=h1grid_np;
    NonParam.(casename).agrid=agrid_np;
    NonParam.(casename).weights=weightshat;
    NonParam.(casename).exitflag=exitflag;
    NonParam.(casename).objective=resnorm;
    Wshat=WsFn(thetahat);
    NonParam.(casename).earningsMean=MeanFromWs(Wshat);
    NonParam.(casename).earningsGini=GiniFromWs(Wshat);
    NonParam.(casename).earningsMedian=MedianFn(thetahat);
    NonParam.(casename).earningsRatioMeanToMedian=NonParam.(casename).earningsMean./NonParam.(casename).earningsMedian;
    % Table 2 (HVY2006 footnote 20)
    NonParam.(casename).Table2=100*sum(abs(ResidualsFn(thetahat)))/(3*nd);
    % Table 3: 100*var(E[PV|a])/var(PV), PV as of age 20
    Mean_PV=sum(weightshat(:).*PV20(:));
    mass_a=sum(weightshat,1);
    CondMean_PV=sum(weightshat.*PV20,1)./max(mass_a,eps);
    NonParam.(casename).Table3=100*sum(mass_a.*(CondMean_PV-Mean_PV).^2)/sum(weightshat(:).*(PV20(:)-Mean_PV).^2);
    % Table 4: ability and human capital at birth; skewness is mean/median, the median is that of the discrete distribution
    mass_h1=sum(weightshat,2);
    Mean_a=sum(mass_a.*agrid_np); Mean_h=sum(mass_h1.*h1grid_np);
    StdDev_a=sqrt(sum(mass_a.*(agrid_np-Mean_a).^2)); StdDev_h=sqrt(sum(mass_h1.*(h1grid_np-Mean_h).^2));
    Median_a=agrid_np(find(cumsum(mass_a)>=0.5,1)); Median_h=h1grid_np(find(cumsum(mass_h1)>=0.5,1));
    Corr_ah=sum(sum(weightshat.*(h1grid_np-Mean_h).*(agrid_np-Mean_a)))/(StdDev_a*StdDev_h);
    NonParam.(casename).Table4=[Mean_a; StdDev_a/Mean_a; Mean_a/Median_a; Mean_h; StdDev_h/Mean_h; Mean_h/Median_h; Corr_ah];
    % Check: mass on the edges of the (h1,a) grid (if it is not small, the grids should be wider)
    NonParam.(casename).massongridedge=sum(weightshat([1,end],:),'all')+sum(weightshat(2:end-1,[1,end]),'all');

    %% Cross-section (HVY2006 Section 5.1), for alpha=0.7, constructed as for the parametric case
    if alpha_c==0.7
        ages_cs=20:1:58;
        popweight=1.01.^(-(ages_cs-20)); popweight=popweight/sum(popweight);
        csvalues=zeros(nbins,length(ages_cs)); csweights=zeros(nbins,length(ages_cs));
        for aa=1:length(ages_cs)
            csvalues(:,aa)=reshape(Enp(:,:,ages_cs(aa)-Params.agejshifter),[],1);
            csweights(:,aa)=popweight(aa)*weightshat(:);
        end
        csvalues=csvalues(:); csweights=csweights(:);
        keep=(csweights>0); csvalues=csvalues(keep); csweights=csweights(keep)/sum(csweights(keep));
        [csvalues,sortindex_cs]=sort(csvalues); csweights=csweights(sortindex_cs);
        cumpop=[0; cumsum(csweights)]; cuminc=[0; cumsum(csweights.*csvalues)]; cuminc=cuminc/cuminc(end);
        csGini=1-sum((cuminc(2:end)+cuminc(1:end-1)).*csweights);
        csMedian=csvalues(find(cumpop(2:end)>=0.5,1));
        [cumpop_u,iu]=unique(cumpop,'last'); % (tiny weights can leave the cumulative population share unchanged in floating point)
        csTopShares=1-interp1(cumpop_u,cuminc(iu),1-[0.2,0.1,0.05,0.01]);
        NonParam.(casename).CrossSection=[csGini, sum(csweights.*csvalues)/csMedian, 100*csTopShares];
    end
end

%% Print the non-parametric results (ours, with the HVY2006 values alongside)
alphastr=sprintf('%8.1f',alpha_vec);
for ss=1:length(startage_vec)
    fprintf('\nTable 2 of HVY2006: mean absolute deviation (%%), NON-parametric, accumulation starts at age %i \n',startage_vec(ss))
    fprintf('alpha               %s \n',alphastr)
    ours=zeros(1,length(alpha_vec));
    for aa=1:length(alpha_vec)
        ours(aa)=NonParam.(sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss))).Table2;
    end
    fprintf('Ours                %s \n',sprintf('%8.1f',ours))
    fprintf('HVY2006             %s \n',sprintf('%8.1f',Table2paperNP.(sprintf('startage%i',startage_vec(ss)))))
end
fprintf('End of outputs relating to Table 2 (non-parametric) \n')

for ss=1:length(startage_vec)
    fprintf('\nTable 3 of HVY2006: percentage of variance in PV of earnings due to learning ability differences, NON-parametric, accumulation starts at age %i \n',startage_vec(ss))
    fprintf('alpha               %s \n',alphastr)
    ours=zeros(1,length(alpha_vec));
    for aa=1:length(alpha_vec)
        ours(aa)=NonParam.(sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss))).Table3;
    end
    fprintf('Ours                %s \n',sprintf('%8.1f',ours))
    fprintf('HVY2006             %s \n',sprintf('%8.1f',Table3paperNP.(sprintf('startage%i',startage_vec(ss)))))
end
fprintf('(PV of earnings as of age 20, for both starting ages, as in HVY2006) \n')
fprintf('End of outputs relating to Table 3 (non-parametric) \n')

Table5rownames={'Mean (a)','Coef. of variation (a)','Skewness (a)','Mean (h1)','Coef. of variation (h1)','Skewness (h1)','Correlation (a,h1)'};
Table5panelnames={'A','B'};
for ss=1:length(startage_vec)
    fprintf('\nTable 4 of HVY2006: ability and human capital at birth, NON-parametric, Panel %s: accumulation starts at age %i (ours, then HVY2006 in brackets) \n',Table5panelnames{ss},startage_vec(ss))
    fprintf('alpha                      %s \n',sprintf('%18.1f',alpha_vec))
    paperT4=Table4paper.(sprintf('startage%i',startage_vec(ss)));
    for rr=1:7
        rowstr='';
        for aa=1:length(alpha_vec)
            rowstr=[rowstr,sprintf('%9.3f (%6.3f)',NonParam.(sprintf('alpha%02i_start%02i',round(10*alpha_vec(aa)),startage_vec(ss))).Table4(rr),paperT4(rr,aa))];
        end
        fprintf('%-26s %s \n',Table5rownames{rr},rowstr)
    end
end
fprintf('End of outputs relating to Table 4 \n')

fprintf('\nCross-section earnings distribution (HVY2006 Section 5.1), NON-parametric, alpha=0.7 \n')
fprintf('%-46s %8s %9s %8s %8s %8s %8s \n','','Gini','Mean/Med','Top20%','Top10%','Top5%','Top1%')
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, ours, non-parametric, from age 10',NonParam.alpha07_start10.CrossSection)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, ours, non-parametric, from age 20',NonParam.alpha07_start20.CrossSection)
fprintf('%-46s %8.3f %9.2f %8.1f %8.1f %8.1f %8.1f \n','Model, HVY2006, non-parametric',[0.327, 1.18, 40.9, 27.0, 17.5, 6.1])
fprintf('End of outputs relating to the cross-section (non-parametric) \n')

fprintf('\nNon-parametric calibrations: objective, lsqnonlin exitflag, and mass on the edges of the (h1,a) grid: \n')
for cc=1:nCases
    casename=sprintf('alpha%02i_start%02i',round(10*CaseList(cc,1)),CaseList(cc,2));
    fprintf('%s: objective %8.5f, exitflag %i, mass on grid edges %8.2e \n',casename,NonParam.(casename).objective,NonParam.(casename).exitflag,NonParam.(casename).massongridedge)
end

save ./SavedOutput/HVY2006_NonParametric.mat NonParam CaseList alpha_vec startage_vec Table2paperNP Table3paperNP Table4paper

if CreateFigures==1
    % HVY2006 Figure 4: as Figure 5, but for the non-parametric case; axes and shape as HVY2006 Fig 4
    figure_c=figure_c+1;
    fig4=figure(figure_c);
    fig4.Position=[100 100 420 1190]; % same shape as HVY2006 Fig 4 (three stacked panels, each about 1.27 times as wide as tall)
    agedata=20:1:58;
    subplot(3,1,1);
    plot(agedata,Fig1_Mean,'ok','MarkerSize',3); hold on;
    plot(agedata,NonParam.alpha07_start10.earningsMean,'-b','LineWidth',1);
    plot(agedata,NonParam.alpha07_start20.earningsMean,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,120]);
    xlabel('Age');
    title('Panel A: Mean Earnings');
    legend({'Data','Model, start age 10','Model, start age 20'},'Location','southeast');
    subplot(3,1,2);
    plot(agedata,Fig1_Gini,'ok','MarkerSize',3); hold on;
    plot(agedata,NonParam.alpha07_start10.earningsGini,'-b','LineWidth',1);
    plot(agedata,NonParam.alpha07_start20.earningsGini,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,0.45]);
    xlabel('Age');
    title('Panel B: Dispersion (Gini Coefficient)');
    subplot(3,1,3);
    plot(agedata,Fig1_Mean./Fig1_Median,'ok','MarkerSize',3); hold on;
    plot(agedata,NonParam.alpha07_start10.earningsRatioMeanToMedian,'-b','LineWidth',1);
    plot(agedata,NonParam.alpha07_start20.earningsRatioMeanToMedian,'--r','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,1.4]);
    xlabel('Age');
    title('Panel C: Skewness (Mean/Median)');
    saveas(fig4,'./SavedOutput/Graphs/HVY2006_Fig4_nonparametric.png');
end
end




%% doPart(6): fixed initial human capital (HVY2006 Section 5.2, Fig 6)
if doPart(6)==1
% Everyone starts (at age 10, alpha=0.7) with the same human capital h1, and only ability differs: a histogram on
% a grid of 20 ability levels (HVY2006 footnote 21). The ability weights and the common h1 are chosen to best fit
% the same targets as before (HVY2006 footnote 18). For a given h1, each ability level gives one earnings path
% (closed-form solution, as in doPart(5)), so the 20 weights are calibrated exactly as in the non-parametric case
% (weights=exp(theta)/sum(exp(theta)), lsqnonlin); the common h1 is then chosen by a grid search.
% The ability grid spans +-2.5 standard deviations of log(a) from the parametric calibration (alpha=0.7, from age 10),
% and the grid for h1 spans +-2.5 standard deviations of log(h1) from the same calibration.
na_f6=20; nh1_f6=81;
optsF6=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','off','MaxFunctionEvaluations',1e6,'MaxIterations',5000,'FunctionTolerance',1e-10,'StepTolerance',1e-12); % FunctionTolerance 1e-10: with 1e-12 one case crawled to the iteration limit with no change in the fit

% The case: alpha=0.7, accumulation from age 10
N_j=49;
Params.J=N_j;
Params.agej=1:1:N_j;
Params.agejshifter=9;
Params.w=(1+Params.g).^(Params.agej-1);
Params.alpha=0.7;
jdata=(20:1:58)-Params.agejshifter; % model periods of the data ages 20-58
Tmean=Fig1_Mean'; Tgini=Fig1_Gini'; Tmeanmedian=Fig1_Mean'./Fig1_Median';
nd=length(Tmean);

cp=CalibParams.alpha07_start10;
logagrid_f6=linspace(cp.mean_logability-2.5*cp.stddev_logability,cp.mean_logability+2.5*cp.stddev_logability,na_f6); % row
agrid_f6=exp(logagrid_f6);
h1candidates=exp(linspace(cp.mean_logh1-2.5*cp.stddev_logh1,cp.mean_logh1+2.5*cp.stddev_logh1,nh1_f6));
% Initial guess for the ability weights: the marginal distribution of ability in the parametric calibration
weights0=exp(-((logagrid_f6-cp.mean_logability)/cp.stddev_logability).^2/2)';
weights0=weights0/sum(weights0);
theta0=log(max(weights0,1e-12));

lambda=zeros(1,N_j+1);
for jj=N_j:-1:1
    lambda(jj)=Params.w(jj)+Params.beta*(1-Params.delta)*lambda(jj+1);
end

F6objective=nan(1,nh1_f6);
F6theta=zeros(na_f6,nh1_f6);
for hh=1:nh1_f6
    %% Earnings path of each ability level, everyone starting with h1=h1candidates(hh)
    hnow=h1candidates(hh)*ones(1,na_f6);
    Ef6=zeros(na_f6,N_j);
    for jj=1:N_j
        xstar=(Params.beta*lambda(jj+1)*Params.alpha*agrid_f6/Params.w(jj)).^(1/(1-Params.alpha));
        xnow=min(hnow,xstar); % l=1 (all time learning) when xstar>=h
        Ef6(:,jj)=(Params.w(jj)*(hnow-xnow))';
        hnow=agrid_f6.*xnow.^Params.alpha+hnow*(1-Params.delta);
    end
    %% Model statistics as functions of theta, as in doPart(5)
    Esorted=zeros(na_f6,nd); sortindex=zeros(na_f6,nd);
    for dd=1:nd
        [Esorted(:,dd),sortindex(:,dd)]=sort(Ef6(:,jdata(dd)));
    end
    SortMat=sparse(1:na_f6*nd,sortindex(:),1,na_f6*nd,na_f6);
    WeightsFn=@(theta) exp(theta-max(theta))/sum(exp(theta-max(theta)));
    WsFn=@(theta) reshape(SortMat*WeightsFn(theta),na_f6,nd);
    MeanFromWs=@(Ws) sum(Ws.*Esorted,1);
    GiniFromWs=@(Ws) 1-sum(Ws.*(2*cumsum(Ws.*Esorted,1)-Ws.*Esorted),1)./sum(Ws.*Esorted,1);
    MedianFromWs=@(Ws,cumWs,ix) Esorted(max(ix-1,1)+(0:nd-1)*na_f6)+(ix>1).*(0.5-cumWs(max(ix-1,1)+(0:nd-1)*na_f6))./Ws(ix+(0:nd-1)*na_f6).*(Esorted(ix+(0:nd-1)*na_f6)-Esorted(max(ix-1,1)+(0:nd-1)*na_f6));
    MedianFn=@(theta) MedianFromWs(WsFn(theta),cumsum(WsFn(theta),1),sum(cumsum(WsFn(theta),1)<0.5,1)+1);
    ResidualsFn=@(theta) [log(Tmean./MeanFromWs(WsFn(theta))), log(Tgini./GiniFromWs(WsFn(theta))), log(Tmeanmedian./(MeanFromWs(WsFn(theta))./MedianFn(theta)))]';
    if all(isfinite(ResidualsFn(theta0)))
        [F6theta(:,hh),F6objective(hh)]=lsqnonlin(ResidualsFn,theta0,[],[],optsF6);
    end
end
[~,hhbest]=min(F6objective);

%% Results at the best common h1
hnow=h1candidates(hhbest)*ones(1,na_f6);
Ef6=zeros(na_f6,N_j);
for jj=1:N_j
    xstar=(Params.beta*lambda(jj+1)*Params.alpha*agrid_f6/Params.w(jj)).^(1/(1-Params.alpha));
    xnow=min(hnow,xstar);
    Ef6(:,jj)=(Params.w(jj)*(hnow-xnow))';
    hnow=agrid_f6.*xnow.^Params.alpha+hnow*(1-Params.delta);
end
weightsF6=exp(F6theta(:,hhbest)-max(F6theta(:,hhbest))); weightsF6=weightsF6/sum(weightsF6);
FixedH1.h1=h1candidates(hhbest);
FixedH1.agrid=agrid_f6;
FixedH1.weights=weightsF6;
FixedH1.objective=F6objective(hhbest);
FixedH1.age=Params.agejshifter+(1:1:N_j);
FixedH1.earningsMean=zeros(1,N_j); FixedH1.earningsGini=zeros(1,N_j); FixedH1.earningsMedian=zeros(1,N_j);
for jj=1:N_j
    [Es,si]=sort(Ef6(:,jj)); Ws=weightsF6(si);
    FixedH1.earningsMean(jj)=sum(Ws.*Es);
    cumWE=cumsum(Ws.*Es)/FixedH1.earningsMean(jj);
    FixedH1.earningsGini(jj)=1-sum(Ws.*(2*cumWE-Ws.*Es/FixedH1.earningsMean(jj)));
    ix=find(cumsum(Ws)>=0.5,1); cumW=cumsum(Ws);
    if ix>1
        FixedH1.earningsMedian(jj)=Es(ix-1)+(0.5-cumW(ix-1))/Ws(ix)*(Es(ix)-Es(ix-1));
    else
        FixedH1.earningsMedian(jj)=Es(1);
    end
end
FixedH1.earningsRatioMeanToMedian=FixedH1.earningsMean./FixedH1.earningsMedian;
FixedH1.Table2=100*sum(abs([log(Tmean./FixedH1.earningsMean(jdata)), log(Tgini./FixedH1.earningsGini(jdata)), log(Tmeanmedian./FixedH1.earningsRatioMeanToMedian(jdata))]))/(3*nd);
[~,jmin]=min(FixedH1.earningsGini(jdata)); FixedH1.GiniMinAge=19+jmin; % age (from 20) at which the Gini is lowest

fprintf('\nHVY2006 Fig 6: fixed initial human capital (common h1), non-parametric ability, alpha=0.7, accumulation from age 10 \n')
if hhbest>1 && hhbest<nh1_f6
    F6gridnote='at an interior grid point';
else
    F6gridnote='ON THE EDGE of the grid, so the grid should be wider';
end
fprintf('Best common h1=%8.3f (grid of %i values from %8.3f to %8.3f; the best is %s) \n',FixedH1.h1,nh1_f6,h1candidates(1),h1candidates(end),F6gridnote)
if ~exist('NonParam','var') % the non-parametric results, for comparison (from doPart(5))
    load('./SavedOutput/HVY2006_NonParametric.mat','NonParam')
end
fprintf('Objective (sum of squared log ratios) %8.5f; mean absolute deviation (%%, as Table 2) %6.2f, against %6.2f with h1 differences (non-parametric, alpha=0.7, from age 10) \n',FixedH1.objective,FixedH1.Table2,NonParam.alpha07_start10.Table2)
fprintf('Gini of earnings at ages 20, 25, 30, 40, 50, 58: model %s; data %s \n',sprintf('%6.3f',FixedH1.earningsGini(jdata([1,6,11,21,31,39]))),sprintf('%6.3f',Fig1_Gini([1,6,11,21,31,39])))
fprintf('The model Gini is lowest at age %i (HVY2006: the bottom of the U-shape is at about age 24) \n',FixedH1.GiniMinAge)
fprintf('End of outputs relating to Fig 6 \n')

save ./SavedOutput/HVY2006_FixedH1.mat FixedH1 h1candidates F6objective

if CreateFigures==1
    % HVY2006 Figure 6: data, and the model with a common initial human capital (from age 10, alpha=0.7); axes and shape as HVY2006 Fig 6
    figure_c=figure_c+1;
    fig6=figure(figure_c);
    fig6.Position=[100 100 420 1190]; % same shape as HVY2006 Fig 6 (three stacked panels, each about 1.27 times as wide as tall)
    agedata=20:1:58;
    subplot(3,1,1);
    plot(agedata,Fig1_Mean,'ok','MarkerSize',3); hold on;
    plot(agedata,FixedH1.earningsMean(jdata),'-b','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,140]);
    xlabel('Age');
    title('Panel A: Mean Earnings');
    legend({'Data','Model, start age 10'},'Location','southeast');
    subplot(3,1,2);
    plot(agedata,Fig1_Gini,'ok','MarkerSize',3); hold on;
    plot(agedata,FixedH1.earningsGini(jdata),'-b','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,0.5]);
    xlabel('Age');
    title('Panel B: Dispersion (Gini Coefficient)');
    subplot(3,1,3);
    plot(agedata,Fig1_Mean./Fig1_Median,'ok','MarkerSize',3); hold on;
    plot(agedata,FixedH1.earningsRatioMeanToMedian(jdata),'-b','LineWidth',1);
    hold off;
    xlim([20,58]); ylim([0,1.4]);
    xlabel('Age');
    title('Panel C: Skewness (Mean/Median)');
    saveas(fig6,'./SavedOutput/Graphs/HVY2006_Fig6.png');
end
end




%% doPart(7): write the LaTeX tables (for LatexDoc/App_HuggettVenturaYaron2006.tex), from the saved results
% Uses only the .mat files saved by doPart(2), (4), (5) and (6), so it needs no GPU and can be run on its own.
if doPart(7)==1
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
LTX4=load('./SavedOutput/HVY2006_Results.mat');
LTX5=load('./SavedOutput/HVY2006_NonParametric.mat');
LTX6=load('./SavedOutput/HVY2006_FixedH1.mat');
LTXdata=load('./SavedOutput/HVY2006_Data.mat','Fig1_Gini');
alphacolheads=sprintf(' & $\\alpha=%1.1f$',LTX4.alpha_vec);
panelnames={'Panel A: accumulation starts at age 10','Panel B: accumulation starts at age 20'};

%% Table 1: parameter values
FID=fopen('./SavedOutput/LatexInputs/HVY2006_Table1.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID,'Model periods & Interest rate & Rental growth & Depreciation rate & Production function \\\\ \\hline \n');
fprintf(FID,'$J=39,49$ & $r=%1.2f$ & $g=%1.4f$ & $\\delta=%1.4f$ & $\\alpha\\in\\{%s\\}$ \\\\ \n',0.04,0.0014,0.0114,strjoin(arrayfun(@(x) sprintf('%1.1f',x),LTX4.alpha_vec,'UniformOutput',false),','));
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

%% Tables 2 and 3: rows non-parametric and parametric (the document shows the HVY2006 original just above, so no HVY2006 values here)
for tt=2:3
    FID=fopen(sprintf('./SavedOutput/LatexInputs/HVY2006_Table%i.tex',tt),'w');
    fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l%s} \n \\hline \\hline \n',repmat('r',1,length(LTX4.alpha_vec)));
    fprintf(FID,'Case%s \\\\ \\hline \n',alphacolheads);
    for ss=1:2
        startage=LTX4.startage_vec(ss);
        fprintf(FID,'\\multicolumn{%i}{l}{\\emph{%s}} \\\\ \n',1+length(LTX4.alpha_vec),panelnames{ss});
        oursNP=zeros(1,length(LTX4.alpha_vec)); oursP=zeros(1,length(LTX4.alpha_vec));
        for aa=1:length(LTX4.alpha_vec)
            casename=sprintf('alpha%02i_start%02i',round(10*LTX4.alpha_vec(aa)),startage);
            if tt==2
                oursNP(aa)=LTX5.NonParam.(casename).Table2; oursP(aa)=LTX4.Results.(casename).Table2;
            else
                oursNP(aa)=LTX5.NonParam.(casename).Table3; oursP(aa)=LTX4.Results.(casename).Table3;
            end
        end
        fprintf(FID,'Non-parametric%s \\\\ \n',sprintf(' & %4.1f',oursNP));
        fprintf(FID,'Parametric%s \\\\ \n',sprintf(' & %4.1f',oursP));
    end
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);
end

%% Tables 4 (non-parametric) and 5 (parametric) (the document shows the HVY2006 original just above, so no HVY2006 values here)
statnames={'Mean ($a$)','Coef. of variation ($a$)','Skewness ($a$)','Mean ($h_1$)','Coef. of variation ($h_1$)','Skewness ($h_1$)','Correlation ($a,h_1$)'};
for tt=4:5
    FID=fopen(sprintf('./SavedOutput/LatexInputs/HVY2006_Table%i.tex',tt),'w');
    fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l%s} \n \\hline \\hline \n',repmat('r',1,length(LTX4.alpha_vec)));
    fprintf(FID,'Statistic%s \\\\ \\hline \n',alphacolheads);
    for ss=1:2
        startage=LTX4.startage_vec(ss);
        fprintf(FID,'\\multicolumn{%i}{l}{\\emph{%s}} \\\\ \n',1+length(LTX4.alpha_vec),panelnames{ss});
        for rr=1:7
            rowstr='';
            for aa=1:length(LTX4.alpha_vec)
                casename=sprintf('alpha%02i_start%02i',round(10*LTX4.alpha_vec(aa)),startage);
                if tt==4
                    ours=LTX5.NonParam.(casename).Table4(rr);
                else
                    ours=LTX4.Results.(casename).Table5(rr);
                end
                if rr==4 % mean h1 is a level, so fewer decimals
                    rowstr=[rowstr,sprintf(' & %5.1f',ours)];
                else
                    rowstr=[rowstr,sprintf(' & %5.3f',ours)];
                end
            end
            fprintf(FID,'%s%s \\\\ \n',statnames{rr},rowstr);
        end
    end
    fprintf(FID,'\\hline \n \\end{tabular*} \n');
    fclose(FID);
end

%% Table 6: persistence, alpha=0.7
corrnames={'Corr$(E_j,E_{j-1})$','Corr$(E_j,E_{j-5})$','Corr$(E_j,E_{j-10})$','Corr$(z_j,z_{j-1})$','Corr$(z_j,z_{j-5})$','Corr$(z_j,z_{j-10})$'};
FID=fopen('./SavedOutput/LatexInputs/HVY2006_Table6.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lrrrr} \n \\hline \\hline \n');
fprintf(FID,' & \\multicolumn{2}{c}{Parametric, from age 10} & \\multicolumn{2}{c}{Parametric, from age 20} \\\\ \n');
fprintf(FID,'Statistic & $j=45$ & $j=40$ & $j=45$ & $j=40$ \\\\ \\hline \n');
for rr=1:6
    if rr==1
        fprintf(FID,'\\multicolumn{5}{l}{\\emph{Panel A: correlation, levels}} \\\\ \n');
    elseif rr==4
        fprintf(FID,'\\multicolumn{5}{l}{\\emph{Panel B: correlation, growth rates}} \\\\ \n');
    end
    fprintf(FID,'%s & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n',corrnames{rr},LTX4.Results.alpha07_start10.Table6(rr,1),LTX4.Results.alpha07_start10.Table6(rr,2),LTX4.Results.alpha07_start20.Table6(rr,1),LTX4.Results.alpha07_start20.Table6(rr,2));
end
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

%% Cross-section earnings distribution (HVY2006 Section 5.1)
FID=fopen('./SavedOutput/LatexInputs/HVY2006_CrossSection.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrr} \n \\hline \\hline \n');
fprintf(FID,' & Gini & Mean/Median & Top 20\\%% & Top 10\\%% & Top 5\\%% & Top 1\\%% \\\\ \\hline \n');
csrows={'Data, replication',LTX4.CrossSectionData; 'Data, HVY2006',LTX4.CrossSectionpaper.data; ...
    'Model, parametric, from age 10',LTX4.Results.alpha07_start10.CrossSection; 'Model, parametric, from age 20',LTX4.Results.alpha07_start20.CrossSection; ...
    'Model, non-parametric, from age 10',LTX5.NonParam.alpha07_start10.CrossSection; 'Model, non-parametric, from age 20',LTX5.NonParam.alpha07_start20.CrossSection; ...
    'Model, non-parametric, HVY2006',LTX4.CrossSectionpaper.model};
for rr=1:size(csrows,1)
    fprintf(FID,'%s & %5.3f & %5.2f & %5.1f & %5.1f & %5.1f & %5.1f \\\\ \n',csrows{rr,1},csrows{rr,2});
end
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

%% Calibrated parameters (parametric case)
FID=fopen('./SavedOutput/LatexInputs/HVY2006_CalibParams.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}rrrrrrr} \n \\hline \\hline \n');
fprintf(FID,'$\\alpha$ & Start age & $\\mu_{\\log h_1}$ & $\\mu_{\\log a}$ & $\\sigma_{\\log h_1}$ & $\\sigma_{\\log a}$ & $\\rho_{\\log h_1,\\log a}$ \\\\ \\hline \n');
for ss=1:2
    for aa=1:length(LTX4.alpha_vec)
        casename=sprintf('alpha%02i_start%02i',round(10*LTX4.alpha_vec(aa)),LTX4.startage_vec(ss));
        cp=LTX4.Results.(casename).CalibParams;
        fprintf(FID,'%1.1f & %i & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',LTX4.alpha_vec(aa),LTX4.startage_vec(ss),cp.mean_logh1,cp.mean_logability,cp.stddev_logh1,cp.stddev_logability,tanh(cp.FTcorr_logh1logability));
    end
end
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

%% Gini of earnings by age: data, baseline model, common ability (Proposition 1), common initial human capital (Fig 6)
ages_g=[20,25,30,40,50,58];
FID=fopen('./SavedOutput/LatexInputs/HVY2006_GiniByAge.tex','w');
fprintf(FID,'\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l%s} \n \\hline \\hline \n',repmat('r',1,length(ages_g)));
fprintf(FID,'Age%s \\\\ \\hline \n',sprintf(' & %i',ages_g));
fprintf(FID,'Data%s \\\\ \n',sprintf(' & %5.3f',LTXdata.Fig1_Gini(ages_g-19)));
fprintf(FID,'Model, baseline (parametric, from age 20)%s \\\\ \n',sprintf(' & %5.3f',LTX4.Results.alpha07_start20.earningsGini(ages_g-19)));
fprintf(FID,'Model, common ability (Proposition 1)%s \\\\ \n',sprintf(' & %5.3f',LTX4.Prop1GiniE(ages_g-19)));
fprintf(FID,'Model, common $h_1$ (Fig 6, from age 10)%s \\\\ \n',sprintf(' & %5.3f',LTX6.FixedH1.earningsGini(ages_g-9)));
fprintf(FID,'\\hline \n \\end{tabular*} \n');
fclose(FID);

fprintf('doPart(7): LaTeX tables written to ./SavedOutput/LatexInputs/ \n')
end


%% Comment:
% What this code does is actually way more computationally complex than the
% exercise of HVY2006. The trick is to notice that none of the parameters
% being calibrated actually matter for the value and policy function. Hence
% HVY2006 solve the value/policy function once, and then do their
% calibration step just looping over the initial dist/stationary
% dist/life-cycle profile steps (ours includes the redundant value/policy
% fn step inside the iterations). Because the value/policy function step is 
% one of the harder parts computationally this makes a big difference. But 
% with 20 odd years of extra computing power we can anyway just
% solve the model even with this extra computation and it shows off
% how the calibration command in VFI Toolkit works.


%% End of run
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','on'); % do not leave an interactive session silently offscreen
end
fprintf('HuggettVenturaYaron2006: run finished (doPart=%s)\n',mat2str(doPart))
diary off






