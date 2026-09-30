%% Replication of Bruggemann (2021) - Higher Taxes at the Top: The Role of Entrepreneurs

% Record the full run (incl. any error) to B2021_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./B2021_Diary.txt','file'), delete('./B2021_Diary.txt'); end
diary ./B2021_Diary.txt

% Subcodes (the ReturnFn, the static-problem/income/tax function family, and the
% two alternative-economy scripts) live in ./B2021subcodes/; MATLAB does not
% search subfolders of the working directory, so put it on the path. Note this
% does NOT change cwd, so the './SavedOutput/...' paths still resolve to here.
% Only this subfolder goes on the path: AlessandroCodes/ and ReplicationMaterials/
% hold the original authors' code, which must not shadow anything.
addpath('./B2021subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
addpath(genpath('./VFIToolkit-matlab/'));

%% Setup
% TO-DO: roughly line 430: % Calculate Entrepreneur entry and exit: SHOULD THIS BE CONDITIONAL ON
    % BEING YOUNG? OR FOR WHOLE POPULATION?

use_jointgrid_fordecisionvariable=1


%% Model size
% Two decision variables
n_l=301; % labor supply (B2021 used 300 for labor supply, which seems pointlessly high but in the interests of replication I'm stuck doing the same)
n_entre=2; % 0=worker, 1=entrepreneur
% One endogenous states
n_asset=501; % assets (B2021 used 480 points, seems too few given she appears to do pure discretization; here we use gridinterplayer)
% Three markov exogenous states
n_age=2; % age (1=young & 2=old)
n_eta=6; % eta (labor productivity) I use eta to denote labor productivity, B2021 called it epsilon
n_theta=4; % theta (entrepreneurial ability)    
% But importantly, we only care about eta & theta when young, we do not
% need them when old. We will therefore set this up as a joint-grid.
% Note: The decision to be an entrepreneur or a worker in B2021 is a static decision
% Hence why entrepreneurship is not an endogenous state.
% Nor does entrepreneurship need to be a decision variable. But then it becomes complex
% to write out all the FnsToEvaluate as we have to constantly figure out if
% choosing entrepreneur or worker.
% Therefore I include entrepreneur-worker as a decision variable, this
% slows the code down (is still fast enough to do everything) but it
% makes all the codes much easier to read and substantially reduces the
% odds of making a mistake.

doPart=[1,1,1,1,1];
% 1: baseline model
% 2: stationary eqm for the different top tax rate
% 3: transition paths for the different top tax rates
% 4: further analysis of 60% top marginal tax rate (the one B2021 finds maximizes the social welfare fn)
% 5: alternative economies
doPart5=[1,1]; % there are two no-entpreneur economies, this allows doing just one
% 5A (NoEntrepreneur) completed 2026-09-21 and is loaded from B2021_doPart5A.mat

headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end
% Note: the root property set here is what actually does the work, and it is not a
% workspace variable, so the 'load ./SavedOutput/B2021_doPartN.mat' calls below cannot undo
% it. The headlessFigures value itself does get overwritten by those loads, which is
% why it rides along in B2021_doPart.mat next to doPart (see the else-branches below).

close all % close any figures, make sure they are all cleanly built from scratch
% Make sure the subfolders to save output exist
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

save ./SavedOutput/B2021_doPart.mat doPart headlessFigures


% "Share of hiring entrepreneurs" in Table 3 appears to be defined as the "Share of Entrepreneurs
% who hire", that is in the model the share for whom their labor demand
% exceeds their own personal labor supply. [This interpretation is based on
% the paragraph describing Table 6 that starts on pg 15 and ends on pg 16.]
%
% Table 6 reports employment by 'number of employees' that the entrepreneur hires, 
% but given that model is just units of time this is not really well defined. 
% Since the average labor supply in producitivity units of a single household 
% in the model is E[eta*l], we count this as the amount that corresponds to
% one employee.
%
% B2021 appears to force the choice of aprime onto the a_grid.
% benchmark.f90 creates 'agrid', and then on line 1678-82 the next period
% assets are chosen on this grid (line 1678: DO j1=1,da ! tomorrow's a').
% Her agrid is only 480 points, with hgrid (for labor supply) being 300
% points. Using pure discretized VFI with this few grid points on assets
% would be enough in Aiyagari, but very unlikely to be in models with 
% substantial wealth inequality, like these entrepreneurial-choice models 
% or like Castaneda, Diaz-Gimenez & Rio-Rull (2003). Note, this pure 
% discretized VFI is also in the codes of Cagetti & De Nardi (2009) which 
% B2021 builds from.
%
% Workers supply l units of labor, but these are times eta the
% productivity units. So labor supply of workers is l*eta. So you would
% expect that entrepreneurs own-labor supply to also be l*eta, but B2021
% just uses l [eqn 2 and 15 tell us that l^e goes into the entrepreneurs
% production, and that l^e=lbar. But there is then a typo in eqn 12, as it
% just uses l in utility, which should presumably be l=l^e. Notice though that this
% means eta 'disappears' when you become an entrepreneur, and it is not
% obvious that this setup shouldn't instead have l^e=lbar*eta in eqn 15.]
% You can see the impact of this in FnsToEvaluate.L below, which instead of
% just being l*eta has to be l*eta*(e==0)+l*(e==1), which seems a bit odd.
%
% When B2021 plots 'deciles', like in Figures 3,4,5 she groups the 20-30
% and the 30-40. Here I plot them separately.
%
% B2021 does 'Social Welfare', here the codes also compute
% 'Behind-the-Veil' welfare just out of interest.

%% Parameters

Params.headlessFigures=headlessFigures; % just so the setting is recorded alongside the run

% Preferences
Params.beta=0.9; % discount factor
Params.sigma1=1.5; % CRRA risk aversion
Params.sigma2=1.7; % Inverse of Frisch elasticity of labor supply
Params.ell=3; % household time endowment (so l=1 represents average hours of working 1/3 of the time; not really used, see creation of l_grid below)
Params.xi=0.716; % weight on disutility of labor supply

% Stochastic ageing
Params.pi_y=0.978; % probability of remaining young
Params.pi_o=0.911; % probability of remaining old

% Entrepreneur's production fn
Params.gamma=0.359; % relative importance of capital
Params.upsilon=0.864; % (decreasing) returns to scale [Lucas span-of-control parameter]
% Collateral constraint
Params.lambda=1.5;
% Entrepreneurs labor input
Params.lbar=1;

% Corporate sector production fn
Params.Z=1; % technology level (normalization; B2021 calls this A)
Params.alpha=0.33; % relative importance of capital

% Depreciation rate
Params.delta=0.06; % deprecation rate of capital (in both sectors)

% Taxes
Params.ybar=1; % average income (initial guess)
Params.tau_c=0.110; % consumption tax rate
Params.d=0.2385; % (times ybar) standard deduction (from income, before income taxes)
Params.tau_i_r1_stat=0.1; % Statutory Income tax rate, bracket 1
Params.tau_i_r2_stat=0.15; % Statutory Income tax rate, bracket 2
Params.tau_i_r3_stat=0.25; % Statutory Income tax rate, bracket 3
Params.tau_i_r4_stat=0.28; % Statutory Income tax rate, bracket 4
Params.tau_i_r5_stat=0.33; % Statutory Income tax rate, bracket 5
Params.tau_i_r6_stat=0.35; % Statutory Income tax rate, bracket 6
Params.tau_i_t1=0; % Income tax treshold 1
Params.tau_i_t2=0.214; % (times ybar) % Income tax treshold 2
Params.tau_i_t3=0.868; % (times ybar) % Income tax treshold 3
Params.tau_i_t4=1.753; % (times ybar) % Income tax treshold 4
Params.tau_i_t5=2.672; % (times ybar) % Income tax treshold 5
Params.tau_i_t6=4.771; % (times ybar) % Income tax treshold 6
Params.tau_i_adj=0.669; % linear scaling factor so that income tax raises 'right' revenue
Params.tau_s=0; % flat-tax on income representing state and local taxes (initial guess)
% B2021 page 7, eqns 4 & 5, says there is a tau_s for state and local taxes. 
% But then in Table 1 there is no tau_s. 
% B2021 pg 14 explains "the linear tax rate tau_s is endogenous and balancing the budget"

% Put the adjustment onto the statutory rates to get the ones the model uses
Params.tau_i_r1=Params.tau_i_adj*Params.tau_i_r1_stat; % Income tax rate, bracket 1
Params.tau_i_r2=Params.tau_i_adj*Params.tau_i_r2_stat; % Income tax rate, bracket 2
Params.tau_i_r3=Params.tau_i_adj*Params.tau_i_r3_stat; % Income tax rate, bracket 3
Params.tau_i_r4=Params.tau_i_adj*Params.tau_i_r4_stat; % Income tax rate, bracket 4
Params.tau_i_r5=Params.tau_i_adj*Params.tau_i_r5_stat; % Income tax rate, bracket 5
Params.tau_i_r6=Params.tau_i_adj*Params.tau_i_r6_stat; % Income tax rate, bracket 6
% Note: I did this because otherwise the 70% tax rate B2021 considers
% becomes a 0.7/0.669>1 tax rate, which is likely problemattic. If you wanted to
% calibrate tau_i_adj you could always set up the tau_i_r* in terms of
% tau_i_adj and tau_i_r*_stat using ParameterizeParamsFn

% Prices
Params.r=0.1; % interest rate (initial guess)
Params.w=0.7; % wage per-producitivity-unit-per-unit-of-time (initial guess)

% Government spending
Params.G=0.15; % government spending (initial guess)
Params.pension=0.4; % pension (initial guess) (B2021 calls this b)
% Calibration targets
Params.GdivYtarget=0.146; % used to determine G
Params.pensionreplacementrate=0.4; % used to determine pension [I assume in B2021 Table 1 where it puts this relative to y it should say ybar; pg 14 suggests it should be ybar]

% Lump-Sum transfers
Params.lumpsum=0; % 0 in baseline, but needed elsewhere to lump-sum transfer any additional tax revenues back to households (or negative if tax revenues are reduced)

%% Grids
% l_grid=linspace(0,Params.ell,n_l)'; % labor supply
% B2021 codes show that actually max l is 1.6 (B2021, benchmark.f90: maxl=1.6 is line 1132, line 1526 uses it to create hgrid, which is the grid on the labor supply)
Params.maxl=1.6;
l_grid=linspace(0,Params.maxl,n_l)'; % labor supply
% make sure lbar is a point in the grid [lbar is the fixed labor supply that entrepreneurs must provide]
[~,lbarindex]=min(abs(l_grid-Params.lbar));
l_grid(lbarindex)=Params.lbar;

entre_grid=[0;1]; % 0=worker, 1=entrepreneur

% Set grid for asset holdings
assetmaxfactor=520; % This is the max assets 
% B2021, benchmark.f90: maxa=520 on line 1126.
% But when I used this, no-one exceeds about 200 anyway, so is just wasting grid points. Regardless, I stick with it.
asset_grid=assetmaxfactor*(linspace(0,1,n_asset).^3)'; % linspace ^3 puts more points near zero, where the curvature of value and policy functions is higher

age_grid=[1;2]; % 1=young, 2=old

pi_age=[Params.pi_y, 1-Params.pi_y;...
    1-Params.pi_o, Params.pi_o]; % transitions between young and old

% B2021 pg 12 "I take the values for the first five levels of the labor ability process 
% from Cagetti & De Nardi (2009)... [and] introdce a high sixth level. [CDN2009 report 
% the grid and transition probabilities in their Appendix A]
eta_grid1to5=[0.2468, 0.4473, 0.7654, 1.3097, 2.3742]'; 
pi_eta1to5=[0.7376, 0.2473, 0.0150, 0.0002, 0.0000;....
    0.1947, 0.5555, 0.2328, 0.0169, 0.0001;...
    0.0113, 0.2221, 0.5333, 0.2221, 0.0113;...
    0.0001 0.0169 0.2328 0.5555 0.1947;...
    0.0000 0.0002 0.0150 0.2473 0.7376];
% Now add in the sixth point
Params.eta6=26; % value of 6th eta point
eta_grid=[eta_grid1to5; Params.eta6];
Params.prob_eta6=0.00160; % probability of going to 6th eta point
% 2026-09-25: was 0.002, which is B2021 Table 2's value ROUNDED to three decimals. Her
% Appendix A transition matrix gives 0.00160 in the sixth column, and her Fortran
% parameters.txt has !pi6 = 1.6005450E-03. The rounded value overstates it by 25%, and
% since eta6=26 is eleven times eta5=2.374 that lands squarely on the top tail: it raised
% the stationary share at eta6 by 24.4% (0.02740 vs 0.02202), mean eta by 8.7% and
% Gini(eta) by 4.6%. That is the likely source of our ybar sitting 6.2% above her
% equilibrium totinc, of Table 3's workers' income Gini at 0.58 against her 0.52, and of
% Figure 6 Panel A turning down too early (the no-entrepreneur economy is workers only).
Params.prob_eta6toeta3=0.071; % from 6th eta point, you either go to 3rd point, or remain in 6th
pi_eta=[pi_eta1to5*(1-Params.prob_eta6), Params.prob_eta6*ones(5,1);...
    0,0,Params.prob_eta6toeta3,0,0,1-Params.prob_eta6toeta3];
pi_eta=pi_eta./sum(pi_eta,2); % renormalize rows (some were 1.0001)

% B2021 Table 2 reports theta grid and transition probabilities
theta_grid=[0,0.682,1.750,2.818]'; % entrepreneurial ability
pi_theta=[0.963, 0.037, 0, 0;...
    0.275, 0.581, 0.144, 0;...
    0, 0.275, 0.581, 0.144;...
    0, 0, 0.304, 0.696]; % B2021 calls this Lambda

% We need the stationary distributions for both eta and theta, as it is
% assumed newborns draw their eta and theta from this stationary dist
[eta_mean,~,~,eta_statdist]=MarkovChainMoments(eta_grid,pi_eta);
[theta_mean,~,~,theta_statdist]=MarkovChainMoments(theta_grid,pi_theta);


%% Get into form for VFI toolkit
if use_jointgrid_fordecisionvariable==0
    n_d = [n_l,n_entre];
else
    n_d = [n_l+1,1]; % hardcodes n_entre=2
    if n_entre~=2
        error('joint grid for decision variable hardcodes n_entre=2')
    end
end
n_a    = n_asset;
n_z    = [n_eta*n_theta+1,1,1]; % joint-grid, n_eta*n_theta points for young, 1 point for old
if use_jointgrid_fordecisionvariable==0
    d_grid = [l_grid; entre_grid];
else
    % d_gridvals
    d_grid=[[l_grid; l_grid(lbarindex)],[entre_grid(1)*ones(n_l,1); entre_grid(2)]];
    % first n_l rows are worker, who can choose any of n_l labor supply values
    % last row is entrepreneur, who can only choose lbar as the sole labor supply value
end
a_grid = asset_grid;
% For z_grid, first 24 rows are young, last row is old
z_grid = [age_grid(1)*ones(n_eta*n_theta,1), repmat(eta_grid,n_theta,1), repelem(theta_grid,n_eta,1); 
          age_grid(2),0,0]; 
% note: use zero for eta and theta when old, so you will get issues if you try to do anything with them

% When agents are 'old' it becomes irrelevant what their eta & theta
% are, so we use a joint-grid on z and do not track them.
% When agents become 'young' the "new born household's two abilities are
% uncorrelated with the abilities of the parent household", so use 
% pi_newborn based on the stationary distributions of eta and theta (which are i.i.d)
% When agents are 'young', the transitions are based on pi_eta and pi_theta we created above.
pi_etatheta=kron(pi_theta, pi_eta); % in reverse order [young-young]

etatheta_statdist=kron(theta_statdist,eta_statdist); % in reverse order
pi_newborn=etatheta_statdist'; % i.i.d. based on etatheta_statdist

pi_z=[pi_age(1,1)*pi_etatheta, pi_age(1,2)*ones(n_eta*n_theta,1);...  % top left is young-young, top right is young-old, 
    pi_age(2,1)*pi_newborn, pi_age(2,2)]; % bottom left is old-young, bottom right is old-old

%% Return fn
DiscountFactorParamNames={'beta'};

ReturnFn=@(l,e,aprime,a,age,eta,theta,r,w,sigma1,sigma2,xi,lbar,lambda,delta,gamma,upsilon,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    Bruggemann2021_ReturnFn(l,e,aprime,a,age,eta,theta,r,w,sigma1,sigma2,xi,lbar,lambda,delta,gamma,upsilon,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);

%%
% Initial guesses for some parameters that will be determined in general eqm
Params.r       = 0.05;
Params.w       = 1.2; % 1.4
Params.tau_s   = 0.1; %0
Params.G       = 0.15;
Params.pension = 0.1; %0.4
Params.ybar    = 1; % 1


%% Try solving value fn
vfoptions.gridinterplayer=1; % set to 0 for debugging
vfoptions.ngridinterp=50;
vfoptions.maxaprimediff=20;
vfoptions.lowmemory=1;

tic;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime=toc

simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

% on transition path, use divide-and-conquer (preGI not a relevant option)
vfoptionstpath.gridinterplayer=vfoptions.gridinterplayer;
vfoptionstpath.ngridinterp=vfoptions.ngridinterp;
vfoptionstpath.divideandconquer=1;
vfoptionstpath.level1n=25;
vfoptionstpath.lowmemory=1; % my desktop ran out of memory without this

% Transition-path solver options. These live out here rather than inside the doPart(3)
% block because doPart(5) passes transpathoptions into the two alternative-economy
% subcodes, so a doPart=[0,0,0,0,1] run needs them to exist without Part 3 having run.
transpathoptions.GEnewprice=2; % 2=Anderson acceleration of the shooting map (3=plain shooting)
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice2.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'CapitalMarket','r',0,0.1;... % CapitalMarket is positive is r is to large, so subtract
    'LaborMarket','w',0,0.1;... % LaborMarket is positive is r is to large, so subtract
    'GovBudget','lumpsum',0,0.1}; % GovBudget is positive if lumpsum is too big, so subtract
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

% Now just run the TransitionPath_InfHorz command (all of the other inputs
% are things we had already had to define to be able to solve for the initial and final equilibria)
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;

% figure(1)
% subplot(3,1,1); plot(asset_grid,l_grid(Policy(1,:,1,3,4))) % labor supply
% subplot(3,1,2); plot(asset_grid,entre_grid(Policy(2,:,1,3,4)))
% subplot(3,1,3); plot(asset_grid,asset_grid(Policy(3,:,1,3,4))) % this is jsut lower grid point, but will do


%% Try out the stationary dist command
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions);

%% Aggregates

% Create functions to be evaluated
FnsToEvaluate.K_noncorp = @(l,e,aprime,a,age,eta,theta, r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_kFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % Assets used in non-corporate sector (=entrepreneurs)
FnsToEvaluate.A = @(l,e,aprime,a,age,eta,theta) a; % Total assets of households (workers and entrepreneurs)
FnsToEvaluate.N_noncorp = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % Labor used in non-corporate sector (=entrepreneurs), excluding own labor
FnsToEvaluate.N_lbar = @(l,e,aprime,a,age,eta,theta,lbar) lbar*e*(age==1); % Entrepreneurs own labor supply
FnsToEvaluate.L = @(l,e,aprime,a,age,eta,theta) l*eta*(e==0)+l*(e==1); % Total labor supply
FnsToEvaluate.Y_noncorp =  @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_YnoncorpFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % output of non-corporate sector
FnsToEvaluate.IncomeTaxRevenue = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    B2021_IncomeTaxRevenueFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6); % Tax Revenue from the Income Tax
FnsToEvaluate.ConsumptionTaxRevenue = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    tau_c*B2021_ConsumptionFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6); % Tax Revenue from the ConsumptionTax Tax
FnsToEvaluate.PensionSpending=@(l,e,aprime,a,age,eta,theta,pension) pension*(age==2);
FnsToEvaluate.Income=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)...
    B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension); % taxable income, before deductions. B2021 does not expicitly say how ybar is calculated so I am guessing this is it
FnsToEvaluate.Entrepreneur =  @(l,e,aprime,a,age,eta,theta) (e==1)*(age==1); % Not needed, just want to see it while finding the general eqm
FnsToEvaluate.Young =  @(l,e,aprime,a,age,eta,theta) (age==1); % Not needed, just want to see it while finding the general eqm

% Test the FnsToEvaluate
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

% IntermediateEqns
heteroagentoptions.intermediateEqns.N_corp=@(L,N_noncorp,N_lbar) L-(N_noncorp+N_lbar);
heteroagentoptions.intermediateEqns.K_corp=@(A,K_noncorp) A-K_noncorp;
heteroagentoptions.intermediateEqns.Y_corp=@(K_corp,N_corp,alpha,Z) Z*(K_corp^alpha)*(N_corp^(1-alpha));
heteroagentoptions.intermediateEqns.Y=@(Y_corp,Y_noncorp) Y_corp+Y_noncorp;


%% General equilbrium
% GE params
% r
% w
% tau_s (balance government budget; G+PensionSpending=TaxRevenue(income+consumption)

% Four calibration targets
% ybar is set to target the mean earnings
% G is set to target G/Y
% pension is to target the pension replacement rate
% tau_i_adj is to target the income tax revenue/GDP
% B2021 reports the calibrated value of tau_i_adj, but does not appear to report the exact target, so I do not do it here.

GEPriceParamNames={'r','w','tau_s','G','pension','ybar'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but not in transition)

GeneralEqmEqns.CapitalMarket = @(r,K_corp,N_corp,alpha,delta,Z) r-(alpha*Z*(K_corp/N_corp)^(alpha-1)-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.LaborMarket = @(w,K_corp,N_corp,alpha,Z) w-(1-alpha)*Z*(K_corp/N_corp)^(alpha); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.GovBudget = @(G,PensionSpending,IncomeTaxRevenue,ConsumptionTaxRevenue,lumpsum) G+PensionSpending+lumpsum-IncomeTaxRevenue-ConsumptionTaxRevenue; %Government runs balanced budget [take advantage of the fact that lumpsum is same for all, so adding up across everyone just gives the lumpsum parameter value; note, lumpsum=0 in baseline]
% Three calibration targets
GeneralEqmEqns.Pensions = @(pension,ybar,pensionreplacementrate) pension/ybar - pensionreplacementrate;
GeneralEqmEqns.AvgIncome = @(ybar,Income) ybar - Income; % get ybar to be the average income
GeneralEqmEqns.SizeOfGovernment = @(G,Y,GdivYtarget) G/Y-GdivYtarget; % government spending as percent of GDP
% Note: you can omit the SizeOfGovernment general eqm condition, since G does nothing in the model you can replace G with GdivYtarget in the
% GovBudget, and then after solving the general eqm you set G=GdivYtarget*Y. I don't do this just so that you can more easily change
% the codes without breaking them, and because they solve fast enough anyway.

heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.fminalgo=[8,1]; % fast but not really high accuracy, then a higher accuracy
heteroagentoptions.toleranceGEcondns=[1e-4,1e-5]; % high accuracy on final solve

% Set initial value for general eqm 
% (changed to decent initial guesses based on solution, the comments afterwards show the initial guess the first time I ran)
Params.r=0.02; %0.05
Params.w=1.3; % 1.2
Params.tau_s=0.1; %0.01
Params.G=0.3; %0.15
Params.pension=0.8; %0.1
Params.ybar=2.1; % 1

%% Solve for the initial/baseline stationary general equilbirium
% Create Tables 3,4,5,6
if doPart(1)==1
    heteroagentoptions.constrainpositive={'r','w','tau_s','pension','ybar'};
    
    heteroagentoptions.verbose=2; % Just while debugging
    
    [p_eqm_initial0,GeneralEqmCondn_initial0]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    
    Params.r=p_eqm_initial0.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_initial0.w;
    Params.tau_s=p_eqm_initial0.tau_s;
    Params.G=p_eqm_initial0.G;
    Params.pension=p_eqm_initial0.pension;
    Params.ybar=p_eqm_initial0.ybar;

    %% Re-solve, without the calibration targets to make sure we get a clean general eqm
    GEPriceParamNames={'r','w','tau_s'};
    %heteroagentoptions.constrainpositive={'r'};
    heteroagentoptions.constrainpositive={'w'};
    GeneralEqmEqns=rmfield(GeneralEqmEqns,'SizeOfGovernment');
    GeneralEqmEqns=rmfield(GeneralEqmEqns,'Pensions');
    GeneralEqmEqns=rmfield(GeneralEqmEqns,'AvgIncome');
    [p_eqm_initial,GeneralEqmCondn_initial]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

    Params.r=p_eqm_initial.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_initial.w;
    Params.tau_s=p_eqm_initial.tau_s;


    %% Now that we have the GE, let's calculate a bunch of related objects
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

    PolicyVals=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_z,d_grid,a_grid, vfoptions); % This will give you the policy in terms of values rather than index
    
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);

    %% keep a copy for later
    StationaryDist_initial=StationaryDist;
    V_initial=V;
    Policy_initial=Policy;
    Params_initial=Params;
    
    
    %% Check that we don't hit the top of asset grid (Fig12 is not one from Bruggemann 2021, just something I want to see)
    Fig13=figure(13);
    assetdist_young_e=cumsum(sum(StationaryDist(:,1:24).*shiftdim(PolicyVals(2,:,1:24),1),2),1); % recall: joint-grid on z, and e is second decision variable
    assetdist_young_note=cumsum(sum(StationaryDist(:,1:24).*shiftdim(1-PolicyVals(2,:,1:24),1),2),1); % recall: joint-grid on z, and e is second decision variable
    assetdist_old=cumsum(StationaryDist(:,25),1);
    plot(asset_grid,assetdist_young_note,asset_grid,assetdist_young_e,asset_grid,assetdist_old)
    title('cdf of HHs over assets')
    legend('worker','entrepreneur','retiree')
    saveas(Fig13,'./SavedOutput/Graphs/Bruggemann2021_FigGridCheck.png')
    

    %% Calculate various statistics for Table 3
    FnsToEvaluate.EntrepreneurHire = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar) (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>0); % n>0 (note: n excludes lbar)
    % I am guessing that this is the percentage of entrepreneurs who choose n>0
    % (hire labor not including their own personal labor supply)
    FnsToEvaluate.HoursWorked=@(l,e,aprime,a,age,eta,theta,lbar) l*(1-e)*(age==1)+lbar*e*(age==1); % l for workers and lbar for entrepreneurs
    
    simoptions.conditionalrestrictions.entrepreneurs=@(l,e,aprime,a,age,eta,theta) e*(age==1); % entrepreneurs
    simoptions.conditionalrestrictions.workers=@(l,e,aprime,a,age,eta,theta) (1-e)*(age==1); % workers
    simoptions.conditionalrestrictions.workerretirees=@(l,e,aprime,a,age,eta,theta) (age==2) + (age==1)*(1-e); % workers and retirees
    
    % To be able to calculate share of entrepreneurs among the top 1 percent incomes, we need the cutoff
    simoptions.nquantiles=100; % so we get the cutoff
    % To be easily able to calculate Top 1 percent income share, we want 100 points for lorenz curve
    simoptions.npoints=100;
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    simoptions=rmfield(simoptions,'nquantiles');
    simoptions=rmfield(simoptions,'npoints');

    % Count up the fraction of entrepreneurs among the top 1 percent [not entirely clear what should be done with the 'just retired entrepreneurs', guessing they are defined as not entrepreneurs]
    Params.toponepercentincomecutoff=AllStats.Income.QuantileCutoffs(100);
    FnsToEvaluate2.EntrepreneurTop1p=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,toponepercentincomecutoff)...
        (e==1)*(age==1)*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)>toponepercentincomecutoff);
    AllStats2=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate2,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

    % Just give some feedback so we can see nothing looks odd
    AllStats.Entrepreneur.Mean
    [AllStats.A.Mean, AllStats.L.Mean]
    [AllStats.entrepreneurs.A.Mean, AllStats.workers.A.Mean]
    [AllStats.A.Mean, AllStats.K_noncorp.Mean, AllStats.A.Mean-AllStats.K_noncorp.Mean] % assets, capital to entrepreneurs, capital to corporate
    [AllStats.L.Mean, AllStats.N_noncorp.Mean,AllStats.L.Mean-AllStats.N_lbar.Mean-AllStats.N_noncorp.Mean] % labor supply, labor to entrepreneurs, labor to corporate
    % Corporate labour nets out BOTH the labour entrepreneurs hire and their own lbar (which is
    % already inside Y_noncorp via their production fn); cf. the GE condition N_corp and line 913.
    Output_corp=Params.Z*((AllStats.A.Mean-AllStats.K_noncorp.Mean)^Params.alpha)*((AllStats.L.Mean-AllStats.N_lbar.Mean-AllStats.N_noncorp.Mean)^(1-Params.alpha));
    Y=Output_corp+AllStats.Y_noncorp.Mean;
    [Y,AllStats.Y_noncorp.Mean, Output_corp]

    % Calculate Entrepreneur entry and exit: SHOULD THIS BE CONDITIONAL ON
    % BEING YOUNG? OR FOR WHOLE POPULATION?
    simoptions.transprobs={'Entrepreneur'};
    CorrelationAndTransitionProbStats=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [], n_d, n_a, n_z, d_grid, a_grid,z_grid,pi_z,simoptions);
    % CorrelationAndTrasitionProbStats.Entrepreneur.TransitionProbs should give same entry/exit as below Panel Data simulation. It does! :)
    
    % Simulate a panel so we can see Entrepreneur entry and exit
    simoptions.simperiods=500;
    simoptions.numbersims=1000;
    SimPanelValues=SimPanelValues_InfHorz(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions);
    % Count entry, exit, and also 'exit of new entrants'
    NonEntrants=0;
    Entrants=0;
    NonExits=0;
    Exits=0;
    NewEntrantStays=0;
    NewEntrantExits=0;
    for ii=1:simoptions.numbersims
        for tt=3:simoptions.simperiods-1
            if SimPanelValues.Entrepreneur(tt-2,ii)==0
                if SimPanelValues.Entrepreneur(tt-1,ii)==1
                    Entrants=Entrants+1;
                    if SimPanelValues.Entrepreneur(tt,ii)==1
                        NewEntrantStays=NewEntrantStays+1;
                    elseif SimPanelValues.Entrepreneur(tt,ii)==0
                        NewEntrantExits=NewEntrantExits+1;
                    end                    
                elseif SimPanelValues.Entrepreneur(tt-1,ii)==0
                    NonEntrants=NonEntrants+1;
                end
            elseif SimPanelValues.Entrepreneur(tt-2,ii)==1
                if SimPanelValues.Entrepreneur(tt-1,ii)==1
                    NonExits=NonExits+1;
                elseif SimPanelValues.Entrepreneur(tt-1,ii)==0
                    Exits=Exits+1;
                end
            end
        end
    end
    % Out of interest, what is the exit rate of new entrants?
    ExitRateOfNewEntrants=NewEntrantExits/(NewEntrantExits+NewEntrantStays);
    fprintf('Exit rate of new entrants (out of interest, not in a table/figure): %1.3f \n', ExitRateOfNewEntrants)
    % Comment to self: with entrepreneurship as a decision rather than an endogenous state this exit rate of new 
    % entrants is only about half that of Kitao (2008), which also means it is about half that of US data as Kitao (2008) was close to data.
    % Entry and exit rates
    EntryRate=Entrants/(Entrants+NonEntrants);
    ExitRate=Exits/(Exits+NonExits);
    
    % Save copies of entry and exit rates for Table 10
    Table10_initial=[CorrelationAndTransitionProbStats.Entrepreneur.TransitionProbs(1,2), CorrelationAndTransitionProbStats.Entrepreneur.TransitionProbs(2,1), AllStats.Entrepreneur.Mean];
    
    % Table 3
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table3.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lll} \n');
    fprintf(FID, '\\multicolumn{3}{l}{Calibration targets and results} \\\\ \\hline  \n');
    fprintf(FID, '  & Target & Model  \\\\ \\hline \n');
    fprintf(FID, 'Overall Economy  &  &  \\\\ \n');
    fprintf(FID, '   Capital-output ratio  & 2.65 & %8.2f \\\\ \n', AllStats.A.Mean/Y);
    fprintf(FID, '   Top 1 percent income share  & 0.17 & %8.2f \\\\ \n', 1-AllStats.Income.LorenzCurve(99));
    % B2021's Fortran (benchmark.f90, subroutine medw_ratio, line 2586) builds the W group as
    % (prgridyw+prgridow)/(1-totentr), i.e. every non-entrepreneur INCLUDING retirees, not just
    % young workers. Her E group is prgridye/totentr, which is our entrepreneurs restriction.
    fprintf(FID, '   Ratio of median net worth E/W  & 7.26 & %8.2f \\\\ \n', AllStats.entrepreneurs.A.Median/AllStats.workerretirees.A.Median);
    fprintf(FID, '   Wealth Gini  & 0.85 & %8.2f \\\\ \n', AllStats.A.Gini);
    fprintf(FID, '   Average Federal Income Tax rate  & 0.12 & %8.2f \\\\ \n', AllStats.IncomeTaxRevenue.Mean/AllStats.Income.Mean);
    fprintf(FID, ' & &  \\\\ \n');
    fprintf(FID, 'Entrepreneurs & &  \\\\ \n');
    fprintf(FID, '   Fraction of entrepreneurs  & 0.07 & %8.2f \\\\ \n', AllStats.Entrepreneur.Mean);
    fprintf(FID, '   Entry rate & 0.02 & %8.2f \\\\ \n', EntryRate);
    fprintf(FID, '   Exit rate & 0.22 & %8.2f \\\\ \n', ExitRate);
    fprintf(FID, '   Entrepreneurs share of total income & 0.17 & %8.2f \\\\ \n', (AllStats.entrepreneurs.Income.Mean*AllStats.Entrepreneur.Mean)/AllStats.Income.Mean);
    fprintf(FID, '   Entrepreneurs income Gini & 0.65 & %8.2f \\\\ \n', AllStats.entrepreneurs.Income.Gini);
    fprintf(FID, '   Share of Entrepreneurs among top 1 percent income & 0.35 & %8.2f \\\\ \n', AllStats2.EntrepreneurTop1p.Mean*100); % the *100 is because mass of Top 1 percent is 0.01 by definition
    fprintf(FID, '   Share of entrepreneurs who hire & 0.66 & %8.2f \\\\ \n', AllStats.entrepreneurs.EntrepreneurHire.Mean); % Guessing this is their labor demand (including own labor) as share of total labor supply. Note, this is share of 'hours worked' rather than 'number of employees'
    fprintf(FID, ' & &  \\\\ \n');
    fprintf(FID, 'Workers & &  \\\\ \n');
    fprintf(FID, '   Average working time  & 1.00 & %8.2f \\\\ \n', AllStats.workers.HoursWorked.Mean);
    fprintf(FID, '   Workers income Gini  & 0.51 & %8.2f \\\\ \n', AllStats.workerretirees.Income.Gini);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: Targets not checked as part of replication. Renamed "Share of hiring entrepreneurs" to "Share of entrepreneurs who hire", as I find this clearer. \n');
    fprintf(FID, 'We include the revenue from $$tau_s$$ in the federal income tax revenue rates, since it is a federal tax and levied on income, original appears not to include it. \n');
    fprintf(FID, '"Ratio of median net worth E/W" is the median assets of entrepreneurs (young households who choose to be entrepreneurs) divided by the median assets of everyone else, \n');
    fprintf(FID, 'namely young workers together with retirees. Grouping the retirees in with the workers follows B2021, whose codes build this denominator as the non-entrepreneur share of the asset \n');
    fprintf(FID, 'distribution (benchmark.f90, subroutine medw\\_ratio), rather than as young workers alone. \n');
    fprintf(FID, 'The E group is subtle. Our entrepreneurs are the young households who \\textit{choose} to be entrepreneurs this period, because in this replication the \n');
    fprintf(FID, 'occupation is a decision rather than an endogenous state. B2021 instead has the occupation in the state, and her E group is the mass sitting in the \n');
    fprintf(FID, 'young-entrepreneur block of the state space (benchmark.f90, prgridye), that is, the households who chose to be entrepreneurs \\textit{last} period. \n');
    fprintf(FID, 'Her group therefore excludes this period entrants and includes this period exiters, while ours does the reverse. With an entry rate of 0.02 on a large \n');
    fprintf(FID, 'young-worker mass, entrants are roughly a fifth of the entrepreneur mass and arrive from the worker asset distribution, whereas her exiters are a \n');
    fprintf(FID, 'comparable fraction of asset-rich incumbents. Her E group is thus wealthier than ours by construction, which is the likely reason our ratio (5.7) \n');
    fprintf(FID, 'falls short of hers (7.7) even though both fall short of the 7.26 target. Reproducing her group would require adding the occupation to the state. \n');
    fprintf(FID, 'Two further details of her calculation push the other way: she takes the median as the last grid point whose cdf is at or below 0.50, one point below \n');
    fprintf(FID, 'the usual median and a larger downward bias for entrepreneurs because her 480-point asset grid is sparse at the top, and she divides prgridye by \n');
    fprintf(FID, 'totentr, which are two different masses, so her entrepreneur cdf overshoots one by about 1.5 percent. \n');
    fprintf(FID, 'The workers'' income Gini is over every household with no business capital, which following B2021 (table3.do classifies on kstar>0) means young \n');
    fprintf(FID, 'workers together with the retirees, not young workers alone. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    % Table 4
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table4.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccccc} \n');
    fprintf(FID, ' \\multicolumn{9}{c}{Distributions of Income and Wealth} \\\\ \\hline  \\hline  \n');
    fprintf(FID, ' & \\multicolumn{4}{c}{Income} & \\multicolumn{4}{c}{Wealth} \\\\ \\hline  \n');
    fprintf(FID, '  & Gini & Top 1\\%% & Top 5\\%% & Top 10\\%% &  Gini & Top 1\\%% & Top 5\\%% & Top 10\\%%  \\\\ \\hline \n');
    fprintf(FID, 'All Households  &  &  \\\\ \n');
    fprintf(FID, '   Data   & 0.549 & 17.2 & 33.7 & 44.4  & 0.846 & 34.1 & 60.9 & 74.4 \\\\ \n');
    fprintf(FID, '   Model  & %8.3f & %8.1f & %8.1f & %8.1f  & %8.3f & %8.1f & %8.1f & %8.1f \\\\ \n', AllStats.Income.Gini, 100*(1-AllStats.Income.LorenzCurve(99)), 100*(1-AllStats.Income.LorenzCurve(95)), 100*(1-AllStats.Income.LorenzCurve(90)),   AllStats.A.Gini, 100*(1-AllStats.A.LorenzCurve(99)), 100*(1-AllStats.A.LorenzCurve(95)), 100*(1-AllStats.A.LorenzCurve(90)));
    fprintf(FID, 'Entrepreneurs   &  &  \\\\ \n');
    fprintf(FID, '   Data   & 0.650 & 21.1 & 42.2 & 55.1  & 0.771 & 25.5 & 50.1 & 64.1 \\\\ \n');
    fprintf(FID, '   Model  & %8.3f & %8.1f & %8.1f & %8.1f  & %8.3f & %8.1f & %8.1f & %8.1f \\\\ \n', AllStats.entrepreneurs.Income.Gini, 100*(1-AllStats.entrepreneurs.Income.LorenzCurve(99)), 100*(1-AllStats.entrepreneurs.Income.LorenzCurve(95)), 100*(1-AllStats.entrepreneurs.Income.LorenzCurve(90)),   AllStats.entrepreneurs.A.Gini, 100*(1-AllStats.entrepreneurs.A.LorenzCurve(99)), 100*(1-AllStats.entrepreneurs.A.LorenzCurve(95)), 100*(1-AllStats.entrepreneurs.A.LorenzCurve(90)));
    fprintf(FID, 'Workers and Retirees  &  &  \\\\ \n');
    fprintf(FID, '   Data   & 0.518 & 14.5 & 29.7 & 40.6  & 0.832 & 31.4 & 57.3 & 71.3 \\\\ \n');
    fprintf(FID, '   Model  & %8.3f & %8.1f & %8.1f & %8.1f  & %8.3f & %8.1f & %8.1f & %8.1f \\\\ \n', AllStats.workerretirees.Income.Gini, 100*(1-AllStats.workerretirees.Income.LorenzCurve(99)), 100*(1-AllStats.workerretirees.Income.LorenzCurve(95)), 100*(1-AllStats.workerretirees.Income.LorenzCurve(90)),   AllStats.workerretirees.A.Gini, 100*(1-AllStats.workerretirees.A.LorenzCurve(99)), 100*(1-AllStats.workerretirees.A.LorenzCurve(95)), 100*(1-AllStats.workerretirees.A.LorenzCurve(90)));
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: Data not checked during replication \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    % Average income taxes for top percentiles of income
    FnsToEvaluate3.IncomeTaxRevenue=FnsToEvaluate.IncomeTaxRevenue;
    FnsToEvaluate3.Income=FnsToEvaluate.Income;
    % Can just get the means by centile and easily calc from these
    simoptions.nquantiles=100;
    AllStats3=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate3,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    simoptions=rmfield(simoptions,'nquantiles');
    
    Top1_AvgIncomeTaxRate=AllStats3.IncomeTaxRevenue.QuantileMeans(100)/AllStats3.Income.QuantileMeans(100);
    Top3_AvgIncomeTaxRate=sum(AllStats3.IncomeTaxRevenue.QuantileMeans(98:100))/sum(AllStats3.Income.QuantileMeans(98:100)); % note, because quantiles are equal mass, we just sum them to get the income of top 3, same for revenue
    Top5_AvgIncomeTaxRate=sum(AllStats3.IncomeTaxRevenue.QuantileMeans(96:100))/sum(AllStats3.Income.QuantileMeans(96:100)); % note, because quantiles are equal mass, we just sum them to get the income of top 5, same for revenue
    Top10_AvgIncomeTaxRate=sum(AllStats3.IncomeTaxRevenue.QuantileMeans(91:100))/sum(AllStats3.Income.QuantileMeans(91:100));  % note, because quantiles are equal mass, we just sum them to get the income of top 10, same for revenue

    % Table 5
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table5.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n');
    fprintf(FID, ' \\multicolumn{5}{c}{Average Tax Rates for Top Percentiles of the Income Distribution (\\%%)} \\\\ \\hline  \\hline  \n');
    fprintf(FID, '  & Top 1\\%% & Top 3\\%% & Top 5\\%% & Top 10\\%% \\\\ \\hline \n');
    fprintf(FID, '   Data   & 23.4 & 21.9 & 20.6 & 18.5  \\\\ \n');
    fprintf(FID, '   Model  & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Top1_AvgIncomeTaxRate, 100*Top3_AvgIncomeTaxRate, 100*Top5_AvgIncomeTaxRate, 100*Top10_AvgIncomeTaxRate);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: Data not checked during replication. \n');
    fprintf(FID, 'We include the revenue from $$tau_s$$ in the average tax rates, since it is a federal tax and levied on income, original appears not to include it. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    % Note that the model does not actualyl have a concept of someone
    % working or not. So any attempt to do 'number of employees' is a bit
    % hand-waving. B2021 resolves it by taking one efficiency unit of hired
    % labor to be one employee, and then binning n as if it were a headcount
    % rounded to the nearest integer, so her cutoffs are 5.5, 10.5 and 20.5
    % (Stata/Dofiles/table6.do). We follow her exactly: raw n, no rescaling by
    % average labor supply, and the same half-integer cutoffs.
    FnsToEvaluate4.Entrepreneur1to5e = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
        (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>0)*(B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)<=5.5); % B2021 bins raw n at 5.5 (table6.do)
    FnsToEvaluate4.Entrepreneur6to10e = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
        (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>5.5)*(B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)<=10.5); % B2021 bins raw n at 10.5 (table6.do)
    FnsToEvaluate4.Entrepreneur11to20e = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
        (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>10.5)*(B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)<=20.5); % B2021 bins raw n at 20.5 (table6.do)
    FnsToEvaluate4.Entrepreneur20eplus = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
        (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>20.5); % B2021 bins raw n at 20.5 (table6.do)
    AllStats4=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate4,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);


    Table6_basic=[AllStats.EntrepreneurHire.Mean/AllStats.Entrepreneur.Mean, AllStats4.Entrepreneur1to5e.Mean/AllStats.Entrepreneur.Mean, AllStats4.Entrepreneur6to10e.Mean/AllStats.Entrepreneur.Mean, AllStats4.Entrepreneur11to20e.Mean/AllStats.Entrepreneur.Mean, AllStats4.Entrepreneur20eplus.Mean/AllStats.Entrepreneur.Mean];

    % Table 6
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table6.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n');
    fprintf(FID, ' \\multicolumn{3}{c}{Firm Size Distribution: Data and Model} \\\\ \\hline  \\hline  \n');
    fprintf(FID, '  & Data & Model \\\\ \\hline \n');
    fprintf(FID, '   Share of entrepreneurs who hire & 0.661 & %8.3f \\\\ \n', Table6_basic(1)); % Guessing this is their labor demand (including own labor) as share of total labor supply. Note, this is share of 'hours worked' rather than 'number of employees'
    fprintf(FID, '   1-5 employees & 0.692 & %8.3f \\\\ \n', Table6_basic(2)/Table6_basic(1));
    fprintf(FID, '   6-10 employees & 0.119 & %8.3f \\\\ \n', Table6_basic(3)/Table6_basic(1));
    fprintf(FID, '   11-20 employees & 0.065 & %8.3f \\\\ \n', Table6_basic(4)/Table6_basic(1));
    fprintf(FID, '   More than 20 employees & 0.125 & %8.3f \\\\ \n', Table6_basic(5)/Table6_basic(1));
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: Data not checked during replication. Note: Rows 2-to-5 are expressed as a fraction of row 1 (so rows 2-5 sum to one). \n');
    fprintf(FID, 'The model has no notion of a headcount of employees, only of efficiency units of hired labor $n$. Following B2021 (Stata/Dofiles/table6.do) we take \n');
    fprintf(FID, 'one efficiency unit to be one employee and bin $n$ as a headcount rounded to the nearest integer, so the cutoffs are 5.5, 10.5 and 20.5 rather than \n');
    fprintf(FID, '5, 10 and 20, and $n$ is used as it is rather than rescaled by average labor supply. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    V_initial_check=ValueFnFromPolicy_InfHorz(Policy,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    max(abs(V_initial_check(:)-V_initial(:)))


    % For the CEV welfare evaluations later, we need
    DistOfYoung_initial=StationaryDist_initial(:,1:end-1)/sum(sum(StationaryDist_initial(:,1:end-1)));
    BehindTheVeil_initial=sum(sum(V_initial(:,1:end-1).*DistOfYoung_initial));
    SocialWelfare_initial=sum(sum(V_initial.*StationaryDist_initial));
    % And we need to calculate the lifetime utility of consumption (the value fn omitting the utility of leisure).
    % Easiest it just to set parameter chi to zero, so there is no utility of leisure
    xi_backup=Params.xi;
    Params.xi=0; % no utility of leisure, so consumption is the only source of utility
    V_ConsOnly_initial=ValueFnFromPolicy_InfHorz(Policy,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    % Put chi back into Params
    Params.xi=xi_backup;
    BehindTheVeil_ConsOnly_initial=sum(sum(V_ConsOnly_initial(:,1:end-1).*DistOfYoung_initial));
    SocialWelfare_ConsOnly_initial=sum(sum(V_ConsOnly_initial.*StationaryDist_initial));

    
    %% keep a copy for later
    AllStats_initial=AllStats;
    AllStats2_initial=AllStats2;
    Y_initial=Y;

    % Save only what Part 1 produces. Parts 2-5 each load their own file in turn, so
    % nothing has to be forwarded through here. Note Params, FnsToEvaluate, simoptions and
    % heteroagentoptions are all MODIFIED by Part 1 (the eqm prices are written back into
    % Params, two FnsToEvaluate are appended, conditionalrestrictions are added to
    % simoptions, and constrainpositive is narrowed to {'w'} for Part 2), so they must be
    % saved even though the setup section above also defines them.
    save ./SavedOutput/B2021_doPart1.mat Params Params_initial FnsToEvaluate simoptions heteroagentoptions ...
        p_eqm_initial GeneralEqmCondn_initial V_initial Policy_initial StationaryDist_initial V_ConsOnly_initial ...
        DistOfYoung_initial BehindTheVeil_initial BehindTheVeil_ConsOnly_initial ...
        SocialWelfare_initial SocialWelfare_ConsOnly_initial ...
        AllStats_initial AllStats2_initial Y_initial Table10_initial ...
        AllStats AllStats2 AllStats3 AllStats4 Table6_basic EntryRate ExitRate ExitRateOfNewEntrants
else
    load ./SavedOutput/B2021_doPart1.mat
    load ./SavedOutput/B2021_doPart.mat doPart headlessFigures
    Params.headlessFigures=headlessFigures; % the load above overwrote Params
end





%% Find the top marginal tax rate
% The Top Marginal Tax Rate in the baseline is tau_i_r6=tau_i_adj*tau_i_r6_stat, where the statutory rate is 35%
tau_TMTR_vec=0.35:0.05:1.00; % Consider tax rates from the baseline 35% up to 100%, doing 5 percentage point intervals
% Part 5 needs a longer grid than Parts 2-4: B2021's Figure 6 runs its x-axis out to an
% effective TMTR of 80 percent (Panel A) and 70 percent (Panel B), while 1.00 statutory is
% only 66.9 percent effective. 1.20 statutory = 80.3 percent effective covers both panels.
% Kept separate from tau_TMTR_vec so that Figure 1 and Parts 2-4 are unchanged.
tau_TMTR_vec_part5=0.35:0.05:1.20;
% Convert the statutory rates to the effective rates used by the model
tau_i_r6_vec=Params.tau_i_adj*tau_TMTR_vec;
% Note: 35% statutory rate is the initial eqm, so one of these reforms is 'nothing happens'
% Note: 100% statutory rate is a 70% effective rate, which is the largest tax rate considered in Figure 1 of B2021

% Change lumpsum to non-zero, as it now will be (is going to be found in general eqm)
Params.lumpsum=-0.1; % First is a lower tax rate, so transfers will be negative

%% Solve the stationary equilibria for different top marginal tax rates
if doPart(2)==1
    Figure1_CEV=zeros(length(tau_i_r6_vec),1);
    Figure1_CEV_BehindTheVeil=zeros(length(tau_i_r6_vec),1);
    
    GEPriceParamNames={'r','w','lumpsum'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but not in transition)
    
    clear GeneralEqmEqns
    GeneralEqmEqns.CapitalMarket = @(r,K_corp,N_corp,alpha,delta,Z) r-(alpha*Z*(K_corp^(alpha-1))*(N_corp^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    GeneralEqmEqns.LaborMarket = @(w,K_corp,N_corp,alpha,Z) w-(1-alpha)*Z*(K_corp^(alpha))*(N_corp^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    GeneralEqmEqns.GovBudget = @(G,PensionSpending,IncomeTaxRevenue,ConsumptionTaxRevenue,lumpsum) G+PensionSpending+lumpsum-IncomeTaxRevenue-ConsumptionTaxRevenue; %Government runs balanced budget [take advantage of the fact that lumpsum is same for all, so adding up across everyone just gives the lumpsum parameter value; note, lumpsum=0 in baseline]
    % B2021, pg 16-17: "When increasing the effective TMTR, I keep the level of government spending (including transfers to retirees) as
    % well as other tax parameters such as the tax brackets and standard deductions at their benchmark level. Any additional tax revenue
    % generated by the tax increase is redistributed through a tax-free lump-sum transfer to all households.
    
    FinalEqmResults=struct();
    
    for ii=1:length(tau_i_r6_vec)
        fprintf('Now doing final eqm %i of %i \n',ii,length(tau_i_r6_vec))
        Params.tau_i_r6=tau_i_r6_vec(ii);
        % In the first iteration the baseline eqm is used for r,w (and we set lumpsum to -0.1)
        % For ii>=2, just uses the previous iteration as the initial guess for the current ii

        [p_eqm_final,GeneralEqmCondn_final]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
        
        Params.r=p_eqm_final.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
        Params.w=p_eqm_final.w;
        Params.lumpsum=p_eqm_final.lumpsum;

        [V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
        StationaryDist_final=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
        
        FinalEqmResults(ii).p_eqm_final=p_eqm_final;
        FinalEqmResults(ii).GeneralEqmCondn_final=GeneralEqmCondn_final;
        FinalEqmResults(ii).V_final=V_final;
        FinalEqmResults(ii).Policy_final=Policy_final;
        FinalEqmResults(ii).StationaryDist_final=StationaryDist_final;
        FinalEqmResults(ii).Params_final=Params;

        % Consumption Equivalent Variation (welfare gain/loss)
        % We calculate both the 'behind-the-veil' expected lifetime utility in the baseline and each final eqm.
        DistOfYoung_final=StationaryDist_final(:,1:end-1)/sum(sum(StationaryDist_final(:,1:end-1)));
        BehindTheVeil_final=sum(sum(V_final(:,1:end-1).*DistOfYoung_final));
        % And the Social Welfare (which is what B2021 uses). The difference is that social welfare is based on the whole agent dist, while Behind-the-Veil is based on just the young.
        SocialWelfare_final=sum(sum(V_final.*StationaryDist_final));
        % B2021 ignores leisure keeping it constant, and this together with the fact
        % that utility is separable and the utility of consumption is a CES utility fn,
        % we get that CEV is given by the following simple formula
        WelfareChange_BehindTheVeil=BehindTheVeil_final-SocialWelfare_initial;
        CEV_BehindTheVeil=(WelfareChange_BehindTheVeil/BehindTheVeil_ConsOnly_initial+1)^(1/(1-Params.sigma1)) -1;
        WelfareChange=SocialWelfare_final-SocialWelfare_initial;
        CEV=(WelfareChange/SocialWelfare_ConsOnly_initial+1)^(1/(1-Params.sigma1)) -1;
        % Key to this formula being so simple is the strong assumption to
        % just ignore leisure; together with separable utility and a CES
        % utility of consumption.
        
        % Note: B2021 defines CEV in eqn (19) on page
        % Shows how to derive this simplified formula in Appendix C, on pg 32.

        Figure1_CEV(ii)=CEV;
        Figure1_CEV_BehindTheVeil(ii)=CEV_BehindTheVeil;
    end

    % V_final/Policy_final/StationaryDist_final/p_eqm_final are just the last loop iteration;
    % Parts 3 and 4 re-extract the ones they want from FinalEqmResults(ii), so they are not saved.
    save ./SavedOutput/B2021_doPart2.mat FinalEqmResults Figure1_CEV Figure1_CEV_BehindTheVeil
    % Check GeneralEqmCondn_final to be sure the transition all worked fine
else
    load ./SavedOutput/B2021_doPart2.mat
    load ./SavedOutput/B2021_doPart.mat doPart headlessFigures
    Params.headlessFigures=headlessFigures; % the load above overwrote Params
end


%% Solve the transition paths for the different top marginal tax rates
% Create Figure 1
if doPart(3)==1 
    Figure1_CEVpath=zeros(length(tau_i_r6_vec),1);
    FigureA1_LumpSumTPathPeriod1=zeros(length(tau_i_r6_vec),1);
    FigureA1_LumpSum=zeros(length(tau_i_r6_vec),1);

    % For this we need the following extra objects: PricePathOld,
    % PriceParamNames, ParamPath, ParamPathNames, T, V_final,
    % StationaryDist_initial (already calculated V_final & StationaryDist_initial above)
    
    % Number of time periods to allow for the transition (if you set T too low
    % it will cause problems, too high just means run-time will be longer).
    T=100; % Bruggemann (2021) graphs suggest she uses 100 periods
        
    FnsToEvaluate_TransPath.K_noncorp=FnsToEvaluate.K_noncorp;
    FnsToEvaluate_TransPath.A=FnsToEvaluate.A;
    FnsToEvaluate_TransPath.N_noncorp=FnsToEvaluate.N_noncorp;
    FnsToEvaluate_TransPath.N_lbar=FnsToEvaluate.N_lbar;
    FnsToEvaluate_TransPath.L=FnsToEvaluate.L;
    FnsToEvaluate_TransPath.IncomeTaxRevenue=FnsToEvaluate.IncomeTaxRevenue;
    FnsToEvaluate_TransPath.ConsumptionTaxRevenue=FnsToEvaluate.ConsumptionTaxRevenue;
    FnsToEvaluate_TransPath.PensionSpending=FnsToEvaluate.PensionSpending;

    % A version that contains additional model stats we want to look at but don't need for solving the transition path
    FnsToEvaluate_TransPath2=FnsToEvaluate_TransPath;
    FnsToEvaluate_TransPath2.Y_noncorp=FnsToEvaluate.Y_noncorp;

    
    TransPathGeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,N_lbar,L,alpha,delta,Z) r-(alpha*Z*((A-K_noncorp)^(alpha-1))*((L-N_lbar-N_noncorp)^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,N_lbar,L,alpha,Z) w-(1-alpha)*Z*((A-K_noncorp)^(alpha))*((L-N_lbar-N_noncorp)^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.GovBudget = @(G,PensionSpending,IncomeTaxRevenue,ConsumptionTaxRevenue,lumpsum) G+PensionSpending+lumpsum-IncomeTaxRevenue-ConsumptionTaxRevenue;  %Government runs balanced budget
    % Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.
    

    TPathResults=struct();
    
    % load ./SavedOutput/B2021_doPart3_ii.mat
    % for ii=14:length(tau_i_r6_vec)
    for ii=1:length(tau_i_r6_vec)
        fprintf('Now doing transition %i of %i \n',ii,length(tau_i_r6_vec))
        Params.tau_i_r6=tau_i_r6_vec(ii);
        Params.r=p_eqm_initial.r; % use benchmark economy
        Params.w=p_eqm_initial.w; % use benchmark economy
        Params.lumpsum=0; % use benchmark economy

        if ii==1 % Be exact rather than approximate in the case where initial and final are the same
            p_eqm_final.r=p_eqm_initial.r;
            p_eqm_final.w=p_eqm_initial.w;
            p_eqm_final.lumpsum=0;
            V_final=V_initial;
            Policy_final=Policy_initial;
        else
            p_eqm_final=FinalEqmResults(ii).p_eqm_final;
            V_final=FinalEqmResults(ii).V_final;
            Policy_final=FinalEqmResults(ii).Policy_final;
        end
        
        % We want to look at a one off unanticipated path of tau_i_r6 (the top marginal tax rate).
        ParamPath.tau_i_r6=Params.tau_i_r6*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
        % (the way ParamPath is set is designed to allow for a series of changes in the parameters)

        % We need to give an initial guess for the price path
        PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
        PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final.w, floor(T/3))'; p_eqm_final.w*ones(T-floor(T/3),1)];
        PricePath0.lumpsum=[linspace(0, p_eqm_final.lumpsum, floor(T/3))'; p_eqm_final.lumpsum*ones(T-floor(T/3),1)];
        
        % Solve for r, w and lumpsum together, straight from the PricePath0 guess above.
        transpathoptions.maxiter=300;
        [PricePath,GeneralEqmCondnPath]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions,vfoptionstpath);
        
        [VPath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);

        if ii==1 % ii=1 is the null reform: tau unchanged, prices and V_final all set to
            % the baseline, so the transition must reproduce the baseline exactly.
            fprintf('NULL-REFORM CHECK (ii=1, nothing changes; all three should be ~0)\n');
            fprintf('   max|VPath(:,:,1)-V_initial|  = %g\n', gather(max(abs(VPath(:,:,1)-V_initial),[],'all')));
            fprintf('   max|PricePath.r - r_initial| = %g\n', gather(max(abs(PricePath.r-p_eqm_initial.r))));
            fprintf('   PricePath.lumpsum(1)         = %g\n', gather(PricePath.lumpsum(1)));
        end
        % NB: transpathoptions goes in the 11th slot and simoptions in the 12th. Passing only 11
        % arguments put simoptions into the transpathoptions slot and left simoptions unset, so
        % the function defaulted simoptions.gridinterplayer=0 and iterated the whole agent
        % distribution with no interpolation layer, dumping mass that belongs between coarse
        % nodes j and j+1 onto j.
        AgentDistPath=AgentDistOnTransPath_InfHorz(StationaryDist_initial, PricePath,ParamPath,PolicyPath,n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
        AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate_TransPath2,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
        
        Y_corp_Path=Params.Z*((AggVarsPath.A.Mean-AggVarsPath.K_noncorp.Mean).^(Params.alpha)).*((AggVarsPath.L.Mean-AggVarsPath.N_lbar.Mean-AggVarsPath.N_noncorp.Mean).^(1-Params.alpha));
        YPath=Y_corp_Path+AggVarsPath.Y_noncorp.Mean;

        TPathResults(ii).ParamPath=ParamPath;
        TPathResults(ii).PricePath=PricePath;
        TPathResults(ii).GeneralEqmCondnPath=GeneralEqmCondnPath;
        TPathResults(ii).VPath=gather(VPath);
        TPathResults(ii).PolicyPath=gather(PolicyPath);
        TPathResults(ii).AgentDistPath=gather(AgentDistPath);
        TPathResults(ii).AggVarsPath=AggVarsPath;
        TPathResults(ii).YPath=YPath;

        % Consumption Equivalent Variation (welfare gain/loss)
        % We calculate both the Social-Welfare (over the agent dist) and the the 'Behind-the-Veil' (based on the 'newborns'/young) in the baseline and each final eqm.
        % This calculation is essentially the same as we did to compare stationary equilibria, just using first period of VPath instead of V_final.
        % See previous for explanation and details.
        BehindTheVeil_TPath=sum(sum(VPath(:,1:end-1,1,1,1).*DistOfYoung_initial));
        SocialWelfare_TPath=sum(sum(VPath(:,1:end,1,1,1).*StationaryDist_initial));

        WelfareChange_BehindTheVeil=BehindTheVeil_TPath-BehindTheVeil_initial;
        CEVpath_BehindTheVeil=(WelfareChange_BehindTheVeil/BehindTheVeil_ConsOnly_initial+1)^(1/(1-Params.sigma1)) -1;
        WelfareChange=SocialWelfare_TPath-SocialWelfare_initial;
        CEVpath=(WelfareChange/SocialWelfare_ConsOnly_initial+1)^(1/(1-Params.sigma1)) -1;

        Figure1_CEVpath(ii)=CEVpath;
        Figure1_CEVpath_BehindTheVeil(ii)=CEVpath_BehindTheVeil;

        % For Figure A1 we keep
        FigureA1_LumpSumTPathPeriod1(ii)=PricePath.lumpsum(1);
        FigureA1_LumpSum(ii)=p_eqm_final.lumpsum;

        % Crash-recovery checkpoint. Saving the whole workspace here wrote up to 337MB
        % fourteen times per run; TPathResults is the only large thing that belongs in it.
        save ./SavedOutput/B2021_doPart3_ii.mat TPathResults T ii ...
            Figure1_CEVpath Figure1_CEVpath_BehindTheVeil FigureA1_LumpSum FigureA1_LumpSumTPathPeriod1
    end
    
    % Figure 1: Welfare Maximizing Top Marginal Tax Rate
    fig1=figure(1);
    plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEVpath,'r--o')
    hold on
    plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEV,'b-o')
    hold off
    legend({'transition','stationary eqm comparison'},'Location','northeast')
    ylabel('CEV (in percent)')
    xlabel('Effective Top Marginal Tax Rate')
    title('Welfare-Maximizing Top Marginal Tax Rate')
    saveas(fig1,'./SavedOutput/Graphs/Bruggemann2021_Fig1.png')
    
    % Figure A1: Revenue Maximizing Top Marginal Tax Rate (Laffer Curve)
    fig7=figure(7);
    plot(Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSumTPathPeriod1/Params.ybar,'r--o')
    hold on
    plot(Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSum/Params.ybar,'b-o')
    hold off
    legend({'1st period of transition','stationary eqm comparison'},'Location','northeast')
    ylabel('Lump-sum transfer (in percent of average income)')
    xlabel('Effective Top Marginal Tax Rate')
    title('Revenue-Maximizing Top Marginal Tax Rate')
    saveas(fig7,'./SavedOutput/Graphs/Bruggemann2021_FigA1.png')


    % Re-Do these figures, but statutory tax rates
    % Figure 1: Welfare Maximizing Top Marginal Tax Rate
    fig10=figure(10);
    plot(tau_TMTR_vec,100*Figure1_CEVpath,'r--o')
    hold on
    plot(tau_TMTR_vec,100*Figure1_CEV,'b-o')
    hold off
    legend({'transition','stationary eqm comparison'},'Location','northeast')
    ylabel('CEV (in percent)')
    xlabel('Statutory Top Marginal Tax Rate')
    title('Welfare-Maximizing Top Marginal Tax Rate')
    saveas(fig10,'./SavedOutput/Graphs/Bruggemann2021_Fig1alt.png')

    % Figure A1: Revenue Maximizing Top Marginal Tax Rate (Laffer Curve)
    fig11=figure(11);
    plot(tau_TMTR_vec,100*FigureA1_LumpSumTPathPeriod1/Params.ybar,'r--o')
    hold on
    plot(tau_TMTR_vec,100*FigureA1_LumpSum/Params.ybar,'b-o')
    hold off
    legend({'1st period of transition','stationary eqm comparison'},'Location','northeast')
    ylabel('Lump-sum transfer (in percent of average income)')
    xlabel('Statutory Top Marginal Tax Rate')
    title('Revenue-Maximizing Top Marginal Tax Rate')
    saveas(fig11,'./SavedOutput/Graphs/Bruggemann2021_FigA1alt.png')


    % Re-Do with Behind-the-Veil Welfare
    % Figure 1: Welfare Maximizing Top Marginal Tax Rate
    fig12=figure(12);
    plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEVpath_BehindTheVeil,'r--o')
    hold on
    plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEV_BehindTheVeil,'b-o')
    hold off
    legend({'transition','stationary eqm comparison'},'Location','northeast')
    ylabel('CEV (in percent)')
    xlabel('Effective Top Marginal Tax Rate')
    title('Behind-the-Veil Welfare-Maximizing Top Marginal Tax Rate')
    saveas(fig12,'./SavedOutput/Graphs/Bruggemann2021_Fig1_BehindTheVeil.png')

    % Check GeneralEqmCondnPath to be sure the transition all worked fine

    % Split: the summary series are tiny and everything downstream reads them, while
    % TPathResults holds 14 full transition paths (VPath/PolicyPath/AgentDistPath at
    % 501-by-25-by-100 each) and is only needed by Part 4.
    save ./SavedOutput/B2021_doPart3.mat T Figure1_CEVpath Figure1_CEVpath_BehindTheVeil ...
        FigureA1_LumpSum FigureA1_LumpSumTPathPeriod1
    save ./SavedOutput/B2021_doPart3_TPathResults.mat TPathResults
else
    load ./SavedOutput/B2021_doPart3.mat
    load ./SavedOutput/B2021_doPart.mat doPart headlessFigures
    Params.headlessFigures=headlessFigures; % the load above overwrote Params
    if doPart(4)==1 % only Part 4 reads TPathResults, and it is large, so skip the load otherwise
        if exist('./SavedOutput/B2021_doPart3_TPathResults.mat','file')
            load ./SavedOutput/B2021_doPart3_TPathResults.mat
        end % a B2021_doPart3.mat saved before the split still holds TPathResults itself, hence the guard
        if ~exist('TPathResults','var')
            error('Part 4 needs TPathResults, but neither B2021_doPart3.mat nor B2021_doPart3_TPathResults.mat supplied it; re-run Part 3')
        end
    end
end


%% Lots more analysis of the welfare-maximizing top marginal tax rate
% Create Figures 2, 3, 4 & 5 and Tables 7, 8, 9 & 10 [Fig 2&3 do not explicitly say they are for the welfare-maximizing TMTR?]
if doPart(4)==1
    % We want to look closer at the impact of the welfare maximizing top marginal tax rate
    [~,ii]=max(Figure1_CEVpath);
    vv=tau_TMTR_vec(ii);
    Params.tau_i_r6_stat=vv;
    Params.tau_i_r6=Params.tau_i_adj*vv;

    % Load the final eqm
    p_eqm_final=FinalEqmResults(ii).p_eqm_final;
    V_final=FinalEqmResults(ii).V_final;
    Policy_final=FinalEqmResults(ii).Policy_final;
    StationaryDist_final=FinalEqmResults(ii).StationaryDist_final;
    Params_final=FinalEqmResults(ii).Params_final;

    % Load the path
    ParamPath=TPathResults(ii).ParamPath;
    PricePath=TPathResults(ii).PricePath;
    VPath=gpuArray(TPathResults(ii).VPath);
    PolicyPath=gpuArray(TPathResults(ii).PolicyPath);
    AgentDistPath=gpuArray(TPathResults(ii).AgentDistPath);
    YPath=TPathResults(ii).YPath;

    % Create a version of AggVarsPath that contains more FnsToEvaluate (things we didnt need when finding the general eqm)
    AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

    % take a look at these, should be essentially zero (say 1e-5)
    GeneralEqmCondn_final=FinalEqmResults(ii).GeneralEqmCondn_final;
    GeneralEqmCondnPath=TPathResults(ii).GeneralEqmCondnPath; 


    % For Figure 2, we need percent changes
    YPath_pch=100*(YPath-Y_initial)/Y_initial;
    APath_pch=100*(AggVarsPath.A.Mean-AllStats_initial.A.Mean)/AllStats_initial.A.Mean;
    LPath_pch=100*(AggVarsPath.L.Mean-AllStats_initial.L.Mean)/AllStats_initial.L.Mean;
    rPath_pch=100*(PricePath.r-Params_initial.r)/Params_initial.r;
    wPath_pch=100*(PricePath.w-Params_initial.w)/Params_initial.w;
    % Note: for lump-sum is as percent of 'average income'
    lumpsumPath_pch=100*PricePath.lumpsum/Params.ybar;

    % Figure 2
    fig2=figure(2);
    subplot(1,2,1); plot(1:1:T,YPath_pch,'k-',1:1:T,APath_pch,'r--',1:1:T,LPath_pch,'b-.')
    legend({'output (Y)','capital (K)','effective labor (N)'},'Location','southoutside')
    ylabel('Changes (in percent)')
    xlabel('Number of period since tax increase')
    subplot(1,2,2); 
    yyaxis left
    plot(1:1:T,rPath_pch,'k-',1:1:T,wPath_pch,'r--')
    ylabel('Price changes (in percent)')
    yyaxis right
    plot(1:1:T,lumpsumPath_pch,'b-.') % I'm just guessing that 'average income' here refers to ybar
    ylabel('Lump-sum transfer // (in percent of average income)')
    legend({'interest rate (r)','wage (w)','lump-sum transfer (on right y-axis)'},'Location','southoutside')
    xlabel('Number of period since tax increase')
    sgtitle('Changes in Aggregate Variables along the Transition')
    saveas(fig2,'./SavedOutput/Graphs/Bruggemann2021_Fig2.png')
    
    %% Entry and Exit
    % Calculate Entrepreneur entry and exit
    FnsToEvaluate5.Entrepreneur=FnsToEvaluate.Entrepreneur;
    simoptions.transprobs={'Entrepreneur'};
    
    % Entry and exit in the final eqm
    CorrelationAndTransitionProbStats=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist_final, Policy_final, FnsToEvaluate5,Params_final, [], n_d, n_a, n_z, d_grid, a_grid,z_grid,pi_z,simoptions);

    % Fraction of entrepreneurs in the final eqm (must be evaluated on the final eqm dist & policy, not the benchmark one)
    AggVars_final=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_final, Policy_final, FnsToEvaluate5, Params_final, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

    % Save copies of entry and exit rates for Table 10
    Table10_final=[CorrelationAndTransitionProbStats.Entrepreneur.TransitionProbs(1,2), CorrelationAndTransitionProbStats.Entrepreneur.TransitionProbs(2,1), AggVars_final.Entrepreneur.Mean];
    
    % Entrepreneur entry/exit rates over the transition path (B2021 replication just needs those from first period of transition, but easier to just calc the whole transition path)
    CorrelationAndTransitionProbStatsTPath=EvalFnOnTransPath_AutoCorrTransProbs_InfHorz(FnsToEvaluate5, AgentDistPath, PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,pi_z,simoptions);
    
    % Save copies of the entry and exit rates for first-to-second period of the transition path
    Table10_tpath=[CorrelationAndTransitionProbStatsTPath.Entrepreneur.TransitionProbs(1,2,1), CorrelationAndTransitionProbStatsTPath.Entrepreneur.TransitionProbs(2,1,1), AggVarsPath.Entrepreneur.Mean(1)];
    
    simoptions=rmfield(simoptions,'transprobs'); % clean up simoptions    
    
    %% What we need from the initial dist
    % Most the the figures we want to do are based on cutoffs for the agent distribution by income
    % So first, we just want to centile cutoffs for income
    simoptions.nquantiles=100;
    FnsToEvaluate_Income.Income=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)...
        B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension);
    AllStats_incomecutoffs=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_initial,Policy_initial,FnsToEvaluate_Income,Params_initial,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    incomecutoffs=AllStats_incomecutoffs.Income.QuantileCutoffs;
    simoptions=rmfield(simoptions,'nquantiles'); % clean up simoptions
    % Set up the cutoffs as parameters
    Params.incomecutoffs_10=incomecutoffs(10); % cutoff for 10th percentile of income dist
    Params.incomecutoffs_20=incomecutoffs(20);
    Params.incomecutoffs_30=incomecutoffs(30);
    Params.incomecutoffs_40=incomecutoffs(40);
    Params.incomecutoffs_50=incomecutoffs(50);
    Params.incomecutoffs_60=incomecutoffs(60);
    Params.incomecutoffs_70=incomecutoffs(70);
    Params.incomecutoffs_80=incomecutoffs(80);
    Params.incomecutoffs_90=incomecutoffs(90);
    Params.incomecutoffs_97=incomecutoffs(97);
    % need to also putthem into Params_initial
    Params_initial.incomecutoffs_10=Params.incomecutoffs_10;
    Params_initial.incomecutoffs_20=Params.incomecutoffs_20;
    Params_initial.incomecutoffs_30=Params.incomecutoffs_30;
    Params_initial.incomecutoffs_40=Params.incomecutoffs_40;
    Params_initial.incomecutoffs_50=Params.incomecutoffs_50;
    Params_initial.incomecutoffs_60=Params.incomecutoffs_60;
    Params_initial.incomecutoffs_70=Params.incomecutoffs_70;
    Params_initial.incomecutoffs_80=Params.incomecutoffs_80;
    Params_initial.incomecutoffs_90=Params.incomecutoffs_90;
    Params_initial.incomecutoffs_97=Params.incomecutoffs_97;

    
    %% Change in welfare by income deciles
    Figure3A_CEV=zeros(11,3); % (decile, worker/entrepreneur/old)
    Figure3B_CEV=zeros(11,3); % (decile, worker/entrepreneur/old)

    % I need indicators for the income deciles in terms of the state-space.
    % Easiest it to use ValuesOnGrid for income, and then use the cut-offs
    % with this.
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy_initial,FnsToEvaluate_Income,Params_initial,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    Indicator_IncomeDecile1=logical(ValuesOnGrid.Income<Params.incomecutoffs_10);
    Indicator_IncomeDecile2=logical((Params.incomecutoffs_10<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_20));
    Indicator_IncomeDecile3=logical((Params.incomecutoffs_20<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_30));
    Indicator_IncomeDecile4=logical((Params.incomecutoffs_30<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_40));
    Indicator_IncomeDecile5=logical((Params.incomecutoffs_40<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_50));
    Indicator_IncomeDecile6=logical((Params.incomecutoffs_50<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_60));
    Indicator_IncomeDecile7=logical((Params.incomecutoffs_60<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_70));
    Indicator_IncomeDecile8=logical((Params.incomecutoffs_70<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_80));
    Indicator_IncomeDecile9=logical((Params.incomecutoffs_80<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_90));
    Indicator_IncomePercentile90to97=logical((Params.incomecutoffs_90<=ValuesOnGrid.Income).*(ValuesOnGrid.Income<Params.incomecutoffs_97));
    Indicator_IncomePercentile97plus=logical(Params.incomecutoffs_97<=ValuesOnGrid.Income);

    % We also use indicators for entrepreneur/worker/old
    Indicator_Old=[zeros(1,n_z(1)-1),1].*ones(n_a,1);
    FnsToEvaluateIndicator.Entrepreneur=@(l,e,aprime,a,age,eta,theta) (e==1)*(age==1);
    FnsToEvaluateIndicator.Worker=@(l,e,aprime,a,age,eta,theta) (e==0)*(age==1);
    ValuesOnGrid2=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy_initial,FnsToEvaluateIndicator,Params_initial,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);


    % Now for each of these, calculate the CEV
    VPath_period1=VPath(:,:,:,:,1);

    for dd=1:11
        if dd==1
            Indicator_IncomeDecile=Indicator_IncomeDecile1;
        elseif dd==2
            Indicator_IncomeDecile=Indicator_IncomeDecile2;
        elseif dd==3
            Indicator_IncomeDecile=Indicator_IncomeDecile3;
        elseif dd==4
            Indicator_IncomeDecile=Indicator_IncomeDecile4;
        elseif dd==5
            Indicator_IncomeDecile=Indicator_IncomeDecile5;
        elseif dd==6
            Indicator_IncomeDecile=Indicator_IncomeDecile6;
        elseif dd==7
            Indicator_IncomeDecile=Indicator_IncomeDecile7;
        elseif dd==8
            Indicator_IncomeDecile=Indicator_IncomeDecile8;
        elseif dd==9
            Indicator_IncomeDecile=Indicator_IncomeDecile9;
        elseif dd==10
            Indicator_IncomeDecile=Indicator_IncomePercentile90to97;
        elseif dd==11
            Indicator_IncomeDecile=Indicator_IncomePercentile97plus;
        end

        % Recall that we have a joint-grid on z, so indexing is not trivial
        for pp=1:3
            if pp==1
                Indicator_pp=ValuesOnGrid2.Worker;
            elseif pp==2
                Indicator_pp=ValuesOnGrid2.Entrepreneur;
            elseif pp==3
                Indicator_pp=Indicator_Old;
            end

            jointIndicator=logical(Indicator_pp.*Indicator_IncomeDecile);

            SocialWelfare_initial=sum(V_initial(jointIndicator).*StationaryDist_initial(jointIndicator),1); 
            SocialWelfare_final=sum(V_final(jointIndicator).*StationaryDist_final(jointIndicator),1);
            SocialWelfare_ConsOnly_initial=sum(V_ConsOnly_initial(jointIndicator).*StationaryDist_initial(jointIndicator),1);
            Figure3B_CEV(dd,pp)=squeeze((SocialWelfare_final-SocialWelfare_initial)./SocialWelfare_ConsOnly_initial+1).^(1/(1-Params.sigma1)) -1;
            SocialWelfare_TPath=sum(VPath_period1(jointIndicator).*StationaryDist_initial(jointIndicator),1);
            Figure3A_CEV(dd,pp)=squeeze((SocialWelfare_TPath-SocialWelfare_initial)./SocialWelfare_ConsOnly_initial+1).^(1/(1-Params.sigma1)) -1;
        end
    end

    xaxissetup=["0 - 10%", "10 - 20%", "20 - 30%", "30 - 40%", "40 - 50%", "50 - 60%", "60 - 70%", "70 - 80%", "80 - 90%", "90 - 97%", "Top 3%"];

    % Figure 3: Change in welfare by income deciles
    fig3=figure(3);
    subplot(1,2,1); bar(xaxissetup, 100*Figure3A_CEV')
    legend({'Workers','Entrepreneurs','Old'},'Location','southoutside')
    ylabel('CEV (in percent)')
    title('Panel A: Transition')
    subplot(1,2,2); bar(xaxissetup, 100*Figure3B_CEV(:,1:2)')
    legend({'Workers','Entrepreneurs'},'Location','southoutside')
    ylabel('CEV (in percent)')
    title('Panel B: Stationary Eqm Comparison')
    saveas(fig3,'./SavedOutput/Graphs/Bruggemann2021_Fig3.png')
    % B2021 Figure 3 appears to be in meaningless utils. Here I do Figure 3 in CEV.
    % This shouldn't change the qualitative results (utils are ordinal) but
    % may change how much bigger or smaller they are (utils are not
    % cardinal).
    
    %% Change in various outcomes by income deciles
    
    % Define conditional restrictions for each of the income groups that form the x-axis of Figures 3,4,5
    % But need to do one version for workers (e=0)
    simoptions.conditionalrestrictions.wIncomeDecile1=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_10) ...
        (e==0)*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_10);
    simoptions.conditionalrestrictions.wIncomeDecile2=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_10,incomecutoffs_20) ...
        (e==0)*(incomecutoffs_10<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_20);
    simoptions.conditionalrestrictions.wIncomeDecile3=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_20,incomecutoffs_30) ...
        (e==0)*(incomecutoffs_20<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_30);
    simoptions.conditionalrestrictions.wIncomeDecile4=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_30,incomecutoffs_40) ...
        (e==0)*(incomecutoffs_30<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_40);
    simoptions.conditionalrestrictions.wIncomeDecile5=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_40,incomecutoffs_50) ...
        (e==0)*(incomecutoffs_40<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_50);
    simoptions.conditionalrestrictions.wIncomeDecile6=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_50,incomecutoffs_60) ...
        (e==0)*(incomecutoffs_50<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_60);
    simoptions.conditionalrestrictions.wIncomeDecile7=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_60,incomecutoffs_70) ...
        (e==0)*(incomecutoffs_60<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_70);
    simoptions.conditionalrestrictions.wIncomeDecile8=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_70,incomecutoffs_80) ...
        (e==0)*(incomecutoffs_70<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_80);
    simoptions.conditionalrestrictions.wIncomeDecile9=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_80,incomecutoffs_90) ...
        (e==0)*(incomecutoffs_80<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_90);
    simoptions.conditionalrestrictions.wIncomep90to97=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_90,incomecutoffs_97) ...
        (e==0)*(incomecutoffs_90<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_97);
    simoptions.conditionalrestrictions.wIncomepTop3=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_97) ...
        (e==0)*(incomecutoffs_97<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension));
    % and one version for entrepreneurs (e=1)
    simoptions.conditionalrestrictions.eIncomeDecile1=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_10) ...
        e*(age==1)*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_10);
    simoptions.conditionalrestrictions.eIncomeDecile2=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_10,incomecutoffs_20) ...
        e*(age==1)*(incomecutoffs_10<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_20);
    simoptions.conditionalrestrictions.eIncomeDecile3=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_20,incomecutoffs_30) ...
        e*(age==1)*(incomecutoffs_20<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_30);
    simoptions.conditionalrestrictions.eIncomeDecile4=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_30,incomecutoffs_40) ...
        e*(age==1)*(incomecutoffs_30<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_40);
    simoptions.conditionalrestrictions.eIncomeDecile5=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_40,incomecutoffs_50) ...
        e*(age==1)*(incomecutoffs_40<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_50);
    simoptions.conditionalrestrictions.eIncomeDecile6=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_50,incomecutoffs_60) ...
        e*(age==1)*(incomecutoffs_50<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_60);
    simoptions.conditionalrestrictions.eIncomeDecile7=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_60,incomecutoffs_70) ...
        e*(age==1)*(incomecutoffs_60<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_70);
    simoptions.conditionalrestrictions.eIncomeDecile8=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_70,incomecutoffs_80) ...
        e*(age==1)*(incomecutoffs_70<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_80);
    simoptions.conditionalrestrictions.eIncomeDecile9=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_80,incomecutoffs_90) ...
        e*(age==1)*(incomecutoffs_80<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_90);
    simoptions.conditionalrestrictions.eIncomep90to97=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_90,incomecutoffs_97) ...
        e*(age==1)*(incomecutoffs_90<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension))*(B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)<incomecutoffs_97);
    simoptions.conditionalrestrictions.eIncomepTop3=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,incomecutoffs_97) ...
        e*(age==1)*(incomecutoffs_97<=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension));
    
    % The functions we want to evaluate in each decile are for the most part already FnsToEvalute
    FnsToEvaluate6.Consumption= @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
        B2021_ConsumptionFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);
    FnsToEvaluate6.A=FnsToEvaluate.A;
    FnsToEvaluate6.L=FnsToEvaluate.L;
    FnsToEvaluate6.Y_noncorp=FnsToEvaluate.Y_noncorp;
    FnsToEvaluate6.K_noncorp=FnsToEvaluate.K_noncorp;
    FnsToEvaluate6.N_noncorp=FnsToEvaluate.N_noncorp;
    FnsToEvaluate6.N_lbar=FnsToEvaluate.N_lbar; % needed by the corporate/entrepreneurial labour split below
    FnsToEvaluate6.IncomeTaxRevenue=FnsToEvaluate.IncomeTaxRevenue;
    FnsToEvaluate6.ConsumptionTaxRevenue=FnsToEvaluate.ConsumptionTaxRevenue;
    FnsToEvaluate6.Entrepreneur=FnsToEvaluate.Entrepreneur;
    AllStats_initial=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_initial,Policy_initial,FnsToEvaluate6,Params_initial,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    % A conditional restriction with zero mass gets no stats fields (the toolkit writes them
    % only under RestrictedSampleMass>0), while the Figure 4/5 assignments below index them
    % directly. Fill NaN so an empty group plots as a gap instead of erroring.
    restrictnames=fieldnames(simoptions.conditionalrestrictions);
    statnames={'A','Consumption','L','K_noncorp','N_noncorp'};
    for rr=1:length(restrictnames)
        for ff=1:length(statnames)
            if AllStats_initial.(restrictnames{rr}).RestrictedSampleMass==0
                AllStats_initial.(restrictnames{rr}).(statnames{ff}).Mean=NaN;
            end
        end
    end
    % Calculate output
    Output_corp_initial=Params_initial.Z*((AllStats_initial.A.Mean-AllStats_initial.K_noncorp.Mean)^Params_initial.alpha)*((AllStats_initial.L.Mean-AllStats_initial.N_lbar.Mean-AllStats_initial.N_noncorp.Mean)^(1-Params_initial.alpha));
    Y_initial=Output_corp_initial+AllStats_initial.Y_noncorp.Mean;

    % Now calculate the same for the final eqm, and the (first period of) transition path.    
    
    %% Calculate some things from final stationary eqm
    Params.r=p_eqm_final.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_final.w;
    Params.lumpsum=p_eqm_final.lumpsum;
    % Params.tau_i_r6 is already set to the top marginal tax rate from reform

    [V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
    AllStats_final=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_final,Policy_final, FnsToEvaluate6,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    % A conditional restriction with zero mass gets no stats fields (the toolkit writes them
    % only under RestrictedSampleMass>0), while the Figure 4/5 assignments below index them
    % directly. Fill NaN so an empty group plots as a gap instead of erroring.
    restrictnames=fieldnames(simoptions.conditionalrestrictions);
    statnames={'A','Consumption','L','K_noncorp','N_noncorp'};
    for rr=1:length(restrictnames)
        for ff=1:length(statnames)
            if AllStats_final.(restrictnames{rr}).RestrictedSampleMass==0
                AllStats_final.(restrictnames{rr}).(statnames{ff}).Mean=NaN;
            end
        end
    end

    % Calculate output
    Output_corp_final=Params.Z*((AllStats_final.A.Mean-AllStats_final.K_noncorp.Mean)^Params.alpha)*((AllStats_final.L.Mean-AllStats_final.N_lbar.Mean-AllStats_final.N_noncorp.Mean)^(1-Params.alpha));
    Y_final=Output_corp_final+AllStats_final.Y_noncorp.Mean;
    
    %% Calculate some things from transition path (we only want first period, but easier to do the whole thing)
    tic;
    AllStatsPath=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate6,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    % A conditional restriction with zero mass gets no stats fields (the toolkit writes them
    % only under RestrictedSampleMass>0), while the Figure 4/5 assignments below index them
    % directly. Fill NaN so an empty group plots as a gap instead of erroring.
    restrictnames=fieldnames(simoptions.conditionalrestrictions);
    statnames={'A','Consumption','L','K_noncorp','N_noncorp'};
    for rr=1:length(restrictnames)
        for ff=1:length(statnames)
            if AllStatsPath.(restrictnames{rr}).RestrictedSampleMass(1)==0
                AllStatsPath.(restrictnames{rr}).(statnames{ff}).Mean(1)=NaN;
            end
        end
    end
    toc
    
    %% Figure 4: Changes in Average Choices by Occupation and Income Level:
    % Consumption (c), Savings (a), Hours worked (l)
    % Note: I report the change in the mean, not the mean of the change
    
    % Consumption (c)
    Figure4A=zeros(11,2,2); % decile, entrepreneur/worker, transition/stationary eqm
    % decile 1
    Figure4A(1,1,1)=(AllStats_final.eIncomeDecile1.Consumption.Mean-AllStats_initial.eIncomeDecile1.Consumption.Mean)/AllStats_initial.eIncomeDecile1.Consumption.Mean;
    Figure4A(1,1,2)=(AllStatsPath.eIncomeDecile1.Consumption.Mean(1)-AllStats_initial.eIncomeDecile1.Consumption.Mean)/AllStats_initial.eIncomeDecile1.Consumption.Mean;
    Figure4A(1,2,1)=(AllStats_final.wIncomeDecile1.Consumption.Mean-AllStats_initial.wIncomeDecile1.Consumption.Mean)/AllStats_initial.wIncomeDecile1.Consumption.Mean;
    Figure4A(1,2,2)=(AllStatsPath.wIncomeDecile1.Consumption.Mean(1)-AllStats_initial.wIncomeDecile1.Consumption.Mean)/AllStats_initial.wIncomeDecile1.Consumption.Mean;
    % decile 2
    Figure4A(2,1,1)=(AllStats_final.eIncomeDecile2.Consumption.Mean-AllStats_initial.eIncomeDecile2.Consumption.Mean)/AllStats_initial.eIncomeDecile2.Consumption.Mean;
    Figure4A(2,1,2)=(AllStatsPath.eIncomeDecile2.Consumption.Mean(1)-AllStats_initial.eIncomeDecile2.Consumption.Mean)/AllStats_initial.eIncomeDecile2.Consumption.Mean;
    Figure4A(2,2,1)=(AllStats_final.wIncomeDecile2.Consumption.Mean-AllStats_initial.wIncomeDecile2.Consumption.Mean)/AllStats_initial.wIncomeDecile2.Consumption.Mean;
    Figure4A(2,2,2)=(AllStatsPath.wIncomeDecile2.Consumption.Mean(1)-AllStats_initial.wIncomeDecile2.Consumption.Mean)/AllStats_initial.wIncomeDecile2.Consumption.Mean;
    % decile 3
    Figure4A(3,1,1)=(AllStats_final.eIncomeDecile3.Consumption.Mean-AllStats_initial.eIncomeDecile3.Consumption.Mean)/AllStats_initial.eIncomeDecile3.Consumption.Mean;
    Figure4A(3,1,2)=(AllStatsPath.eIncomeDecile3.Consumption.Mean(1)-AllStats_initial.eIncomeDecile3.Consumption.Mean)/AllStats_initial.eIncomeDecile3.Consumption.Mean;
    Figure4A(3,2,1)=(AllStats_final.wIncomeDecile3.Consumption.Mean-AllStats_initial.wIncomeDecile3.Consumption.Mean)/AllStats_initial.wIncomeDecile3.Consumption.Mean;
    Figure4A(3,2,2)=(AllStatsPath.wIncomeDecile3.Consumption.Mean(1)-AllStats_initial.wIncomeDecile3.Consumption.Mean)/AllStats_initial.wIncomeDecile3.Consumption.Mean;
    % decile 4
    Figure4A(4,1,1)=(AllStats_final.eIncomeDecile4.Consumption.Mean-AllStats_initial.eIncomeDecile4.Consumption.Mean)/AllStats_initial.eIncomeDecile4.Consumption.Mean;
    Figure4A(4,1,2)=(AllStatsPath.eIncomeDecile4.Consumption.Mean(1)-AllStats_initial.eIncomeDecile4.Consumption.Mean)/AllStats_initial.eIncomeDecile4.Consumption.Mean;
    Figure4A(4,2,1)=(AllStats_final.wIncomeDecile4.Consumption.Mean-AllStats_initial.wIncomeDecile4.Consumption.Mean)/AllStats_initial.wIncomeDecile4.Consumption.Mean;
    Figure4A(4,2,2)=(AllStatsPath.wIncomeDecile4.Consumption.Mean(1)-AllStats_initial.wIncomeDecile4.Consumption.Mean)/AllStats_initial.wIncomeDecile4.Consumption.Mean;
    % decile 5
    Figure4A(5,1,1)=(AllStats_final.eIncomeDecile5.Consumption.Mean-AllStats_initial.eIncomeDecile5.Consumption.Mean)/AllStats_initial.eIncomeDecile5.Consumption.Mean;
    Figure4A(5,1,2)=(AllStatsPath.eIncomeDecile5.Consumption.Mean(1)-AllStats_initial.eIncomeDecile5.Consumption.Mean)/AllStats_initial.eIncomeDecile5.Consumption.Mean;
    Figure4A(5,2,1)=(AllStats_final.wIncomeDecile5.Consumption.Mean-AllStats_initial.wIncomeDecile5.Consumption.Mean)/AllStats_initial.wIncomeDecile5.Consumption.Mean;
    Figure4A(5,2,2)=(AllStatsPath.wIncomeDecile5.Consumption.Mean(1)-AllStats_initial.wIncomeDecile5.Consumption.Mean)/AllStats_initial.wIncomeDecile5.Consumption.Mean;
    % decile 6
    Figure4A(6,1,1)=(AllStats_final.eIncomeDecile6.Consumption.Mean-AllStats_initial.eIncomeDecile6.Consumption.Mean)/AllStats_initial.eIncomeDecile6.Consumption.Mean;
    Figure4A(6,1,2)=(AllStatsPath.eIncomeDecile6.Consumption.Mean(1)-AllStats_initial.eIncomeDecile6.Consumption.Mean)/AllStats_initial.eIncomeDecile6.Consumption.Mean;
    Figure4A(6,2,1)=(AllStats_final.wIncomeDecile6.Consumption.Mean-AllStats_initial.wIncomeDecile6.Consumption.Mean)/AllStats_initial.wIncomeDecile6.Consumption.Mean;
    Figure4A(6,2,2)=(AllStatsPath.wIncomeDecile6.Consumption.Mean(1)-AllStats_initial.wIncomeDecile6.Consumption.Mean)/AllStats_initial.wIncomeDecile6.Consumption.Mean;
    % decile 7
    Figure4A(7,1,1)=(AllStats_final.eIncomeDecile7.Consumption.Mean-AllStats_initial.eIncomeDecile7.Consumption.Mean)/AllStats_initial.eIncomeDecile7.Consumption.Mean;
    Figure4A(7,1,2)=(AllStatsPath.eIncomeDecile7.Consumption.Mean(1)-AllStats_initial.eIncomeDecile7.Consumption.Mean)/AllStats_initial.eIncomeDecile7.Consumption.Mean;
    Figure4A(7,2,1)=(AllStats_final.wIncomeDecile7.Consumption.Mean-AllStats_initial.wIncomeDecile7.Consumption.Mean)/AllStats_initial.wIncomeDecile7.Consumption.Mean;
    Figure4A(7,2,2)=(AllStatsPath.wIncomeDecile7.Consumption.Mean(1)-AllStats_initial.wIncomeDecile7.Consumption.Mean)/AllStats_initial.wIncomeDecile7.Consumption.Mean;
    % decile 8
    Figure4A(8,1,1)=(AllStats_final.eIncomeDecile8.Consumption.Mean-AllStats_initial.eIncomeDecile8.Consumption.Mean)/AllStats_initial.eIncomeDecile8.Consumption.Mean;
    Figure4A(8,1,2)=(AllStatsPath.eIncomeDecile8.Consumption.Mean(1)-AllStats_initial.eIncomeDecile8.Consumption.Mean)/AllStats_initial.eIncomeDecile8.Consumption.Mean;
    Figure4A(8,2,1)=(AllStats_final.wIncomeDecile8.Consumption.Mean-AllStats_initial.wIncomeDecile8.Consumption.Mean)/AllStats_initial.wIncomeDecile8.Consumption.Mean;
    Figure4A(8,2,2)=(AllStatsPath.wIncomeDecile8.Consumption.Mean(1)-AllStats_initial.wIncomeDecile8.Consumption.Mean)/AllStats_initial.wIncomeDecile8.Consumption.Mean;
    % decile 9
    Figure4A(9,1,1)=(AllStats_final.eIncomeDecile9.Consumption.Mean-AllStats_initial.eIncomeDecile9.Consumption.Mean)/AllStats_initial.eIncomeDecile9.Consumption.Mean;
    Figure4A(9,1,2)=(AllStatsPath.eIncomeDecile9.Consumption.Mean(1)-AllStats_initial.eIncomeDecile9.Consumption.Mean)/AllStats_initial.eIncomeDecile9.Consumption.Mean;
    Figure4A(9,2,1)=(AllStats_final.wIncomeDecile9.Consumption.Mean-AllStats_initial.wIncomeDecile9.Consumption.Mean)/AllStats_initial.wIncomeDecile9.Consumption.Mean;
    Figure4A(9,2,2)=(AllStatsPath.wIncomeDecile9.Consumption.Mean(1)-AllStats_initial.wIncomeDecile9.Consumption.Mean)/AllStats_initial.wIncomeDecile9.Consumption.Mean;
    % percentile 90 to 97
    Figure4A(10,1,1)=(AllStats_final.eIncomep90to97.Consumption.Mean-AllStats_initial.eIncomep90to97.Consumption.Mean)/AllStats_initial.eIncomep90to97.Consumption.Mean;
    Figure4A(10,1,2)=(AllStatsPath.eIncomep90to97.Consumption.Mean(1)-AllStats_initial.eIncomep90to97.Consumption.Mean)/AllStats_initial.eIncomep90to97.Consumption.Mean;
    Figure4A(10,2,1)=(AllStats_final.wIncomep90to97.Consumption.Mean-AllStats_initial.wIncomep90to97.Consumption.Mean)/AllStats_initial.wIncomep90to97.Consumption.Mean;
    Figure4A(10,2,2)=(AllStatsPath.wIncomep90to97.Consumption.Mean(1)-AllStats_initial.wIncomep90to97.Consumption.Mean)/AllStats_initial.wIncomep90to97.Consumption.Mean;
    % percentile Top 3
    Figure4A(11,1,1)=(AllStats_final.eIncomepTop3.Consumption.Mean-AllStats_initial.eIncomepTop3.Consumption.Mean)/AllStats_initial.eIncomepTop3.Consumption.Mean;
    Figure4A(11,1,2)=(AllStatsPath.eIncomepTop3.Consumption.Mean(1)-AllStats_initial.eIncomepTop3.Consumption.Mean)/AllStats_initial.eIncomepTop3.Consumption.Mean;
    Figure4A(11,2,1)=(AllStats_final.wIncomepTop3.Consumption.Mean-AllStats_initial.wIncomepTop3.Consumption.Mean)/AllStats_initial.wIncomepTop3.Consumption.Mean;
    Figure4A(11,2,2)=(AllStatsPath.wIncomepTop3.Consumption.Mean(1)-AllStats_initial.wIncomepTop3.Consumption.Mean)/AllStats_initial.wIncomepTop3.Consumption.Mean;

    % Savings (a)
    Figure4B=zeros(11,2,2); % decile, entrepreneur/worker, transition/stationary eqm
    % decile 1
    Figure4B(1,1,1)=(AllStats_final.eIncomeDecile1.A.Mean-AllStats_initial.eIncomeDecile1.A.Mean)/AllStats_initial.eIncomeDecile1.A.Mean;
    Figure4B(1,1,2)=(AllStatsPath.eIncomeDecile1.A.Mean(1)-AllStats_initial.eIncomeDecile1.A.Mean)/AllStats_initial.eIncomeDecile1.A.Mean;
    Figure4B(1,2,1)=(AllStats_final.wIncomeDecile1.A.Mean-AllStats_initial.wIncomeDecile1.A.Mean)/AllStats_initial.wIncomeDecile1.A.Mean;
    Figure4B(1,2,2)=(AllStatsPath.wIncomeDecile1.A.Mean(1)-AllStats_initial.wIncomeDecile1.A.Mean)/AllStats_initial.wIncomeDecile1.A.Mean;
    % decile 2
    Figure4B(2,1,1)=(AllStats_final.eIncomeDecile2.A.Mean-AllStats_initial.eIncomeDecile2.A.Mean)/AllStats_initial.eIncomeDecile2.A.Mean;
    Figure4B(2,1,2)=(AllStatsPath.eIncomeDecile2.A.Mean(1)-AllStats_initial.eIncomeDecile2.A.Mean)/AllStats_initial.eIncomeDecile2.A.Mean;
    Figure4B(2,2,1)=(AllStats_final.wIncomeDecile2.A.Mean-AllStats_initial.wIncomeDecile2.A.Mean)/AllStats_initial.wIncomeDecile2.A.Mean;
    Figure4B(2,2,2)=(AllStatsPath.wIncomeDecile2.A.Mean(1)-AllStats_initial.wIncomeDecile2.A.Mean)/AllStats_initial.wIncomeDecile2.A.Mean;
    % decile 3
    Figure4B(3,1,1)=(AllStats_final.eIncomeDecile3.A.Mean-AllStats_initial.eIncomeDecile3.A.Mean)/AllStats_initial.eIncomeDecile3.A.Mean;
    Figure4B(3,1,2)=(AllStatsPath.eIncomeDecile3.A.Mean(1)-AllStats_initial.eIncomeDecile3.A.Mean)/AllStats_initial.eIncomeDecile3.A.Mean;
    Figure4B(3,2,1)=(AllStats_final.wIncomeDecile3.A.Mean-AllStats_initial.wIncomeDecile3.A.Mean)/AllStats_initial.wIncomeDecile3.A.Mean;
    Figure4B(3,2,2)=(AllStatsPath.wIncomeDecile3.A.Mean(1)-AllStats_initial.wIncomeDecile3.A.Mean)/AllStats_initial.wIncomeDecile3.A.Mean;
    % decile 4
    Figure4B(4,1,1)=(AllStats_final.eIncomeDecile4.A.Mean-AllStats_initial.eIncomeDecile4.A.Mean)/AllStats_initial.eIncomeDecile4.A.Mean;
    Figure4B(4,1,2)=(AllStatsPath.eIncomeDecile4.A.Mean(1)-AllStats_initial.eIncomeDecile4.A.Mean)/AllStats_initial.eIncomeDecile4.A.Mean;
    Figure4B(4,2,1)=(AllStats_final.wIncomeDecile4.A.Mean-AllStats_initial.wIncomeDecile4.A.Mean)/AllStats_initial.wIncomeDecile4.A.Mean;
    Figure4B(4,2,2)=(AllStatsPath.wIncomeDecile4.A.Mean(1)-AllStats_initial.wIncomeDecile4.A.Mean)/AllStats_initial.wIncomeDecile4.A.Mean;
    % decile 5
    Figure4B(5,1,1)=(AllStats_final.eIncomeDecile5.A.Mean-AllStats_initial.eIncomeDecile5.A.Mean)/AllStats_initial.eIncomeDecile5.A.Mean;
    Figure4B(5,1,2)=(AllStatsPath.eIncomeDecile5.A.Mean(1)-AllStats_initial.eIncomeDecile5.A.Mean)/AllStats_initial.eIncomeDecile5.A.Mean;
    Figure4B(5,2,1)=(AllStats_final.wIncomeDecile5.A.Mean-AllStats_initial.wIncomeDecile5.A.Mean)/AllStats_initial.wIncomeDecile5.A.Mean;
    Figure4B(5,2,2)=(AllStatsPath.wIncomeDecile5.A.Mean(1)-AllStats_initial.wIncomeDecile5.A.Mean)/AllStats_initial.wIncomeDecile5.A.Mean;
    % decile 6
    Figure4B(6,1,1)=(AllStats_final.eIncomeDecile6.A.Mean-AllStats_initial.eIncomeDecile6.A.Mean)/AllStats_initial.eIncomeDecile6.A.Mean;
    Figure4B(6,1,2)=(AllStatsPath.eIncomeDecile6.A.Mean(1)-AllStats_initial.eIncomeDecile6.A.Mean)/AllStats_initial.eIncomeDecile6.A.Mean;
    Figure4B(6,2,1)=(AllStats_final.wIncomeDecile6.A.Mean-AllStats_initial.wIncomeDecile6.A.Mean)/AllStats_initial.wIncomeDecile6.A.Mean;
    Figure4B(6,2,2)=(AllStatsPath.wIncomeDecile6.A.Mean(1)-AllStats_initial.wIncomeDecile6.A.Mean)/AllStats_initial.wIncomeDecile6.A.Mean;
    % decile 7
    Figure4B(7,1,1)=(AllStats_final.eIncomeDecile7.A.Mean-AllStats_initial.eIncomeDecile7.A.Mean)/AllStats_initial.eIncomeDecile7.A.Mean;
    Figure4B(7,1,2)=(AllStatsPath.eIncomeDecile7.A.Mean(1)-AllStats_initial.eIncomeDecile7.A.Mean)/AllStats_initial.eIncomeDecile7.A.Mean;
    Figure4B(7,2,1)=(AllStats_final.wIncomeDecile7.A.Mean-AllStats_initial.wIncomeDecile7.A.Mean)/AllStats_initial.wIncomeDecile7.A.Mean;
    Figure4B(7,2,2)=(AllStatsPath.wIncomeDecile7.A.Mean(1)-AllStats_initial.wIncomeDecile7.A.Mean)/AllStats_initial.wIncomeDecile7.A.Mean;
    % decile 8
    Figure4B(8,1,1)=(AllStats_final.eIncomeDecile8.A.Mean-AllStats_initial.eIncomeDecile8.A.Mean)/AllStats_initial.eIncomeDecile8.A.Mean;
    Figure4B(8,1,2)=(AllStatsPath.eIncomeDecile8.A.Mean(1)-AllStats_initial.eIncomeDecile8.A.Mean)/AllStats_initial.eIncomeDecile8.A.Mean;
    Figure4B(8,2,1)=(AllStats_final.wIncomeDecile8.A.Mean-AllStats_initial.wIncomeDecile8.A.Mean)/AllStats_initial.wIncomeDecile8.A.Mean;
    Figure4B(8,2,2)=(AllStatsPath.wIncomeDecile8.A.Mean(1)-AllStats_initial.wIncomeDecile8.A.Mean)/AllStats_initial.wIncomeDecile8.A.Mean;
    % decile 9
    Figure4B(9,1,1)=(AllStats_final.eIncomeDecile9.A.Mean-AllStats_initial.eIncomeDecile9.A.Mean)/AllStats_initial.eIncomeDecile9.A.Mean;
    Figure4B(9,1,2)=(AllStatsPath.eIncomeDecile9.A.Mean(1)-AllStats_initial.eIncomeDecile9.A.Mean)/AllStats_initial.eIncomeDecile9.A.Mean;
    Figure4B(9,2,1)=(AllStats_final.wIncomeDecile9.A.Mean-AllStats_initial.wIncomeDecile9.A.Mean)/AllStats_initial.wIncomeDecile9.A.Mean;
    Figure4B(9,2,2)=(AllStatsPath.wIncomeDecile9.A.Mean(1)-AllStats_initial.wIncomeDecile9.A.Mean)/AllStats_initial.wIncomeDecile9.A.Mean;
    % percentile 90 to 97
    Figure4B(10,1,1)=(AllStats_final.eIncomep90to97.A.Mean-AllStats_initial.eIncomep90to97.A.Mean)/AllStats_initial.eIncomep90to97.A.Mean;
    Figure4B(10,1,2)=(AllStatsPath.eIncomep90to97.A.Mean(1)-AllStats_initial.eIncomep90to97.A.Mean)/AllStats_initial.eIncomep90to97.A.Mean;
    Figure4B(10,2,1)=(AllStats_final.wIncomep90to97.A.Mean-AllStats_initial.wIncomep90to97.A.Mean)/AllStats_initial.wIncomep90to97.A.Mean;
    Figure4B(10,2,2)=(AllStatsPath.wIncomep90to97.A.Mean(1)-AllStats_initial.wIncomep90to97.A.Mean)/AllStats_initial.wIncomep90to97.A.Mean;
    % percentile Top 3
    Figure4B(11,1,1)=(AllStats_final.eIncomepTop3.A.Mean-AllStats_initial.eIncomepTop3.A.Mean)/AllStats_initial.eIncomepTop3.A.Mean;
    Figure4B(11,1,2)=(AllStatsPath.eIncomepTop3.A.Mean(1)-AllStats_initial.eIncomepTop3.A.Mean)/AllStats_initial.eIncomepTop3.A.Mean;
    Figure4B(11,2,1)=(AllStats_final.wIncomepTop3.A.Mean-AllStats_initial.wIncomepTop3.A.Mean)/AllStats_initial.wIncomepTop3.A.Mean;
    Figure4B(11,2,2)=(AllStatsPath.wIncomepTop3.A.Mean(1)-AllStats_initial.wIncomepTop3.A.Mean)/AllStats_initial.wIncomepTop3.A.Mean;

    
    % Hours worked (l) [B2021 only plots workers, as those for entrepreneurs will have 0% change, but I calculate them anyway]
    Figure4C=zeros(11,2,2); % decile, entrepreneur/worker, transition/stationary eqm
    % decile 1
    Figure4C(1,1,1)=(AllStats_final.eIncomeDecile1.L.Mean-AllStats_initial.eIncomeDecile1.L.Mean)/AllStats_initial.eIncomeDecile1.L.Mean;
    Figure4C(1,1,2)=(AllStatsPath.eIncomeDecile1.L.Mean(1)-AllStats_initial.eIncomeDecile1.L.Mean)/AllStats_initial.eIncomeDecile1.L.Mean;
    Figure4C(1,2,1)=(AllStats_final.wIncomeDecile1.L.Mean-AllStats_initial.wIncomeDecile1.L.Mean)/AllStats_initial.wIncomeDecile1.L.Mean;
    Figure4C(1,2,2)=(AllStatsPath.wIncomeDecile1.L.Mean(1)-AllStats_initial.wIncomeDecile1.L.Mean)/AllStats_initial.wIncomeDecile1.L.Mean;
    % decile 2
    Figure4C(2,1,1)=(AllStats_final.eIncomeDecile2.L.Mean-AllStats_initial.eIncomeDecile2.L.Mean)/AllStats_initial.eIncomeDecile2.L.Mean;
    Figure4C(2,1,2)=(AllStatsPath.eIncomeDecile2.L.Mean(1)-AllStats_initial.eIncomeDecile2.L.Mean)/AllStats_initial.eIncomeDecile2.L.Mean;
    Figure4C(2,2,1)=(AllStats_final.wIncomeDecile2.L.Mean-AllStats_initial.wIncomeDecile2.L.Mean)/AllStats_initial.wIncomeDecile2.L.Mean;
    Figure4C(2,2,2)=(AllStatsPath.wIncomeDecile2.L.Mean(1)-AllStats_initial.wIncomeDecile2.L.Mean)/AllStats_initial.wIncomeDecile2.L.Mean;
    % decile 3
    Figure4C(3,1,1)=(AllStats_final.eIncomeDecile3.L.Mean-AllStats_initial.eIncomeDecile3.L.Mean)/AllStats_initial.eIncomeDecile3.L.Mean;
    Figure4C(3,1,2)=(AllStatsPath.eIncomeDecile3.L.Mean(1)-AllStats_initial.eIncomeDecile3.L.Mean)/AllStats_initial.eIncomeDecile3.L.Mean;
    Figure4C(3,2,1)=(AllStats_final.wIncomeDecile3.L.Mean-AllStats_initial.wIncomeDecile3.L.Mean)/AllStats_initial.wIncomeDecile3.L.Mean;
    Figure4C(3,2,2)=(AllStatsPath.wIncomeDecile3.L.Mean(1)-AllStats_initial.wIncomeDecile3.L.Mean)/AllStats_initial.wIncomeDecile3.L.Mean;
    % decile 4
    Figure4C(4,1,1)=(AllStats_final.eIncomeDecile4.L.Mean-AllStats_initial.eIncomeDecile4.L.Mean)/AllStats_initial.eIncomeDecile4.L.Mean;
    Figure4C(4,1,2)=(AllStatsPath.eIncomeDecile4.L.Mean(1)-AllStats_initial.eIncomeDecile4.L.Mean)/AllStats_initial.eIncomeDecile4.L.Mean;
    Figure4C(4,2,1)=(AllStats_final.wIncomeDecile4.L.Mean-AllStats_initial.wIncomeDecile4.L.Mean)/AllStats_initial.wIncomeDecile4.L.Mean;
    Figure4C(4,2,2)=(AllStatsPath.wIncomeDecile4.L.Mean(1)-AllStats_initial.wIncomeDecile4.L.Mean)/AllStats_initial.wIncomeDecile4.L.Mean;
    % decile 5
    Figure4C(5,1,1)=(AllStats_final.eIncomeDecile5.L.Mean-AllStats_initial.eIncomeDecile5.L.Mean)/AllStats_initial.eIncomeDecile5.L.Mean;
    Figure4C(5,1,2)=(AllStatsPath.eIncomeDecile5.L.Mean(1)-AllStats_initial.eIncomeDecile5.L.Mean)/AllStats_initial.eIncomeDecile5.L.Mean;
    Figure4C(5,2,1)=(AllStats_final.wIncomeDecile5.L.Mean-AllStats_initial.wIncomeDecile5.L.Mean)/AllStats_initial.wIncomeDecile5.L.Mean;
    Figure4C(5,2,2)=(AllStatsPath.wIncomeDecile5.L.Mean(1)-AllStats_initial.wIncomeDecile5.L.Mean)/AllStats_initial.wIncomeDecile5.L.Mean;
    % decile 6
    Figure4C(6,1,1)=(AllStats_final.eIncomeDecile6.L.Mean-AllStats_initial.eIncomeDecile6.L.Mean)/AllStats_initial.eIncomeDecile6.L.Mean;
    Figure4C(6,1,2)=(AllStatsPath.eIncomeDecile6.L.Mean(1)-AllStats_initial.eIncomeDecile6.L.Mean)/AllStats_initial.eIncomeDecile6.L.Mean;
    Figure4C(6,2,1)=(AllStats_final.wIncomeDecile6.L.Mean-AllStats_initial.wIncomeDecile6.L.Mean)/AllStats_initial.wIncomeDecile6.L.Mean;
    Figure4C(6,2,2)=(AllStatsPath.wIncomeDecile6.L.Mean(1)-AllStats_initial.wIncomeDecile6.L.Mean)/AllStats_initial.wIncomeDecile6.L.Mean;
    % decile 7
    Figure4C(7,1,1)=(AllStats_final.eIncomeDecile7.L.Mean-AllStats_initial.eIncomeDecile7.L.Mean)/AllStats_initial.eIncomeDecile7.L.Mean;
    Figure4C(7,1,2)=(AllStatsPath.eIncomeDecile7.L.Mean(1)-AllStats_initial.eIncomeDecile7.L.Mean)/AllStats_initial.eIncomeDecile7.L.Mean;
    Figure4C(7,2,1)=(AllStats_final.wIncomeDecile7.L.Mean-AllStats_initial.wIncomeDecile7.L.Mean)/AllStats_initial.wIncomeDecile7.L.Mean;
    Figure4C(7,2,2)=(AllStatsPath.wIncomeDecile7.L.Mean(1)-AllStats_initial.wIncomeDecile7.L.Mean)/AllStats_initial.wIncomeDecile7.L.Mean;
    % decile 8
    Figure4C(8,1,1)=(AllStats_final.eIncomeDecile8.L.Mean-AllStats_initial.eIncomeDecile8.L.Mean)/AllStats_initial.eIncomeDecile8.L.Mean;
    Figure4C(8,1,2)=(AllStatsPath.eIncomeDecile8.L.Mean(1)-AllStats_initial.eIncomeDecile8.L.Mean)/AllStats_initial.eIncomeDecile8.L.Mean;
    Figure4C(8,2,1)=(AllStats_final.wIncomeDecile8.L.Mean-AllStats_initial.wIncomeDecile8.L.Mean)/AllStats_initial.wIncomeDecile8.L.Mean;
    Figure4C(8,2,2)=(AllStatsPath.wIncomeDecile8.L.Mean(1)-AllStats_initial.wIncomeDecile8.L.Mean)/AllStats_initial.wIncomeDecile8.L.Mean;
    % decile 9
    Figure4C(9,1,1)=(AllStats_final.eIncomeDecile9.L.Mean-AllStats_initial.eIncomeDecile9.L.Mean)/AllStats_initial.eIncomeDecile9.L.Mean;
    Figure4C(9,1,2)=(AllStatsPath.eIncomeDecile9.L.Mean(1)-AllStats_initial.eIncomeDecile9.L.Mean)/AllStats_initial.eIncomeDecile9.L.Mean;
    Figure4C(9,2,1)=(AllStats_final.wIncomeDecile9.L.Mean-AllStats_initial.wIncomeDecile9.L.Mean)/AllStats_initial.wIncomeDecile9.L.Mean;
    Figure4C(9,2,2)=(AllStatsPath.wIncomeDecile9.L.Mean(1)-AllStats_initial.wIncomeDecile9.L.Mean)/AllStats_initial.wIncomeDecile9.L.Mean;
    % percentile 90 to 97
    Figure4C(10,1,1)=(AllStats_final.eIncomep90to97.L.Mean-AllStats_initial.eIncomep90to97.L.Mean)/AllStats_initial.eIncomep90to97.L.Mean;
    Figure4C(10,1,2)=(AllStatsPath.eIncomep90to97.L.Mean(1)-AllStats_initial.eIncomep90to97.L.Mean)/AllStats_initial.eIncomep90to97.L.Mean;
    Figure4C(10,2,1)=(AllStats_final.wIncomep90to97.L.Mean-AllStats_initial.wIncomep90to97.L.Mean)/AllStats_initial.wIncomep90to97.L.Mean;
    Figure4C(10,2,2)=(AllStatsPath.wIncomep90to97.L.Mean(1)-AllStats_initial.wIncomep90to97.L.Mean)/AllStats_initial.wIncomep90to97.L.Mean;
    % percentile Top 3
    Figure4C(11,1,1)=(AllStats_final.eIncomepTop3.L.Mean-AllStats_initial.eIncomepTop3.L.Mean)/AllStats_initial.eIncomepTop3.L.Mean;
    Figure4C(11,1,2)=(AllStatsPath.eIncomepTop3.L.Mean(1)-AllStats_initial.eIncomepTop3.L.Mean)/AllStats_initial.eIncomepTop3.L.Mean;
    Figure4C(11,2,1)=(AllStats_final.wIncomepTop3.L.Mean-AllStats_initial.wIncomepTop3.L.Mean)/AllStats_initial.wIncomepTop3.L.Mean;
    Figure4C(11,2,2)=(AllStatsPath.wIncomepTop3.L.Mean(1)-AllStats_initial.wIncomepTop3.L.Mean)/AllStats_initial.wIncomepTop3.L.Mean;

    % Figure 4
    fig4A=figure(4);
    subplot(1,2,1); bar(xaxissetup, Figure4A(:,:,2)')
    ylabel('Percent change in average choices')
    title('First Period of Transition')
    subplot(1,2,2); bar(xaxissetup, Figure4A(:,:,1)')
    ylabel('Percent change in average choices')
    title('Stationary Eqm comparison')
    sgtitle('Consumption (c)') 
    saveas(fig4A,'./SavedOutput/Graphs/Bruggemann2021_Fig4A.png')

    fig4B=figure(15);
    subplot(1,2,1); bar(xaxissetup, Figure4B(:,:,2)')
    ylabel('Percent change in average choices')
    title('First Period of Transition')
    subplot(1,2,2); bar(xaxissetup, Figure4B(:,:,1)')
    ylabel('Percent change in average choices')
    title('Stationary Eqm comparison')
    sgtitle('Savings (a)') 
    saveas(fig4B,'./SavedOutput/Graphs/Bruggemann2021_Fig4B.png')

    fig4C=figure(16);
    subplot(1,2,1); bar(xaxissetup, Figure4C(:,:,2)')
    ylabel('Percent change in average choices')
    title('First Period of Transition')
    subplot(1,2,2); bar(xaxissetup, Figure4C(:,:,1)')
    ylabel('Percent change in average choices')
    title('Stationary Eqm comparison')
    sgtitle('Hours Worked (l)')
    legend({'Workers','Entrepreneurs'},'Location','southoutside')
    saveas(fig4C,'./SavedOutput/Graphs/Bruggemann2021_Fig4C.png')
    
    %% Figure 5: Changes to Entrepreneurs' Average Choices by Income Level:
    % Investment (k), Hiring (n)
    
    % Investment (k)
    Figure5A=zeros(11,1,2); % decile, entrepreneur, transition/stationary eqm
    % decile 1
    Figure5A(1,1,1)=(AllStats_final.eIncomeDecile1.K_noncorp.Mean-AllStats_initial.eIncomeDecile1.K_noncorp.Mean)/AllStats_initial.eIncomeDecile1.K_noncorp.Mean;
    Figure5A(1,1,2)=(AllStatsPath.eIncomeDecile1.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile1.K_noncorp.Mean)/AllStats_initial.eIncomeDecile1.K_noncorp.Mean;
    % decile 2
    Figure5A(2,1,1)=(AllStats_final.eIncomeDecile2.K_noncorp.Mean-AllStats_initial.eIncomeDecile2.K_noncorp.Mean)/AllStats_initial.eIncomeDecile2.K_noncorp.Mean;
    Figure5A(2,1,2)=(AllStatsPath.eIncomeDecile2.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile2.K_noncorp.Mean)/AllStats_initial.eIncomeDecile2.K_noncorp.Mean;
    % decile 3
    Figure5A(3,1,1)=(AllStats_final.eIncomeDecile3.K_noncorp.Mean-AllStats_initial.eIncomeDecile3.K_noncorp.Mean)/AllStats_initial.eIncomeDecile3.K_noncorp.Mean;
    Figure5A(3,1,2)=(AllStatsPath.eIncomeDecile3.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile3.K_noncorp.Mean)/AllStats_initial.eIncomeDecile3.K_noncorp.Mean;
    % decile 4
    Figure5A(4,1,1)=(AllStats_final.eIncomeDecile4.K_noncorp.Mean-AllStats_initial.eIncomeDecile4.K_noncorp.Mean)/AllStats_initial.eIncomeDecile4.K_noncorp.Mean;
    Figure5A(4,1,2)=(AllStatsPath.eIncomeDecile4.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile4.K_noncorp.Mean)/AllStats_initial.eIncomeDecile4.K_noncorp.Mean;
    % decile 5
    Figure5A(5,1,1)=(AllStats_final.eIncomeDecile5.K_noncorp.Mean-AllStats_initial.eIncomeDecile5.K_noncorp.Mean)/AllStats_initial.eIncomeDecile5.K_noncorp.Mean;
    Figure5A(5,1,2)=(AllStatsPath.eIncomeDecile5.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile5.K_noncorp.Mean)/AllStats_initial.eIncomeDecile5.K_noncorp.Mean;
    % decile 6
    Figure5A(6,1,1)=(AllStats_final.eIncomeDecile6.K_noncorp.Mean-AllStats_initial.eIncomeDecile6.K_noncorp.Mean)/AllStats_initial.eIncomeDecile6.K_noncorp.Mean;
    Figure5A(6,1,2)=(AllStatsPath.eIncomeDecile6.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile6.K_noncorp.Mean)/AllStats_initial.eIncomeDecile6.K_noncorp.Mean;
    % decile 7
    Figure5A(7,1,1)=(AllStats_final.eIncomeDecile7.K_noncorp.Mean-AllStats_initial.eIncomeDecile7.K_noncorp.Mean)/AllStats_initial.eIncomeDecile7.K_noncorp.Mean;
    Figure5A(7,1,2)=(AllStatsPath.eIncomeDecile7.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile7.K_noncorp.Mean)/AllStats_initial.eIncomeDecile7.K_noncorp.Mean;
    % decile 8
    Figure5A(8,1,1)=(AllStats_final.eIncomeDecile8.K_noncorp.Mean-AllStats_initial.eIncomeDecile8.K_noncorp.Mean)/AllStats_initial.eIncomeDecile8.K_noncorp.Mean;
    Figure5A(8,1,2)=(AllStatsPath.eIncomeDecile8.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile8.K_noncorp.Mean)/AllStats_initial.eIncomeDecile8.K_noncorp.Mean;
    % decile 9
    Figure5A(9,1,1)=(AllStats_final.eIncomeDecile9.K_noncorp.Mean-AllStats_initial.eIncomeDecile9.K_noncorp.Mean)/AllStats_initial.eIncomeDecile9.K_noncorp.Mean;
    Figure5A(9,1,2)=(AllStatsPath.eIncomeDecile9.K_noncorp.Mean(1)-AllStats_initial.eIncomeDecile9.K_noncorp.Mean)/AllStats_initial.eIncomeDecile9.K_noncorp.Mean;
    % percentile 90 to 97
    Figure5A(10,1,1)=(AllStats_final.eIncomep90to97.K_noncorp.Mean-AllStats_initial.eIncomep90to97.K_noncorp.Mean)/AllStats_initial.eIncomep90to97.K_noncorp.Mean;
    Figure5A(10,1,2)=(AllStatsPath.eIncomep90to97.K_noncorp.Mean(1)-AllStats_initial.eIncomep90to97.K_noncorp.Mean)/AllStats_initial.eIncomep90to97.K_noncorp.Mean;
    % percentile Top 3
    Figure5A(11,1,1)=(AllStats_final.eIncomepTop3.K_noncorp.Mean-AllStats_initial.eIncomepTop3.K_noncorp.Mean)/AllStats_initial.eIncomepTop3.K_noncorp.Mean;
    Figure5A(11,1,2)=(AllStatsPath.eIncomepTop3.K_noncorp.Mean(1)-AllStats_initial.eIncomepTop3.K_noncorp.Mean)/AllStats_initial.eIncomepTop3.K_noncorp.Mean;

    % Hiring (n)
    Figure5B=zeros(11,2,2); % decile, entrepreneur/worker, transition/stationary eqm
    % decile 1
    Figure5B(1,1,1)=(AllStats_final.eIncomeDecile1.N_noncorp.Mean-AllStats_initial.eIncomeDecile1.N_noncorp.Mean)/AllStats_initial.eIncomeDecile1.N_noncorp.Mean;
    Figure5B(1,1,2)=(AllStatsPath.eIncomeDecile1.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile1.N_noncorp.Mean)/AllStats_initial.eIncomeDecile1.N_noncorp.Mean;
    % decile 2
    Figure5B(2,1,1)=(AllStats_final.eIncomeDecile2.N_noncorp.Mean-AllStats_initial.eIncomeDecile2.N_noncorp.Mean)/AllStats_initial.eIncomeDecile2.N_noncorp.Mean;
    Figure5B(2,1,2)=(AllStatsPath.eIncomeDecile2.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile2.N_noncorp.Mean)/AllStats_initial.eIncomeDecile2.N_noncorp.Mean;
    % decile 3
    Figure5B(3,1,1)=(AllStats_final.eIncomeDecile3.N_noncorp.Mean-AllStats_initial.eIncomeDecile3.N_noncorp.Mean)/AllStats_initial.eIncomeDecile3.N_noncorp.Mean;
    Figure5B(3,1,2)=(AllStatsPath.eIncomeDecile3.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile3.N_noncorp.Mean)/AllStats_initial.eIncomeDecile3.N_noncorp.Mean;
    % decile 4
    Figure5B(4,1,1)=(AllStats_final.eIncomeDecile4.N_noncorp.Mean-AllStats_initial.eIncomeDecile4.N_noncorp.Mean)/AllStats_initial.eIncomeDecile4.N_noncorp.Mean;
    Figure5B(4,1,2)=(AllStatsPath.eIncomeDecile4.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile4.N_noncorp.Mean)/AllStats_initial.eIncomeDecile4.N_noncorp.Mean;
    % decile 5
    Figure5B(5,1,1)=(AllStats_final.eIncomeDecile5.N_noncorp.Mean-AllStats_initial.eIncomeDecile5.N_noncorp.Mean)/AllStats_initial.eIncomeDecile5.N_noncorp.Mean;
    Figure5B(5,1,2)=(AllStatsPath.eIncomeDecile5.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile5.N_noncorp.Mean)/AllStats_initial.eIncomeDecile5.N_noncorp.Mean;
    % decile 6
    Figure5B(6,1,1)=(AllStats_final.eIncomeDecile6.N_noncorp.Mean-AllStats_initial.eIncomeDecile6.N_noncorp.Mean)/AllStats_initial.eIncomeDecile6.N_noncorp.Mean;
    Figure5B(6,1,2)=(AllStatsPath.eIncomeDecile6.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile6.N_noncorp.Mean)/AllStats_initial.eIncomeDecile6.N_noncorp.Mean;
    % decile 7
    Figure5B(7,1,1)=(AllStats_final.eIncomeDecile7.N_noncorp.Mean-AllStats_initial.eIncomeDecile7.N_noncorp.Mean)/AllStats_initial.eIncomeDecile7.N_noncorp.Mean;
    Figure5B(7,1,2)=(AllStatsPath.eIncomeDecile7.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile7.N_noncorp.Mean)/AllStats_initial.eIncomeDecile7.N_noncorp.Mean;
    % decile 8
    Figure5B(8,1,1)=(AllStats_final.eIncomeDecile8.N_noncorp.Mean-AllStats_initial.eIncomeDecile8.N_noncorp.Mean)/AllStats_initial.eIncomeDecile8.N_noncorp.Mean;
    Figure5B(8,1,2)=(AllStatsPath.eIncomeDecile8.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile8.N_noncorp.Mean)/AllStats_initial.eIncomeDecile8.N_noncorp.Mean;
    % decile 9
    Figure5B(9,1,1)=(AllStats_final.eIncomeDecile9.N_noncorp.Mean-AllStats_initial.eIncomeDecile9.N_noncorp.Mean)/AllStats_initial.eIncomeDecile9.N_noncorp.Mean;
    Figure5B(9,1,2)=(AllStatsPath.eIncomeDecile9.N_noncorp.Mean(1)-AllStats_initial.eIncomeDecile9.N_noncorp.Mean)/AllStats_initial.eIncomeDecile9.N_noncorp.Mean;
    % percentile 90 to 97
    Figure5B(10,1,1)=(AllStats_final.eIncomep90to97.N_noncorp.Mean-AllStats_initial.eIncomep90to97.N_noncorp.Mean)/AllStats_initial.eIncomep90to97.N_noncorp.Mean;
    Figure5B(10,1,2)=(AllStatsPath.eIncomep90to97.N_noncorp.Mean(1)-AllStats_initial.eIncomep90to97.N_noncorp.Mean)/AllStats_initial.eIncomep90to97.N_noncorp.Mean;
    % percentile Top 3
    Figure5B(11,1,1)=(AllStats_final.eIncomepTop3.N_noncorp.Mean-AllStats_initial.eIncomepTop3.N_noncorp.Mean)/AllStats_initial.eIncomepTop3.N_noncorp.Mean;
    Figure5B(11,1,2)=(AllStatsPath.eIncomepTop3.N_noncorp.Mean(1)-AllStats_initial.eIncomepTop3.N_noncorp.Mean)/AllStats_initial.eIncomepTop3.N_noncorp.Mean;

    % Figure 5
    fig5=figure(5);
    subplot(1,2,1); bar(xaxissetup, squeeze(Figure5A(:,1,:))')
    ylabel('Percent change in average choices')
    title('Investment (k)')
    subplot(1,2,2); bar(xaxissetup, squeeze(Figure5A(:,1,:))')
    ylabel('Percent change in average choices')
    title('Hiring (n)')
    sgtitle('Changes to Entrepreneurs Average Choices by Income Level')
    legend({'First period of Transition','Stationary Eqm Comparison'},'Location','southoutside')
    saveas(fig5,'./SavedOutput/Graphs/Bruggemann2021_Fig5.png')
    

    %% Solve the partial eqm for the top marginal tax rate reform
    % We just set all the parameters to what they were in the initial
    % stationary general eqm, and then use the tax rate from the reform.
    Params.r=p_eqm_initial.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_initial.w;
    Params.tau_s=p_eqm_initial.tau_s;
    % Params.tau_i_r6 is already set to the top marginal tax rate from reform

    [V_final_PE,Policy_final_PE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final_PE=StationaryDist_InfHorz(Policy_final_PE,n_d,n_a,n_z,pi_z, simoptions);
    AllStats_final_PE=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_final_PE,Policy_final_PE, FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    % Calculate output
    Output_corp_final_PE=Params.Z*((AllStats_final_PE.A.Mean-AllStats_final_PE.K_noncorp.Mean)^Params.alpha)*((AllStats_final_PE.L.Mean-AllStats_final_PE.N_lbar.Mean-AllStats_final_PE.N_noncorp.Mean)^(1-Params.alpha));
    Y_final_PE=Output_corp_final_PE+AllStats_final_PE.Y_noncorp.Mean;

    %% Following are the levels, then in Table 7 we report them as percentages
    Table7_GE=[Y_final,AllStats_final.A.Mean,AllStats_final.L.Mean,AllStats_final.IncomeTaxRevenue.Mean+AllStats_final.ConsumptionTaxRevenue.Mean,p_eqm_final.r,p_eqm_final.w];
    Table7_PE=[Y_final_PE,AllStats_final_PE.A.Mean,AllStats_final_PE.L.Mean,AllStats_final_PE.IncomeTaxRevenue.Mean+AllStats_final_PE.ConsumptionTaxRevenue.Mean,p_eqm_initial.r,p_eqm_initial.w];
    Table7_initial=[Y_initial,AllStats_initial.A.Mean,AllStats_initial.L.Mean,AllStats_initial.IncomeTaxRevenue.Mean+AllStats_initial.ConsumptionTaxRevenue.Mean,p_eqm_initial.r,p_eqm_initial.w];

    % Table 7
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table7.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \n');
    fprintf(FID, ' \\multicolumn{7}{c}{Aggregate Changes Between Stationary Equilibria after Increasing the Effective TMTR to %4.1f percent} \\\\ \\hline  \\hline  \n', 100*Params.tau_i_r6);
    fprintf(FID, '  & Y & K & N & T & r & w \\\\ \\hline \n');
    fprintf(FID, '  $\\tau^{max}$= %4.1f \\%%, GE (percent) & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table7_GE-Table7_initial)./Table7_initial);
    fprintf(FID, '  $\\tau^{max}$= %4.1f \\%%, PE (percent) & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table7_PE-Table7_initial)./Table7_initial);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: PE stands for partial equilibrium (prices at benchmark level), GE for general equilibrium (prices adjust). \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    % Following are the levels, then in Table 8 we report them as percentages
    Table8_GE_corp=[Y_final-AllStats_final.Y_noncorp.Mean,AllStats_final.A.Mean-AllStats_final.K_noncorp.Mean,AllStats_final.L.Mean-AllStats_final.N_lbar.Mean-AllStats_final.N_noncorp.Mean];
    Table8_PE_corp=[Y_final_PE-AllStats_final_PE.Y_noncorp.Mean,AllStats_final_PE.A.Mean-AllStats_final_PE.K_noncorp.Mean,AllStats_final_PE.L.Mean-AllStats_final_PE.N_lbar.Mean-AllStats_final_PE.N_noncorp.Mean];
    Table8_initial_corp=[Y_initial-AllStats_initial.Y_noncorp.Mean,AllStats_initial.A.Mean-AllStats_initial.K_noncorp.Mean,AllStats_initial.L.Mean-AllStats_initial.N_lbar.Mean-AllStats_initial.N_noncorp.Mean];

    Table8_GE_noncorp=[AllStats_final.Y_noncorp.Mean,AllStats_final.K_noncorp.Mean,AllStats_final.N_noncorp.Mean,AllStats_final.Entrepreneur.Mean];
    Table8_PE_noncorp=[AllStats_final_PE.Y_noncorp.Mean,AllStats_final_PE.K_noncorp.Mean,AllStats_final_PE.N_noncorp.Mean,AllStats_final_PE.Entrepreneur.Mean];
    Table8_initial_noncorp=[AllStats_initial.Y_noncorp.Mean,AllStats_initial.K_noncorp.Mean,AllStats_initial.N_noncorp.Mean,AllStats_initial.Entrepreneur.Mean];

    % Table 8
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table8.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n');
    fprintf(FID, ' \\multicolumn{5}{c}{Aggregate Changes Between Stationary Equilibria after Increasing the Effective TMTR to %4.1f percent: By Sector} \\\\ \\hline  \\hline  \n', 100*Params.tau_i_r6);
    fprintf(FID, '  & Y & K & N & \\#E \\\\ \\hline \n');
    fprintf(FID, 'Corporate Sector & & & & \\\\ \n');
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, GE (percent) & %8.1f & %8.1f & %8.1f &  \\\\ \n', 100*Params.tau_i_r6, 100*(Table8_GE_corp-Table8_initial_corp)./Table8_initial_corp);
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, PE (percent) & %8.1f & %8.1f & %8.1f &  \\\\ \n', 100*Params.tau_i_r6, 100*(Table8_PE_corp-Table8_initial_corp)./Table8_initial_corp);
    fprintf(FID, 'Entrepreneurial Sector & & & & \\\\ \n');
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, GE (percent) & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table8_GE_noncorp-Table8_initial_noncorp)./Table8_initial_noncorp);
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, PE (percent) & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table8_PE_noncorp-Table8_initial_noncorp)./Table8_initial_noncorp);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: PE stands for partial equilibrium (prices at benchmark level), GE for general equilibrium (prices adjust). \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    Table9_final=[(AllStats_final.A.Mean-AllStats_final.K_noncorp.Mean)/AllStats_final.A.Mean, AllStats_final.K_noncorp.Mean/AllStats_final.A.Mean, (AllStats_final.L.Mean-AllStats_final.N_lbar.Mean-AllStats_final.N_noncorp.Mean)/(AllStats_final.L.Mean-AllStats_final.N_lbar.Mean), AllStats_final.N_noncorp.Mean/(AllStats_final.L.Mean-AllStats_final.N_lbar.Mean)];
    Table9_initial=[(AllStats_initial.A.Mean-AllStats_initial.K_noncorp.Mean)/AllStats_initial.A.Mean, AllStats_initial.K_noncorp.Mean/AllStats_initial.A.Mean, (AllStats_initial.L.Mean-AllStats_initial.N_lbar.Mean-AllStats_initial.N_noncorp.Mean)/(AllStats_initial.L.Mean-AllStats_initial.N_lbar.Mean), AllStats_initial.N_noncorp.Mean/(AllStats_initial.L.Mean-AllStats_initial.N_lbar.Mean)];
    
    % Table 9
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table9.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n');
    fprintf(FID, ' \\multicolumn{5}{c}{Relative Factor Allocation in Corporate and Entrepreneurial Sector in the Benchmark Stationary} \\\\ \n');
    fprintf(FID, ' \\multicolumn{5}{c}{Equilibrium and After Increasing the Effective TMTR to %4.1f percent} \\\\ \\hline  \\hline  \n', 100*Params.tau_i_r6);
    fprintf(FID, '  & \\multicolumn{2}{c}{Capital} & \\multicolumn{2}{c}{Labor} \\\\ \\hline \n');
    fprintf(FID, '  & Corporate & Entrepreneurial & Corporate & Entrepreneurial \\\\ \n');
    fprintf(FID, '  Benchmark (percent) & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Table9_initial);
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, GE (percent) & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*Table9_final);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: GE means general equilibrium (prices adjust); only the general equilibrium is reported here. \n');
    fprintf(FID, 'The labor shares are shares of \\textit{market} labor, that is, of total efficiency units of labor hired by the two sectors; \n');
    fprintf(FID, 'entrepreneurs own hours, which are supplied inelastically at $\\bar{l}$ and are not hired on the labor market, are excluded from \n');
    fprintf(FID, 'the denominator. This follows the original Fortran codes, in which corporate labor is totlcorp=toteffl-hiredlabe and \n');
    fprintf(FID, 'entrepreneurial labor is hiredlabe. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    % Table 10
    FID = fopen('./SavedOutput/LatexInputs/Bruggemann2021_Table10.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n');
    fprintf(FID, ' \\multicolumn{4}{c}{Effects on Entrepreneurial Sector} \\\\ \\hline  \\hline  \n');
    fprintf(FID, '  & Entry & Exit & Fraction \\\\ \n');
    fprintf(FID, '  & (in percent) & (in percent) & (in percent) \\\\ \\hline \n');
    fprintf(FID, '  Benchmark & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Table10_initial);
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, stationary eqm & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table10_final-Table10_initial)./Table10_initial);
    fprintf(FID, '  $\\tau^{max}$=%4.1f \\%%, transition path & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Params.tau_i_r6, 100*(Table10_tpath-Table10_initial)./Table10_initial);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: The benchmark row reports levels, the other two rows report percent changes relative to the benchmark. \n');
    fprintf(FID, 'The stationary eqm row compares the two stationary equilibria; the transition path row uses the first period of the transition. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);


    % Part 5 needs nothing from Part 4, so this is purely so the tables and figures can be
    % rebuilt without re-solving. All of it is small: the whole workspace was 346MB only
    % because TPathResults and the gathered VPath/PolicyPath/AgentDistPath were still in it.
    save ./SavedOutput/B2021_doPart4.mat Table7_initial Table7_GE Table7_PE ...
        Table8_initial_corp Table8_GE_corp Table8_PE_corp ...
        Table8_initial_noncorp Table8_GE_noncorp Table8_PE_noncorp ...
        Table9_initial Table9_final Table10_final Table10_tpath ...
        Figure3A_CEV Figure3B_CEV Figure4A Figure4B Figure4C Figure5A Figure5B ...
        AllStats_initial AllStats_final AllStats_final_PE AllStatsPath AllStats_incomecutoffs ...
        Y_initial Y_final Y_final_PE Params_final
else
    load ./SavedOutput/B2021_doPart4.mat
    load ./SavedOutput/B2021_doPart.mat doPart headlessFigures
    Params.headlessFigures=headlessFigures; % the load above overwrote Params
end



%% Solve the alternative equilbria
% Create Figure 6
if doPart(5)==1

    if doPart5(1)==1
        % Solve No Entrepreneur Model 1: theta=0 is only value of theta
        output1=Bruggemann2021_NoEntrepreneur(tau_TMTR_vec_part5,Params_initial,n_d,n_a,n_eta,d_grid,a_grid,z_grid,pi_age,pi_eta,eta_statdist,T,ReturnFn,FnsToEvaluate,DiscountFactorParamNames,vfoptions,simoptions,heteroagentoptions,transpathoptions,vfoptionstpath);
        save ./SavedOutput/B2021_doPart5A.mat output1
    else
        load ./SavedOutput/B2021_doPart5A.mat
    end
    
    if doPart5(2)==1
        % Solve No Entrepreneur Model 2: eliminate eta_6, the largest value of eta
        output2=Bruggemann2021_NoHighestAbility(tau_TMTR_vec_part5,Params_initial,n_d,n_a,n_eta,n_theta,d_grid,a_grid,age_grid,eta_grid,theta_grid,pi_age,pi_eta,pi_theta,eta_statdist,theta_statdist,T,ReturnFn,FnsToEvaluate,DiscountFactorParamNames,vfoptions,simoptions,heteroagentoptions,transpathoptions,vfoptionstpath);
        save ./SavedOutput/B2021_doPart5B.mat output2
    else
        load ./SavedOutput/B2021_doPart5B.mat
    end
    
    % Figure 6: I do transition paths as well, B2021 only has stationary eqm comparisons
    fig6=figure(6);
    subplot(1,2,1); plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEVpath, Params.tau_i_adj*tau_TMTR_vec_part5,100*output1.CEV_TPath, Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEV, Params.tau_i_adj*tau_TMTR_vec_part5,100*output1.CEV_StationaryEqm)
    ylabel('CEV (in percent)')
    xlabel('Effective Top Marginal Tax Rate')
    xlim([0.20 0.80]) % same x-axis as B2021 Figure 6 No Entrepreneur
    title('No Entrepreneur')
    legend({'Baseline: First period of Transition','No Entrepreneurs: First period of Transition','Baseline: Stationary Eqm Comparison','No Entrepreneurs: Stationary Eqm Comparison'},'Location','southoutside')
    subplot(1,2,2); plot(Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEVpath, Params.tau_i_adj*tau_TMTR_vec_part5,100*output2.CEV_TPath, Params.tau_i_adj*tau_TMTR_vec,100*Figure1_CEV, Params.tau_i_adj*tau_TMTR_vec_part5,100*output2.CEV_StationaryEqm)
    ylabel('CEV (in percent)')
    xlabel('Effective Top Marginal Tax Rate')
    xlim([0.20 0.70]) % same x-axis as B2021 Figure 6 No Highest Labor Ability
    title('No Highest Labor Ability')
    legend({'Baseline: First period of Transition','No Highest Labor Ability: First period of Transition','Baseline: Stationary Eqm Comparison','No Highest Labor Ability: Stationary Eqm Comparison'},'Location','southoutside')
    sgtitle('Optimal Tax Rates under Alternative Specifications')
    saveas(fig6,'./SavedOutput/Graphs/Bruggemann2021_Fig6.png')
    
    
    % Figure 7: Revenue Maximizing Top Marginal Tax Rate (Laffer Curve)
    % This is not a figure of B2021. It is a version of Figure 6 showing the alternative
    % no entrepreneur models, but for Figure A1 which shows Laffer curves
    fig14=figure(14);
    subplot(1,2,1); plot(Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSumTPathPeriod1/Params.ybar, Params.tau_i_adj*tau_TMTR_vec_part5,100*output1.LafferCurve_LumpSumTPathPeriod1/Params.ybar, Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSum/Params.ybar, Params.tau_i_adj*tau_TMTR_vec_part5,100*output1.LafferCurve_LumpSum/Params.ybar)
    ylabel('Lump-sum transfer (in percent of average income)')
    xlabel('Effective Top Marginal Tax Rate')
    xlim([0.20 0.80])
    title('No Entrepreneur')
    legend({'Baseline: First period of Transition','No Entrepreneurs: First period of Transition','Baseline: Stationary Eqm Comparison','No Entrepreneurs: Stationary Eqm Comparison'},'Location','southoutside')
    subplot(1,2,2); plot(Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSumTPathPeriod1/Params.ybar, Params.tau_i_adj*tau_TMTR_vec_part5,100*output2.LafferCurve_LumpSumTPathPeriod1/Params.ybar, Params.tau_i_adj*tau_TMTR_vec,100*FigureA1_LumpSum/Params.ybar, Params.tau_i_adj*tau_TMTR_vec_part5,100*output2.LafferCurve_LumpSum/Params.ybar)
    ylabel('Lump-sum transfer (in percent of average income)')
    xlabel('Effective Top Marginal Tax Rate')
    xlim([0.20 0.70])
    title('No Highest Labor Ability')
    legend({'Baseline: First period of Transition','No Highest Labor Ability: First period of Transition','Baseline: Stationary Eqm Comparison','No Highest Labor Ability: Stationary Eqm Comparison'},'Location','southoutside')
    sgtitle('Revenue-Maximizing Top Marginal Tax Rate under Alternative Specifications')
    saveas(fig14,'./SavedOutput/Graphs/Bruggemann2021_Fig7.png')


else
    load ./SavedOutput/B2021_doPart5A.mat
    load ./SavedOutput/B2021_doPart5B.mat
    load ./SavedOutput/B2021_doPart.mat doPart headlessFigures
    Params.headlessFigures=headlessFigures; % the load above overwrote Params
end


%% End of run
diary off
