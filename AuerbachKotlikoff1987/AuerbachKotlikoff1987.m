% Replication of: Auerbach & Kotlikoff (1987) - Dynamic Fiscal Policy,
% Chapter 5.

% INCOMPLETE:Calculate all the equilibria, and have transition paths but
% they do not yet incorporate the welfare-compenstating transfers.

% A few lines needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
try % Server has 16 cores, but is shared with other users, so use max of 8.
    parpool(8)
    gpuDevice(1)
catch % Desktop has less than 8, so will give error, on desktop it is fine to use all available cores.
    parpool
end
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;

%%

n_l=51;
n_a=251;
N_j=55;

n_z=1; z_grid=0; pi_z=1; % Deterministic model, so no exogenous shocks z.

%% Parameters

Params.J=N_j;
Params.jj=1:1:Params.J;

% Discount rate
Params.delta=0.015;
Params.beta=1/(1+Params.delta);

% Constant population growth rate
Params.n=0.01;

% Utility function parameters
Params.rho=0.8;
Params.gamma=0.25;
Params.alpha=1.5;
% Params.ej=4.47+0.033*Params.jj-0.00067*Params.jj.^2; % Footnote on pg 52.
% Given Figure 5.2 I suspected may be a typo (that or noone is working more than 0.1 of time)
% Lines 47-50, 62 & 1393 of aktax3.for in codes https://ideas.repec.org/c/dge/qmrbcd/90.html confirm that it should actually be the following
Params.ej=exp(0.47+0.033*Params.jj-0.00067*Params.jj.^2)/(exp(0.47+0.033-0.00067));
% That is, they normalize ej at age 1 to take value of 1.

% Tax parameters for baseline
Params.tau_phi=0.15; % Proportional income tax. This is the only tax used as part of baseline.
Params.tau_pi=0; % Progressivity of income tax.
Params.tau_c=0; % Consumption tax.
Params.tau_w=0; % Earnings (wage income) tax.
Params.tau_k=0; % Capital income tax
Params.z=0; % Fraction of investment which can be deducted from corporate income tax.

% Firm production function
Params.epsilon=0.25; % Capital share in the Cobb-Douglas production fn
Params.sigma=1; % Not actually used (have hard-coded the Cobb-Douglas production function rather than more general CES prodn fn)
Params.A=0.892657593; % Technology level
% Investment adjustment-costs
b=0; % Set to b=10 for Chapter 9
% Note: there is no capital depreciation.


%% Grids
a_grid=gpuArray(linspace(0,7,n_a(1))'); % Is impossition of NO BORROWING CORRECT HERE OR DO THEY ALLOW IT TO BE NEGATIVE???
l_grid=gpuArray(linspace(0,1,n_l)');

%% General eqm variables: give some initial values
GEPriceParamNames={'wt','rt','Gt'};
Params.wt=1; % wage rate
Params.rt=0.067; % interest rate (after corporate taxes, but before income taxes)
Params.Gt=0.7; % Government consumption

%% Initial distribution of agents at age 20 (jj=1)
jequaloneDist=zeros(n_a,1); 
jequaloneDist(1,1)=1; % All agents born with zero assets

%% Return Function
DiscountFactorParamNames={'beta'};
ReturnFn=@(l,aprime,a,z,jj,wt,rt,rho,gamma,alpha,ej,tau_phi,tau_pi,tau_c,tau_w) AuerbachKotlikoff1987_ReturnFn(l,aprime,a,z,jj,wt,rt,rho,gamma,alpha,ej,tau_phi,tau_pi,tau_c,tau_w) 
ReturnFnParamNames={'jj','wt','rt','rho','gamma','alpha','ej','tau_phi','tau_pi','tau_c','tau_w'}; %It is important that these are in same order as they appear in 'AuerbachKotlikoff1987_ReturnFn'

%% Misc
n_d=n_l;
d_grid=l_grid;

sprintf('Grid sizes:')
n_a
n_z

vfoptions.verbose=1;
vfoptions.policy_forceintegertype=1; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)

%% Now solve the value function iteration problem (this is just for trying out codes before jumping to general eqm)

tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
toc

%% Agents stationary distribution
% Based on the constant population growth rate n, we get that the
% population weights will be
Params.mewj=ones(1,Params.J);
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj);

simoptions.nsims=4*10^5;
simoptions.ncores=NCores;
simoptions.iterate=1;
AgeWeights={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeights,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);


%% General Equilibrium Conditions (aka. marketclearance)

% Steady State Aggregates (important that ordering of Names and Functions is the same)
SSvalueParamNames=struct();
SSvalueParamNames(1).Names={};
SSvalueParamNames(2).Names={'ej'};
SSvalueParamNames(3).Names={'wt','rt','ej','tau_phi','tau_pi','tau_c','tau_w','tau_k'};
SSvaluesFn_1 = @(d_val,aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
SSvaluesFn_2 = @(d_val,aprime_val,a_val,z_val,ej) (1-d_val)*ej; % Aggregate labour supply (in efficiency units)
SSvaluesFn_3 = @(d_val,aprime_val,a_val,z_val, wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w,tau_k) AuerbachKotlikoff1987_TotalTaxRevenue(d_val,aprime_val,a_val,z_val,wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w,tau_k); % Total tax revenues
SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_3};

% General Equilibrium Equations
GeneralEqmEqnParamNames=struct();
GeneralEqmEqnParamNames(1).Names={'epsilon'};
GeneralEqmEqn_1 = @(AggVars,p,epsilon) p(1)-(1-epsilon)*(AggVars(1)^epsilon)*(AggVars(2)^(-epsilon)); % The requirement that the wage rate equals the marginal product of (efficiency units of) labor.
GeneralEqmEqnParamNames(2).Names={'epsilon'};
GeneralEqmEqn_2 = @(AggVars,p,epsilon) p(2)-(epsilon)*(AggVars(1)^(epsilon-1))*(AggVars(2)^(1-epsilon)); % Rate of return on assets is related to Marginal Product of Capital and Tobin's q
%^PRESENTLY IGNORES TOBINS q
GeneralEqmEqnParamNames(3).Names={};
GeneralEqmEqn_3 = @(AggVars,p) p(3)-AggVars(3); % Government balanced budget constraint
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2,GeneralEqmEqn_3};


%% Test
SSvalues_AggVars=SSvalues_AggVars_FHorz_Case1(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)

%% Solve for the General Equilibrium
% Use the toolkit to find the equilibrium price index

% GEPriceParamNames={'wt','rt','Gt'}; % Already delared above.
wt_grid=linspace(0.5,2,21)'*Params.wt; %bad
rt_grid=linspace(0.5,2,21)'*Params.rt; %okay
Gt_grid=linspace(0.5,2,21)'*Params.Gt; %good
p_grid=[wt_grid,rt_grid,Gt_grid];

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
n_p=[length(wt_grid),length(rt_grid),length(Gt_grid)];
heteroagentoptions.pgrid=p_grid;
heteroagentoptions.verbose=1;
[p_eqm_init,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeights,n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% findeqmtime=toc
Params.wt=p_eqm_init(1);
Params.rt=p_eqm_init(2);
Params.Gt=p_eqm_init(3);

save ./SavedOutput/AuerbachKotlikoff1987_GE.mat Params
 
%% Get Value fn, policy fn, agents dist in GE.
%Params.wt=1; % wage rate
%Params.rt=0.067; % interest rate (after corporate taxes, but before income taxes)

[V_init, Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeights,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

% Rather than have the population normalized to 1 (I prefer this as then
% the model can be thought of in terms of probability distributions)
% Auerbach & Kotlikoff (1987) normalize their population so that the mass
% of agents of 'age' j=1 is mass one. The following is a correcting factor:
Params.PopnCorrection=gather(1/sum(sum(sum(StationaryDist_init(:,:,1)))));

%% For baseline (initial) steady-state, replicate Table 5.1 and Figure 5.2

SSvalueParamNames(4).Names={'wt','rt','ej','tau_phi','tau_pi','tau_c','tau_w'};
SSvaluesFn_4 = @(d_val,aprime_val,a_val,z_val,wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w) AuerbachKotlikoff1987_ConsumptionFn(d_val,aprime_val,a_val,z_val,wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w); % Private Consumption
SSvalueParamNames(5).Names={'delta'};
SSvaluesFn_5 = @(d_val,aprime_val,a_val,z_val,delta) aprime_val-(1-delta)*a_val; % National Savings rate (gross rate)
SSvalueParamNames(6).Names={'wt','ej'};
SSvaluesFn_6 = @(d_val,aprime_val,a_val,z_val,wt,ej) wt*(1-d_val)*ej; % Earnings (pre-tax)
SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_3,SSvaluesFn_4,SSvaluesFn_5,SSvaluesFn_6};
SSvalues_AggVars=SSvalues_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)
SSvalues_AggVars=Params.PopnCorrection*SSvalues_AggVars; % Renormalize to AK1987 population size

if Params.sigma==1 % Cobb-Douglas
    Y=Params.A*(SSvalues_AggVars(1)^Params.epsilon)*(SSvalues_AggVars(2)^(1-Params.epsilon));
end

% Table 5_1
FID = fopen('./SavedOutput/LatexInputs/AuerbachKotlikoff_Table5_1.tex', 'w');
fprintf(FID, 'The Base Case Steady-State \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcclc} \n \\hline \\hline \n');
fprintf(FID, 'Capital Stock & %8.1f  & \\quad & Private Consumption & %8.2f \\\\ \n \\hline \n', SSvalues_AggVars(1),SSvalues_AggVars(4));
fprintf(FID, 'Labor Supply  & %8.1f  & & Capital-Output ratio & %8.2f \\\\ \n \\hline \n', SSvalues_AggVars(2),SSvalues_AggVars(1)/Y);
fprintf(FID, 'Wage & %8.3f  &  & National Savings rate & %8.2f \\%% \\\\ \n \\hline \n',          Params.wt,(Params.delta*SSvalues_AggVars(1))/Y);
fprintf(FID, 'Pre-tax interest rate & %8.2f \\%%  &   & Income Tax rate & %8.2f \\%% \\\\ \n \\hline \n', 100*Params.rt,100*Params.tau_phi);
fprintf(FID, 'National Income & %8.2f  &   & Social Security Tax rate & %8.2f \\%% \\\\ \n \\hline \n', Y,0);
fprintf(FID, 'Government Consumption & %8.2f  &   & Social Security replacement rate & %8.2f \\%% \\\\ \n \\hline \n', Params.PopnCorrection*Params.Gt,0);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on grid sizes of $n_l=%d$ for leisure, and $n_a=%d$ for assets. \\\\ \n', n_d, n_a);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

SimLifeCycleProfiles=SimLifeCycleProfiles_FHorz_Case1(jequaloneDist,Policy_init, SSvaluesFn,SSvalueParamNames,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z, simoptions);

% Figure 5.2
figure(1)
plot(1:1:N_j, SimLifeCycleProfiles(4,:,1),1:1:N_j, SimLifeCycleProfiles(6,:,1))
legend('Consumption','Earnings')
saveas(gcf,'./SavedOutput/Graphs/Figure_5_2.pdf')

% plot(Params.ej)


%% Calculate some transition paths

% Auerbach & Kotlikoff (1987) define the tax reforms in Chapter 5
% based on keeping government spending the same and setting tax rates to 
% raise this revenue.
% To implement this we have to switch our definition of the GE prices, so
% that instead of finding G we find the tax rate (for the relevant reform).

%% Consumption Tax

%% Calculate the final equilibrium
Params.tau_phi=0;
GEPriceParamNames={'wt','rt','tau_c'};

% First, calculate the final general equilibrium.
GeneralEqmEqnParamNames(3).Names={'Gt'};
GeneralEqmEqn_3 = @(AggVars,p, Gt) Gt-AggVars(3); % Government balanced budget constraint (role of GEPriceParamNames already in AggVars(3) which is total tax revenues)
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2,GeneralEqmEqn_3};

tauc_grid=linspace(0,0.2,21)';
p_grid=[wt_grid,rt_grid,tauc_grid];

disp('Calculating price vector corresponding to the final stationary eqm')
n_p=[length(wt_grid),length(rt_grid),length(tauc_grid)];
heteroagentoptions.pgrid=p_grid;
heteroagentoptions.verbose=1;
[p_eqm_final,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeights,n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

Params.wt=p_eqm_final(1);
Params.rt=p_eqm_final(2);
Params.tau_c=p_eqm_final(3);
[V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);


%%
save ./SavedOutput/AuerbachKotlikoff1987_TransitionPathSetup.mat V_final StationaryDist_init n_d n_a n_z N_j pi_z d_grid a_grid z_grid ReturnFn SSvaluesFn GeneralEqmEqns Params DiscountFactorParamNames ReturnFnParamNames SSvalueParamNames GeneralEqmEqnParamNames p_eqm_init p_eqm_final

% load ./SavedOutput/AuerbachKotlikoff1987_TransitionPathSetup.mat V_final StationaryDist_init n_d n_a n_z N_j pi_z d_grid a_grid z_grid ReturnFn SSvaluesFn GeneralEqmEqns Params DiscountFactorParamNames ReturnFnParamNames SSvalueParamNames GeneralEqmEqnParamNames p_eqm_init p_eqm_final

%% Now, the transition path
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)
Params.tau_phi=0.15;
Params.wt=p_eqm_init(1);
Params.rt=p_eqm_init(2);
Params.tau_c=0;

p_init=[p_eqm_init(1);p_eqm_init(2);0];
p_final=p_eqm_final;

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=150

% We want to look at a one off unanticipated change of tax rate tau_phi. ParamPath & PathParamNames are thus given by
ParamPath=zeros(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)
ParamPathNames={'tau_phi'}; % This is the parameter that gets changed 'away' from it's initial value.
% Consumption tax is considered a GEPriceParam, hence why it is not here.

% We need to give an initial guess for the price path on interest rates
PricePath0_1=[linspace(p_init(1), p_final(1), floor(T/2))'; p_final(1)*ones(T-floor(T/2),1)]; 
PricePath0_2=[linspace(p_init(2), p_final(2), floor(T/2))'; p_final(2)*ones(T-floor(T/2),1)]; 
PricePath0_3=[linspace(p_init(3), p_final(3), floor(T/2))'; p_final(3)*ones(T-floor(T/2),1)];
PricePath0=[PricePath0_1,PricePath0_2,PricePath0_3];% PricePath0 is matrix of size T-by-'number of prices'
PricePathNames={'rt','wt','tau_c'};

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;

%%
vfoptions.policy_forceintegertype=1
PricePathNew=TransitionPath_Case1_Fhorz(PricePath0, PricePathNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, AgeWeights, SSvalueParamNames, GeneralEqmEqnParamNames,transpathoptions);




