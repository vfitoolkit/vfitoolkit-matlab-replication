% Replication of Restuccia & Urrutia (2004) - Intergenerational Persistence of Earnings: The Role of Early and College Education

% FOR FINAL VERSION, INCREASE SOME GRIDS, INCLUDING n_p

% Open questions relating to replication:
% What is the value of tauchenq used by RU2004? (Based on first column of
% Table 2 I appear to have same value as RU2004; q=3).
% Am I correct in assuming pi (their notation, here denoted b) remains constant between j=1 and j=2, and
% then only changes following AR(1) between j=2 and j=1?

% Note: I redefine theta to have a different meaning here to in the paper.
% RU2004: theta=0 for partial uni, =1 for completed uni. Here: theta=0 for
% no uni, =1 for partial uni, =2 for completed uni. (I also introduce a
% bhatprime, which in RU2004 notation is just called bhat)

% Solving this model requires a substantial rewrite of the value function
% problem into the form recognised by the VFI Toolkit (essentially, into the
% notation of Stokey, Lucas & Prescott, 1989). 

% The model is a general equilibrium 2-period OLG
% The value function itself is a Case 2, 2-period, finite horizon problem.
% Case 2, because of both the shocks to (the endogenous state) human
% capital between periods 1 and 2, and the shocks to (the endogenous state)
% of whether a person who goes to university completes their degree or only
% partially completes their degree (theta in notation of RU2004).

% The VFI Toolkit does not allow the number of endogenous states, exogenous
% states, and decision variables to be different for different ages but it can
% do so implicitly. It does allow for the number of points and the grid values to change
% depending on age (so the actually interpretation of, e.g., an endogenous
% state could change completely from one age to another. By setting the
% number of points on for a given state/decision to 1, the VFI Toolkit can
% essentially react as if the number of state/decision variables is
% changing for different ages. This trick will be used in this replication,
% it is not needed per-se to solve the model, but does make it much faster.

% Two decision variables: acquired ability of child (bhatprime, only j==1)
%                         schooling decision (s, only j==1)
% Three endogenous state variables: human capital (h, both j==1 & 2), 
%                                   acquired ability of child (bhat, only j==2),
%                                   schooling (s, only j==2)
% Two stochastic exogenous variables: innate ability of child (b, both j==1 & 2)
%                                     college completion shock (thetahat, only j==2, determines the value of theta)

% Age
% Note: Since the return function is different for different ages (not just
% that some parameter values differ by age), a parameter for age will need
% to be passed to the return function (this is not done automatically as
% many models do not require it).

Params.J=2; % Ages 1 and 2

% A few lines needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
% addpath(genpath('./subcodes/'))
try % Server has 20 cores, but is shared with other users, so use max of 12.
    parpool(12)
    gpuDevice(1)
catch % Desktop has less than 12, so will give error, on desktop it is fine to use all available cores.
    parpool
end
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;
simoptions.ncores=NCores;

% simoptions.parallel=3;

%% Declare the model parameters
% Note that w and g will be determined in General equilbrium, so these are really just initial guesses.

% Age
Params.agej=[1,2];

% Preference
Params.beta=0.52; % rate at which agents discount future
Params.sigma=1.5; % Risk-aversion

% Technology
Params.A=1;

% Education
Params.gamma=0.24;
Params.f=0.64;
Params.nupperbar=0.25;
Params.nlowerbar=0.125;
Params.pupperbar=1.48;
Params.plowerbar=0.86;
Params.psi0=0.27;
Params.psi1=1.02;

% Life cycle
Params.zeta=1.12;

% Income tax and university subsidy
Params.tau=0.039; % =(g+kappa*F)/y, see bottom of pg 1361 of UR2004; This eqn is one of the GE constraints of the model and will determine g
Params.kappa0=0.36;
Params.kappa1=1; % RU2004 don't mention this in baseline model, but is required for the policy experiment. See their description of the second policy experiment (their pg 1373).

% Stochastic process for innate abilities
Params.sigma_b=0.48; % RU2004 contains a typo, their model description has sigma_b as the standard deviation of epsilon in first eqn on pg 1359 (std dev of innovation to the AR(1) in logs process on innate ability b). But it is clear from Table 2 that the reported calibrated value in Table 1 is actually the standard deviation of log(b) itself, and not of the innovations.
Params.rho_b=0.2; % UR2004 just call this rho

% The earnings shock: only used for Table 7, not in baseline model.
Params.earningsshocksize=0.2; % 0.2 for 'small shock', 0.5 for 'large shock'

%% General eqm variables: give some initial values
GEPriceParamNames={'g'};
Params.g=0.021; % Public spending on education
Params.w=1; % Wage (no need to actually determine this as a GE price, as we know from the equations that it must be equal to 1.)
            % Firms production function is constant returns to scale in H, and so as long as w=A (=1) the firm is indifferent about H^f.


%% The size of the grids depends on age
vfoptions.agedependentgrids=[1,1,1];
% These three binary numbers indicate age dependence for d (decision variables), a (endogenous
% state variables), and z (exogenous state) variables repectively. When a
% grid depends on age it must be passed as either a function or a structure. When a grid does not
% depend on age it can be passed as usual, or as a function/structure. It is always
% assumed that when z is age dependent that both the grid and the
% transition matrix are age dependent; these should be passed as two
% outputs of a single function in place of z_grid, and the pi_z input
% should be used to pass the relevant parameter names 'AgeDependentGridParamNames'.

% % % Grid sizes to use
% RU2004, Notes:  We make the state space discrete. We construct log-spaced grids for human capital with 60 points, 
% innate ability with 15 points, and acquired ability with 60 points.
n_bhatprime=150;
n_s=2;

n_h=150;
n_bhat=150;
%n_s=2

n_b=15; % 15 is the value used by RU2004
n_thetahat=51; % A uniform [0,1] random variable used to determine theta according to whether thetahat is greater or less than q(bhat)
n_earnings=1; % The earnings shock required for Table 7 (no earnings shock in baseline model)
              % n_earnings can only take value of 1 or 2 (baseline, and Table 7 respectively)

N_j=Params.J; % Number of periods in finite horizon

% When sizes depend on age, you express the dependence on age as different column for each age.
n_d=[n_bhatprime,1; n_s, 1]; % bhatprime,s
n_a=[n_h, n_h; 1, n_bhat;1,n_s]; % h,bhat,s
n_z=[n_b,n_b; 1,n_thetahat; n_earnings, n_earnings]; % b, thetahat % Independent of age, although the transition probabilities are dependent on age (by my reading b is unchanged from age 1 to 2, and then changes from 2 to 1; I am not 100% certain of this interpretation of RU2004)


% % % Grids
% RU2004, Notes:  "We make the state space discrete. We construct log-spaced grids for human capital with 60 points, innate ability with 15 points, and acquired ability with 60 points."
% RU2004 Tech Appendix: 
% "We choose the grid for acquired ability bhat with a lower bound equal to bhat(1)=b(1)*g^gamma and upper bound bhat(60)= b(15)*(5g)^gamma. 
%  The human capital grid has lower bound h(1) = bhat(1) and upper bound h(60) = pupperbar*bhat(60)."
% Note: the choice of lower and upper grid limits of bhat are based on the equation in bottom left hand side of pg 1358 that determines bhat.
% Note: the upper grid limit of h is then based on eqn at bottom right of pg 1359 about human capital and college completion
% Remark: allowing a grid to depend on a parameter to be determined in
% general equilibrium is not good practice, so rather than 'g' itself, I use a fixed 'g_forgrid' to give same value in all cases.

% A few parameters the grids will need
Params.tauchenq=3; % Based on email communication with Carlos Urrutia it appears they set q=3
Params.n_b=n_b;
Params.n_bhat=n_bhat;
Params.n_h=n_h;
Params.n_thetahat=n_thetahat;
Params.n_earnings=n_earnings;
Params.g_forgrid=0.021; % This is the initial guess for g, and close the the eqm value of RU2004.
Params.maxgridptparam=5; % RU2004 use 5, I use 8. It is a hyperparameter that determines the max values of the bhat and h grids.
% Function for the age dependent d_grid
AgeDependentGridParamNames.d_grid={'agej','g_forgrid','gamma', 'sigma_b','rho_b', 'n_b','n_bhat' 'tauchenq','maxgridptparam'};
d_gridfn=@(agej,g,gamma, sigma_b,rho_b, n_b,n_bhat, tauchenq,maxgridptparam) RestucciaUrrutia2004_d_gridfn(agej,g,gamma, sigma_b,rho_b, n_b,n_bhat, tauchenq,maxgridptparam); 
% Function for the age dependent a_grid
AgeDependentGridParamNames.a_grid={'agej','g_forgrid','gamma','pupperbar', 'sigma_b','rho_b', 'n_b','n_bhat','n_h' 'tauchenq','maxgridptparam'};
a_gridfn=@(agej,g,gamma,pupperbar, sigma_b,rho_b, n_b, n_bhat, n_h, tauchenq,maxgridptparam) RestucciaUrrutia2004_a_gridfn(agej,g,gamma, pupperbar, sigma_b,rho_b, n_b,n_bhat,n_h, tauchenq,maxgridptparam);
% The function for the age dependent z_grid, also returns age dependent pi_z
AgeDependentGridParamNames.z_grid={'agej', 'sigma_b','rho_b', 'earningsshocksize', 'n_b', 'n_thetahat','n_earnings', 'tauchenq'};
z_gridfn=@(agej, sigma_b,rho_b, earningsshocksize, n_b, n_thetahat,n_earnings, tauchenq) RestucciaUrrutia2004_z_gridfn(agej, sigma_b,rho_b, earningsshocksize, n_b, n_thetahat,n_earnings, tauchenq);

% We need the lower and upper bounds of h_grid and bhat_grid as they are some of
% the inputs to the Phi_aprime function (see later). So create and get those.
agej=2;
a_grid=RestucciaUrrutia2004_a_gridfn(agej,Params.g_forgrid,Params.gamma, Params.pupperbar,Params.sigma_b,Params.rho_b, Params.n_b,Params.n_bhat,Params.n_h, Params.tauchenq, Params.maxgridptparam);
Params.log_lb_h=log(a_grid(1));
Params.log_ub_h=log(a_grid(n_h));
Params.loghgridspacing=log(a_grid(2))-log(a_grid(1)); % The h grid is uniform in logs
Params.log_lb_bhat=log(a_grid(n_h+1));
Params.log_ub_bhat=log(a_grid(n_h+n_bhat));
Params.logbhatgridspacing=log(a_grid(n_h+2))-log(a_grid(n_h+1)); % The bhat grid is uniform in logs

% To take a look at all the grids
vfoptions.parallel=2; % Elsewhere this would be done automatically inside the various commands as it is the default
daz_gridstructure=AgeDependentGrids_Create_daz_gridstructure(n_d,n_a,n_z,N_j,d_gridfn, a_gridfn, z_gridfn, AgeDependentGridParamNames, Params, vfoptions);


%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(bhatprime,sprime,h,bhat,s,b,thetahat,earningsshock,agej,w,g,gamma,tau,sigma,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_ReturnFn(bhatprime,sprime,h,bhat,s,b,thetahat,earningsshock,agej,w,g,gamma,tau,sigma,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1)
ReturnFnParamNames={'agej','w','g','gamma','tau','sigma','kappa0','kappa1','plowerbar','pupperbar','nlowerbar','nupperbar','f','psi0','psi1'}; %It is important that these are in same order as they appear in 'RestucciaUrrutia2004_ReturnFn'

%% Phi, used by Case2 to figure out which next-period-endogenous-state
% Set it up as a function (could get user to set it up as a sparse matrix,
% but don't want sparse matrices to be something user needs to understand,
% also want different state space for each age so would need a structure of
% sparse matrices (one for each age). Function is just much easier.
Case2_Type=12; %  phi_aprime(d,a,z)
vfoptions.phiaprimedependsonage=1;
Phi_aprime=@(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar,pupperbar,zeta,psi0,psi1,log_lb_h,log_ub_h,loghgridspacing,n_h,log_lb_bhat,log_ub_bhat,logbhatgridspacing,n_bhat) RestucciaUrrutia2004_PhiaprimeFn(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar,pupperbar,zeta,psi0,psi1,log_lb_h,log_ub_h,loghgridspacing,n_h,log_lb_bhat,log_ub_bhat,logbhatgridspacing,n_bhat)
PhiaprimeParamNames={'agej','plowerbar','pupperbar','zeta','psi0','psi1','log_lb_h','log_ub_h','loghgridspacing','n_h','log_lb_bhat','log_ub_bhat','logbhatgridspacing','n_bhat'};

%%
vfoptions.dynasty=1;
vfoptions.lowmemory=1
vfoptions.verbose=1;
vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)


%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

disp('Test ValueFnIter')
tic;
[V, Policy]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_gridfn, a_gridfn, z_gridfn, AgeDependentGridParamNames, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, vfoptions);
toc
% CHALLENGES: Need both for each generation to care about the next, and to
% allow the dimensionality of the problem to depend on age. These should
% both be set up within ValeFnIter_Case2_Fhorz. Caring about next
% generation can be set using vfoptions, dimensionality depending on age
% can be detected automatically based on inputs and will cause V and Policy
% to be structures (rather than matrices), with each age being a different field of the
% structure and containing a matrix which is the value (or policy) function for that age.

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.
% sum(sum(sum(sum(Policy==n_a))))

surf(permute(V.j001,[1,4,2,3]))

%surf(permute(V.j001,[1,4,2,3])-permute(V2.j001,[1,4,2,3]))

%% Initial distribution of agents at birth (j=1)
% When using 'dynasty' this is less an 'initial distribution', and instead
% just an initial guess, as the actual distribution will anyway be
% independent of this.
jequaloneDist=zeros([n_a(:,1)',n_z(:,1)'],'gpuArray');
jequaloneDist(sub2ind_homemade([n_a(:,1)',n_z(:,1)'],[30,1,1,8,1,1]))=1; % This will do, is anyway fairly unimportant

%% Agents age distribution
% This model has neither stochastic probability of death, nor population growth, so just.
Params.mewj=ones(1,Params.J)/2;

simoptions.numbersims=10^4;
simoptions.iterate=1;
AgeWeightParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

simoptions.agedependentgrids=vfoptions.agedependentgrids;
simoptions.dynasty=vfoptions.dynasty;
simoptions.phiaprimedependsonage=vfoptions.phiaprimedependsonage;

simoptions.verbose=1;

simoptions.parallel=3

%% Test
disp('Test StationaryDist')

% simoptions.dynasty=1
tic;
StationaryDist=StationaryDist_FHorz_Case2(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,d_gridfn, a_gridfn, z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,Params,PhiaprimeParamNames,simoptions);
toc

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% StationaryDist moments (important that ordering of Names and Functions is
% the same). Things that are needed as part of evaluating the general equilibrium conditions.
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'agej','plowerbar', 'pupperbar','nlowerbar','nupperbar','psi0','psi1'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1) RestucciaUrrutia2004_HFn(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1); % H, Human capital
FnsToEvaluateParamNames(2).Names={'agej','w','kappa0','kappa1','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_kappaF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1); % kappaF (a3_val==1, is that the state s=1 for old)
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2};
% Note that the aggregate labour supply is actually entirely exogenous and so I could just precompute it, but am feeling lazy.

% General Equilibrium Equations
% Recall that GEPriceParamNames={'g'}; In following lines p is the vector of these and so, e.g., p(1) is g.
GeneralEqmEqnParamNames=struct();
GeneralEqmEqnParamNames(1).Names={'tau','A'};
GeneralEqmEqn_1 = @(AggVars,p,tau,A) p(1)+AggVars(2)-tau*A*AggVars(1); % g+kappaF-tau*Y=0 (using Y=A*H): Government budget balance
GeneralEqmEqns={GeneralEqmEqn_1};

% Age 1: 
% RestucciaUrrutia2004_H(0.5,1,0.5,1,1,1,1,1,0.86,1.48, 0.125,0.25,0.27,1.02)

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_gridfn, a_gridfn, z_gridfn, simoptions, AgeDependentGridParamNames);
disp('Test AggVars: Done')
GeneralEqmConditionsVec=real(GeneralEqmConditions_Case2(AggVars,[Params.g], GeneralEqmEqns, Params,GeneralEqmEqnParamNames));

AggVars

%% Solve for the General Equilibrium
% Use the toolkit to find the equilibrium price index. There are two ways
% to do this. In what follows I use the 'search' approach to calculate the
% (initial) General equilibrium. But the commented-out lines that follow it
% show how to set up the grid approach.

% Things seem to work, so lets cut back on the printed feedback
vfoptions.verbose=0;
simoptions.verbose=0;

% Without p_grid, just searching. Use n_p=0. (Setting the actual algorithm
% used to 'search' can be done with heteroagentoptions.fminalgo)
heteroagentoptions.verbose=1;
n_p=51; % MAKE THIS MORE LIKE 101 for final replication codes
if n_p>0
    heteroagentoptions.p_grid=linspace(0.0001,0.05,n_p);
end
[p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case2_FHorz(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.g=p_eqm(1);
save ./SavedOutput/RestucciaUrrutia2004.mat Params p_eqm

% plot(heteroagentoptions.p_grid,GeneralEqmEqnsValues)

% Using p_grid. This can be helpful if you want to, e.g., look for
% possibility of multiple equilibria.
% % GEPriceParamNames={'r','b','T'}; % Already delared above.
% r_grid=linspace(0.5,2,21)'*Params.r;
% b_grid=linspace(0.5,2,21)'*Params.b;
% T_grid=linspace(0.5,2,21)'*Params.T;
% p_grid=[r_grid,b_grid, T_grid];
% 
% disp('Calculating price vector corresponding to the stationary eqm')
% n_p=[length(r_grid),length(b_grid),length(T_grid)];
% heteroagentoptions.pgrid=p_grid;
% heteroagentoptions.verbose=1;
% [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeights,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% Params.r=p_eqm(1);
% Params.b=p_eqm(2);
% Params.T=p_eqm(3);
% save ./SavedOutput/Huggett1996grid.mat Params

%% Compute a few things about the equilibrium
load ./SavedOutput/RestucciaUrrutia2004.mat Params p_eqm

[V, Policy]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_gridfn, a_gridfn, z_gridfn, AgeDependentGridParamNames, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, vfoptions);
StationaryDist=StationaryDist_FHorz_Case2(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,d_gridfn, a_gridfn, z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,Params,PhiaprimeParamNames,simoptions);

% surf(permute(StationaryDist.j001,[1,4,2,3]))
save ./SavedOutput/RestucciaUrrutia2004_2.mat V Policy StationaryDist

% %% Take a look at a few things as part of sanity checking results
% 
% % Policies
% Policy.j001(1,:,1,1,1) % acquired ability decisions (note: do people try to leave the grid?)
% Policy.j001(1,:,1,1,8)
% Policy.j001(1,:,1,1,15)
% Policy.j001(2,:,1,1,1) % Enrollment decisions 
% Policy.j001(2,:,1,1,8)
% Policy.j001(2,:,1,1,15)
% % Where in distributions of h, and h,bhat are everyone.
% plot(cumsum(sum(sum(sum(StationaryDist.j001,4),3),2))) % human capital, h
% plot(shiftdim(cumsum(sum(sum(sum(StationaryDist.j001,1),2),3)),3)) % innate ability b
% plot(cumsum(sum(sum(sum(sum(StationaryDist.j002,5),4),3),2))) % human capital, h
% plot(cumsum(sum(sum(sum(sum(StationaryDist.j002,5),4),3),1))) % acquired ability bhat
% surf(sum(sum(sum(StationaryDist.j002,5),4),3))

%%

load ./SavedOutput/RestucciaUrrutia2004_2.mat V Policy StationaryDist
% The following are all just based off the baseline model
RestucciaUrrutia2004_Tables1to4 % Actually includes Figures 1 and 2 as well
simoptions.parallel=3; % This is what it was prior to RestucciaUrrutia2004_Table1to4, but it gets changed by RestucciaUrrutia2004_Table1to4.


%%
% No longer need these
clear V Policy StationaryDist

% For Tables 5 to 8, for the baseline model
disp('Now doing Table Column 1')
TableColumn(1)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% RestucciaUrrutia2004_TableColumnFn() creates all the statistics reported in the rows of any of Tables 5 to 8.

%% Tables 5 and 6 require doing a bunch of robustness checks to alternative parameter values

% Table 5
rho_vec=[0.1,0.2,0.3];
for rho_c=1:length(rho_vec)
    Params.rho_b=rho_vec(rho_c);
    TableColumn(rho_c+1)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.rho_b=0.2; % Return to baseline value
sigmab_vec=[0.4,0.48,0.6];
for sigmab_c=1:length(sigmab_vec)
    Params.sigma_b=sigmab_vec(sigmab_c);
    TableColumn(sigmab_c+4)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.sigma_b=0.48; % Return to baseline value

save ./SavedOutput/RestucciaUrrutia2004_3.mat TableColumn

FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table5.tex', 'w');
fprintf(FID, 'Sensitivity Analysis with Respect to $\\rho$ and $\\sigma_{\\pi}$ \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccc} \n \\hline \\hline \n');
fprintf(FID, ' & Benchmark & & & &  \\\\ \n');
fprintf(FID, ' & $\\rho=0.2$ & $\\rho=0.1$ & $\\rho=0.3$ & $\\sigma_{\\pi}=0.4$ & $\\sigma_{\\pi}=0.6$ \\\\ \n');
fprintf(FID, ' & $\\sigma_{\\pi}=0.48$ & & & & \\\\ \\hline \n');
fprintf(FID, 'Intergenerational correlation & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfInnateAbility, TableColumn(2).IntergenerationalCorrelationOfInnateAbility, TableColumn(4).IntergenerationalCorrelationOfInnateAbility, TableColumn(5).IntergenerationalCorrelationOfInnateAbility, TableColumn(7).IntergenerationalCorrelationOfInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(2).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(4).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(5).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(7).IntergenerationalCorrelationOfAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfEarnings, TableColumn(2).IntergenerationalCorrelationOfEarnings, TableColumn(4).IntergenerationalCorrelationOfEarnings, TableColumn(5).IntergenerationalCorrelationOfEarnings, TableColumn(7).IntergenerationalCorrelationOfEarnings);
fprintf(FID, 'Disparity: std(log x)     & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogInnateAbility, TableColumn(2).StdDevLogInnateAbility, TableColumn(4).StdDevLogInnateAbility, TableColumn(5).StdDevLogInnateAbility, TableColumn(7).StdDevLogInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogAcquiredAbility, TableColumn(2).StdDevLogAcquiredAbility, TableColumn(4).StdDevLogAcquiredAbility, TableColumn(5).StdDevLogAcquiredAbility, TableColumn(7).StdDevLogAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogEarnings, TableColumn(2).StdDevLogEarnings, TableColumn(4).StdDevLogEarnings, TableColumn(5).StdDevLogEarnings, TableColumn(7).StdDevLogEarnings);
fprintf(FID, 'Other aggregate statistics    & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ College enrollment   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).CollegeEnrollmentRate, TableColumn(2).CollegeEnrollmentRate, TableColumn(4).CollegeEnrollmentRate, TableColumn(5).CollegeEnrollmentRate, TableColumn(7).CollegeEnrollmentRate);
fprintf(FID, '$\\quad$ Private early/GDP (as \\%)  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', 100*TableColumn(1).PrivateEarlyEducDivGDP, 100*TableColumn(2).PrivateEarlyEducDivGDP, 100*TableColumn(4).PrivateEarlyEducDivGDP, 100*TableColumn(5).PrivateEarlyEducDivGDP, 100*TableColumn(7).PrivateEarlyEducDivGDP);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: All other parameters are the same as in the benchmark economy; public expenditures in early education adjusted to balance the government budget. \\\\ \n');
fprintf(FID, 'Note: In model notation these columns are: $b$, $\\hat{b}$, and $wh$. I follow Restuccia \\& Urrutia (2004) in reporting as cross-sectional the number conditional on being an elderly household, not cross-sectional over the whole model economy. \\\\ \n');
fprintf(FID, 'Restuccia \\& Urrutia (2004) explain calculation of intergenerational correlation of earnings at bottom of pg 1363. I assume the intergeneration correlations of (log) innate and acquired ability are calculated by the analagous regressions (with modification for acquired ability as is only observed for old). \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% TableColumn(1).CollegeEnrollmentRate=TableColumn(1).CollegeEnrollmentRate.Mean(1)

% Table 6
gamma_vec=[0.1,0.24,0.4];
for gamma_c=1:length(gamma_vec)
    Params.gamma=gamma_vec(gamma_c);
    TableColumn(gamma_c+7)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.gamma=0.24; % Return to baseline value
pupperbar_vec=[1.4,1.48,1.6];
for pupperbar_c=1:length(pupperbar_vec)
    Params.pupperbar=pupperbar_vec(pupperbar_c);
    TableColumn(pupperbar_c+10)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.pupperbar=1.48; % Return to baseline value

save ./SavedOutput/RestucciaUrrutia2004_3.mat TableColumn

FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table6.tex', 'w');
fprintf(FID, 'Sensitivity Analysis with Respect to $\\gamma$ and $\\bar{p}$ \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccc} \n \\hline \\hline \n');
fprintf(FID, ' & Benchmark & & & &  \\\\ \n');
fprintf(FID, ' & $\\gamma=0.24$ & $\\gamma=0.1$ & $\\gamma=0.4$ & $\\bar{p}=1.4$ & $\\bar{p}=1.6$ \\\\ \n');
fprintf(FID, ' & $\\bar{p}=1.48$ & & & & \\\\ \\hline \n');
fprintf(FID, 'Intergenerational correlation & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfInnateAbility, TableColumn(8).IntergenerationalCorrelationOfInnateAbility, TableColumn(10).IntergenerationalCorrelationOfInnateAbility, TableColumn(11).IntergenerationalCorrelationOfInnateAbility, TableColumn(13).IntergenerationalCorrelationOfInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(8).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(10).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(11).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(13).IntergenerationalCorrelationOfAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfEarnings, TableColumn(8).IntergenerationalCorrelationOfEarnings, TableColumn(10).IntergenerationalCorrelationOfEarnings, TableColumn(11).IntergenerationalCorrelationOfEarnings, TableColumn(13).IntergenerationalCorrelationOfEarnings);
fprintf(FID, 'Disparity: std(log x)     & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogInnateAbility, TableColumn(8).StdDevLogInnateAbility, TableColumn(10).StdDevLogInnateAbility, TableColumn(11).StdDevLogInnateAbility, TableColumn(13).StdDevLogInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogAcquiredAbility, TableColumn(8).StdDevLogAcquiredAbility, TableColumn(10).StdDevLogAcquiredAbility, TableColumn(11).StdDevLogAcquiredAbility, TableColumn(13).StdDevLogAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogEarnings, TableColumn(8).StdDevLogEarnings, TableColumn(10).StdDevLogEarnings, TableColumn(11).StdDevLogEarnings, TableColumn(13).StdDevLogEarnings);
fprintf(FID, 'Other aggregate statistics    & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ College enrollment         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).CollegeEnrollmentRate, TableColumn(8).CollegeEnrollmentRate, TableColumn(10).CollegeEnrollmentRate, TableColumn(11).CollegeEnrollmentRate, TableColumn(13).CollegeEnrollmentRate);
fprintf(FID, '$\\quad$ Private early/GDP (as \\%%) & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', 100*TableColumn(1).PrivateEarlyEducDivGDP, 100*TableColumn(8).PrivateEarlyEducDivGDP, 100*TableColumn(10).PrivateEarlyEducDivGDP, 100*TableColumn(11).PrivateEarlyEducDivGDP, 100*TableColumn(13).PrivateEarlyEducDivGDP);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: All other parameters are the same as in the benchmark economy; public expenditures in early education adjusted to balance the government budget. \\\\ \n');
fprintf(FID, 'Note: In model notation these columns are: $b$, $\\hat{b}$, and $wh$. I follow Restuccia \\& Urrutia (2004) in reporting as cross-sectional the number conditional on being an elderly household, not cross-sectional over the whole model economy. \\\\ \n');
fprintf(FID, 'Restuccia \\& Urrutia (2004) explain calculation of intergenerational correlation of earnings at bottom of pg 1363. I assume the intergeneration correlations of (log) innate and acquired ability are calculated by the analagous regressions (with modification for acquired ability as is only observed for old). \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Table 7 uses the 'earnings shock' (RU2004 denote this by 'z', in these codes it is referred to as 'earningsshock'
n_earnings=2;
n_z=[n_b,n_b; 1,n_thetahat; n_earnings, n_earnings]; % b, thetahat % Independent of age, although the transition probabilities are dependent on age (by my reading b is unchanged from age 1 to 2, and then changes from 2 to 1; I am not 100% certain of this interpretation of RU2004)
Params.n_earnings=n_earnings;
Params.earningsshocksize=0.2;
% When using 'dynasty' this is less an 'initial distribution', and instead just an initial guess, as 
% the actual distribution will anyway beindependent of this.
jequaloneDist=zeros([n_a(:,1)',n_z(:,1)'],'gpuArray');
jequaloneDist(sub2ind_homemade([n_a(:,1)',n_z(:,1)'],[30,1,1,8,1,1]))=0.5; % This will do, is anyway fairly unimportant
jequaloneDist(sub2ind_homemade([n_a(:,1)',n_z(:,1)'],[30,1,1,8,1,2]))=0.5; % This will do, is anyway fairly unimportant
% Note: initial dist is 50-50 over earnings shock
TableColumn(14)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.earningsshocksize=0.5;
TableColumn(15)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

n_earnings=1;Params.n_earnings=n_earnings; % Return to baseline value
save ./SavedOutput/RestucciaUrrutia2004_3.mat TableColumn

% Table 7
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table7.tex', 'w');
fprintf(FID, 'Adding Post-College Permanent Earnings Shocks \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & Benchmark & Small shock & Large shock \\\\ \n');
fprintf(FID, ' & Economy & $z=\\{0.8,1.2\\}$ & $z=\\{0.5,1.5\\}$ \\\\ \\hline \n');
fprintf(FID, 'Intergenerational correlation & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfInnateAbility, TableColumn(14).IntergenerationalCorrelationOfInnateAbility, TableColumn(15).IntergenerationalCorrelationOfInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(14).IntergenerationalCorrelationOfAcquiredAbility, TableColumn(15).IntergenerationalCorrelationOfAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfEarnings, TableColumn(14).IntergenerationalCorrelationOfEarnings, TableColumn(15).IntergenerationalCorrelationOfEarnings);
fprintf(FID, 'Disparity: std(log x)     & & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Innate ability   & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogInnateAbility, TableColumn(14).StdDevLogInnateAbility, TableColumn(15).StdDevLogInnateAbility);
fprintf(FID, '$\\quad$ Acquired ability & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogAcquiredAbility, TableColumn(14).StdDevLogAcquiredAbility, TableColumn(15).StdDevLogAcquiredAbility);
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).StdDevLogEarnings, TableColumn(14).StdDevLogEarnings, TableColumn(15).StdDevLogEarnings);
fprintf(FID, 'Correlation with log acquired ability & & & \\\\ \n ');
fprintf(FID, '$\\quad$ College enrollment     & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).Correlationwithlogacqabilt_collegeenroll, TableColumn(14).Correlationwithlogacqabilt_collegeenroll, TableColumn(15).Correlationwithlogacqabilt_collegeenroll);
fprintf(FID, '$\\quad$ Educational attainment & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).Correlationwithlogacqabilt_educattain, TableColumn(14).Correlationwithlogacqabilt_educattain, TableColumn(15).Correlationwithlogacqabilt_educattain);
fprintf(FID, '$\\quad$ Log earnings           & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).Correlationwithlogacqabilt_logearnings, TableColumn(14).Correlationwithlogacqabilt_logearnings, TableColumn(15).Correlationwithlogacqabilt_logearnings);
fprintf(FID, 'Other aggregate statistics       & & & \\\\ \n ');
fprintf(FID, '$\\quad$ College enrollment         & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).CollegeEnrollmentRate, TableColumn(14).CollegeEnrollmentRate, TableColumn(15).CollegeEnrollmentRate);
fprintf(FID, '$\\quad$ Private early/GDP (as \\%) & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).PrivateEarlyEducDivGDP, TableColumn(14).PrivateEarlyEducDivGDP, TableColumn(15).PrivateEarlyEducDivGDP);
fprintf(FID, '$\\quad$ Average College Premium    & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).AverageCollegePremium, TableColumn(14).AverageCollegePremium, TableColumn(15).AverageCollegePremium);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Same parameters as in the benchmark economy; public expenditures in early education adjusted to balance the government budget. \\\\ \n');
fprintf(FID, 'Note: In model notation these columns are: $b$, $\\hat{b}$, and $wh$. I follow Restuccia \\& Urrutia (2004) in reporting as cross-sectional the number conditional on being an elderly household, not cross-sectional over the whole model economy. \\\\ \n');
fprintf(FID, 'Restuccia \\& Urrutia (2004) explain calculation of intergenerational correlation of earnings at bottom of pg 1363. I assume the intergeneration correlations of (log) innate and acquired ability are calculated by the analagous regressions (with modification for acquired ability as is only observed for old). \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Turn the earnings shocks back off
n_earnings=1;
Params.n_earnings=n_earnings;
Params.earningsshocksize=0;
n_z=[n_b,n_b; 1,n_thetahat; n_earnings, n_earnings]; % b, thetahat % Independent of age, although the transition probabilities are dependent on age (by my reading b is unchanged from age 1 to 2, and then changes from 2 to 1; I am not 100% certain of this interpretation of RU2004)
% When using 'dynasty' this is less an 'initial distribution', and instead just an initial guess, as 
% the actual distribution will anyway beindependent of this.
jequaloneDist=zeros([n_a(:,1)',n_z(:,1)'],'gpuArray');
jequaloneDist(sub2ind_homemade([n_a(:,1)',n_z(:,1)'],[30,1,1,8,1,1]))=1; % This will do, is anyway fairly unimportant


%% Table 8: involves three policy experiments
% First policy experiment: Increase in Early Education Expenditures
Params.tau=0.047;
% Params.g % g is anyway a general equilibrium variable and so determined there.
TableColumn(16)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.tau=0.039; % Return to baseline value

% Second policy experiment: Increase in College Subsidies
Params.tau=0.047;
Params.kappa1=1.025;
TableColumn(17)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.tau=0.039; % Return to baseline value
Params.kappa1=1; % Return to baseline value

% Third policy experiment: Flat College Subsidy
Params.kappa1=Params.kappa0; % RU2004 describe this as setting kappa(wh)=kappa0, but to actually implement this we simply use the following kappa0 and kappa1 values
Params.kappa0=0; 
TableColumn(18)=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.kappa0=0.36; % Return to baseline value
Params.kappa1=1; % Return to baseline value

save ./SavedOutput/RestucciaUrrutia2004_3.mat TableColumn

% Table 8
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table8.tex', 'w');
fprintf(FID, 'Policy Experiments \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, ' & & Increase & Increase & Flat \\\\ \n');
fprintf(FID, ' & & in early & in college & college \\\\ \n');
fprintf(FID, ' & Benchmark & expenditures & subsidy & subsidy \\\\ \\hline \n');
fprintf(FID, 'Intergenerational correlation & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Earnings         & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfEarnings, TableColumn(16).IntergenerationalCorrelationOfEarnings, TableColumn(17).IntergenerationalCorrelationOfEarnings, TableColumn(18).IntergenerationalCorrelationOfEarnings);
fprintf(FID, '$\\quad$ Educational attainment   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfEducationalAttainment, TableColumn(16).IntergenerationalCorrelationOfEducationalAttainment, TableColumn(17).IntergenerationalCorrelationOfEducationalAttainment, TableColumn(18).IntergenerationalCorrelationOfEducationalAttainment);
fprintf(FID, '$\\quad$ Consumption & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).IntergenerationalCorrelationOfConsumption, TableColumn(16).IntergenerationalCorrelationOfConsumption, TableColumn(17).IntergenerationalCorrelationOfConsumption, TableColumn(18).IntergenerationalCorrelationOfConsumption);
fprintf(FID, 'Expenditures (percent of GDP)  & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ Private early education  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).ExpendAsPercentOfGDP_PrivateEarlyEducation, TableColumn(16).ExpendAsPercentOfGDP_PrivateEarlyEducation, TableColumn(17).ExpendAsPercentOfGDP_PrivateEarlyEducation, TableColumn(18).ExpendAsPercentOfGDP_PrivateEarlyEducation);
fprintf(FID, '$\\quad$ Public early eduction    & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).ExpendAsPercentOfGDP_PublicEarlyEducation, TableColumn(16).ExpendAsPercentOfGDP_PublicEarlyEducation, TableColumn(17).ExpendAsPercentOfGDP_PublicEarlyEducation, TableColumn(18).ExpendAsPercentOfGDP_PublicEarlyEducation);
fprintf(FID, '$\\quad$ Private college education & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).ExpendAsPercentOfGDP_PrivateCollegeEducation, TableColumn(16).ExpendAsPercentOfGDP_PrivateCollegeEducation, TableColumn(17).ExpendAsPercentOfGDP_PrivateCollegeEducation, TableColumn(18).ExpendAsPercentOfGDP_PrivateCollegeEducation);
fprintf(FID, '$\\quad$ Public college education & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).ExpendAsPercentOfGDP_PublicCollegeEducation, TableColumn(16).ExpendAsPercentOfGDP_PublicCollegeEducation, TableColumn(17).ExpendAsPercentOfGDP_PublicCollegeEducation, TableColumn(18).ExpendAsPercentOfGDP_PublicCollegeEducation);
fprintf(FID, 'Other aggregate statistics    & & & & \\\\ \n ');
fprintf(FID, '$\\quad$ College enrollment rate & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).CollegeEnrollmentRate, TableColumn(16).CollegeEnrollmentRate, TableColumn(17).CollegeEnrollmentRate, TableColumn(18).CollegeEnrollmentRate);
fprintf(FID, '$\\quad$ College dropout rate    & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).CollegeDropoutRate, TableColumn(16).CollegeDropoutRate, TableColumn(17).CollegeDropoutRate, TableColumn(18).CollegeDropoutRate);
fprintf(FID, '$\\quad$ Aggregate human capital & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).AggregateHumanCapital, TableColumn(16).AggregateHumanCapital, TableColumn(17).AggregateHumanCapital, TableColumn(18).AggregateHumanCapital);
fprintf(FID, '$\\quad$ Aggregate consumption   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n ', TableColumn(1).AggregateConsumption, TableColumn(16).AggregateConsumption, TableColumn(17).AggregateConsumption, TableColumn(18).AggregateConsumption);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: All other parameters are the same as in the benchmark economy; public expenditures in early education adjusted to balance the government budget. \\\\ \n');
fprintf(FID, 'Note: In model notation these columns are: $b$, $\\hat{b}$, and $wh$. I follow Restuccia \\& Urrutia (2004) in reporting as cross-sectional the number conditional on being an elderly household for $b$ and $wh$ (young for $\\hat{b}$), not cross-sectional over the whole model economy. \\\\ \n');
fprintf(FID, 'Restuccia \\& Urrutia (2004) explain calculation of intergenerational correlation of earnings at bottom of pg 1363. I assume the intergeneration correlations of (log) innate and acquired ability are calculated by the analagous regressions (with modification for acquired ability as is only observed for old). \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Figure 3
figure(3)
subplot(2,2,1); bar(TableColumn(1).FigureData.CollegeEnrollment)
title('Benchmark')
ylim([0,100])
ylabel('Percent of Quintile')
subplot(2,2,2); bar(TableColumn(18).FigureData.CollegeEnrollment)
title('Flat College Subsidy')
ylim([0,100])
subplot(2,2,3); bar(TableColumn(16).FigureData.CollegeEnrollment)
title('Increase Early Expenditure')
ylim([0,100])
ylabel('Percent of Quintile')
subplot(2,2,4); bar(TableColumn(17).FigureData.CollegeEnrollment)
title('Increase College Subsidy')
ylim([0,100])
saveas(gcf,'./SavedOutput/Graphs/RestucciaUrrutia2004_Figure3.png')

% Figure 4
figure(4)
subplot(2,2,1); bar(TableColumn(1).FigureData.CollegeDropoutRate)
title('Benchmark')
ylim([0,60])
ylabel('Percent of Quintile')
subplot(2,2,2); bar(TableColumn(18).FigureData.CollegeDropoutRate)
title('Flat College Subsidy')
ylim([0,60])
subplot(2,2,3); bar(TableColumn(16).FigureData.CollegeDropoutRate)
title('Increase Early Expenditure')
ylim([0,60])
ylabel('Percent of Quintile')
subplot(2,2,4); bar(TableColumn(17).FigureData.CollegeDropoutRate)
title('Increase College Subsidy')
ylim([0,60])
saveas(gcf,'./SavedOutput/Graphs/RestucciaUrrutia2004_Figure4.png')

% TableColumn.FigureData.CollegeEnrollment
% TableColumn.FigureData.CollegeCompletionRate