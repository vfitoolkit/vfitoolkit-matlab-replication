% Replicates the results of Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap

% QUESTIONS OUTSTANDING
% Do GL2017 change T for Fiscal Policy transitions? (I do, increase to T=50)
% Is my interpretation of the Fischer deflation experiment correct? (that B increases and then stays at this new higher level)
% Do GL2017 change tauchen q when reducing the number of grid points to 5
% in durable goods model? (probably, they appear to use Tauchen-Hussey, not Tauchen as reported in paper)
% Codes provided by GL2017 on website of Lorenzoni only cover baseline model for flexible and NK transitions. (so do not answer the above questions)


% A few lines that get 'NCores', the number of CPUs that can be used by codes when running on parallel CPUs.
% addpath(genpath('./MatlabToolkits/'))
parpool
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;

SkipInitialFinal=0
OverwriteWithGL2017Grid=1
OverwriteWithGL2017CodeParams=1

%% To translate Guerrieri & Lorenzoni (2017) into the standard setup of VFI Toolkit I use following:
% d variables: n_it
% aprime variables: b_it+1
% a variables: b_it
% z variables: theta_it, e_it

simoptions.parallel=2 % 4: Sparse matrix, but then put result on gpu
vfoptions.lowmemory=0
transpathoptions.lowmemory=1 % Essentially vfoptions.lowmemory=1 for the transition path.

%% Set some basic variables

n_d=41 
n_a=2^10 % Guerrieri & Lorenzoni (2017) use 200 points for agent distribution; VFI is done the same but they use the EGM (endo grid method) and then use ConesaKrueger style probabilistic weights to nearest grid point for agent dist simulation
n_theta=13; % Guerrieri & Lorenzoni (2017), pg 1438, states that they use a "12 state markov chain, following the approach in Tauchen (1986)". This description can be misinterpreted as theta is in fact the combination of a two-state markov on employed/unemployed with a 12-state Tauchen approx of AR(1) applied to employment. (This is clear from their codes) [The precise wording of GL2017 is correct, just easily misread.]
n_r=1001;

%Parameters (mostly from G&L2017, Table 1)
Params.beta=0.9711; %Model period is one-quarter of a year.

Params.gamma=4; % Coefficient of relative risk aversion
Params.eta=1.5; % Curvature of utility from leisure
Params.psi=12.48; % Coefficient on leisure in utility

Params.pi_eu=0.0573; % Transition to unemployment (the last three comes from their code, not paper)
Params.pi_ue=0.882; % Transition to employment
% Note: you cannot change pi_eu or pi_ue without changing the return
% function because the implied uncondtional probability of being unemployed
% is hard coded there.
Params.rho=0.967; % Persistence of productivity shock (G&L2017 call this rho)
Params.sigmasq_epsilon=0.017; % Variance of the shock to the log-AR(1) process on labour productivity.
Params.tauchenq=2.1; % I have reverse engineered this value from the grid of GL2017 (copy of their grid is included in their codes). They themselves reverse engineered choice of roughly 2.1 so that the variance of the resulting process (variance of 12-state markov logz) is as close as possible to what it should be (variance of true AR(1) logz).
% This choice of tauchenq is crucial to the results/replication of GL2017. It means
% that the min and max productivity shocks in the model, while they have the correct variance, have a range which
% is roughly just +-standard deviation. Because this range of shocks is so (empirically unrealistically) small
% lots of agents in the model end up near the borrowing constraint (which is not true using, e.g., tauchenq=3).
% Without lots of agents near the borrowing constraint the credit crisis (change in phi) has very little effect on the model.

Params.v=0.1; % Unemployment benefit
Params.B=1.6; % Bond supply
Params.Bprime=Params.B; % Bond supply is unchanging (for most of the paper); this is needed as it is part of government budget constraint that determines lump-sum tax tau_t. (Obviously B=Bprime must hold in any stationary general eqm.)
% For creating graphs later it will be useful to store both the initial and final values of phi
Params.phi_initial=0.959;
Params.phi_final=0.525;
Params.phi=Params.phi_initial; % Borrowing limit

% Params.r % General equilibrium interest rate (I use this instead of q)
% Params.tau % Lump-sum taxes, determined in general equilibrium (can implement it directly inside the ReturnFn)

Params.omega=0; % This is not actually needed for anything until we get to the 'Sticky Wage' model (Section 4 of GL2017, pg 1450)
% In the New Keynesian Sticky Wages model this appears as a wedge that
% increases the utility of leisure (intended to capture a fall in labor
% demand as a result of wages not falling). See GL2017 for explanation.


%% GL2017 codes param values
% The codes 'compute_steady_states.m' provided on Lorenzoni's website
% actually use the following calibrated values (contents of par.mat that
% comes with the codes).
if OverwriteWithGL2017CodeParams==1
    Params.beta=0.9774;
    Params.v=0.1670;
    Params.B=2.6712;
    Params.Bprime=Params.B;
    Params.psi=15.8819;
    % For creating graphs later it will be useful to store both the initial and final values of phi
    Params.phi_initial=1.6005;
    Params.phi_final=0.8767;
    Params.phi=Params.phi_initial; % Borrowing limit
end
% Note that the values of psi, B, phi_initial, and phi_final in particular are
% substantially different from those reported in Table 1 of GL2017. The
% differences for B, phi_initial and phi_final are simply because GL2017 do
% not actually ever report these, they report these divided by annual GDP
% (roughly 1.64; e.g. in paper phi_initial=1.6/1.64=0.925). However this is
% never made clear in the paper (the targets are described as being in these 'as fraction of
% annual GDP' terms, but there is no indication the actual reported
% parameter values are).


%% Some Toolkit options
simoptions.ncores=NCores; % Number of CPU cores
tauchenoptions.parallel=1

% Create markov process for the exogenous income (based on idea of employment and unemployment states, following Imrohoroglu, 1989).
[theta1_grid, pi_theta1]=TauchenMethod(0, Params.sigmasq_epsilon, Params.rho, n_theta-1, Params.tauchenq,tauchenoptions);
z_grid=[0; exp(theta1_grid)];
% G&L2017, pg 1438 "when first employed, workers draw theta from its unconditional distribution"; so here compute the unconditional distribution
pistar_theta1=ones(n_theta-1,1)/(n_theta-1);
for ii=1:10^4 % There is a more efficient form to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_theta1=pi_theta1'*pistar_theta1; 
end

% Floden & Linde (2001) report annual values for the AR(1) on log(z) as
% rho_FL=0.9136, sigmasq_epsilon_FL=0.0426; can calculate sigmasq_z_FL=sigmasq_epsilon_FL/(1-rho_FL^2)=0.2577;
% GL2017, footnote 11 give the formulae for annual rho and sigmasqz in terms of the quarterly:
rho=Params.rho;
sigmasq_epsilon=Params.sigmasq_epsilon;
sigmasq_epsilon_annual=(1/(4^2))*(4+6*rho+4*rho^2+2*rho^3)*(sigmasq_epsilon/(1-rho^2)); % Gives 0.2572 which is very close to FL
autocovariance_annual=(1/(4^2))*(rho+2*rho^2+3*rho^3+4*rho^4+3*rho^5+2*rho^6+rho^7)*(sigmasq_epsilon/(1-rho^2));
rho_annual=autocovariance_annual/sigmasq_epsilon_annual; % Gives 0.9127, which is close to FL. Note that this depends only on quarterly rho (the quarterly sigmasq_epsilon disappears when doing the division of annual autocovariance by annual variance)
% These check out, so the values in GL2017 for rho and sigmasq_epsilon are correct.

pi_z=[(1-Params.pi_ue), Params.pi_ue*pistar_theta1'; Params.pi_eu*ones(n_theta-1,1),(1-Params.pi_eu)*pi_theta1];
% Rows did not sum to one due to rounding errors at order of 10^(-11), fix this
pi_z=pi_z./sum(pi_z,2);
pistar_z=ones(n_theta,1)/n_theta;
for ii=1:10^4 %  % There is a more efficient way to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_z=pi_z'*pistar_z; % Formula could be used to find stationary dist of the employment unemployment process, then just combine with stationary dist of theta1, which is already calculated
end
% "The average level of theta is chosen so that yearly output in the initial steady state is normalized to 1"
z_grid=z_grid/sum(z_grid.*pistar_z);
% Double-check that this is 1
% sum(z_grid.*pistar_z)

% That the "normalized to 1" refers to E[theta] and not E[n*theta] is clear from setting
% v=0.1 to satisfy "For the unemployment benefit, we also follow Shimer
% (2005) and set it to 40% of average labor income." (pg 1438)
% Note, they do not actually ever normalize to 1 in the codes, GL2017 has E[theta]=1.0718

% Either their reported parameter values are slightly incorrect or the
% Tauchen method code they used is not accurate, as they get slightly
% different grid, so I will just load theirs directly
if OverwriteWithGL2017Grid==1
    load ./PaperMaterials/replication-codes-for-Credit-Crises-2017-1f7cb32/inc_process.mat
    theta = [0; exp(x)];
    z_grid=theta;
    
    tol_dist = 1e-10; % tolerance lvl for distribution
    S     = length(theta);
    fin   = 0.8820;        % job-finding probability
    sep   = 0.0573;        % separation probability
    % new transition matrix
    Pr = [1-fin, fin*pr; sep*ones(S-1, 1), (1-sep)*Pr];
    % find new invariate distribution
    pr  = [0, pr];
    dif = 1;
    while dif > tol_dist
      pri = pr*Pr;
      dif = max(abs(pri-pr));
      pr  = pri;
    end
    
    pi_z=Pr;
    pistar_z=pr';
    
    % Note: GL2017 codes do not do the normalization of z_grid, they just have
    % E[theta]=1.07 rather than 1 as stated in the paper.
end

%% Grids
% Set grid for asset holdings
Params.alowerbar=-1.25*Params.phi; % This seems reasonable (No-one can go below -Params.phi in any case). Note that Fischer deflation experiment (second last part of paper) won't work unless this is below -1.1*phi
Params.aupperbar=20; % Not clear exactly what value is appropriate, have gone with this, and checked that increasing it makes no difference.
a_grid=(Params.aupperbar-Params.alowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a)')-1)+Params.alowerbar;
% GL2017 codes use alowerbar of -2 and aupperbar of 50, but since no-one
% goes anywhere near 50 this seems excessive (even for them 0.9995 of population ends up below 20.8690). 
% Of their 200 points, only 183 are actually above -phi=-0.925 (only 171 are above final -phi=-0.525)
% Because GL2017 use a weighted probablility to approach to put agents on
% grid (policy fn allows off grid choices), and do not enforce the budget constraint on this, they in fact end
% up with 0.19% of the initial stationary distribution violating the budget constraint! (0.9981 mass is above phi1)
% For final distributions this becomes 0.96% violating budget constraint (0.9904 mass is above phi1)
% Of their 200 points, only 115 are between -phi1 and 20.86.
% The actual magnitude of the numerical error induced by this is small enough not to matter in their final results.

% Set grid for interest rate, r
r_grid=linspace(-0.02,0.04,n_r)'; % This grid is substantially wider than the actual likely equilibrium values and so is somewhat overkill.
% % Set grid for tax rate
% tau_grid=linspace(); % Can calculate this from the gov budget constraint

%Bring model into the notational conventions used by the toolkit
d_grid=linspace(0,1,n_d)'; % Labor supply
p_grid=r_grid;

n_z=n_theta;
n_p=n_r;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
%FnsToEvaluateParamNames(2).Names={'v'};
%FnsToEvaluateFn_2 = @(d_val, aprime_val,a_val,z_val,v) v*(z_val==0); % Total unemployment benefits
FnsToEvaluate={FnsToEvaluateFn_1}; %, FnsToEvaluateFn_2};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqnParamNames(1).Names={'B'};
GeneralEqmEqn_1 = @(AggVars,p,B) AggVars(1)-B; %The requirement that the aggregate assets (lending and borrowing; by government and private) equal zero
GeneralEqmEqns={GeneralEqmEqn_1};

%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v,B,Bprime,omega) GuerrieriLorenzoni2017_ReturnFn(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v,B,Bprime,omega);
ReturnFnParamNames={'r', 'gamma', 'psi', 'eta', 'phi', 'v','B','Bprime','omega'}; %It is important that these are in same order as they appear in 'GuerrieriLorenzoni2017_ReturnFn'

%% Solve the initial stationary equilibrium
% GL2007 misleadingly refer to this as the initial steady-state equilibrium, which it
% is not. It is the inital stationary equilibrium. (there are plenty of shocks at the idiosyncratic level, hence not steady-state which means the absence of shocks)

V0=ones(n_a,n_z,'gpuArray');
%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'}; %,'tau'

heteroagentoptions.verbose=1;
heteroagentoptions.pgrid=p_grid;
if SkipInitialFinal==0
    disp('Calculating price vector corresponding to the stationary eqm')
    % tic;
    [p_eqm_initial,p_eqm_index_initial, MarketClearance_initial]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    % findeqmtime=toc
    Params.r=p_eqm_initial;
    save ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params p_eqm_initial p_eqm_index_initial MarketClearance_initial n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat
end


%% Now that we know what the equilibrium price is, lets calculate a bunch of

if SkipInitialFinal==0
    % %other things associated with the equilibrium
    p_eqm_index_initial
    disp('Calculating various equilibrium objects')
    [~,Policy_initial]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
    
    % PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);
    
    StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
    
%     sum(sum(StationaryDist_initial))-1
%     max(max(abs(StationaryDist_initial-StationaryDist_initialOLD)))
%     sum(sum(abs(StationaryDist_initial-StationaryDist_initialOLD)))
%     max(abs(cumsum(sum(StationaryDist_initial,2))-cumsum(sum(StationaryDist_initialOLD,2))))
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
    
%     AggVars_initial
    
    % eqm_MC=real(GeneralEqmConditions_Case1(AggVars,Params.q, GeneralEqmEqns, Params, GeneralEqmEqnParamNames));
    
    save ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params p_eqm_initial p_eqm_index_initial MarketClearance_initial Policy_initial StationaryDist_initial AggVars_initial n_d n_a n_z
    % load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params p_eqm_initial p_eqm_index_initial MarketClearance_initial Policy_initial StationaryDist_initial
end
% Note, these things will have been loaded above anyway if skipping

% Some things you might want to take a look at just to see what is going on.
% p_eqm_initial
% AggVars_initial
% [Params.B, Params.Bprime]
% plot(p_grid, MarketClearance_initial)
% 
% Policy_initial(2,1000:end,11:13)
% plot(shiftdim(Policy_initial(2,:,:),1))
% plot(shiftdim(Policy_initial(2,:,:),1)-(1:1:n_a)'.*ones(n_a,n_z))
% plot(sum(StationaryDist_initial,2))
% plot(cumsum(sum(StationaryDist_initial,2)))

%% Table 1
FID = fopen('./SavedOutput/LatexInputs/GuerrieriLorenzoni2017_Table1.tex', 'w');
fprintf(FID, 'Parameters Values \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llll} \n \\hline \\hline \n');
fprintf(FID, 'Parameter & Explanation & Value & Target/Source \\\\ \n \\hline \n');
fprintf(FID, '$\\beta$ &  Discount Factor  & %8.4f  & Interest rate r=2.5\\%%  \\\\ \n', Params.beta);
fprintf(FID, '$\\gamma$ & Coefficient of relative risk aversion  & %i  &   \\\\ \n', Params.gamma);
fprintf(FID, '$\\eta$ & Curvature of utility of leisure  & %8.1f  & Average Frisch elasticity=1  \\\\ \n', Params.eta);
fprintf(FID, '$\\psi$ & Coefficient on leisure in utility & %8.2f  & Average hours worked 0.4 of endowment \\\\ \n', Params.psi);
fprintf(FID, '$\\rho$ & Persistence of productivity shock & %8.3f  & Persistence of wage process \\\\ \n', Params.rho);
fprintf(FID, '$\\sigma_{\\epsilon}$ & Variance of productivity shock & %8.3f  & Variance of wage process \\\\ \n', Params.sigmasq_epsilon);
fprintf(FID, '$\\pi_{e,u}$ & Transition to unemployment & %8.3f  & Shimer (2005)  \\\\ \n', Params.pi_eu);
fprintf(FID, '$\\pi_{u,e}$ & Transition to employment & %8.3f  & Shimer (2005)  \\\\ \n', Params.pi_ue);
fprintf(FID, '$\\upsilon$ & Unemployment benefit & %8.2f  & 40\\%% of average labor income \\\\ \n', Params.v);
fprintf(FID, '$B$ & Bond supply & %8.1f  & Liquid assets (flow of funds) \\\\ \n', Params.B);
fprintf(FID, '$phi$ & Borrowing limit & %8.3f  & Total gross debt (flow of funds) \\\\ \n', Params.phi);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: values of $B$ and $phi$ differ from those in Guerrieri \& Lorenzoni (2017). The actual model parameters are reported here, while those in paper are the parameter divided by annual output (annual output equals 4 times quarterly output; model is quarterly). \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Final stationary equilibrium
Params.phi=Params.phi_final;

if SkipInitialFinal==0
    [p_eqm_final,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    save ./SavedOutput/GuerrieriLorenzoni2017_final.mat Params p_eqm_final p_eqm_index_final MarketClearance_final
    
    Params.r=p_eqm_final;
    [V_final,Policy_final]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
    
    % PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);
    
    StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
    
    AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
    
    save ./SavedOutput/GuerrieriLorenzoni2017_final.mat Params p_eqm_final p_eqm_index_final MarketClearance_final V_final Policy_final StationaryDist_final AggVars_final n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_final.mat
end

%% Compute Annual GDP

% GL2017 describe Figure 1 as "Figure I shows the optimal values of
% consumption and labor supply as a function of the initial level of bond
% holdings". This is incorrect. The x-axis is not the level of bond holdings, it is
% bond-holdings as a fraction of annual output. I follow GL2017 in what
% they plot, adding the footnote to bottom of figure to explain what is
% actually graphed.
% In fact this is true of many of the Figures of GL2017, which label the
% x-axis as either b or B, but are in fact reporting b (or B) divided by
% annual output.
% To do this I need to compute annual GDP: quarterly output is y=theta*n.
% Following few lines do this (together with multiplication by 4 to make it annual)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % Output
FnsToEvaluateExtra={FnsToEvaluateFn_output};

QuarterlyOutput_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluateExtra,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
QuarterlyOutput_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluateExtra,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

AnnualOutput_initial=4*QuarterlyOutput_initial;
AnnualOutput_final=4*QuarterlyOutput_final;

%% Figure 1
figure(1)
l_a=length(n_a);
l_z=length(n_z);
Fig1FnsToEvaluateParamNames(1).Names={'r', 'v', 'B', 'Bprime'};
Fig1FnsToEvaluateFn_consumption = @(d_val, aprime_val,a_val,z_val,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d_val, aprime_val, a_val, z_val,r, v, B, Bprime); % Consumption
PolicyValues=PolicyInd2Val_Case1(Policy_initial,n_d,n_a,n_z,d_grid,a_grid, 2);
permuteindexes=[1+(1:1:(l_a+l_z)),1];
PolicyValuesPermute=permute(PolicyValues,permuteindexes);

% Note: 'EvalFnOnAgentDist_Grid_Case1' is not really indended for use by the end user. It
% is an internal function of the VFI Toolkit and hence it requires the
% parameters in the form of a vector, rather than working with the
% parameter structure and the names of the parameters.
ConsumptionDecision=EvalFnOnAgentDist_Grid_Case1(Fig1FnsToEvaluateFn_consumption,[Params.r,Params.v, Params.B, Params.Bprime],PolicyValuesPermute,n_d,n_a,n_z,a_grid,z_grid,2);
subplot(2,1,1); plot(a_grid,ConsumptionDecision(:,2),a_grid,ConsumptionDecision(:,8))
% legend('Mean','10th','25th','Median')
title({'Consumption'})
xlabel('Bond holdings')
% Labour supply is the d variable
subplot(2,1,2); plot(a_grid,d_grid(Policy_initial(1,:,2)),a_grid,d_grid(Policy_initial(1,:,8)))
legend('\theta^2','\theta^8');
title({'Labor Supply'})
xlabel('Bond holdings')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure1.pdf'])

% Now a version that reproduces GL2017 exactly
subplot(2,1,1); plot(a_grid/AnnualOutput_initial,ConsumptionDecision(:,2),a_grid/AnnualOutput_initial,ConsumptionDecision(:,8))
% legend('Mean','10th','25th','Median')
title({'Consumption'})
xlabel('Bond holdings as fraction of annual output')
ylabel('(Quarterly)')
% Labour supply is the d variable
subplot(2,1,2); plot(a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,2)),a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,8)))
legend('\theta^2','\theta^8');
title({'Labor Supply'})
xlabel('Bond holdings as fraction of annual output')
ylabel('(Quarterly)')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure1_GL2017.pdf'])


%% Figure 2
% This is a bit awkward to plot, but notice that MarketClearance is just
% measuring how large a change in B would be required to make that interest rate an equilibrium. 
figure(2)
% Need to change to annual interest rates (model period is quarterly)
annualinterestrate=100*((1+heteroagentoptions.pgrid).^4-1); % Note: GL2017 multiply by 4 rather than doing the explicit compounding.
plot(Params.B+MarketClearance_initial, annualinterestrate, Params.B+MarketClearance_final, annualinterestrate)
legend('initial stationary eqm','final stationary eqm');
% title({'Bond Market Equilibrum in Stationary Eqm'})
xlabel('aggregate bond holdings')
ylabel('(annual) interest rate')
xlim([0.5,3])
ylim([0,4])
xline(Params.B,'-','Bond supply') % Requires at least Matlab R2018B
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure2.pdf'])

% Now a version that reproduces GL2017 exactly
plot((Params.B+MarketClearance_initial)/AnnualOutput_initial, annualinterestrate, (Params.B+MarketClearance_final)/AnnualOutput_final, annualinterestrate)
legend('initial stationary eqm','final stationary eqm');
% title({'Bond Market Equilibrum in Stationary Eqm'})
xlabel('aggregate bond holdings')
ylabel('(annual) interest rate')
xlim([0.5,3])
ylim([0,4])
xline(Params.B/AnnualOutput_initial,'-','Bond supply') % Requires at least Matlab R2018B
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure2_GL2017.pdf'])

% I deliberately change the description to 'stationary equilibrium' instead
% of 'steady state'. I consider the later a slightly misleading
% description, albeit common in the literature, as the model contains
% plenty of shocks that are active (at agent level) and steady-state is used in 
% representative agent models to designate the equilibrium in the model with no 
% shocks (subtly, this is different to stationary equilibrium, the model in which shocks 
% exists, but there simply haven't been any for an infinite period of time). It is only the
% cross-section that is stationary, agents are moving around. But based on
% the existence of shocks in the model I consider stationary a more accurate description/definition.


%% Figure 4
figure(4)
% Unclear from GL2017 paper what is plotted in the top panel.
% From their codes it is clear that it is 'mean based on the unconditional
% distribution of the exogenous shock' (note that this is not the same
% as the 'mean based on the distribution of agents at that asset level').
subplot(2,1,1); plot(a_grid,sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid,a_grid,sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)
subplot(2,1,2); plot(a_grid,sum(StationaryDist_initial,2), a_grid,sum(StationaryDist_final,2))
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure4.pdf'])

% Now a version that reproduces GL2017 exactly
% Note: they are reporting x-axes that are 'as a fraction of initial annual
% output' for the final equilibrium results, not the final annual output.
% (Although in practice initial and final annual output are almost equal)
% For the top panel in fact everything (both axes) is as a fraction of
% initial annual output.
subplot(2,1,1); plot(a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial,a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial)
subplot(2,1,2); plot(a_grid/AnnualOutput_initial,sum(StationaryDist_initial,2), a_grid/AnnualOutput_initial,sum(StationaryDist_final,2))
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure4_GL2017.pdf'])


%%
% Free up space on gpu
clear ConsumptionDecision
clear Policy_initial
clear PolicyValues PolicyValuesPermute
clear StationaryDist_final

% load ./SavedOutput/GuerrieriLorenzoni2017_final.mat V_final
% load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat StationaryDist_initial

%% Compute the transition path
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=25

% We want to look at a one off unanticipated path of phi. ParamPath & PathParamNames are thus given by
ParamPath=Params.phi_final*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath(1:6)=temp(2:7); % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)
ParamPathNames={'phi'};

% We need to give an initial guess for the price path on interest rates
% PricePath0=[linspace(p_eqm_initial, p_eqm_final, floor(T/2))'; p_eqm_final*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0=[linspace(-0.01, p_eqm_final, floor(T/3))'; p_eqm_final*ones(T-floor(T/3),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0=p_eqm_final*ones(T,1);
PricePathNames={'r'};


% Rewrite the aggregate variable to be next period bonds rather than
% current bonds as this is the actual timing of the decision which the
% interest rate (r) effects
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) aprime_val; % Aggregate assets decisions
FnsToEvaluate={FnsToEvaluateFn_1};
% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqnParamNames(1).Names={'Bprime'};
GeneralEqmEqn_1 = @(AggVars,p,Bprime) p-0.1*(AggVars(1)-Bprime); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
GeneralEqmEqns={GeneralEqmEqn_1};

% [transpathoptions.GEnewprice=1 means that the GeneralEqmEqns should be
% expressed as how to generate a new guess for the price based on the
% current guess; transpathoptions.GEnewprice=0 means the GeneralEqmEqns
% should be expressed as for the standard general eqm conditions, namely
% equations that take the value of 0 in general eqm.]

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1
transpathoptions.verbose=1
PricePath=TransitionPath_Case1(PricePath0, PricePathNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat PricePath n_d n_a n_z

% For later we will keep another copy
PricePath_Flex=PricePath;

%% Figure 3
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat
load ./SavedOutput/GuerrieriLorenzoni2017_final.mat
load ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat 

Fig3FnsToEvaluateParamNames(1).Names={};
Fig3FnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig3FnsToEvaluateParamNames(2).Names={};
Fig3FnsToEvaluateFn_debt = @(d_val, aprime_val,a_val,z_val) -a_val*(a_val<0); % debt is (minus of) negative assets
Fig3FnsToEvaluate={Fig3FnsToEvaluateFn_output, Fig3FnsToEvaluateFn_debt};

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig3FnsToEvaluate,Params, Fig3FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1(Fig3FnsToEvaluate, Fig3FnsToEvaluateParamNames,PricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

Output_pch=([AggVars_initial(1); AggVarsPath(:,1)]-AggVars_initial(1))/AggVars_initial(1);

figure(3)
% Borrowing limit
subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath])
title('borrowing constraint')
% household debt-to-GDP ratio
subplot(2,2,2); plot(0:1:T,[AggVars_initial(2)./AggVars_initial(1); AggVarsPath(:,2)./AggVarsPath(:,1)])
title('household debt-to-annual-GDP ratio')
% interest rate
subplot(2,2,3); plot(0:1:T,100*(((1+[p_eqm_initial; PricePath]).^4)-1)) % converts to annual rate by compounding (not just multiplying by 4)
title('annual interest rate')
% output
subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
title('output')
% ylabel('percent deviation from inital output in stationary eqm')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure3.pdf'])

% Now a version that reproduces GL2017 exactly.
% Borrowing limit
subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath]./AnnualOutput_initial)
title('borrowing constraint as fraction-of-initial-annual-output')
% household debt-to-GDP ratio
subplot(2,2,2); plot(0:1:T,[AggVars_initial(2); AggVarsPath(1:end,2)]./AggVars_initial(1))
title('household debt-to-initial-annual-GDP ratio')
% interest rate
subplot(2,2,3); plot(0:1:T,100*4*[p_eqm_initial; PricePath]) % converts to annual rate by compounding (not just multiplying by 4)
title('annual interest rate')
% output
subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
% ylabel('percent deviation from inital output in stationary eqm')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure3_GL2017.pdf'])


%% Figure 5
Fig5FnsToEvaluateParamNames(1).Names={'r', 'v', 'B', 'Bprime'};
Fig5FnsToEvaluateFn_consumption = @(d_val, aprime_val,a_val,z_val,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d_val, aprime_val, a_val, z_val,r, v, B, Bprime); % Consumption
Fig5FnsToEvaluate={Fig5FnsToEvaluateFn_consumption};

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig5FnsToEvaluate,Params, Fig5FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, Fig5FnsToEvaluateParamNames,PricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);
% Only debt limit reduction path
UnchangedPricePath=p_eqm_initial*ones(T,1);
AggVarsPath_partial_onlydebtlimit=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, Fig5FnsToEvaluateParamNames,UnchangedPricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);
% Only interest rate path
UnchangedParamPath=Params.phi_initial*ones(T,1);
AggVarsPath_partial_onlyinterestrate=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, Fig5FnsToEvaluateParamNames,PricePath,PricePathNames, UnchangedParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

AggVarsPath_GE_pch=100*(AggVarsPath_GE-AggVars_initial)./AggVars_initial;
AggVarsPath_partial_onlydebtlimit_pch=100*(AggVarsPath_partial_onlydebtlimit-AggVars_initial)./AggVars_initial;
AggVarsPath_partial_onlyinterestrate_pch=100*(AggVarsPath_partial_onlyinterestrate-AggVars_initial)./AggVars_initial;

figure(5)
plot(0:1:T,[0;AggVarsPath_GE_pch], 0:1:T,[0;AggVarsPath_partial_onlydebtlimit_pch], 0:1:T,[0;AggVarsPath_partial_onlyinterestrate_pch])
title('Consumption Response Deviation')
legend('General eqm response','Partial eqm response to debt limit reduction','Partial eqm response to interest rate changes')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure5.pdf'])


%% Figure 6
% Unlike the other figures relating to the transition, which just require
% the distribution of agents at each time period, this figure requires
% following individuals along the transition path based on where they
% started. Hence will use SimPanelValues_TransPath_Case1() rather than
% EvalFnOnTransPath_AggVars_Case1(); actually the later could/should also be used
% here, but want to show the different options available as part of VFI Toolkit.

% For each of the four transition paths simulate 100 paths drawing from the relevant initial percentile, then take mean.
% (Alternative would be a much bigger panel based on drawing from the actual inital stationary distribution, then 
% just restrict panel based on percentiles and take means, but this would involve much more computation time)
simoptions.numbersims=100;

% We are going to want the values for consumption
Fig6FnsToEvaluateParamNames(1).Names={'r', 'v', 'B', 'Bprime'};
Fig6FnsToEvaluateFn_consumption = @(d_val, aprime_val,a_val,z_val,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d_val, aprime_val, a_val, z_val,r, v, B, Bprime); % Consumption
Fig6FnsToEvaluate={Fig6FnsToEvaluateFn_consumption};

% First, figure out the asset values that correspond to the percentiles
assetdist=cumsum(sum(StationaryDist_initial,2));
[~, prctileindexes]=min(abs(assetdist-(1:1:100)/100)); % prctilevalues_doublecheck should be approx equal to prctilevalues
% prctilevalues_doublecheck=assetdist(prctileindexes); % should give [0.01,0.02, etc up to 1]

% 1st percentile (note, the 1st percentile are going to be those closest to the borrowing constraint)
% Set up the appropriate inital distribution and simulate a panel data set (over the transition) from this.
InitialDist_1stpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_1stpercentile(prctileindexes(1),:)=StationaryDist_initial(prctileindexes(1),:)./sum(StationaryDist_initial(prctileindexes(1),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
% Everything else is just completely standard
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_1stpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, PricePathNames, ParamPath, ParamPathNames, T, V_final, InitialDist_1stpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
% The line in the figure is just the mean for each time period of these (I am guessing), but expressed as 
% percent deviation from steady state. [Not obvious if I should take mean and then percent deviation, or 
% take percent deviation and then mean; have gone with the former.]
SimPanelValues=shiftdim(SimPanelValues,1);
% Fig6_1stPercentileTrace=(mean(SimPanelValues,2)-mean(SimPanelValues(1,:)))/mean(SimPanelValues(1,:));
Fig6_1stPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial)/AggVars_initial;

% Now just repeat for 10th, 20th and 50th percentiles
InitialDist_10thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_10thpercentile(prctileindexes(10),:)=StationaryDist_initial(prctileindexes(10),:)./sum(StationaryDist_initial(prctileindexes(10),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_10thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, PricePathNames, ParamPath, ParamPathNames, T, V_final, InitialDist_10thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_10thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial)/AggVars_initial;
InitialDist_20thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_20thpercentile(prctileindexes(20),:)=StationaryDist_initial(prctileindexes(20),:)./sum(StationaryDist_initial(prctileindexes(20),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_20thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, PricePathNames, ParamPath, ParamPathNames, T, V_final, InitialDist_20thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_20thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial)/AggVars_initial;
InitialDist_50thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_50thpercentile(prctileindexes(50),:)=StationaryDist_initial(prctileindexes(50),:)./sum(StationaryDist_initial(prctileindexes(50),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_50thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, PricePathNames, ParamPath, ParamPathNames, T, V_final, InitialDist_50thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_50thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial)/AggVars_initial;

figure(6)
plot(0:1:T-1, [0; Fig6_1stPercentileTrace], 0:1:T-1, [0;Fig6_10thPercentileTrace], 0:1:T-1, [0;Fig6_20thPercentileTrace], 0:1:T-1, [0;Fig6_50thPercentileTrace])
title('Consumption Response by Percentile in Initial Wealth Distribution')
legend('1st percentile','10th percentile','20th percentile','50th percentile')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure6.pdf'])


%% Figure 7
Fig7FnsToEvaluateParamNames(1).Names={};
Fig7FnsToEvaluateFn_employment = @(d_val, aprime_val,a_val,z_val) d_val; %n_it in notation of GL2017
Fig7FnsToEvaluate={Fig7FnsToEvaluateFn_employment};

Employment_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig7FnsToEvaluate,Params, Fig7FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
EmploymentPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig7FnsToEvaluate, Fig7FnsToEvaluateParamNames,PricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

EmploymentPath_GE_pch=100*(EmploymentPath_GE-Employment_initial)./Employment_initial;

figure(7)
plot(0:1:T,[0; EmploymentPath_GE_pch])
title('Employment Response')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure7.pdf'])


%% Sticky Wages
% What GL2017 refer to as sticky wages actually enters the model as an increase in
% the utility of leisure. (They justify this as decreasing hours worked
% without decreasing the actual real-age-per-unit-of-time-worked which is
% intended to capture a decrease in labour demand by (unmodelled) firms.)

% Only two changes: first is that the return function needs to be modified for this.
% Second is general equilibrium conditions, these now need to enforce r>=0
% and determine omega>0 as that which is needed to acheive that r=0 as a general
% eqm outcome, whenever omega=0 would otherwise mean r<0.
% Then compute transition path (omega=0 in initial and final states, so these remain unchanged).
% The key is all in writing the general equilibrium conditions in right way
% so that omega=0, except when this would lead to r<0, and in these cases
% need to pick omega so that r=0.

% We need to give an initial guess for the price path on interest rates and
% omega. Lets just start with the flexible prices general eqm path for
% interest rates, and with omega=0 for all periods as our initial guess.
PricePath0=[PricePath, zeros(T,1)];
PricePathNames_NK={'r','omega'};

% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqnParamNames(1).Names={'Bprime'};
% Only change to this is to enforce that will only try to decrease interest rate when r>=0 is satisfied, never when it is less than zero. When r<0
% instead as 0.01 (note that transition path codes will use 0.9*old+0.1*new anyway, so essentially adding 0.001, a tenth of one percentage point, 
% to interest rate)
GeneralEqmEqn_1 = @(AggVars,p,Bprime) (p(1)>=0)*(p(1)-0.1*(AggVars(1)-Bprime))+(p(1)<0)*(p(2)>0)*(p(1)+2*abs(p(1))+0.0001); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
GeneralEqmEqnParamNames(2).Names={};
GeneralEqmEqn_2 = @(AggVars,p,Bprime) (p(1)<0)*(p(2)+0.003)+(p(1)>=0)*(p(2)>0.003)*(p(2)-0.003); % If r>=0 then send omega towards zero (as in this case it returns max{0,omega-0.005}). r<0 then increase omega (in this case it returns omega+0.005). (Note: the choice of +0.03 in principle will be good or bad depending on the update rule for the transition path; given that weight on old price is 0.9 this will shift omega up by (1-0.9)*0.02 whenever the interest rate is negative)
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2};
% Remark: This approach to how to update omega is different from that used
% by GL2017. They instead force the ZLB on the interest rate r, and then use
% omega to ensure that goods markets clear (Y=C; remember this is an
% endowment economy). omega is updated based on how far the model is from
% goods market clearance.
% Remark: The above updating rules took a few attempts to come up with (the r>=0
% constraint is a 'knife-edge' and so it was not trivial to find something
% that settles down rather than 'jump' back and forth between r<0 and r>0).

% [transpathoptions.GEnewprice=1 means that the GeneralEqmEqns should be
% expressed as how to generate a new guess for the price based on the
% current guess; transpathoptions.GEnewprice=0 means the GeneralEqmEqns
% should be expressed as for the standard general eqm conditions, namely
% equations that take the value of 0 in general eqm.]

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1
transpathoptions.verbose=1
transpathoptions.tolerance=10^(-4) % will run until r and omega settle to four digits
PricePath_NK=TransitionPath_Case1(PricePath0, PricePathNames_NK, ParamPath, ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

%% Figure 8 (I do an additional Figure 17 that shows the 'wedges' and employment)
Fig8FnsToEvaluateParamNames(1).Names={};
Fig8FnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig8FnsToEvaluateParamNames(2).Names={};
Fig8FnsToEvaluateFn_employment = @(d_val, aprime_val,a_val,z_val) d_val; %n_it in notation of GL2017
Fig8FnsToEvaluate={Fig8FnsToEvaluateFn_output, Fig8FnsToEvaluateFn_employment};

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig8FnsToEvaluate,Params, Fig8FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
% Difference between following two lines is PricePath vs PricePath_NK
AggVarsPath_Flex=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, Fig8FnsToEvaluateParamNames,PricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);
AggVarsPath_NK=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, Fig8FnsToEvaluateParamNames,PricePath_NK,PricePathNames_NK, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

OutputPath_pch_Flex=(AggVarsPath_Flex(:,1)-AggVars_initial(1))./AggVars_initial(1);
OutputPath_pch_NK=(AggVarsPath_NK(:,1)-AggVars_initial(1))./AggVars_initial(1);
% My extras
EmploymentPath_pch_Flex=(AggVarsPath_Flex(:,2)-AggVars_initial(2))./AggVars_initial(2); % This has already been calculated.
EmploymentPath_pch_NK=(AggVarsPath_NK(:,2)-AggVars_initial(2))./AggVars_initial(2);

figure(8)
% interest rate
subplot(2,1,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_Flex], 0:1:T,4*100*[p_eqm_initial; PricePath_NK(:,1)])
title('annual interest rate')
legend('flex price', 'NK fix price (ZLB)')
% output
subplot(2,1,2); plot(0:1:T,100*[0; OutputPath_pch_Flex],0:1:T,100*[0; OutputPath_pch_NK])
title('output')
% ylabel('percent deviation from inital output in stationary eqm')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure8.pdf'])

% Not actually part of GL2017, but let's take a look at the path for omega (the wedge on labour supply that 'enforces' the zero lower bound on interest rates)
figure(17)
% wedges
subplot(2,1,1); plot(0:1:T,[0; PricePath_NK(:,2)])
title('omega (NK wedge on wages)')
% employment
subplot(2,1,2); plot(0:1:T,100*[0; EmploymentPath_pch_Flex],0:1:T,100*[0; EmploymentPath_pch_NK])
title('employment')
legend('flex price', 'NK fix price (ZLB)')
% ylabel('percent deviation from inital output in stationary eqm')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_ExtraFigure1.pdf'])

%%
ParamPath_baseline=ParamPath;

%% Do the alternative calibrations. There are four.
% The relevant parameter values are all given in the online appendix of GL2017 (pg 4 and 5).
% Note that some of the parameter values being set below are actually just the same as their baseline values.

% All of them graph subset of the interest rate, output, employment
% For each, need to calculate initial and final general eqm, then transition path for flex and NK
% The script 'GuerrieriLorezoni2017_altcalibration' is mostly just a copy
% paste of the above that generates the relevant outputs for whatever
% calibration is currently implicit in Ṕarams.

AnnualOutput_initial=gather(AnnualOutput_initial); % Note: model is quarterly, so this is four times the model period output

% Need to plot the baseline for some
Full_OutputPath_pch_Flex_baseline=100*[0; OutputPath_pch_Flex];
Full_PricePath_Flex_baseline=4*100*[p_eqm_initial; PricePath_Flex];

% Is not quite clear from GL2017 if the following are exactly what they do,
% in the sense that I just use the existing AnnualOutput_initial to figure
% out the parameter values for B and phi from the ones they report in the
% online appendix. They probably use the AnnualOutput_initial that would be
% relevant in each of these alternative calibrations. This could be backed
% out from use of their codes; if anyone wants to do so and send me the
% 'corrected' values I would be grateful, but I cannot be bothered right
% now.

% Figures 9 and 11 will 'fail' to replicate. This is 'intentional'. They
% calculate transition paths under a zero lower bound on interest rates in
% which the final equilibrium has a negative interest rate.
% Mathematically this can be accomadated using a wedge, as GL2017 do for
% the all the zero lower bound transitions, in the final equilibrium
% definition as well, but the economic justification underpinning this
% would become that fixed wages can never fall, even after an infinite
% amount of time passes. If you want to mathematically replicate these just
% add the wedge 'omega' and the 'second' general eqm condition (from
% transition paths) to the codes that determine the final eqm. I choose
% instead to 'fail' to replicate to highlight that these alternative calibrations
% actually involve modifications to the definition of the final eqm (or
% equivalently to the definition of an NK (fixed wage) transition path.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

%% First alternative: Median Wealth Calibration
Params.beta=0.9615;
Params.eta=1.5;
Params.psi=19.03;
Params.v=0.1;
Params.B=0.5241*AnnualOutput_initial;
Params.Bprime=Params.B;
Params.phi_initial=0.6031*AnnualOutput_initial;

altcalib_figurenumber=9;

GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

%% Second alternative: Low phi Calibration
Params.beta=0.9755;
Params.eta=0.25;
Params.psi=1.389;
Params.v=0.1;
Params.B=1.6*AnnualOutput_initial;
Params.Bprime=Params.B;
Params.phi_initial=0.7077*AnnualOutput_initial;

altcalib_figurenumber=10;

GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

%% Third alternative: Median wealth, low phi Calibration
Params.beta=0.9544;
Params.eta=0.25;
Params.psi=1.427;
Params.v=0.1;
Params.B=0.429*AnnualOutput_initial;
Params.Bprime=Params.B;
Params.phi_initial=0.7184*AnnualOutput_initial;

altcalib_figurenumber=11;

GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

%% Fourth alternative: Calibration with gamma=6
% This one is kind of pointless. It turns up the risk aversion, but because
% Tauchen q=2 (hyperparameter when doing Tauchen approx of the AR(1) log-income process)
% the min and max values of the log-income process are only +/-1 standard
% deviation anyway. So model economy has had all the actual risk removed
% from it before starting.
Params.beta=0.9756;
Params.eta=1.5;
Params.psi=101.22;
Params.v=0.1;
Params.B=1.6*AnnualOutput_initial;
Params.Bprime=Params.B;
Params.phi_initial=1.1*AnnualOutput_initial;
% And presumably also (but not explicitly mentioned in GL2017 appendix)
Params.gamma=6;

altcalib_figurenumber=12;

GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

%%
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat StationaryDist_initial

%% Fiscal Policy
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
% Pg 1457 of GL2017: since taxes are lump-sum it is in principle
% possible to fully offset the tightening of the borrowing limit be
% transfering the 'same amount' to agents.
% The experiment that will be performed here is to look at policies that
% only partially offest the long run change in borrowing limit phi.

% Note: After coding the PathForB below it becomes clear that T=25 will not
% be enough to reach a new stationary recursive equilibrium (B_t has not yet reached Bprime). So have
% increased to T=50.
T=50;
% Not perfect, but B_t now reaches very close to Bprime and is almost constant.
% Based on Figures, it appears GL2017 kept T=25.

% Two transfer policies are assessed, the first where government
% temporarily reduces the lump-sum tax tau for all households, and the
% second where the government temporarily raises the unemployment benefit.

% Both transfer policies assume a fixed sequence of
% government deficits implied by the bond supply following the path,
%   B_t = rho_b^t Binitial + (1-rho_b^t)*Bfinal   
% where Bfinal is the new long-run level of B.
% Set Bfinal=1.2*B, and rho_b=0.95.
Params.rho_b=0.95;
Params.Bfinal=1.2*Params.B;
ParamPathForB=(Params.rho_b.^(0:1:T-1)*Params.B+(1-Params.rho_b.^(0:1:T-1))*Params.Bfinal)';
% Because of how the current lump-sum tax rate 'tau' is calculated we need
% B_t over the whole path. We will actually also need B_tplus1 as well as
% this is key to the general equilibrium condition. (See above discussion
% when the general equilibrium condition was first described and introduced
% when computing the initial general equilibrium.)
ParamPathForBprime=[ParamPathForB(2:end); Params.rho_b^(T)*Params.B+(1-Params.rho_b^(T))*Params.Bfinal];
% Force that Bprime=B in the final period.
ParamPathForBprime(end)=ParamPathForB(end);

% Before we get into the transition paths themselves we need to calculate
% the final general equilibrium.
Params.phi=Params.phi_final;
% Params.Binitial=Params.B; % So can restore after calculation. % Just load the inital params again instead; I do this roughly 10-15 lines below
Params.B=ParamPathForBprime(end); % Use ParamPathForBprime(end) rather than Params.Bfinal
Params.Bprime=Params.B;
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};
GeneralEqmEqnParamNames(1).Names={'B'};
GeneralEqmEqn_1 = @(AggVars,p,B) AggVars(1)-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
GeneralEqmEqns={GeneralEqmEqn_1}; %, GeneralEqmEqn_2};
[p_eqm_final_Fiscal,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_final_Fiscal;
[V_final_Fiscal,~]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
% Params.phi=Params.phi_initial;  
% Params.B=Params.Binitial;
% Params.Bprime=Params.B;
% Done, now back to the transition path.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

% The two transfer policies are then different ways of spending the associated deficit implied by the bond supply path.

% First, we consider a policy where the tax tau_t adjusts to balance government budget to 
% generate the relevant deficit (period-by-period). tau_t is determined inside the return function, 
% and is already determined in this way.

% The second, transition considers a policy in which unemployment benefits are increased by 50% 
% for the first two years after the shock. Afterward the unemployment
% benefit reverts to its initial value and throughout the transition path
% tau_t adjusts to satisfy the government budget constraint.
% This will simply involve adding the unemployment benefit path for v to
% ParamPath, and will be done below after we finish the first case.

% Set the new path for parameters to be modelled. In addition to the same
% path on phi as before, we now have a path for B.
ParamPath_Fiscal1=[Params.phi_final*ones(T,1),ParamPathForB,ParamPathForBprime]; % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath_Fiscal1(1:6,1)=temp(2:7)'; % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)
ParamPathNames_Fiscal1={'phi','B','Bprime'};
% The rest of the setup is just same thing as the baseline case with
% flexible prices. So most of the following lines are just copy-paste from there.

% We need to give an initial guess for the price path on interest rates
PricePath0_Fiscal=[linspace(p_eqm_initial, p_eqm_final_Fiscal, floor(T/2))'; p_eqm_final_Fiscal*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0_Fiscal=[linspace(-0.01, p_eqm_final_Fiscal, floor(T/3))'; p_eqm_final_Fiscal*ones(T-floor(T/3),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0=p_eqm_final*ones(T,1);
PricePathNames_Fiscal={'r'}; % Note: tau is also being determined, but this is done inside the ReturnFn

% Rewrite the aggregate variable to be next period bonds rather than
% current bonds as this is the actual timing of the decision which the
% interest rate (r) effects
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) aprime_val; % Aggregate assets decisions
FnsToEvaluate={FnsToEvaluateFn_1};
% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqnParamNames(1).Names={'Bprime'};
GeneralEqmEqn_1 = @(AggVars,p,Bprime) p-0.1*(AggVars(1)-Bprime); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
GeneralEqmEqns={GeneralEqmEqn_1};

% [transpathoptions.GEnewprice=1 means that the GeneralEqmEqns should be
% expressed as how to generate a new guess for the price based on the
% current guess; transpathoptions.GEnewprice=0 means the GeneralEqmEqns
% should be expressed as for the standard general eqm conditions, namely
% equations that take the value of 0 in general eqm.]

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1
transpathoptions.verbose=1
PricePath_Fiscal1=TransitionPath_Case1(PricePath0_Fiscal, PricePathNames_Fiscal, ParamPath_Fiscal1, ParamPathNames_Fiscal1, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fiscalpolicy.mat PricePath_Fiscal1

% Now for the second transition.
% The second, transition considers a policy in which unemployment benefits are increased by 50% 
% for the first two years after the shock. Afterward the unemployment
% benefit reverts to its initial value and throughout the transition path
% tau_t adjusts to satisfy the government budget constraint.
% This will simply involve adding the unemployment benefit path for v to
% ParamPath, and will be done below after we finish the first case.

% Seems to be having problems converging because it goes wrong near end of transition path during first few iterations and never recovers.
% transpathoptions.weightscheme=4 % 4: a gradually opening window, combined with expontially decreasing weighting of update
% Only change is to ParamPath,
% ParamPath_Fiscal2=[ParamPath_Fiscal1(:,1),Params.v*[1*ones(8,1); ones(T-8,1)],ParamPath_Fiscal1(:,2:3)];
% ParamPathNames_Fiscal2={'phi','v','B','Bprime'};
ParamPath_Fiscal2=[ParamPath_Fiscal1,Params.v*[1.5*ones(8,1); ones(T-8,1)]];
ParamPathNames_Fiscal2={'phi','B','Bprime','v'};
% Final equilibrium is same as for the other Fiscal Policy transition path.
PricePath_Fiscal2=TransitionPath_Case1(PricePath0_Fiscal, PricePathNames_Fiscal, ParamPath_Fiscal2, ParamPathNames_Fiscal2, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fiscalpolicy.mat PricePath_Fiscal1 PricePath_Fiscal2

% Now we have done the Fiscal Policy, just need to create the relevant Figure.
FiscPolFnsToEvaluateParamNames(1).Names={};
FiscPolFnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
FiscPolFnsToEvaluate={FiscPolFnsToEvaluateFn_output};

% Note that Params is currently all the initial values.

% Get the 
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FiscPolFnsToEvaluate,Params, FiscPolFnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath_Fiscal1=EvalFnOnTransPath_AggVars_Case1(FiscPolFnsToEvaluate, FiscPolFnsToEvaluateParamNames,PricePath_Fiscal1,PricePathNames_Fiscal, ParamPath_Fiscal1, ParamPathNames_Fiscal1, Params, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);
AggVarsPath_Fiscal2=EvalFnOnTransPath_AggVars_Case1(FiscPolFnsToEvaluate, FiscPolFnsToEvaluateParamNames,PricePath_Fiscal2,PricePathNames_Fiscal, ParamPath_Fiscal2, ParamPathNames_Fiscal2, Params, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

Output_pch_Fiscal1=([AggVars_initial(1); AggVarsPath_Fiscal1(:,1)]-AggVars_initial(1))/AggVars_initial(1);
Output_pch_Fiscal2=([AggVars_initial(1); AggVarsPath_Fiscal2(:,1)]-AggVars_initial(1))/AggVars_initial(1);

% Figure 13
figure(13)
% interest rate
subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_Flex; ones(T-length(PricePath_Flex),1)*PricePath_Flex(end)],0:1:T,4*100*[p_eqm_initial; PricePath_Fiscal1],0:1:T,4*100*[p_eqm_initial; PricePath_Fiscal2])
title('annual interest rate')
% output
subplot(1,2,2); plot(0:1:T,100*[0; OutputPath_pch_Flex; ones(T-length(OutputPath_pch_Flex),1)*OutputPath_pch_Flex(end)], 0:1:T,100*Output_pch_Fiscal1, 0:1:T,100*Output_pch_Fiscal2)
title('output')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure13.pdf'])


%% Fischer Debt Deflation
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
% Page 1459, "A simple experiment in which we assume that all debt (private and public)
% is in nominal terms and there is an unexpected 10% reduction in the price
% level at date t=0. After that, the debt level is constant."
% Not entirely clear but this seems to mean they simply adjust the initial
% stationary distribution of agents, and change 'B'. Note that it seems
% like this means that B will remain fixed at this new level as Figure XIV does not show the 
% economy returning to its initial steady state, nor the baseline final steady state but instead 
% to a new higher level, but this is not clear from the text alone.

% Inflation factor
InflationFactor=0.9; % An unexpected 10% reduction in the price level (=1 means no change; 1.1 would be a 10% inflation, 0.9 would be a 10% deflation)

% Start with the public debt.
Params.B=(1/InflationFactor)*Params.B;
Params.Bprime=(1/InflationFactor)*Params.Bprime;
% Now the private debt, figure out where on a_grid people get shifted to.
% (Note, the way I do this will only work if no-one is in the very ends of the grid (deflation would push them off the grid, but this is not possible so they would instead be 'distorted' back onto the grid).)
newa_grid=(1/InflationFactor)*a_grid; % (absolute value of) debts decrease in real terms, savings decrease in real terms. [Opposite if it is a deflation, as in this example.]
[~,indextoshiftto]=min(abs(a_grid-newa_grid'));
checkthislookslikeagrid=a_grid(indextoshiftto); % It does
newStationaryDist_initial=zeros(size(StationaryDist_initial),'gpuArray');
for ii=1:length(a_grid)
    newStationaryDist_initial(indextoshiftto(ii),:)=newStationaryDist_initial(indextoshiftto(ii),:)+StationaryDist_initial(ii,:);
end
% The rest is just to recalculate the final distribution (based on new B & Bprime)
Params.phi=Params.phi_final;
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};
GeneralEqmEqnParamNames(1).Names={'B'};
GeneralEqmEqn_1 = @(AggVars,p,B) AggVars(1)-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
GeneralEqmEqns={GeneralEqmEqn_1}; %, GeneralEqmEqn_2};
[p_eqm_final_Fisher,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_final_Fisher;
[V_final,~]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.phi=Params.phi_initial;
% Then the transition path. The param path is trivial, as it is just the
% new parameters, together with the usual path on borrowing limit (phi).
ParamPath=[Params.phi_final*ones(T,1),[Params.B,Params.Bprime].*ones(T,2)];
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath(1:6,1)=temp(2:7)'; % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there

ParamPathNames={'phi','B','Bprime'};
PricePath0=p_eqm_final_Fisher*ones(T,1);
PricePathNames={'r'};
% Rewrite the aggregate variable to be next period bonds rather than
% current bonds as this is the actual timing of the decision which the
% interest rate (r) effects
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) aprime_val; % Aggregate assets decisions
FnsToEvaluate={FnsToEvaluateFn_1};
% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqnParamNames(1).Names={'Bprime'};
GeneralEqmEqn_1 = @(AggVars,p,Bprime) p-0.1*(AggVars(1)-Bprime); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
GeneralEqmEqns={GeneralEqmEqn_1};

% Now just run the TransitionPath_Case1 command
PricePath_FisherDef=TransitionPath_Case1(PricePath0, PricePathNames, ParamPath, ParamPathNames, T, V_final, newStationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fischerdeflation.mat PricePath_FischerDeflation V_final p_eqm_final

% Now we have done the Fischer Deflation, just need to create the relevant Figure.
FisherDefFnsToEvaluateParamNames(1).Names={};
FisherDefFnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
FisherDefFnsToEvaluate={FisherDefFnsToEvaluateFn_output};

% Already have the AggVars_initial from previously, so just use those.
AggVarsPath_FisherDef=EvalFnOnTransPath_AggVars_Case1(FisherDefFnsToEvaluate, FisherDefFnsToEvaluateParamNames,PricePath_FisherDef,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

Output_pch_FisherDef=([AggVars_initial(1); AggVarsPath_FisherDef(:,1)]-AggVars_initial(1))/AggVars_initial(1);

% Figure 14
figure(14)
% interest rate
subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_Flex;ones(T-length(PricePath_Flex),1)*PricePath_Flex(end)], 0:1:T,4*100*[p_eqm_initial; PricePath_FisherDef])
legend('Baseline','Fisher Deflation')
title('annual interest rate')
% output
subplot(1,2,2); plot(0:1:T,100*[0; OutputPath_pch_Flex; ones(T-length(OutputPath_pch_Flex),1)*OutputPath_pch_Flex(end)], 0:1:T,100*Output_pch_FisherDef)
title('output')
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure14.pdf'])

%% Durable Goods
% Extend the model to introduce durable goods. These will be modelled as a
% second asset. They are less liquid (captured in model by an adjustment
% cost) and can be used as collateral. Consumers will want durable goods as 
% the durable goods are simply added directly into the utility function.
% The extended model all includes a spread between borrowing and lending
% interest rates.
% Notice that because there is no supply side to the durable goods (the
% supply of durable goods is perfectly elastic) the addition of durable
% goods does not involve adding any complication whatsoever to the general
% equilibrium aspects of the model.

load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

% Based on graphs it appears GL2017 use T=50.
T=50;

% Because this involves solving a whole other model with different
% dimensions, return function, etc., I have done it in a seperate function
% that just returns the few things that are needed for the graphs.
% As inputs it takes the parameters, and the grid sizes for those parts of 
% the model that are not changed.

% The following lines are the options relating to the durable goods model
% that can be set, and which are inputs to the function.
n_k=2^7; % Number of grid points for the durable asset.
Params.kupperbar=10; % (max value for grid, the min is hard-coded as zero)

ReduceTauchenqForDurables=1; % Reduces Tauchen approximation of theta to 5-states. This is done by GL2017 to ease computation (footnote 40 of online appendix).

%Parameters that are changed or newly introduced for the durable goods model (G&L2017, Table 1)
Params.beta=0.9713; %Model period is one-quarter of a year.
Params.psi=2.54; % Coefficient on leisure in utility
Params.alpha=0.7; % Coefficient on non-durables (Newly introduced) 
Params.delta=0.0129; % Durables depreciation rate (Newly introduced)
Params.zeta=0.15; % Proportional loss on durable sales (Newly introduced)
Params.chi=0.01; % Intermediation cost (Newly introduced)
Params.v=0.16; % Unemployment benefit
Params.phi_k=0.8; % Borrowing limit (Replaces phi: is now a collateral constraint on borrowing against durable goods)
% NOTE: These come from Table A.4 of online appendix. It contains footnote
% "The quantities v and B are expressed in terms of yearly aggregate
% output." No such footnote appears on Table 1 which gives their
% (identical) values in baseline calibration. It is unclear what to make of
% this (codes appear to treat v & B values as those for model period,
% namely quarterly).

% Do two experiments.
% The first it to change phi_k to 0.56 gradually over a linear path that lasts six quarters.
% The second is a transitory increase in chi of 0.06 (to 0.07), this is done
% instantly and then the transitory increase decays back to zero at a rate of decay of 0.6.
ParamPath_phi_k=0.56*ones(T,1);
ParamPath_phi_k(1:6)=linspace(Params.phi_k,0.56,6)';
ParamPath_chi=(Params.chi+0.06.^(0:1:T-1))';

ParamPath_Durables1=ParamPath_phi_k;
ParamPath_Durables2=ParamPath_chi;
ParamPathNames_Durables1={'phi_k'};
ParamPathNames_Durables2={'chi'};

% We need to give an initial guess for the price path on interest rates.
% This is done inside the durable goods function once initial and final
% stationary general eqm interest rates are known. (Both paths have same
% initial general eqm, but I am going to be lazy and just calculate it
% twice, this allows the durable goods fn to be more generalized.)

% The general equilibrium conditions are unchanged. These are declared
% inside the durable good function.
figurenumber=15;
Output_DurableGoods1=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a,n_theta,n_k,n_p,d_grid,a_grid,z_grid,pi_z,pistar_z, ReduceTauchenqForDurables, T, ParamPath_Durables1, ParamPathNames_Durables1, vfoptions, simoptions, heteroagentoptions, transpathoptions, tauchenoptions);
figurenumber=16;
Output_DurableGoods2=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a,n_theta,n_k,n_p,d_grid,a_grid,z_grid,pi_z,pistar_z, ReduceTauchenqForDurables, T, ParamPath_Durables2, ParamPathNames_Durables2, vfoptions, simoptions, heteroagentoptions, transpathoptions, tauchenoptions);




