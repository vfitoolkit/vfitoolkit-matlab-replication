function Output=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a, T, ParamPath)
% Replicates the results of Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap
% For the Durable Goods model of Section 7.

%% To translate durable goods model in Section 7 of Guerrieri & Lorenzoni (2017) into the standard setup of VFI Toolkit I use following:
% d variables: n_it
% aprime variables: b_it+1, k_it+1
% a variables: b_it, k_it
% z variables: theta_it, e_it

% Keep the initial value of whatever the ParamPath parameter is
% There is only one, but am just going to allow for possibility of multiple anyway
ParamPathNames=fieldnames(ParamPath);
for ii=1:length(ParamPathNames) % Note: it is actually just length 1, but whatever
    initialvalforParamPathParams(ii)=Params.(ParamPathNames{ii});
end

%% Set some basic variables

% For durable goods model GL2017 reduce the grid size to 5 for the productivity shock on theta, and set tauchenq=1.
% Following is really just copy-paste from main script, just redoing with different n_theta value and tauchenq value. 
n_theta=6;
% Params.tauchenq=1; % Done in main script
[theta1_grid,pi_theta1]=discretizeAR1_Tauchen(0,Params.rho,sqrt(Params.sigmasq_epsilon),n_theta-1,Params.tauchenq);
z_grid=[0; exp(theta1_grid)];
pistar_theta1=ones(n_theta-1,1)/(n_theta-1);
for ii=1:10^4 % G&L2017, pg 1438 "when first employed, workers draw theta from its unconditional distribution"
    pistar_theta1=pi_theta1'*pistar_theta1; % There is a more efficient form to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
end
pi_z=[(1-Params.pi_ue), Params.pi_ue*pistar_theta1'; Params.pi_eu*ones(n_theta-1,1),(1-Params.pi_eu)*pi_theta1];
pi_z=pi_z./sum(pi_z,2);
pistar_z=ones(n_theta,1)/n_theta;
for ii=1:10^4 %  % There is a more efficient way to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_z=pi_z'*pistar_z; % Formula could be used to find stationary dist of the employment unemployment process, then just combine with stationary dist of theta1, which is already calculated
end
z_grid=z_grid/sum(z_grid.*pistar_z);

%% Grids
% Need to update the grids for d and a based on new number of grid points
d1_grid=linspace(0,1,n_d(1))'; % Labor supply

a1_grid=(Params.aupperbar-Params.alowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a(1))')-1)+Params.alowerbar; % assets

% Set grid for durable good holdings
Params.klowerbar=0;
a2_grid=(Params.kupperbar-Params.klowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a(2))')-1)+Params.klowerbar; % durable goods

a2range=Params.kupperbar-Params.klowerbar;
d2_grid=linspace(-a2range/6,a2range/6,n_d(2))'; % up to +- 1/6 of the range of a2

% Switch endogenous states into form required by VFI Toolkit
d_grid=[d1_grid; d2_grid];
a_grid=[a1_grid; a2_grid];
n_z=n_theta;

%% Durable goods as an experience asset
% To use an experienceasset, we need to define aprime(d,a) 
% [in notation of the current model, a2prime(d2,a2)]

vfoptions.experienceasset=1; % Using an experience asset
% Note: by default, assumes it is the last d variable that controls the
% evolution of the experience asset (and that the last a variable is
% the experience asset).

% aprimeFn gives the value of a2prime
vfoptions.aprimeFn=@(d2,a2) d2+a2;
% The first two inputs must be (d,a) [in the sense of aprime(d,a)], then any parameters

% We also need to tell simoptions about the experience asset
simoptions.experienceasset=1;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid; % Needed to handle aprimeFn 
simoptions.a_grid=a_grid; % Needed to handle aprimeFn


%%
FnsToEvaluate.Aprime = @(d1,d2,aprime,a,k,z) aprime; % Aggregate assets (which is this periods state)
GeneralEqmEqns.BondMarketClearence = @(Aprime,Bprime) Aprime-Bprime; %The requirement that the aggregate assets (lending and borrowing) equal zero

%% 
DiscountFactorParamNames={'beta'};

% Change the durable goods return function.
ReturnFn=@(d1, d2, bprime, b, k, z, r, alpha, gamma, psi, eta, phi_k, delta, zeta, chi, v, B, Bprime, omega) ...
    GuerrieriLorenzoni2017_DurableGoods_ReturnFn(d1, d2, bprime, b, k, z, r, alpha, gamma, psi, eta, phi_k, delta, zeta, chi, v, B, Bprime, omega);

%% TEST
Params.beta=0.8;
vfoptions.verbose=1;

% vfoptions=rmfield(vfoptions,'maxiter');

% vfoptions.divideandconquer=0; % for transition path, turn on divide-and-conquer
% vfoptions.maxhowards=0;
% tic;
% [V0,Policy0]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
% time0=toc;
% % THIS IS JUST AS A COMPARISON

vfoptions.divideandconquer=1; % for transition path, turn on divide-and-conquer
vfoptions.maxhowards=0;
% vfoptions.maxiter=2;
tic;
[V1,Policy1]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time1=toc;
% This works

vfoptions.divideandconquer=1; % for transition path, turn on divide-and-conquer
vfoptions.maxhowards=1;
% vfoptions.maxiter=2;
tic;
[V2,Policy2]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time2=toc;
% HOWARDS IS NOT YET WORKING

tic;
StationaryDist=StationaryDist_Case1(Policy1,n_d,n_a,n_z,pi_z, simoptions, Params);
disttime=toc;

tic;
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy1, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2,simoptions);
aggvarstime=toc;

max(abs(V1(:)-V2(:)))
max(abs(Policy1(:)-Policy2(:)))
[time1,time2]
disttime
aggvarstime

% Get same answer but Howards is broken



%% Solve the initial stationary equilibrium
vfoptions.divideandconquer=1; % Using an experience assetç
heteroagentoptions.verbose=1;

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};

fprintf('Calculating initial eqm (for durable goods) \n')
[p_eqm_initial, ~,MarketClearance_initial]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_initial.r;
[~,Policy_initial]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_initial.mat Params p_eqm_initial StationaryDist_initial AggVars_initial Policy_initial MarketClearance_initial

%% Final stationary equilibrium
% Only change is
% ParamPathNames=fieldnames(ParamPath); % Done above
for ii=1:length(ParamPathNames)  % Note: it is actually just length 1, but whatever
    temp=ParamPath.(ParamPathNames{ii});
    Params.(ParamPathNames{ii})=temp(end);
end

fprintf('Calculating final eqm (for durable goods) \n')
[p_eqm_final, ~,MarketClearance_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_final.r;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_final.mat Params p_eqm_final StationaryDist_final AggVars_final Policy_final MarketClearance_final

%%
% Free up space on gpu
clear StationaryDist_final

%% Compute the transition path

fprintf('Calculating general eqm transitiona (for durable goods) \n')

% We need to give an initial guess for the price path on interest rates
PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final.r, floor(T/2))'; p_eqm_final.r*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'

% Rewrite the aggregate variable to be next period bonds rather than current bonds as this is the actual 
% timing of the decision which the interest rate (r) effects
TransPathFnsToEvaluate.Aprime = @(d, aprime,kprime,a,k,z) aprime; % Aggregate assets decisions
% Rewrite the General Eqm conditions as rules for updating the price
TransPathGeneralEqmEqns.BondMarket = @(Aprime,Bprime) Aprime-Bprime; % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)

transpathoptions.GEnewprice=3; % 3 allows you to say which general eqm eqn is used to update which general eqm price/param, and how
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'BondMarket','r',0,0.01};... % first GE condition will be positive if r is too big (as this will lead to too much Aprime), so subtract
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.



PricePath=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  TransPathFnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames,transpathoptions,vfoptions,simoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_transition.mat PricePath

DurableGoodsFigFnsToEvaluate.output = @(d, aprime, kprime,a,k,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
DurableGoodsFigFnsToEvaluate.durablespurchases = @(d, aprime, kprime,a,k,z,delta) kprime-delta*k; % am not certain of this definition for "purchases"
DurableGoodsFigFnsToEvaluate.nondurablespurchases = @(d, aprime, kprime,a,k,z,r, delta, zeta, chi, v, B, Bprime)  GuerrieriLorenzoni2017_DurableGoods_ConsumptionFn(d, aprime, kprime,a,k,z,r, delta, zeta, chi, v, B, Bprime);

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, DurableGoodsFigFnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1(DurableGoodsFigFnsToEvaluate,PricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn,transpathoptions);

Output_pch=([AggVars_initial.output.Mean; AggVarsPath.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
durablespurchases_pch=([AggVars_initial.durablespurchases.Mean; AggVarsPath.durablespurchases.Mean]-AggVars_initial.durablespurchases.Mean)/AggVars_initial.durablespurchases.Mean;
nondurablespurchases_pch=([AggVars_initial.nondurablespurchases.Mean; AggVarsPath.nondurablespurchases.Mean]-AggVars_initial.nondurablespurchases.Mean)./AggVars_initial.nondurablespurchases.Mean;

% figure(figurenumber)
% % First parameter in ParamPath
% subplot(2,2,1); plot(0:1:T,[initialvalforParamPathParams(1); ParamPath.(ParamPathNames{1})])
% % interest rate
% subplot(2,2,2); plot(0:1:T,[p_eqm_initial.r; PricePath.r])
% % output
% subplot(2,2,3); plot(0:1:T,4*100*Output_pch)
% % durable and nondurable purchases
% subplot(2,2,4); plot(0:1:T,durablespurchases_pch,0:1:T,nondurablespurchases_pch)
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure',num2str(figurenumber),'.pdf'])
save(['./SavedOutput/GL2017_fig',num2str(figurenumber),'.mat'],'initialvalforParamPathParams','ParamPath','p_eqm_initial','PricePath','Output_pch','durablespurchases_pch','nondurablespurchases_pch','T','AggVars_initial','AggVarsPath')

Output.p_eqm_initial=p_eqm_initial;
Output.PricePath=PricePath;
Output.initialvalforParamPathParams=initialvalforParamPathParams;
Output.ParamPath=ParamPath;
Output.Output_pch=Output_pch;
Output.durablespurchases_pch=durablespurchases_pch;
Output.nondurablespurchases_pch=nondurablespurchases_pch;






