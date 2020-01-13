function Output=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a,n_theta,n_k,n_p,d_grid,a_grid,z_grid,pi_z,pistar_z, ReduceTauchenqForDurables, T, ParamPath, ParamPathNames, vfoptions, simoptions, heteroagentoptions, transpathoptions, tauchenoptions)
% Replicates the results of Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap
% For the Durable Goods model of Section 7.


%% To translate durable goods model in Section 7 of Guerrieri & Lorenzoni (2017) into the standard setup of VFI Toolkit I use following:
% d variables: n_it
% aprime variables: b_it+1, k_it+1
% a variables: b_it, k_it
% z variables: theta_it, e_it

% Keep the initial value of whatever the ParamPath parameter is
% There is only one, but am just going to allow for possibility of multiple anyway
for ii=1:length(ParamPathNames)
    initialvalforParamPathParams(ii)=Params.(ParamPathNames{ii});
end

%% Set some basic variables

% For durable goods model GL2017 reduce the grid size to 5 for the productivity shock on theta.
if ReduceTauchenqForDurables==1
    % Following is really just copy-paste from main script, but now with different n_theta value. See there for comments/explanation.
    n_theta=6;
    [theta1_grid, pi_theta1]=TauchenMethod(0, Params.sigmasq_epsilon, Params.rho, n_theta-1, Params.tauchenq,tauchenoptions);
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
end

%% Grids
% Set grid for durable good holdings
Params.klowerbar=0;
k_grid=(Params.kupperbar-Params.klowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_k)')-1)+Params.klowerbar;

% Switch endogenous states into form required by VFI Toolkit
n_a=[n_a,n_k];
a_grid=[a_grid;k_grid];
n_z=n_theta;

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};
GeneralEqmEqnParamNames(1).Names={'Bprime'};
GeneralEqmEqn_1 = @(AggVars,p,Bprime) AggVars(1)-Bprime; %The requirement that the aggregate assets (lending and borrowing) equal zero
GeneralEqmEqns={GeneralEqmEqn_1};

%% 
DiscountFactorParamNames={'beta'};

% Change the durable goods return function.
ReturnFn=@(d_val, a1prime_val, a2prime_val, a1_val, a2_val, z_val,r, alpha, gamma, psi, eta, phi_k,delta, zeta, chi, v, B, Bprime,omega) GuerrieriLorenzoni2017_DurableGoods_ReturnFn(d_val, a1prime_val, a2prime_val, a1_val, a2_val, z_val,r, alpha, gamma, psi, eta, phi_k,delta, zeta, chi, v, B, Bprime, omega);
ReturnFnParamNames={'r', 'alpha', 'gamma', 'psi', 'eta', 'phi_k','delta', 'zeta', 'chi', 'v', 'B', 'Bprime','omega'}; %It is important that these are in same order as they appear in 'GuerrieriLorenzoni2017_ReturnFn'

%% Solve the initial stationary equilibrium

V0=ones([n_a,n_z],'gpuArray');
%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};

fprintf('Calculating initial eqm (for durable goods) \n')
[p_eqm_initial,p_eqm_index_initial, GeneralEqmCondns_initial]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_initial;
[~,Policy_initial]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_initial.mat Params p_eqm_initial p_eqm_index_initial GeneralEqmCondns_initial StationaryDist_initial AggVars_initial

%% Final stationary equilibrium
% Only change is
for ii=1:size(ParamsPath,2)
    Params.(ParamsPathParamNames{ii})=ParamsPath(end,ii);
end

fprintf('Calculating final eqm (for durable goods) \n')
[p_eqm_final,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_final;
[V_final,Policy_final]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_final.mat Params p_eqm_final p_eqm_index_final MarketClearance_final StationaryDist_final AggVars_final

%%
% Free up space on gpu
clear ConsumptionDecision
clear Policy_initial
clear PolicyValues PolicyValuesPermute
clear StationaryDist_final

%% Compute the transition path

fprintf('Calculating general eqm transitiona (for durable goods) \n')

% We need to give an initial guess for the price path on interest rates
PricePath0=[linspace(p_eqm_initial, p_eqm_final, floor(T/2))'; p_eqm_final*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
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
PricePath=TransitionPath_Case1(PricePath0, PricePathNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

save ./SavedOutput/GuerrieriLorenzoni2017_durablegoods_transition.mat PricePath

DurableGoodsFigFnsToEvaluateParamNames(1).Names={};
DurableGoodsFigFnsToEvaluateFn_output = @(d_val, a1prime_val, a2prime_val,a1_val,a2_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
DurableGoodsFigFnsToEvaluateParamNames(2).Names={};
DurableGoodsFigFnsToEvaluateFn_durablespurchases = @(d_val, a1prime_val, a2prime_val,a1_val,a2_val,z_val) a2prime_val; % am not certain of this definition for "purchases"
DurableGoodsFigFnsToEvaluateParamNames(3).Names={'r', 'delta', 'zeta', 'chi', 'v', 'B', 'Bprime'};
DurableGoodsFigFnsToEvaluateFn_nondurablespurchases = @(d_val, a1prime_val, a2prime_val,a1_val,a2_val,z_val,r, delta, zeta, chi, v, B, Bprime)  GuerrieriLorenzoni2017_DurableGoods_ConsumptionFn(d_val, a1prime_val, a2prime_val, a1_val, a2_val, z_val,r, delta, zeta, chi, v, B, Bprime);
DurableGoodsFigFnsToEvaluate={DurableGoodsFigFnsToEvaluateFn_output, DurableGoodsFigFnsToEvaluateFn_durablespurchases,DurableGoodsFigFnsToEvaluateFn_nondurablespurchases};

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, DurableGoodsFigFnsToEvaluate,Params, DurableGoodsFigFnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1(DurableGoodsFigFnsToEvaluate, DurableGoodsFigFnsToEvaluateParamNames,PricePath,PricePathNames, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

Output_pch=([AggVars_initial(1); AggVarsPath(:,1)]-AggVars_initial(1))/AggVars_initial(1);
durablespurchases_pch=([AggVars_initial(2); AggVarsPath(:,2)]-AggVars_initial(2))/AggVars_initial(2);
nondurablespurchases_pch=([AggVars_initial(3); AggVarsPath(:,3)]-AggVars_initial(3))/AggVars_initial(3);


figure(figurenumber)
% First parameter in ParamPath
subplot(2,2,1); plot(0:1:T,[initialvalforParamPathParams(1); ParamPath])
% interest rate
subplot(2,2,2); plot(0:1:T,[p_eqm_initial; PricePath])
% output
subplot(2,2,3); plot(0:1:T,4*100*Output_pch)
% durable and nondurable purchases
subplot(2,2,4); plot(0:1:T,durablespurchases_pch,0:1:T,nondurablespurchases_pch)
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure',num2str(figurenumber),'.pdf'])

Output.p_eqm_initial=p_eqm_initial;
Output.PricePath=PricePath;
Output.initialvalforParamPathParams=initialvalforParamPathParams;
Output.ParamPath=ParamPath;
Output.Output_pch=Output_pch;
Output.durablespurchases_pch=durablespurchases_pch;
Output.nondurablespurchases_pch=nondurablespurchases_pch;






