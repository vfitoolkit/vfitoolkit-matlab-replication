% This script is called by GuerrieriLorenzoni2017.m
% Following is largely just a copy of a chunk of GuerrieriLorenzoni2017.m file.

% Just use the grids, parameter values, etc that have already been set up by GuerrieriLorenzoni2017.m

Params.phi=Params.phi_initial;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqnParamNames(1).Names={'B'};
GeneralEqmEqn_1 = @(AggVars,p,B) AggVars(1)-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
GeneralEqmEqns={GeneralEqmEqn_1}; %, GeneralEqmEqn_2};

% % Following are already set by GuerrieriLorenzoni2017.m
% DiscountFactorParamNames={'beta'};
% 
% ReturnFn=@(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v,B,Bprime,omega) GuerrieriLorenzoni2017_ReturnFn(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v,B,Bprime,omega);
% ReturnFnParamNames={'r', 'gamma', 'psi', 'eta', 'phi', 'v','B','Bprime','omega'}; %It is important that these are in same order as they appear in 'GuerrieriLorenzoni2017_ReturnFn'

%% Solve the initial stationary equilibrium
% GL2007 misleadingly refer to this as the initial steady-state equilibrium, which it
% is not. It is the inital stationary equilibrium. (there are plenty of shocks at the idiosyncratic level, hence not steady-state which means the absence of shocks)

V0=ones(n_a,n_z,'gpuArray');
%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'}; %,'tau'

heteroagentoptions.verbose=1;
heteroagentoptions.pgrid=p_grid;

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
[p_eqm_initial,p_eqm_index_initial, MarketClearance_initial]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% findeqmtime=toc
Params.r=p_eqm_initial;

% Now that we know what the equilibrium price is, lets calculate a bunch of other things associated with the equilibrium
p_eqm_index_initial
disp('Calculating various equilibrium objects')
[~,Policy_initial]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
% AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

% Final stationary equilibrium
%  Only change is
Params.phi=Params.phi_final;

[p_eqm_final,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm_final;
[V_final,Policy_final]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

% Free up space on gpu
clear ConsumptionDecision
clear PolicyValues PolicyValuesPermute
clear StationaryDist_final

%% Compute the flexible price transition path
% We want to look at a one off unanticipated path of phi. ParamPath & PathParamNames are thus given by
ParamPath=Params.phi_final*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath(1:6)=temp(2:7); % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)
ParamPathNames={'phi'};

% We need to give an initial guess for the price path on interest rates
% PricePath0=[linspace(p_eqm_initial, p_eqm_final, floor(T/2))'; p_eqm_final*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0=[linspace(-0.01, p_eqm_final, floor(T/3))'; p_eqm_final*ones(T-floor(T/3),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0=p_eqm_final*ones(T,1);
PricePathNames_flex={'r'};

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

transpathoptions.weightscheme=1
transpathoptions.verbose=1
PricePath_flex=TransitionPath_Case1(PricePath0, PricePathNames_flex, ParamPath, ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);

%% Sticky Wages
if altcalib_figurenumber==9 || altcalib_figurenumber==11
    PricePath0=[PricePath_flex, zeros(T,1)];
    PricePathNames_NK={'r','omega'};
    
    % Rewrite the General Eqm conditions as rules for updating the price
    transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
    GeneralEqmEqnParamNames(1).Names={'Bprime'};
    % Only change to this is to enforce that will only try to decrease interest
    % rate when r>=0 is satisfied, never when it is less than zero. When r<0
    % instead as 0.01 (note that transition path codes will use 0.9*old+0.1*new anyway, so essentially adding 0.001, a tenth of one percentage point, to interest rate)
    GeneralEqmEqn_1 = @(AggVars,p,Bprime) (p(1)>=0)*(p(1)-0.1*(AggVars(1)-Bprime))+(p(1)<0)*(p(2)>0)*(p(1)+2*abs(p(1))+0.0001); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
    GeneralEqmEqnParamNames(2).Names={};
    GeneralEqmEqn_2 = @(AggVars,p,Bprime) (p(1)<0)*(p(2)+0.003)+(p(1)>=0)*(p(2)>0.003)*(p(2)-0.003); % If r>=0 then send omega towards zero (as in this case it returns max{0,omega-0.005}). r<0 then increase omega (in this case it returns omega+0.005). (Note: the choice of +0.03 in principle will be good or bad depending on the update rule for the transition path; given that weight on old price is 0.9 this will shift omega up by (1-0.9)*0.02 whenever the interest rate is negative)
    GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2};
    
    PricePath_NK=TransitionPath_Case1(PricePath0, PricePathNames_NK, ParamPath, ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,transpathoptions);
end
%% Figure for Alternative Calibration
AltCalibFnsToEvaluateParamNames(1).Names={};
AltCalibFnsToEvaluateFn_output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
AltCalibFnsToEvaluateParamNames(2).Names={};
AltCalibFnsToEvaluateFn_employment = @(d_val, aprime_val,a_val,z_val) d_val; %n_it in notation of GL2017
AltCalibFnsToEvaluate={AltCalibFnsToEvaluateFn_output, AltCalibFnsToEvaluateFn_employment};

% For the initial dist will need
Params.phi=Params.phi_initial;
Params.r=p_eqm_initial;

% Get the 
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, AltCalibFnsToEvaluate,Params, AltCalibFnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath_flex=EvalFnOnTransPath_AggVars_Case1(AltCalibFnsToEvaluate, AltCalibFnsToEvaluateParamNames,PricePath_flex,PricePathNames_flex, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

Output_pch_flex=([AggVars_initial(1); AggVarsPath_flex(:,1)]-AggVars_initial(1))/AggVars_initial(1);
Employment_pch_flex=([AggVars_initial(2); AggVarsPath_flex(:,2)]-AggVars_initial(2))/AggVars_initial(2);

if altcalib_figurenumber==9 || altcalib_figurenumber==11
    AggVarsPath_NK=EvalFnOnTransPath_AggVars_Case1(AltCalibFnsToEvaluate, AltCalibFnsToEvaluateParamNames,PricePath_NK,PricePathNames_NK, ParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);
    Output_pch_NK=([AggVars_initial(1); AggVarsPath_NK(:,1)]-AggVars_initial(1))/AggVars_initial(1);
    Employment_pch_NK=([AggVars_initial(2); AggVarsPath_NK(:,2)]-AggVars_initial(2))/AggVars_initial(2);
end

% For some of the alternative calibrations employment is plotted as well.
figure(altcalib_figurenumber)
if altcalib_figurenumber==9
    % interest rate
    subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_flex],0:1:T,4*100*[p_eqm_initial; PricePath_NK(:,1)])
    title('annual interest rate')
    % output
    subplot(1,2,2); plot(0:1:T,Output_pch_flex, 0:1:T,Output_pch_NK)
    title('output (pct deviation)')
    % ylabel('percent deviation from inital output in stationary eqm')
elseif altcalib_figurenumber==10
    % interest rate
    subplot(1,3,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_flex])
    title('annual interest rate')
    % output
    subplot(1,3,2); plot(0:1:T,Output_pch_flex)
    title('output (pct deviation)')
    % ylabel('percent deviation from inital output in stationary eqm')
    % employment
    subplot(1,3,3); plot(0:1:T,Employment_pch_flex,0:1:T)
    title('employment')
    % ylabel('percent deviation from inital employment in stationary eqm')
elseif altcalib_figurenumber==11
    % interest rate
    subplot(1,3,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_flex],0:1:T,4*100*[p_eqm_initial; PricePath_NK(:,1)])
    title('annual interest rate')
    % output
    subplot(1,3,2); plot(0:1:T,Output_pch_flex, 0:1:T,Output_pch_NK)
    title('output (pct deviation)')
    % ylabel('percent deviation from inital output in stationary eqm')
    % employment
    subplot(1,3,3); plot(0:1:T,Employment_pch_flex,0:1:T,Employment_pch_NK)
    title('employment')
    % ylabel('percent deviation from inital employment in stationary eqm')
elseif altcalib_figurenumber==12
    % interest rate
    subplot(1,2,1); plot(0:1:T,Full_PricePath_Flex_baseline, 0:1:T,4*100*[p_eqm_initial; PricePath_flex])
    title('annual interest rate')
    % output
    subplot(1,2,2); plot(0:1:T,Full_OutputPath_pch_Flex_baseline, 0:1:T,Output_pch_flex)
    title('output (pct deviation)')
    legend('baseline', 'higher risk aversion')
    % ylabel('percent deviation from inital output in stationary eqm')
end
saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.pdf'])





