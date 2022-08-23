% This script is called by GuerrieriLorenzoni2017.m
% Following is largely just a copy of a chunk of GuerrieriLorenzoni2017.m file.

% Just use the grids, parameter values, etc that have already been set up by GuerrieriLorenzoni2017.m

Params.phi=Params.phi_initial;

%% Solve the initial stationary equilibrium
% GL2007 misleadingly refer to this as the initial steady-state equilibrium, which it
% is not. It is the inital stationary equilibrium. (there are plenty of shocks at the idiosyncratic level, hence not steady-state which means the absence of shocks)

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'}; %,'tau'

heteroagentoptions.verbose=1;
heteroagentoptions.pgrid=p_grid;

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% [p_eqm_initial,p_eqm_index_initial, MarketClearance_initial]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);
[p_eqm_initial,p_eqm_index_initial, MarketClearance_initial]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);
% findeqmtime=toc
Params.r=p_eqm_initial.r;

% Now that we know what the equilibrium price is, lets calculate a bunch of other things associated with the equilibrium
p_eqm_index_initial
disp('Calculating various equilibrium objects')
[~,Policy_initial]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
% AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% Final stationary equilibrium
%  Only change is
Params.phi=Params.phi_final;
% [p_eqm_final,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);
[p_eqm_final,p_eqm_index_final, MarketClearance_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);
Params.r=p_eqm_final.r;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% Free up space on gpu
clear ConsumptionDecision Policy_final
clear PolicyValues PolicyValuesPermute
clear StationaryDist_final

%% Compute the flexible price transition path
clear PricePath0 ParamPath
% We want to look at a one off unanticipated path of phi. ParamPath & PathParamNames are thus given by
ParamPath.phi=Params.phi_final*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath.phi(1:6)=temp(2:7); % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
PricePath0.r=[linspace(-0.01, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % PricePath0 is matrix of size T-by-'number of prices'

fprintf('Starting flex prices transition in alt calibration \n')
whos
PricePath_flex=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  TransPathFnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions);

%% Sticky Wages
if altcalib_figurenumber==9 || altcalib_figurenumber==11
    PricePath0.r=PricePath_flex.r;
    PricePath0.omega=zeros(T,1);

    fprintf('Starting sticky wages transition in alt calibration \n')
    whos
    PricePath_NK=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  TransPathFnsToEvaluate, NKTransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_NK);
end

%% Figure for Alternative Calibration
AltCalibFnsToEvaluate.output = @(d_val, aprime_val,a_val,z_val) d_val*z_val; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
AltCalibFnsToEvaluate.employment = @(d_val, aprime_val,a_val,z_val) d_val; %n_it in notation of GL2017

% For the initial dist will need
Params.phi=Params.phi_initial;
Params.r=p_eqm_initial.r;

% Get the 
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, AltCalibFnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
AggVarsPath_flex=EvalFnOnTransPath_AggVars_Case1(AltCalibFnsToEvaluate, PricePath_flex, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn,transpathoptions);

Output_pch_flex=([AggVars_initial.output.Mean; AggVarsPath_flex.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
Employment_pch_flex=([AggVars_initial.employment.Mean; AggVarsPath_flex.employment.Mean]-AggVars_initial.employment.Mean)/AggVars_initial.employment.Mean;

if altcalib_figurenumber==9 || altcalib_figurenumber==11
    AggVarsPath_NK=EvalFnOnTransPath_AggVars_Case1(AltCalibFnsToEvaluate,PricePath_NK, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn,transpathoptions_NK);
    Output_pch_NK=([AggVars_initial.output.Mean; AggVarsPath_NK.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
    Employment_pch_NK=([AggVars_initial.employment.Mean; AggVarsPath_NK.employment.Mean]-AggVars_initial.employment.Mean)/AggVars_initial.employment.Mean;
end


if CreateFigures==1
    % For some of the alternative calibrations employment is plotted as well.
    figure(altcalib_figurenumber)
    if altcalib_figurenumber==9
        % interest rate
        subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_flex.r],0:1:T,4*100*[p_eqm_initial.r; PricePath_NK.r])
        title('annual interest rate')
        % output
        subplot(1,2,2); plot(0:1:T,Output_pch_flex, 0:1:T,Output_pch_NK)
        title('output (pct deviation)')
        % ylabel('percent deviation from inital output in stationary eqm')
    elseif altcalib_figurenumber==10
        % interest rate
        subplot(1,3,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_flex.r])
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
        subplot(1,3,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_flex.r],0:1:T,4*100*[p_eqm_initial.r; PricePath_NK.r])
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
        subplot(1,2,1); plot(0:1:T,Full_PricePath_Flex_baseline, 0:1:T,4*100*[p_eqm_initial.r; PricePath_flex.r])
        title('annual interest rate')
        % output
        subplot(1,2,2); plot(0:1:T,Full_OutputPath_pch_Flex_baseline, 0:1:T,Output_pch_flex)
        title('output (pct deviation)')
        legend('baseline', 'higher risk aversion')
        % ylabel('percent deviation from inital output in stationary eqm')
    end
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.pdf'])
else
    if altcalib_figurenumber==9
        save(['./SavedOutput/GL2017Figs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.mat'], 'T', 'p_eqm_initial', 'PricePath_flex', 'PricePath_NK', 'Output_pch_flex', 'Output_pch_NK')
    elseif altcalib_figurenumber==10
        save(['./SavedOutput/GL2017Figs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.mat'], 'T', 'p_eqm_initial', 'PricePath_flex', 'Output_pch_flex', 'Employment_pch_flex')
    elseif altcalib_figurenumber==11
        save(['./SavedOutput/GL2017Figs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.mat'], 'T', 'p_eqm_initial', 'PricePath_flex', 'PricePath_NK', 'Output_pch_flex', 'Output_pch_NK', 'Employment_pch_flex', 'Employment_pch_NK')
    elseif altcalib_figurenumber==12
        save(['./SavedOutput/GL2017Figs/GuerrieriLorenzoni2017_Figure',num2str(altcalib_figurenumber),'.mat'], 'T', 'p_eqm_initial', 'PricePath_flex', 'Output_pch_flex', 'Full_PricePath_Flex_baseline', 'Full_OutputPath_pch_Flex_baseline')
    end
end

%% Clean up after itself
clear Policy_initial StationaryDist_initial V_final
disp('What is leftover from altcalib')
whos

