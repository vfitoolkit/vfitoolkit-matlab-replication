function Output=RestucciaRogerson2008_Fn(fixedsubsidy_c, normalize_employment, n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,ReturnFn,DiscountFactorParamNames, FnsToEvaluate, GEPriceParamNames, GeneralEqmEqns, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions)
% Only substantial change from the baseline codes is the need to determine
% Params.subsidyrate in general eqm so as to get K to equal Kbaseline.
% Note that "if Params.upsilon(:,1)>0" means that this is only done when
% subsidies are being used.

vfoptions.parallel=2; % use gpu

if sum(Params.upsilon(:,1))>0 && fixedsubsidy_c==0 % If using subsidies
    
    % Obviously we have to add subsidies. Less obvious is that to 'enforce'
    % Kbaseline we not only need to add it as a general eqm condition, we
    % also now need to include labour market clearance (and so include Ne
    % as a general eqm price parameter) because our earlier trick of just
    % 'renormalizing' afterwards will break the Kbaseline condition.
    GEPriceParamNames={'w','subsidyrate', 'Ne'};
        
    heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0,0};
    
    % A 'condlentry' general equilibrium condition will take values of greater
    % than zero for firms that decide to enter, less than zero for first that
    % decide not to enter (or more accurately, after entry decision they draw
    % their state, and then decide to cancel/abort their entry).
    
    GeneralEqmEqns.CondlEntry = @(ValueFn,beta) beta*ValueFn-0; % Conditional entry condition.
    GeneralEqmEqns.Entry = @(EValueFn,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
    GeneralEqmEqns.K = @(kbar,Kbaseline) kbar-Kbaseline; % K is unchanged (we choose subsidies to satisfy this, see RR2008 for explanation of why they wish to do this)
    GeneralEqmEqns.LabourMarket = @(nbar) nbar-1; % Labour market clearance (labour supply is perfectly inelastic and equal to one; we need labour demand to be equal to this; this condition determines mass of entering agents, Ne)
end 

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,~, ~]=HeteroAgentStationaryEqm_Case1(0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
Params.w=p_eqm.w;
Params.ebar=p_eqm.ebar;
if sum(Params.upsilon(:,1))>0 && fixedsubsidy_c==0  % If using subsidies
    Params.subsidyrate=p_eqm.subsidyrate;
    Params.Ne=p_eqm.Ne;
end

% Calculate some things in the general eqm
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);


AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

if sum(Params.upsilon(:,1))==0 || fixedsubsidy_c==1 % If no subsidies (or subsidy is fixed)
    % Without subsides we are not attempting to keep K unchanged (equal to
    % Kbaseline).
    % So we can do the labour market clearance condition here
    % as a renormalization without having to worry about this breaking 
    % the K equals Kbaseline condition.
    % Impose the labour market clearance, which involves calculating Ne.
    % AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
    InitialNe=Params.Ne;
    Params.Ne=1/AggVars.nbar.Aggregate; % AggValues is presently equal to Nbar. This line is imposing/satisfying the labour market clearance condition.
    StationaryDist.mass=StationaryDist.mass*(Params.Ne/InitialNe); % Take advantage of linearity of the stationary distribution in new entrants distribution.
end

Output.Params=Params;

%% Calculate various Aggregates, mostly related to Table 3

% 
% FnsToEvaluateFn_kbar=FnsToEvaluate{1};
% FnsToEvaluateFn_nbar=FnsToEvaluate{2};
% FnsToEvaluateFn_output=FnsToEvaluate{3};
% FnsToEvaluateFn_subsidy=FnsToEvaluate{4};
% FnsToEvaluateFn_outputofsubsidised=FnsToEvaluate{5};

Output.Y=AggVars.output.Aggregate;
Output.N=AggVars.nbar.Aggregate;
Output.K=AggVars.kbar.Aggregate;
Output.KdivY=Output.K/Output.Y;
Output.S=AggVars.subsidy.Aggregate;
Output.Ys=AggVars.outputofsubsidised.Aggregate;
% Params.w
Output.mass=StationaryDist.mass; % E in notation of RR2008
Output.TFP=(Output.Y/Output.N)./((Output.K/Output.N)^Params.alpha); % RR2008 call this 'A'.
Output.Ys_divY=AggVars.outputofsubsidised.Aggregate/AggVars.output.Aggregate; % The variable Ys/Y represents the output share of establishments that are receiving a subsidy
Output.SdivY=AggVars.subsidy.Aggregate/AggVars.output.Aggregate; % The variable S /Y is the total subsidies paid out to establishments receiving subsidies as a fraction of output
Output.tau_s=Params.subsidyrate; % The variable tau_s is the size of the subsidy required to generate a steady-state capital stock equal to that in the distortion-free economy

%% Calculate outputs related to Table 4
% Is just the TFP number, but for various different cases. So nothing further to calculate.

%% Tables 5 to 9 are all based on the outputs already generated. 
% They just report the same statistics for different economies/cases.

end