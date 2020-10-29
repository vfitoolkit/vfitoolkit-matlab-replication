function Output=RestucciaRogerson2008_Fn(fixedsubsidy_c, normalize_employment, n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,ReturnFn,DiscountFactorParamNames, ReturnFnParamNames, FullFnsToEvaluate, GEPriceParamNames, GeneralEqmEqnParamNames, GeneralEqmEqns, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions)
% Only substantial change from the baseline codes is the need to determine
% Params.subsidyrate in general eqm so as to get K to equal Kbaseline.
% Note that "if Params.upsilon(:,1)>0" means that this is only done when
% subsidies are being used.

Parallel=simoptions.parallel;

FnsToEvaluateFn_kbar=FullFnsToEvaluate{1};
FnsToEvaluateFn_nbar=FullFnsToEvaluate{2};
FnsToEvaluateFn_output=FullFnsToEvaluate{3};
FnsToEvaluateFn_subsidy=FullFnsToEvaluate{4};
FnsToEvaluateFn_outputofsubsidised=FullFnsToEvaluate{5};

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate={};

if sum(Params.upsilon(:,1))>0 && fixedsubsidy_c==0 % If using subsidies
    
    % Obviously we have to add subsidies. Less obvious is that to 'enforce'
    % Kbaseline we not only need to add it as a general eqm condition, we
    % also now need to include labour market clearance (and so include Ne
    % as a general eqm price parameter) because our earlier trick of just
    % 'renormalizing' afterwards will break the Kbaseline condition.
    GEPriceParamNames={'w','subsidyrate', 'Ne'};
    
    FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
%     FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
    FnsToEvaluateParamNames(2).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
%     FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (((1-taurate*z2_val)*z1_val*gamma)/w)^(1/(1-gamma)) *((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^(alpha/(1-gamma)); % which evaluates to Nbar in the aggregate
    FnsToEvaluate={FnsToEvaluateFn_kbar,FnsToEvaluateFn_nbar};
    
    heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0,0};
    
    % A 'condlentry' general equilibrium condition will take values of greater
    % than zero for firms that decide to enter, less than zero for first that
    % decide not to enter (or more accurately, after entry decision they draw
    % their state, and then decide to cancel/abort their entry).
    
    GeneralEqmEqnParamNames(1).Names={'beta'};
    GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0; % Conditional entry condition.
    GeneralEqmEqnParamNames(2).Names={'beta','ce'};
    GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
    GeneralEqmEqnParamNames(3).Names={'Kbaseline'};
    GeneralEqmEqn_K = @(AggVars,GEprices,Kbaseline) AggVars(1)-Kbaseline; % K is unchanged (we choose subsidies to satisfy this, see RR2008 for explanation of why they wish to do this)
    GeneralEqmEqnParamNames(4).Names={};
    GeneralEqmEqn_LabourMarket = @(AggVars,GEprices) AggVars(2)-1; % Labour market clearance (labour supply is perfectly inelastic and equal to one; we need labour demand to be equal to this; this condition determines mass of entering agents, Ne)
    GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry,GeneralEqmEqn_K, GeneralEqmEqn_LabourMarket};
end 

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,~, ~]=HeteroAgentStationaryEqm_Case1(0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
Params.w=p_eqm.w;
Params.ebar=p_eqm.ebar;
if sum(Params.upsilon(:,1))>0 && fixedsubsidy_c==0  % If using subsidies
    Params.subsidyrate=p_eqm.subsidyrate;
    Params.Ne=p_eqm.Ne;
end

% Calculate some things in the general eqm
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

if sum(Params.upsilon(:,1))==0 || fixedsubsidy_c==1 % If no subsidies (or subsidy is fixed)
    % Without subsides we are not attempting to keep K unchanged (equal to
    % Kbaseline).
    % So we can do the labour market clearance condition here
    % as a renormalization without having to worry about this breaking 
    % the K equals Kbaseline condition.
    % Impose the labour market clearance, which involves calculating Ne.
    FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
    % FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (((1-taurate*z2_val)*z1_val*gamma)/w)^(1/(1-gamma)) *((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^(alpha/(1-gamma)); % which evaluates to Nbar in the aggregate
    FnsToEvaluate={FnsToEvaluateFn_nbar};
    AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
    InitialNe=Params.Ne;
    Params.Ne=1/AggValues; % AggValues is presently equal to Nbar. This line is imposing/satisfying the labour market clearance condition.
    StationaryDist.mass=StationaryDist.mass*(Params.Ne/InitialNe); % Take advantage of linearity of the stationary distribution in new entrants distribution.
end

Output.Params=Params;

%% Calculate various Aggregates, mostly related to Table 3
FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
% FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
FnsToEvaluateParamNames(2).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
% FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*z1_val*gamma)/w)^(1/(1-gamma)) *((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^(alpha/(1-gamma));
FnsToEvaluateParamNames(3).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
% FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) z1_val*(((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^alpha)*(((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^gamma);
FnsToEvaluateParamNames(4).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
% FnsToEvaluateFn_subsidy = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) subsidyrate*(z2_val<0)*(z1_val*(((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^alpha)*(((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma)))^gamma)); % The subsidy rate * receives subsidy * output
FnsToEvaluateParamNames(5).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
% FnsToEvaluateFn_outputofsubsidised
FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_nbar, FnsToEvaluateFn_output, FnsToEvaluateFn_subsidy, FnsToEvaluateFn_outputofsubsidised};

AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

Output.Y=AggValues(3);
Output.N=AggValues(2);
Output.K=AggValues(1);
Output.KdivY=Output.K/Output.Y;
Output.S=AggValues(4);
Output.Ys=AggValues(5);
% Params.w
Output.mass=StationaryDist.mass; % E in notation of RR2008
Output.TFP=(Output.Y/Output.N)./((Output.K/Output.N)^Params.alpha); % RR2008 call this 'A'.
Output.Ys_divY=AggValues(5)/AggValues(3); % The variable Ys/Y represents the output share of establishments that are receiving a subsidy
Output.SdivY=AggValues(4)/AggValues(3); % The variable S /Y is the total subsidies paid out to establishments receiving subsidies as a fraction of output
Output.tau_s=Params.subsidyrate; % The variable tau_s is the size of the subsidy required to generate a steady-state capital stock equal to that in the distortion-free economy

%% Calculate outputs related to Table 4
% Is just the TFP number, but for various different cases. So nothing further to calculate.

%% Tables 5 to 9 are all based on the outputs already generated. 
% They just report the same statistics for different economies/cases.

end