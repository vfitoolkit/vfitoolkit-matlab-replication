function Output=ConesaKrueger1999_Fn(IndivProdShock,Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid,pi_z,ReturnFn, ReturnFnParamNames, DiscountFactorParamNames, jequaloneDist, GEPriceParamNames, AgeWeightsParamNames, T, ParamPath, PricePath0, GeneralEqmEqns_Transition, vfoptions, simoptions, heteroagentoptions,transpathoptions, Names_i, PTypeDistParamNames)
% Conesa & Krueger (1999) - Social Security Reform with Heterogeneous Agents

% Note that the initial and final eqm do depend on 'IndivProductivityShock'
% but not on 'PolicyReform'. I just anyway recompute them for 'PolicyReform'.

vfoptionslowmemtemp=vfoptions.lowmemory;
vfoptions.lowmemory=0; % Use this for stationary GE, then switch back for transition path
%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Steady State Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate_K = @(d_val,aprime_val,a_val,z_val) a_val; % Aggregate assets K
if IndivProdShock<=5
    FnsToEvaluateParamNames(2).Names={'I_j','epsilon_j'};
    FnsToEvaluate_N = @(d_val,aprime_val,a_val,z_val,I_j,epsilon_j) I_j*epsilon_j*z_val*d_val; % Aggregate labour supply (in efficiency units), CK1999 call this N
elseif IndivProdShock>=6
    FnsToEvaluateParamNames(2).Names={'I_j','epsilon_j','zeta'};
    FnsToEvaluate_N = @(d_val,aprime_val,a_val,z_val,I_j,epsilon_j,zeta) I_j*epsilon_j*exp(z_val+zeta)*d_val; % Aggregate labour supply (in efficiency units), CK1999 call this N
end
FnsToEvaluateParamNames(3).Names={'sj','n'};
FnsToEvaluate_L = @(d_val,aprime_val,a_val,z_val,sj,n) (1-sj)*aprime_val/(1+n); % Tr, accidental bequest transfers % The divided by (1+n) is due to population growth and because this is Tr_{t+1}
FnsToEvaluateParamNames(4).Names={'I_j'};
FnsToEvaluate_4 = @(d_val,aprime_val,a_val,z_val,I_j) I_j; % Fraction of population of working age (note that in principle this can calculated directly, but am computing it here as a way to double-check)
FnsToEvaluate={FnsToEvaluate_K,FnsToEvaluate_N,FnsToEvaluate_L,FnsToEvaluate_4};

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','tau','SS','Tr'}; In following lines p is the vector of these and so, e.g., p(2) is G.
GeneralEqmEqnParamNames=struct();
GeneralEqmEqnParamNames(1).Names={'theta','alpha','delta'};
GeneralEqmEqn_1 = @(AggVars,GEprices,theta,alpha,delta) GEprices(1)-(theta*(alpha)*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqnParamNames(2).Names={'b'};
GeneralEqmEqn_2 = @(AggVars,GEprices,b) GEprices(2)-b*(1-AggVars(4))/AggVars(4); % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
GeneralEqmEqnParamNames(3).Names={'b','theta','alpha'};
GeneralEqmEqn_3 = @(AggVars,GEprices,b,theta,alpha) GEprices(3)-b*(theta*(1-alpha)*((AggVars(1)/AggVars(2))^alpha))*AggVars(2)/AggVars(4); % Social Security adds up, based on eqn (2.12) b*w*N/(fraction working age) [this is eqn 2.12]
GeneralEqmEqnParamNames(4).Names={};
GeneralEqmEqn_4 = @(AggVars,GEprices) GEprices(4)-AggVars(3); % Accidental bequests (adjusted for population growth) are equal to transfers received (this is essentially eqn (2.14))
GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_2, GeneralEqmEqn_3, GeneralEqmEqn_4};

%% Solve for the initial General Equilibrium
% Use the toolkit to find the equilibrium price index. 
% In what follows I use the 'search' approach to calculate the
% (initial) General equilibrium.

if IndivProdShock<=5
    [p_eqm_init,p_eqm_index, GeneralEqmEqnsValues_init]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
elseif IndivProdShock>=6
    % There are different permanent types of agents
    [p_eqm_init,p_eqm_index, GeneralEqmEqnsValues_init]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, 0, pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, PTypeDistParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.r=p_eqm_init.r;
Params.tau=p_eqm_init.tau;
Params.SS=p_eqm_init.SS;
Params.Tr=p_eqm_init.Tr;
save ./SavedOutput/ConesaKrueger1999.mat Params p_eqm_init

% Some things about the initial general equilibrium  we will need for the transition path.
if IndivProdShock<=5
    [V_init, Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
    StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
    AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)
elseif IndivProdShock>=6
    [V_init, Policy_init]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
    StationaryDist_init=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
    AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid, simoptions); % The 2 is for Parallel (use GPU)
end

% Calculate the relevant entries for Table 5.
% Add a calculation of average time worked for working age
FnsToEvaluateParamNames(5).Names={};
FnsToEvaluate_5 = @(d_val,aprime_val,a_val,z_val) d_val; % Hours worked
FnsToEvaluate={FnsToEvaluate_K,FnsToEvaluate_N,FnsToEvaluate_L,FnsToEvaluate_4,FnsToEvaluate_5};

% Calculate the relevant entries for Table 5.
if IndivProdShock<=5
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)
    AggVars_MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid,2); % The 2 is for Parallel (use GPU)
elseif IndivProdShock>=6
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, Names_i,d_grid, a_grid, z_grid, simoptions);
    AggVars_MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid, z_grid,simoptions);
end

StationaryEqmStats(1).b=Params.b;
StationaryEqmStats(1).r=Params.r;
StationaryEqmStats(1).K=AggVars(1);
StationaryEqmStats(1).N1=AggVars(2);
KdivN1=AggVars(1)/AggVars(2);
KdivN2=((Params.r+Params.delta)/(Params.theta*Params.alpha))^(1/(Params.alpha-1)); % This is just to double check
StationaryEqmStats(1).w=Params.theta*(1-Params.alpha)*(KdivN1^Params.alpha);
StationaryEqmStats(1).h=100*AggVars(5)/AggVars(4); % total hours by working age divided by fraction working age. (multiply by 100 to express as percent rather than as fraction)
y=Params.theta*(AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha));
StationaryEqmStats(1).y=y;
StationaryEqmStats(1).Kdivy=AggVars(1)/y;
StationaryEqmStats(1).SSdivy=Params.SS/y;
StationaryEqmStats(1).HoursWorked=AggVars(5);

StationaryEqmStats(1).CoeffOfVariance_a=AggVars_MeanMedianStdDev(1,3)/AggVars_MeanMedianStdDev(1,1); % Standard deviation/Mean


% Statistics just for working age
simoptions.agegroupings=[1,Params.Jr]; % Working age, Retired (not actually interested in the numbers for retired)
if IndivProdShock<=5
    AllEmployedStats=LifeCycleProfiles_FHorz_Case1(StationaryDist_init,Policy_init,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
elseif IndivProdShock>=6
    simoptions.groupptypesforstats=1; % The statistics need to be calculated for the 'grouped PTypes', not seperately for each PType which is the default.
    AllEmployedStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_init,Policy_init,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
end

% Not sure from reading CK1999, but seems to be labour hours?
StationaryEqmStats(1).CoeffOfVariance_l=sqrt(AllEmployedStats(5).Variance(1))/AllEmployedStats(5).Mean(1);

% For Figure 1 & 2, need life cycle profile of average hours worked and assets.
% [Actual CK1999 figure looks like it may in fact use 5 year age bins, but I do individual ages here.]
simoptions.agegroupings=1:1:Params.J; % This is actually the default, but am setting it here explicitly anyway.
if IndivProdShock<=5
    LifeCycleProfiles_init=LifeCycleProfiles_FHorz_Case1(StationaryDist_init,Policy_init,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    % Note: initial eqm is labelled 'Pay as you go' in Figure 1.
    Figures_LifeCycleProfiles.Assets_init=LifeCycleProfiles_init(1).Mean;
    Figures_LifeCycleProfiles.HoursWorked_init=LifeCycleProfiles_init(5).Mean;
elseif IndivProdShock>=6
    % For the figures I get the different life-cycle profiles for each PType
    simoptions.groupptypesforstats=0;
    LifeCycleProfiles_init=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_init,Policy_init,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);    
    % Note: initial eqm is labelled 'Pay as you go' in Figure 1.
    for ii=1:length(Names_i)
        temp=LifeCycleProfiles_init.(Names_i{ii});
        Figures_LifeCycleProfiles.Assets_init.(Names_i{ii})=temp(1).Mean;
        Figures_LifeCycleProfiles.HoursWorked_init.(Names_i{ii})=temp(5).Mean;
    end
end


%% Solve for the final General Equilibrium
% Use the toolkit to find the equilibrium price index. There are two ways
% to do this. In what follows I use the 'search' approach to calculate the
% (initial) General equilibrium. But the commented out lines that follow it
% show how to set up the grid approach.

% All three policy reforms have the same end point, namely completely terminate the social security systen:
Params.b=Params.b_final;
% We know what the General Eqm value for SS will be, so may as well set this as our initial guess.
Params.tau=Params.tau_final;
Params.SS=Params.SS_final;

if IndivProdShock<=5
    [p_eqm_final,p_eqm_index, GeneralEqmEqnsValues_final]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
elseif IndivProdShock>=6
    % There are different permanent types of agents
    [p_eqm_final,p_eqm_index, GeneralEqmEqnsValues_final]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, 0, pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, PTypeDistParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end
Params.r=p_eqm_final.r;
Params.tau=p_eqm_final.tau;
Params.SS=p_eqm_final.SS;
Params.Tr=p_eqm_final.Tr;
save ./SavedOutput/ConesaKrueger1999.mat Params p_eqm_init p_eqm_final
% tau and ss should be zero, so if they are close to it then just set them
% to exactly zero (I only do this if close so that otherwise it is clear
% there is a problem)
if abs(Params.tau)<10^(-3)
    Params.tau=0;
end
if abs(Params.SS)<10^(-3)
    Params.SS=0;
end
    

% Some things about the final general equilibrium  we will need for the transition path.
if IndivProdShock<=5
    [V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
    StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
    AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)
elseif IndivProdShock>=6
    [V_final, Policy_final]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
    StationaryDist_final=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
    AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid, simoptions); % The 2 is for Parallel (use GPU)
end

% Calculate the relevant entries for Table 5.
% Add a calculation of average time worked for working age
FnsToEvaluateParamNames(5).Names={};
FnsToEvaluate_5 = @(d_val,aprime_val,a_val,z_val) d_val; % Hours worked
FnsToEvaluate={FnsToEvaluate_K,FnsToEvaluate_N,FnsToEvaluate_L,FnsToEvaluate_4,FnsToEvaluate_5};

if IndivProdShock<=5
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid, 2); % The 2 is for Parallel (use GPU)
    AggVars_MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid,2); % The 2 is for Parallel (use GPU)
elseif IndivProdShock>=6
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, Names_i,d_grid, a_grid, z_grid, simoptions); % The 2 is for Parallel (use GPU)
    AggVars_MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid, z_grid,simoptions); % The 2 is for Parallel (use GPU)    
end

StationaryEqmStats(2).b=Params.b;
StationaryEqmStats(2).r=Params.r;
StationaryEqmStats(2).K=AggVars(1);
KdivN1=AggVars(1)/AggVars(2);
KdivN2=((Params.r+Params.delta)/(Params.theta*Params.alpha))^(1/(Params.alpha-1)); % This is just to double check
StationaryEqmStats(2).w=Params.theta*(1-Params.alpha)*(KdivN1^Params.alpha);
StationaryEqmStats(2).h=100*AggVars(5)/AggVars(4); % total hours by working age divided by fraction working age. (multiply by 100 to express as percent rather than as fraction)
y=Params.theta*(AggVars(1)^Params.alpha)*(AggVars(2)^(1-Params.alpha));
StationaryEqmStats(2).y=y;
StationaryEqmStats(2).Kdivy=AggVars(1)/y;
StationaryEqmStats(2).SSdivy=Params.SS/y;
StationaryEqmStats(2).HoursWorked=AggVars(5);

StationaryEqmStats(2).CoeffOfVariance_a=AggVars_MeanMedianStdDev(1,3)/AggVars_MeanMedianStdDev(1,1); % Standard deviation/Mean

% Statistics just for working age
simoptions.agegroupings=[1,Params.Jr]; % Working age, Retired (not actually interested in the numbers for retired)
if IndivProdShock<=5
    AllEmployedStats=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
elseif IndivProdShock>=6
    simoptions.groupptypesforstats=1; % The statistics need to be calculated for the 'grouped PTypes', not seperately for each PType which is the default.
    AllEmployedStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_final,Policy_final,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
end

% Not sure from reading CK1999, but seems to be labour hours?
StationaryEqmStats(2).CoeffOfVariance_l=sqrt(AllEmployedStats(5).Variance(1))/AllEmployedStats(5).Mean(1);

if IndivProdShock<=5
    EV_SS_numerator=sum(V_final(1,:,1).*StationaryDist_final(1,:,1)); % 0 assets and age 1
    EV_SS_denominator=sum(V_init(1,:,1).*StationaryDist_init(1,:,1)); % 0 assets and age 1
elseif IndivProdShock>=6
    EV_SS_numerator=0;
    EV_SS_denominator=0;
    for ii=1:length(Names_i)
        EV_SS_numerator=sum(V_final.(Names_i{ii})(1,:,1).*StationaryDist_final.(Names_i{ii})(1,:,1)*StationaryDist_final.ptweights(ii)); % 0 assets and age 1
        EV_SS_denominator=sum(V_init.(Names_i{ii})(1,:,1).*StationaryDist_init.(Names_i{ii})(1,:,1)*StationaryDist_init.ptweights(ii)); % 0 assets and age 1
        % Note that ptweights are anyway the same for final and init.
    end
end

StationaryEqmStats(2).EV_SS=(EV_SS_numerator/EV_SS_denominator)^(1/(Params.gamma*(1-Params.sigma)))-1;

% For Figure 1 & 2, need life cycle profile of average hours worked and assets.
% [Actual CK1999 figure looks like it may in fact use 5 year age bins, but
% I do individual ages here.]
simoptions.agegroupings=1:1:Params.J; % This is actually the default, but am setting it here explicitly anyway.
if IndivProdShock<=5
    LifeCycleProfiles_final=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    % Note: final eqm is labelled 'Fully funded' in Figure 1.
    Figures_LifeCycleProfiles.Assets_final=LifeCycleProfiles_final(1).Mean;
    Figures_LifeCycleProfiles.HoursWorked_final=LifeCycleProfiles_final(5).Mean;
elseif IndivProdShock>=6
    % For the figures I get the different life-cycle profiles for each PType
    simoptions.groupptypesforstats=0;
    LifeCycleProfiles_final=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_final,Policy_final,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    for ii=1:length(Names_i)
        temp=LifeCycleProfiles_final.(Names_i{ii});
        Figures_LifeCycleProfiles.Assets_final.(Names_i{ii})=temp(1).Mean;
        Figures_LifeCycleProfiles.HoursWorked_final.(Names_i{ii})=temp(5).Mean;
    end
end

%%
save ./SavedOutput/ConesaKrueger1999_TransitionInputs.mat V_final Policy_final StationaryDist_init n_d n_a n_z N_j pi_z d_grid a_grid z_grid ReturnFn FnsToEvaluate Params DiscountFactorParamNames ReturnFnParamNames AgeWeightsParamNames FnsToEvaluateParamNames p_eqm_init p_eqm_final ParamPath

% load ./SavedOutput/ConesaKrueger1999_TransitionInputs.mat

%% Now, the transition path itself
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)
Params.r=p_eqm_init.r;
Params.tau=p_eqm_init.tau;
Params.SS=p_eqm_init.SS;
Params.Tr=p_eqm_init.Tr;

p_init=p_eqm_init; %0.75*p_eqm_init(2)]; % Ignoring the behavioural responses we would expect revenues (and thus spending) to fall to 0.75 of what they were, so lets use this as a 'starting guess'
p_final=p_eqm_final;

if vfoptions.fastOLG==1
    % We need to give an initial guess for the price path on interest rates
    PricePath0.r=[linspace(p_init.r, p_final.r, floor(T/2))'; p_final.r*ones(T-floor(T/2),1)];
    PricePath0.tau=p_init.tau*(2*ParamPath.b); % Because tau is so closely related to b this seems a reasonable guess (the 2* is so the init value of b of 0.5 is instead 1)
    PricePath0.SS=p_init.SS*(2*ParamPath.b); % Because SS is so closely related to b this seems a reasonable guess (the 2* is so the init value of b of 0.5 is instead 1)PricePath0.SS=[linspace(p_init.SS, p_final.SS, floor(T/2))'; p_final.SS*ones(T-floor(T/2),1)];
    PricePath0.Tr=[linspace(p_init.Tr, p_final.Tr, floor(T/2))'; p_final.Tr*ones(T-floor(T/2),1)];
    % PricePathNames={'r','tau','SS','Tr'};
else
    % Make sure the final is the correct one. (Have the fastOLG==1 solution
    % as initial guess, just need to make sure the end point is fully accurate)
    PricePath0.r(end-25:end)=p_final.r;
    PricePath0.tau(end-25:end)=p_final.tau;
    PricePath0.SS(end-25:end)=p_final.SS;
    PricePath0.Tr(end-25:end)=p_final.Tr;
end

%%
fprintf('Now computing the Transition path itself \n')
vfoptions.lowmemory=vfoptionslowmemtemp

transpathoptions.verbose=1

vfoptions.verbose=0;
if IndivProdShock<=5
    PricePath=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, transpathoptions, simoptions, vfoptions);
elseif IndivProdShock>=6
    PricePath=TransitionPath_Case1_FHorz_PType(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, Names_i, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, transpathoptions, simoptions, vfoptions);
end

save ./SavedOutput/ConesaKrueger1999_Transition.mat PricePath
% load ./SavedOutput/ConesaKrueger1999_Transition.mat PricePath

%% For Figures 3, 9, and 12
% Create graphs of output per capita, interest rates, capital-output ratio,
% and labor supply over the transition path.

Figures_TransitionAggregates.interestrate=[p_eqm_init.r; PricePath.r];

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate_K = @(d_val,aprime_val,a_val,z_val) a_val; % Aggregate assets K
if IndivProdShock<=5
    FnsToEvaluateParamNames(2).Names={'I_j','epsilon_j'};
    FnsToEvaluate_N = @(d_val,aprime_val,a_val,z_val,I_j,epsilon_j) I_j*epsilon_j*z_val*d_val; % Aggregate labour supply (in efficiency units), CK1999 call this N
elseif IndivProdShock>=6
    FnsToEvaluateParamNames(2).Names={'I_j','epsilon_j','zeta'};
    FnsToEvaluate_N = @(d_val,aprime_val,a_val,z_val,I_j,epsilon_j,zeta) I_j*epsilon_j*exp(z_val+zeta)*d_val; % Aggregate labour supply (in efficiency units), CK1999 call this N
end
FnsToEvaluateParamNames(3).Names={'sj','n'};
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluate_L = @(d_val,aprime_val,a_val,z_val) d_val; % Aggregate labour supply (in efficiency units), call this L
FnsToEvaluate={FnsToEvaluate_K,FnsToEvaluate_N,FnsToEvaluate_L};

if IndivProdShock<=5
    % FOLLOWING TWO LINES WERE ABOUT DEBUGGING
    fprintf('Check if this is zero %8.8f \n',max(max(max(max(max(abs(Policy_final-round(Policy_final))))))))
    fprintf('Check if this is nonzero %8.8f \n',min(min(min(min(Policy_final)))))
end

transpathoptions.vfoptions.policy_forceintegertype=1 % FOR SOME REASON THIS OTHERWISE ERRORED WHEN DOING IndivProdShock=0 with Reform 1
if IndivProdShock<=5
    AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, FnsToEvaluateParamNames, PricePath, ParamPath, Params, T, V_final, Policy_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions);
elseif IndivProdShock>=6
    AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, FnsToEvaluateParamNames, PricePath, ParamPath, Params, T, V_final, Policy_final, StationaryDist_init, n_d, n_a, n_z, N_j, Names_i, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);
end

K_Path=[StationaryEqmStats(1).K; AggVarsPath(:,1)];
N1_Path=[StationaryEqmStats(1).N1; AggVarsPath(:,2)];
L_Path=[StationaryEqmStats(1).HoursWorked; AggVarsPath(:,3)];
y_Path=Params.theta*(K_Path.^Params.alpha).*(N1_Path.^(1-Params.alpha));

Figures_TransitionAggregates.outputpercapita=y_Path; % Output per capita is just equal to output as mass of agents is normalized to one.
Figures_TransitionAggregates.capitaloutputratio=K_Path./y_Path;
Figures_TransitionAggregates.laborsupply=N1_Path;
Figures_TransitionAggregates.hoursworked=L_Path;

%% For Figures 4, 5, 6 and 7 we just need to calculate the EV
% (consumption equivalent variation) for each point on the grid between the
% initial period and the first period of the transition (CK1999, bottom pg
% 764). This is eqn (3.1) of CK1999.
% We already have V_init, which in the notation of eqn (3.1) is v_1
% So we just need to compute v_2, the value fn in the first period of the transition.

if IndivProdShock<=5
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, StationaryDist_init, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions);
    EquivVariation=(VPath(:,:,:,1)./V_init).^(1/(Params.gamma*(1-Params.sigma)))-1;
elseif IndivProdShock>=6
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final, Policy_final, StationaryDist_init, Params, n_d, n_a, n_z, N_j, Names_i, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);
    for ii=1:length(Names_i)
        EquivVariation.(Names_i{ii})=(VPath.(Names_i{ii})(:,:,:,1)./V_init.(Names_i{ii})).^(1/(Params.gamma*(1-Params.sigma)))-1; % Note, gamma & sigma are same across fixed types
    end
end


%% For Figure 8 
% Is only actually needed for one single combination of 'IndivProductivityShock'
% and 'PolicyReform'. But just computing it for all of them anyway.

% Partial equilibrium transition path, with prices fixed at their inital values.
PricePath_Fixed.r=p_init.r*ones(T,1);
PricePath_Fixed.tau=p_final.tau*ones(T,1); % Note that this is will be zero
PricePath_Fixed.SS=p_final.SS*ones(T,1); % Note that this is will be zero
PricePath_Fixed.Tr=p_init.Tr*ones(T,1);
% Note that will uses the ParamPath the reflects the policy changes.

if IndivProdShock<=5
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath_Fixed, ParamPath, T, V_final, Policy_final, StationaryDist_init, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions);
    EquivVariation_FixedPrices=(VPath(:,:,:,1)./V_init).^(1/(Params.gamma*(1-Params.sigma)))-1;
elseif IndivProdShock>=6
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath_Fixed, ParamPath, T, V_final, Policy_final, StationaryDist_init, Params, n_d, n_a, n_z, N_j, Names_i, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);
    for ii=1:length(Names_i)
        EquivVariation_FixedPrices.(Names_i{ii})=(VPath.(Names_i{ii})(:,:,:,1)./V_init.(Names_i{ii})).^(1/(Params.gamma*(1-Params.sigma)))-1; % Note, gamma & sigma are same across fixed types
    end
end

%% Things to output
Output.StationaryEqmStats=StationaryEqmStats; % For Table 5
Output.Figures_LifeCycleProfiles=Figures_LifeCycleProfiles; % For Figures 1&2
Output.Figures_TransitionAggregates=Figures_TransitionAggregates; % For Figures 3, 9 & 12
Output.EquivVariation=EquivVariation; % For Figures 4,5,6,7,10,11,13,14
Output.EquivVariation_FixedPrices=EquivVariation_FixedPrices; % For Figure 8

Output.Params=Params;
Output.StationaryDist_init=StationaryDist_init;
Output.PricePath=PricePath;

end