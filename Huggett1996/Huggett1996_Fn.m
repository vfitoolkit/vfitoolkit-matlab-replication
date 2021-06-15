function OutputResults=Huggett1996_Fn(Params, n_a,n_z,N_j, a_grid, ReturnFn, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, SSvaluesFn, SSvalueParamNames, GEPriceParamNames, GeneralEqmEqns, GeneralEqmEqnParamNames, simoptions, vfoptions,CheckUniquenessOfGE)
% Replication of Huggett (1996) - Wealth Distribution in Life Cycle Economies

%% Need to create appropriate z_grid and pi_z
Params.sigmasqy=Params.sigmasqepsilon./(1-Params.gamma.^2); 
% Initial distribution of income
Params.sigmasqy1=0.38;

% Huggett use 17 states equally spaced between +-4*sigmasqy1, with an 18th
% state at 6*sigmasqy1, "The transition probabilities between states are
% calculated by integrating the area under the normal distribution conditional on
% the current value of the state."
z_grid=[linspace(-4*sqrt(Params.sigmasqy1),4*sqrt(Params.sigmasqy1),17),6*sqrt(Params.sigmasqy1)]';
pi_z=nan(18,18);
% Following lines implement the transition matrix, they are largely just a copy of some code from the TauchenMethod() command.
sigma=sqrt(Params.sigmasqepsilon); %stddev of e
for ii=1:length(z_grid)
    pi_z(ii,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2-Params.gamma*z_grid(ii),0,sigma);
    for jj=2:(length(z_grid)-1)
        pi_z(ii,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2-Params.gamma*z_grid(ii),0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2-Params.gamma*z_grid(ii),0,sigma);
    end
    pi_z(ii,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2-Params.gamma*z_grid(ii),0,sigma);
end
z_grid=gpuArray(z_grid);
pi_z=gpuArray(pi_z);

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros(n_a,n_z,'gpuArray');
[trash,zeroassets_index]=min(abs(a_grid));
% Place them on the existing z_grid based on Params.sigmasqy1 (the variance
% of earnings at age 1) under assumption of normal distribution.
% Following lines implement this, they are largely just a copy of some code from the TauchenMethod() command.
sigma=sqrt(Params.sigmasqy1); %stddev of e
for ii=1:length(z_grid)
    jequaloneDist(zeroassets_index,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2,0,sigma);
    for jj=2:(length(z_grid)-1)
        jequaloneDist(zeroassets_index,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2,0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2,0,sigma);
    end
    jequaloneDist(zeroassets_index,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2,0,sigma);
end

%% Solve for the General Equilibrium
% Use the toolkit to find the equilibrium price index. There are two ways
% to do this. In what follows I use the 'search' approach to calculate the
% (initial) General equilibrium. Having used the search I then do a bunch
% of looking around a p_grid to be 'certain' that this is the only
% equilibrium, or at least the only one in the vicinity.

% Without p_grid, just searching. Use n_p=0. (Setting the actual algorithm
% used to 'search' can be done with heteroagentoptions.fminalgo)
heteroagentoptions.verbose=1;
[p_eqm,~, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, 0, pi_z, 0, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm.r;
Params.b=p_eqm.b;
Params.T=p_eqm.T;
fprintf('search finds equilibrium r=%8.2f, b=%8.2f, T=%8.2f \n')
save ./SavedOutput/Huggett1996.mat Params
% load ./SavedOutput/Huggett1996.mat Params

if CheckUniquenessOfGE==1
    % Using p_grid. This can be helpful if you want to, e.g., look for possibility of multiple equilibria.
    % We will use it as an informal double-check of the equilibrium we found previously
    % GEPriceParamNames={'r','b','T'}; % Already delared above.
    r_grid=linspace(0.5,2,21)'*Params.r; %okay
    b_grid=linspace(0.5,2,21)'*Params.b; %okay
    T_grid=linspace(0.5,2,21)'*Params.T; %good
    p_grid=[r_grid,b_grid, T_grid];
    %
    disp('Calculating price vector corresponding to the stationary eqm')
    n_p=[length(r_grid),length(b_grid),length(T_grid)];
    heteroagentoptions.pgrid=p_grid;
    heteroagentoptions.verbose=1;
    [p_eqm,p_eqm_index1, GeneralEqmEqnsValues1]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
    Params.r=p_eqm.r;
    Params.b=p_eqm.b;
    Params.T=p_eqm.T;
    fprintf('course grid search finds equilibrium r=%8.2f, b=%8.2f, T=%8.2f \n')
    
    % Grid search again, just to be super sure and accurate
    r_grid=linspace(0.95,1.05,21)'*Params.r; %okay
    b_grid=linspace(0.95,1.05,21)'*Params.b; %okay
    T_grid=linspace(0.95,1.05,21)'*Params.T; %good
    p_grid=[r_grid,b_grid, T_grid];
    [p_eqm,p_eqm_index2, GeneralEqmEqnsValues2]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
    Params.r=p_eqm.r;
    Params.b=p_eqm.b;
    Params.T=p_eqm.T;
    fprintf('fine grid search finds equilibrium r=%8.2f, b=%8.2f, T=%8.2f \n')
end

%% Compute a few things about the equilibrium
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,0,n_a,n_z,N_j,pi_z,Params,simoptions);

% Start with some basics (these were previously being done inside return function):
% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% From K/Y, substituting the production function gives
KdivY=(KdivL^(1-Params.alpha))/Params.A;
% We know w=MPL (MPL is marginal product of labour)
Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)
% Huggett (1996) calibrates tau to the following (see pg 478 for explanation)
Params.tau=0.195*(1-Params.delta*KdivY);

% Aggregate wealth transfers.% Start with some basics (these were previously being done inside return function):
% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% From K/Y, substituting the production function gives
KdivY=(KdivL^(1-Params.alpha))/Params.A;
% We know w=MPL (MPL is marginal product of labour)
Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)
% Huggett (1996) calibrates tau to the following (see pg 478 for explanation)
Params.tau=0.195*(1-Params.delta*KdivY);

% Aggregate wealth transfers.
AggregateWealthTransers=zeros(1,N_j);
for jj=1:Params.J
    for ii=1:jj
        AggregateWealthTransers(jj)=AggregateWealthTransers(jj)+Params.T*(1+Params.r*(1-Params.tau))^(ii-1);
    end
end
AggregateWealthTransers=sum(Params.mewj.*AggregateWealthTransers);
% Total wealth
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={}; 
FnsToEvaluate_TotalWealth = @(aprime_val,a_val,z_val) a_val;
FnsToEvaluate={FnsToEvaluate_TotalWealth};
TotalWealth=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, 0, n_a, n_z,N_j, 0, a_grid, z_grid);
% Transfer Wealth Ratio
TransferWealthRatio=AggregateWealthTransers/TotalWealth;


% Calculate fraction of population with zero or negative wealth
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={}; % ZeroOrNegAssets
FnsToEvaluate_ZeroOrNegAssets = @(aprime_val,a_val,z_val) (a_val<=0); % Indicator for ZeroOrNegAssets
FnsToEvaluate={FnsToEvaluate_ZeroOrNegAssets};
FractionWithZeroOrNegAssets=100*EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, 0, n_a, n_z,N_j, 0, a_grid, z_grid);

% Calculate wealth lorenz curve (and thus all percentile shares) and also
% the that for earnings (in text at bottom of pg 480, top of pg 481, there
% are a bunch of descriptions of model earnings, conditional on working age)
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate_1 = @(aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluate_1};
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, 0, n_a, n_z,N_j, 0, a_grid, z_grid);
TopWealthShares=100*(1-LorenzCurves([80,95,99],1)); % Need the 20,5, and 1 top shares for Tables of Huggett (1996)
% Calculate the wealth gini
WealthGini=Gini_from_LorenzCurve(LorenzCurves(:,1));

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,0,n_a,n_z,N_j,0,a_grid,z_grid);

%% Put together a bunch of things that I want to return as output

OutputResults.Params=Params;

OutputResults.AgeConditionalStats=AgeConditionalStats;

OutputResults.WealthGini=WealthGini;
OutputResults.TopWealthShares=TopWealthShares;
OutputResults.FractionWithZeroOrNegAssets=FractionWithZeroOrNegAssets;
OutputResults.TransferWealthRatio=TransferWealthRatio;
OutputResults.KdivY=KdivY;

OutputResults.z_grid=gather(z_grid);
OutputResults.pi_z=gather(pi_z);


if CheckUniquenessOfGE==1
    % Some things to allow looking at whether the equilibrium was unique
    OutputResults.p_eqm_index1=p_eqm_index1;
    OutputResults.GeneralEqmEqnsValues1=GeneralEqmEqnsValues1;
    OutputResults.p_eqm_index2=p_eqm_index2;
    OutputResults.GeneralEqmEqnsValues2=GeneralEqmEqnsValues2;
end

OutputResults=gather(OutputResults); % Just to be sure everything is on cpu, not gpu



end
























