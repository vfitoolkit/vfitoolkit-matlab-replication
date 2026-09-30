function Output=Kitao2008_NoEntrepreneurs2(Params,n_asset, asset_grid, KdivYtarget, GdivYtarget, GiniWealthTarget, vfoptions,simoptions,vfoptionstpath)
% Note: without entreprenuers we are essentially solving the Aiyagari (1994) model, just with a different process for eta.

% K2008, Section 3.3 Economy without Entrepreneurs
% "In the second, the process is calibrated so that the distribution of wealth matches the one in the benchmark model with entrepreneurs. In
% particular, I take the approach as in Castaneda, Diaz-Gimenenz & Rio-Rull (2003), where agents face a highly persistent process of labor productivity
% with more than 60% of the population belonging to the lowest grid and a very small probability of possessing superb productivity that is more
% than 1000 times the media in the economy [Footnote: The exact productivity process of Castaneda, Diaz-Gimenenz & Rio-Rull (2003) would
% generate even more wealth concentration in our model and we adjust the value of the highest productivity shock so that we achieve the wealth
% Gini coefficient that is equivalent to the value in our benchmark model.]"

% In understand this to mean I copy-paste eta_grid and pi_eta from Castaneda, Diaz-Gimenenz & Rio-Rull (2003). And then recalibrate the value of
% eta_grid(end) to hit the Gini of wealth.


%% Create eta following Castaneda, Diaz-Gimenenz & Rio-Rull (2003)
n_eta=4; % following is just setting eta_grid and pi_eta
% Much of this is a lightly edited copy-paste of the VFI Toolkit replication of Castaneda, Diaz-Gimenenz & Rio-Rull (2003)
% https://github.com/vfitoolkit/vfitoolkit-matlab-replication/tree/master/CastanedaDiazGimenezRiosRull2003
% CDG2003 has n_eta as it has four points for working age and 4 for retirement.
% Here we are only using the four points for working age (because of how
% transitions work in CDGRR2003 there is a 'working age transitions' that
% is a submatrix of the full transtions so we can easily just take that and
% use it here)

% Process on exogenous shocks
% From Table 5 of CDGRR2003
e1=1; e2=3.15; e3=9.78; e4=1061.00; % Params.e1=1 is a normalization.
% From Javiers codes they do the normalization on e(2), so divide all four through by e2.
% NOTE ON ORDERING: this used to read "e1=e1/e2; e2=1; e3=e3/e2; e4=e4/e2;", which sets e2=1 before
% e3 and e4 are divided through, so e3 and e4 were left at their raw values of 9.78 and 1061 instead
% of 3.105 and 336.8. That mattered, because e4 is the starting value for the calibration below and
% Kitao2008_NoEntrepreneurs2_ExogShockFn hard-codes e3=9.78/3.15, so the two disagreed. Divide first,
% then set e2=1.
e1=e1/e2; e3=e3/e2; e4=e4/e2; e2=1; %#ok<NASGU> (e1, e2 and e3 are hard-coded inside the ExogShockFn, only e4 is used from here)
% Note: eta_grid is then rescaled inside Kitao2008_NoEntrepreneurs2_ExogShockFn so that the
% unconditional mean of eta is unity, which is the normalization Kitao (2008), Section 3.1, uses for
% the benchmark and which economy 'no entrepreneurs 1' therefore inherits. It has to be applied there
% rather than here, because e4 is calibrated to hit the Gini of wealth and so the mean of eta moves.
% (Note that with 61% of the population in the lowest grid point the MEDIAN of eta is e1, not e2, so
% normalizing on e2 does not in fact line the median worker up with the benchmark. Nothing does: the
% Castaneda, Diaz-Gimenez & Rios-Rull (2003) process is far more skewed than the benchmark process,
% so no single scale factor can match both the mean and the median. Kitao (2008) states the mean
% normalization, so that is the one used.)

% The transition matrix (Gamma_ee, from Table 4 of CDGRR2003) is set inside
% Kitao2008_NoEntrepreneurs2_ExogShockFn, so it is not repeated here.

% We have to set up this grid and transition matrix as a function so that
% we can set e4 as part of calibration (K2008 chooses e4, the max grid
% point in eta_grid, to hit the Gini of wealth).
vfoptions.ExogShockFn=@(e4) Kitao2008_NoEntrepreneurs2_ExogShockFn(e4);
simoptions.ExogShockFn=vfoptions.ExogShockFn;
% and put e4 into Params
Params.e4=e4;
% Get the (mean-normalized) initial grid and transition matrix from that same function
[eta_grid,pi_eta]=Kitao2008_NoEntrepreneurs2_ExogShockFn(Params.e4);

%% To be able to have the Gini of wealth as a calibration target we have to set it up as a CustomModelStat
heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions)... 
    Kitao2008_CustomModelStats_WealthGini(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions);


%%
n_d=0;
n_a=n_asset;
n_z=n_eta;
d_grid=0;
a_grid=asset_grid;
z_grid=eta_grid;
pi_z=pi_eta;

DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,eta,sigma,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k)... 
    Kitao2008_NoEntrepreneurs_ReturnFn(aprime,a,eta,sigma,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k);
% The first inputs must be: decision variables, next period endogenous state, endogenous state, exogenous state. Followed by any parameters

% Create functions to be evaluated
FnsToEvaluate.A = @(aprime,a,eta) a; % Total assets of households (workers and entrepreneurs)
FnsToEvaluate.L = @(aprime,a,eta) eta; % Total labor supply
FnsToEvaluate.TaxRevenue = @(aprime,a,eta,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k)... 
    Kitao2008_NoEntrepreneurs_TaxFn(aprime,a,eta,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k); % Tax Revenue

Params.KdivYtarget=KdivYtarget; % Add KdivYtarget target to Params so it can be used as part of the GeneralEqmEqns
Params.GdivYtarget=GdivYtarget; % Add GdivYtarget target to Params so it can be used as part of the GeneralEqmEqns
Params.GiniWealthTarget=GiniWealthTarget; % Add GiniWealthTarget target to Params so it can be used as part of the GeneralEqmEqns

%% General equilbrium
% Add beta and e4 to GE params as we want to recalibrate them to hit the capital-output ratio target and the Gini coeff of Wealth
GEPriceParamNames={'r','w','tau_I','beta','e4'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but not in transition)
% Note: Uses G from Params, which is the value from the main model

GeneralEqmEqns.CapitalMarket = @(r,A,L,alpha,delta) r-(alpha*(A^(alpha-1))*(L^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.LaborMarket = @(w,A,L,alpha) w-(1-alpha)*(A^(alpha))*(L^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.GovBudget = @(A,L,alpha,GdivYtarget,TaxRevenue) GdivYtarget*((A^alpha)*(L^(1-alpha)))-TaxRevenue; %Government runs balanced budget, with G set as a fixed fraction of GDP
% Note: Kitao (2008), Section 3.2, sets G "exogenously given as a fixed fraction of GDP", at 18%, and
% Section 3.3 says that in the economies without entrepreneurs "we maintain the values of baseline
% parameters at the benchmark levels ... and the tax parameters to satisfy the government budget
% constraint". So it is the RATIO G/Y, and not the level of G, that is held at the benchmark value here.
% This matters: these economies lose the entire non-corporate sector, so their GDP is roughly 17-25%
% below the benchmark, and holding G at the benchmark LEVEL would silently give them a government
% worth 21-24% of their own GDP rather than 18%, which is what drove tau_I up to 12.4% (economy 1)
% and 4.2% (economy 2) where Kitao (2008), Table 7, reports 3.03% for both.
% (For the tau_k policy experiments further below it is the LEVEL of G that is held fixed, at this
% economy's own initial steady state value; see Section 5 of Kitao (2008). That is swapped back in
% once the initial stationary general eqm has been solved.)
GeneralEqmEqns.CapitalOutputRatio = @(A,L,KdivYtarget,alpha) KdivYtarget-A/((A^(alpha))*(L^(1-alpha))); % Capital-output ratio is at target level
% Add Gini of wealth as a target [use CustomModelStats to target it]
GeneralEqmEqns.CalibGiniOfWealth = @(GiniWealth,GiniWealthTarget) GiniWealth-GiniWealthTarget; % Gini of wealth is at target level


%% Solve for the stationary general equilbirium
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.fminalgo=[8,1]; % fast but not really high accuracy, then a higher accuracy
heteroagentoptions.toleranceGEcondns=[1e-4,1e-5]; % high accuracy on final solve
% Note: see the note in Kitao2008.m on which MATLAB versions fminalgo=8 (lsqnonlin) needs

heteroagentoptions.multiGEweights=[1,1,1,2,5]; % was not really getting the Gini of wealth accurate, so putting a bigger weight on it to be sure we hit it as it is important to the exercise of 'no entrepreneur 2 economy'

[p_eqm_initial0,GeneralEqmCondn_initial0]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

Params.r=p_eqm_initial0.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_initial0.w;
Params.tau_I=p_eqm_initial0.tau_I;
Params.beta=p_eqm_initial0.beta;
Params.e4=p_eqm_initial0.e4;

% Hardcode the eta grid and get rid of ExogShockFn and CustomModelStats as they are no longer needed
[eta_grid,pi_eta]=Kitao2008_NoEntrepreneurs2_ExogShockFn(Params.e4);
z_grid=eta_grid; 
pi_z=pi_eta;
vfoptions=rmfield(vfoptions,'ExogShockFn');
simoptions=rmfield(simoptions,'ExogShockFn');
heteroagentoptions=rmfield(heteroagentoptions,'CustomModelStats');

% There is a small risk that the calibration target distracts from general eqm, so now drop them and clean up the general eqm
GEPriceParamNames={'r','w','tau_I'};
GeneralEqmEqns=rmfield(GeneralEqmEqns,'CapitalOutputRatio');
GeneralEqmEqns=rmfield(GeneralEqmEqns,'CalibGiniOfWealth');
heteroagentoptions=rmfield(heteroagentoptions,'multiGEweights');
[p_eqm_initial,GeneralEqmCondn_initial]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_initial % The equilibrium values of the GE prices
GeneralEqmCondn_initial
% We will want this eqm later for transition paths

Params.r=p_eqm_initial.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_initial.w;
Params.tau_I=p_eqm_initial.tau_I;



%% Now that we have the GE, let's calculate a bunch of related objects

[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% PolicyValues=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_s,d_grid,a_grid); % This will give you the policy in terms of values rather than index

StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);

StationaryDist_init=StationaryDist; % need this later for the transition path
V_init=V; % need this later for the CEV calculations

AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
% This AggVars is needed for Table 7
Table7.NoEntrepreneurs2.initial.r=p_eqm_initial.r;
Table7.NoEntrepreneurs2.initial.w=p_eqm_initial.w;
Table7.NoEntrepreneurs2.initial.tau_I=p_eqm_initial.tau_I;
Table7.NoEntrepreneurs2.initial.A=AggVars.A.Mean;
Y=(AggVars.A.Mean^(Params.alpha))*(AggVars.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs2.initial.Y=Y;


%% Calculate wealth inequality for Table 3
FnsToEvaluateIneq.Wealth=FnsToEvaluate.A;
simoptions.npoints=100; % Use 100 points for the lorenz cruve
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluateIneq, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

Table3.NoEntrepreneurs2.WealthGini=AllStats.Wealth.Gini;
Table3.NoEntrepreneurs2.WealthShareTop1=1-AllStats.Wealth.LorenzCurve(99);
Table3.NoEntrepreneurs2.WealthShareTop5=1-AllStats.Wealth.LorenzCurve(95);
Table3.NoEntrepreneurs2.WealthShareTop10=1-AllStats.Wealth.LorenzCurve(90);
Table3.NoEntrepreneurs2.WealthShareTop20=1-AllStats.Wealth.LorenzCurve(80);
Table3.NoEntrepreneurs2.WealthShareTop40=1-AllStats.Wealth.LorenzCurve(60);
Table3.NoEntrepreneurs2.WealthShareTop60=1-AllStats.Wealth.LorenzCurve(40);

% From here on we hold the LEVEL of G fixed, at the level implied by this economy's own initial
% stationary general eqm. This is what Kitao (2008), Section 5, does for the policy experiments:
% "All the policy experiments are revenue neutral tax reforms, i.e. we fix the government
% expenditures at the benchmark economy level."
Params.G=Params.GdivYtarget*Y;
GeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; %Government runs balanced budget

%% Now calculate the two 'final' equilibria: tau_k=0 and tau_k=0.4
% and the transition paths associated with them.

Params.taxincome=2; % 1 is tax income, 2 is to tax capital income and labor income seperately
Params.tau_k=0;
Params.tau_I=p_eqm_initial.tau_I; % initial guess

GEPriceParamNames={'r','w','tau_I'}; % tau_I instead of G


% tau_k=0
Params.tau_k=0;

[p_eqm_final1,GeneralEqmCondn_final1]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_final1 % The equilibrium values of the GE prices
GeneralEqmCondn_final1

% Now that we have the GE, let's calculate a bunch of related objects
Params.r=p_eqm_final1.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_final1.w;
Params.tau_I=p_eqm_final1.tau_I;

[V_final,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);
StationaryDist_final1=StationaryDist;

AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

Table7.NoEntrepreneurs2.tau_k_final1.r=p_eqm_final1.r;
Table7.NoEntrepreneurs2.tau_k_final1.w=p_eqm_final1.w;
Table7.NoEntrepreneurs2.tau_k_final1.tau_I=p_eqm_final1.tau_I;
Table7.NoEntrepreneurs2.tau_k_final1.A=AggVars.A.Mean;
Y=(AggVars.A.Mean^(Params.alpha))*(AggVars.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs2.tau_k_final1.Y=Y;

% Setup for transition path
T=50 % Kitao (2008) graphs suggest she uses 50 periods
ParamPath.tau_k=Params.tau_k*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'

% We need to give an initial guess for the price path
PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final1.r, floor(T/3))'; p_eqm_final1.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final1.w, floor(T/3))'; p_eqm_final1.w*ones(T-floor(T/3),1)];
PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final1.tau_I, floor(T/3))'; p_eqm_final1.tau_I*ones(T-floor(T/3),1)];

FnsToEvaluate_TransPath.A = @(aprime,a,eta) a; % Total assets of households (workers and entrepreneurs)
FnsToEvaluate_TransPath.L = @(aprime,a,eta) eta; % Total labor supply
FnsToEvaluate_TransPath.TaxRevenue = @(aprime,a,eta,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k)... 
    Kitao2008_NoEntrepreneurs_TaxFn(aprime,a,eta,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k); % Tax Revenue

TransPathGeneralEqmEqns.CapitalMarket = @(r,A,L,alpha,delta) r-(alpha*(A^(alpha-1))*(L^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
TransPathGeneralEqmEqns.LaborMarket = @(w,A,L,alpha) w-(1-alpha)*(A^(alpha))*(L^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
TransPathGeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; %Government runs balanced budget
% Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.

transpathoptions.GEnewprice=2; % Anderson acceleration of the shooting update
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice2.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'CapitalMarket','r',0,0.1;... % CapitalMarket is postivie is r is to large, so subtract
    'LaborMarket','w',0,0.1;... % LaborMarket is postivie is r is to large, so subtract
    'GovBudget','tau_I',1,0.1}; % GovBudget is positive if tau_I is too small, so add
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.
    % Note: GEnewprice=2 is Anderson acceleration, which accelerates exactly this shooting update, so it
    % needs the same howtoupdate instructions (just under GEnewprice2 rather than GEnewprice3).

% Now just run the TransitionPath_InfHorz command (all of the other inputs
% are things we had already had to define to be able to solve for the initial and final equilibria)
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;


PricePath=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions,simoptions,vfoptionstpath);

[Vpath,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);

% CEV calculations (CES utility fn makes this easy)
CEV=(Vpath(:,:,1)./V_init).^(1/(1-Params.sigma))-1;

gainCEV=sum(sum(StationaryDist_init((CEV>0)))); % based on initial agent dist

Table7.NoEntrepreneurs2.tau_k_final1.gainCEV=gainCEV;

%% Now for tau_k=0.4

Params.tau_k=0.4;
% GE was giving errors, so attempt a better initial guess
Params.r=0.2;
Params.w=1.2;
Params.tau_I=0.1;
% This better initial guess seems to have been enough, GE now solves fine
% (this is not a very good guess, but is enough to converge from)

[p_eqm_final2,GeneralEqmCondn_final2]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_final2 % The equilibrium values of the GE prices
GeneralEqmCondn_final2

% Now that we have the GE, let's calculate a bunch of related objects
Params.r=p_eqm_final2.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_final2.w;
Params.tau_I=p_eqm_final2.tau_I;

[V_final,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);
StationaryDist_final2=StationaryDist;

AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

Table7.NoEntrepreneurs2.tau_k_final2.r=p_eqm_final2.r;
Table7.NoEntrepreneurs2.tau_k_final2.w=p_eqm_final2.w;
Table7.NoEntrepreneurs2.tau_k_final2.tau_I=p_eqm_final2.tau_I;
Table7.NoEntrepreneurs2.tau_k_final2.A=AggVars.A.Mean;
Y=(AggVars.A.Mean^(Params.alpha))*(AggVars.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs2.tau_k_final2.Y=Y;

% Setup for transition path
ParamPath.tau_k=Params.tau_k*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'

% We need to give an initial guess for the price path
PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final2.r, floor(T/3))'; p_eqm_final2.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final2.w, floor(T/3))'; p_eqm_final2.w*ones(T-floor(T/3),1)];
PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final2.tau_I, floor(T/3))'; p_eqm_final2.tau_I*ones(T-floor(T/3),1)];

[PricePath,GECondPath]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions,simoptions,vfoptionstpath);

[Vpath,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);

% CEV calculations (CES utility fn makes this easy)
CEV=(Vpath(:,:,1)./V_init).^(1/(1-Params.sigma))-1;

gainCEV=sum(sum(StationaryDist_init((CEV>0)))); % based on initial agent dist

Table7.NoEntrepreneurs2.tau_k_final2.gainCEV=gainCEV;



%% Check that we don't hit the top of asset grid (Fig14 is not one from Kitao 2008, just something I want to see)
Fig14=figure(14);
assetdist_init=cumsum(sum(StationaryDist_init,2),1); % Note: dist is (n_a,n_z) here, so sum over dim 2 (eta)
assetdist_final1=cumsum(sum(StationaryDist_final1,2),1);
assetdist_final2=cumsum(sum(StationaryDist_final2,2),1);
semilogx(asset_grid,assetdist_init,asset_grid,assetdist_final1,asset_grid,assetdist_final2)
title('cdf of HHs over assets (No Entrepreneur 2 Model)')
xlabel('assets (model units)')
legend('init','final tauk=0','final tauk=0.4','Location','southeast')
saveas(Fig14,'./SavedOutput/Graphs/Kitao2008_Fig14.png')
% The decisive check is numerical rather than visual: if any noticeable mass has piled up at the top
% gridpoint then the grid is too small and the top of the wealth distribution is being truncated.
fprintf('No Entrepreneur 2: mass at the top asset gridpoint (init/tauk=0/tauk=0.4) = %e / %e / %e \n', 1-assetdist_init(end-1), 1-assetdist_final1(end-1), 1-assetdist_final2(end-1))
fprintf('No Entrepreneur 2: mass above one tenth of the top gridpoint            = %e / %e / %e \n', 1-assetdist_init(find(asset_grid>=asset_grid(end)/10,1)), 1-assetdist_final1(find(asset_grid>=asset_grid(end)/10,1)), 1-assetdist_final2(find(asset_grid>=asset_grid(end)/10,1)))


%%
Output.Table3=Table3;
Output.Table7=Table7;
Output.checkGE.p_eqm_initial=p_eqm_initial;
Output.checkGE.p_eqm_final1=p_eqm_final1;
Output.checkGE.p_eqm_final2=p_eqm_final2;
Output.checkGE.GeneralEqmCondn_initial0=GeneralEqmCondn_initial0;
Output.checkGE.GeneralEqmCondn_initial=GeneralEqmCondn_initial;
Output.checkGE.GeneralEqmCondn_final1=GeneralEqmCondn_final1;
Output.checkGE.GeneralEqmCondn_final2=GeneralEqmCondn_final2;
Output.checkGE.GECondPath=GECondPath;




























end