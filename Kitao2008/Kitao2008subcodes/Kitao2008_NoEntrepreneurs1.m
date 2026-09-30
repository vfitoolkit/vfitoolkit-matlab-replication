function Output=Kitao2008_NoEntrepreneurs1(Params,n_asset, n_eta, asset_grid, eta_grid, pi_eta,KdivYtarget,GdivYtarget,vfoptions,simoptions,vfoptionstpath)
% Note: without entreprenuers we are essentially solving the Aiyagari (1994) model, just with a different process for eta.

% K2008, Section 3.3 Economy without Entrepreneurs
% "In both we maintain the values of baseline parameters at the benchmark levels, except for the subjective discount factor to acheive the same capital-output 
% ratio of 2.65 and the tax parameters to satisfy the government budget constraint. The difference between the two worker economies is in the calibration of the 
% idiosyncratic labor productivity process. In the first, we use the same specification as in the benchmark."

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

%% General equilbrium
% Add beta to GE params as we want to recalibrate it to hit the capital-output ratio target
GEPriceParamNames={'r','w','tau_I','beta'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but not in transition)
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

% Set initial value for general eqm
% No need here as they all already have values in Params
Params.r=0.05;
Params.w=1;
Params.tau_I=0.01;
% Params.beta

%% Solve for the stationary general equilbirium
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.fminalgo=[8,1]; % fast but not really high accuracy, then a higher accuracy
heteroagentoptions.toleranceGEcondns=[1e-4,1e-5]; % high accuracy on final solve
% Note: see the note in Kitao2008.m on which MATLAB versions fminalgo=8 (lsqnonlin) needs

[p_eqm_initial0,GeneralEqmCondn_initial0]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

Params.r=p_eqm_initial0.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_initial0.w;
Params.tau_I=p_eqm_initial0.tau_I;
Params.beta=p_eqm_initial0.beta;

% There is a small risk that the calibration target distracts from general eqm, so now drop it and clean up the general eqm
GEPriceParamNames={'r','w','tau_I'};
GeneralEqmEqns=rmfield(GeneralEqmEqns,'CapitalOutputRatio');
[p_eqm_initial,GeneralEqmCondn_initial]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_initial % The equilibrium values of the GE prices
GeneralEqmCondn_initial
% We will want this eqm later for transition paths

Params.r=p_eqm_initial.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
Params.w=p_eqm_initial.w;
Params.tau_I=p_eqm_initial.tau_I;

%% Now that we have the GE, let's calculate a bunch of related objects

[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% PolicyValues=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_s,d_grid,a_grid,vfoptions); % This will give you the policy in terms of values rather than index

StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);

StationaryDist_init=StationaryDist; % need this later for the transition path
V_init=V; % need this later for the CEV calculations

AggVars_init=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
% This AggVars is needed for Table 7
Table7.NoEntrepreneurs1.initial.r=p_eqm_initial.r;
Table7.NoEntrepreneurs1.initial.w=p_eqm_initial.w;
Table7.NoEntrepreneurs1.initial.tau_I=p_eqm_initial.tau_I;
Table7.NoEntrepreneurs1.initial.A=AggVars_init.A.Mean;
Y=(AggVars_init.A.Mean^(Params.alpha))*(AggVars_init.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs1.initial.Y=Y;

%% Calculate wealth inequality for Table 3
FnsToEvaluateIneq.Wealth=FnsToEvaluate.A;
simoptions.npoints=100; % Use 100 points for the lorenz cruve

simoptions

AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluateIneq, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

Table3.NoEntrepreneurs1.WealthGini=AllStats.Wealth.Gini;
Table3.NoEntrepreneurs1.WealthShareTop1=1-AllStats.Wealth.LorenzCurve(99);
Table3.NoEntrepreneurs1.WealthShareTop5=1-AllStats.Wealth.LorenzCurve(95);
Table3.NoEntrepreneurs1.WealthShareTop10=1-AllStats.Wealth.LorenzCurve(90);
Table3.NoEntrepreneurs1.WealthShareTop20=1-AllStats.Wealth.LorenzCurve(80);
Table3.NoEntrepreneurs1.WealthShareTop40=1-AllStats.Wealth.LorenzCurve(60);
Table3.NoEntrepreneurs1.WealthShareTop60=1-AllStats.Wealth.LorenzCurve(40);

% From here on we hold the LEVEL of G fixed, at the level implied by this economy's own initial
% stationary general eqm. This is what Kitao (2008), Section 5, does for the policy experiments:
% "All the policy experiments are revenue neutral tax reforms, i.e. we fix the government
% expenditures at the benchmark economy level."
Params.G=Params.GdivYtarget*Y;
GeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; %Government runs balanced budget

%% Now calculate the two 'final' equilibria: tau_k=0 and tau_k=0.4
% and the transition paths associated with them.

Params.taxincome=2; % 1 is tax income, 2 is to tax capital income and labor income seperately
% Note: Params.tau_I is now irrelevant
Params.tau_k=0;
Params.tau_I=0;

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

Table7.NoEntrepreneurs1.tau_k_final1.r=p_eqm_final1.r;
Table7.NoEntrepreneurs1.tau_k_final1.w=p_eqm_final1.w;
Table7.NoEntrepreneurs1.tau_k_final1.tau_I=p_eqm_final1.tau_I;
Table7.NoEntrepreneurs1.tau_k_final1.A=AggVars.A.Mean;
Y=(AggVars.A.Mean^(Params.alpha))*(AggVars.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs1.tau_k_final1.Y=Y;

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

Table7.NoEntrepreneurs1.tau_k_final1.gainCEV=gainCEV;

%% Now for tau_k=0.4

Params.tau_k=0.4;
% GE was giving errors, so attempt a better initial guess
Params.r=0.2;
Params.w=1.2;
Params.tau_I=0.1;
% This better initial guess seems to have been enough, GE now solves fine
% (this is not a very good guess, but is enough to converge from)

heteroagentoptions.constrainpositive={'tau_I'}; % Got errors because return fn was negative, pretty sure it was because tau_I went negative, so trying this.

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

Table7.NoEntrepreneurs1.tau_k_final2.r=p_eqm_final2.r;
Table7.NoEntrepreneurs1.tau_k_final2.w=p_eqm_final2.w;
Table7.NoEntrepreneurs1.tau_k_final2.tau_I=p_eqm_final2.tau_I;
Table7.NoEntrepreneurs1.tau_k_final2.A=AggVars.A.Mean;
Y=(AggVars.A.Mean^(Params.alpha))*(AggVars.L.Mean^(1-Params.alpha));
Table7.NoEntrepreneurs1.tau_k_final2.Y=Y;

% Setup for transition path
ParamPath.tau_k=Params.tau_k*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'

% We need to give an initial guess for the price path
PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final2.r, floor(T/3))'; p_eqm_final2.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final2.w, floor(T/3))'; p_eqm_final2.w*ones(T-floor(T/3),1)];
PricePath0.tau_I=p_eqm_final2.tau_I*ones(T,1);

[PricePath,GECondPath]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions,simoptions,vfoptionstpath);

[Vpath,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);

% CEV calculations (CES utility fn makes this easy)
CEV=(Vpath(:,:,1)./V_init).^(1/(1-Params.sigma))-1;

gainCEV=sum(sum(StationaryDist_init((CEV>0)))); % based on initial agent dist

Table7.NoEntrepreneurs1.tau_k_final2.gainCEV=gainCEV;



%% Check that we don't hit the top of asset grid (Fig13 is not one from Kitao 2008, just something I want to see)
Fig13=figure(13);
assetdist_init=cumsum(sum(StationaryDist_init,2),1); % Note: dist is (n_a,n_z) here, so sum over dim 2 (eta)
assetdist_final1=cumsum(sum(StationaryDist_final1,2),1);
assetdist_final2=cumsum(sum(StationaryDist_final2,2),1);
plot(asset_grid,assetdist_init,asset_grid,assetdist_final1,asset_grid,assetdist_final2)
title('cdf of HHs over assets (No Entrepreneur 1 Model)')
xlabel('assets (model units)')
legend('init','final tauk=0','final tauk=0.4','Location','southeast')
saveas(Fig13,'./SavedOutput/Graphs/Kitao2008_Fig13.png')
% The decisive check is numerical rather than visual: if any noticeable mass has piled up at the top
% gridpoint then the grid is too small and the top of the wealth distribution is being truncated.
fprintf('No Entrepreneur 1: mass at the top asset gridpoint (init/tauk=0/tauk=0.4) = %e / %e / %e \n', 1-assetdist_init(end-1), 1-assetdist_final1(end-1), 1-assetdist_final2(end-1))
fprintf('No Entrepreneur 1: mass above one tenth of the top gridpoint            = %e / %e / %e \n', 1-assetdist_init(find(asset_grid>=asset_grid(end)/10,1)), 1-assetdist_final1(find(asset_grid>=asset_grid(end)/10,1)), 1-assetdist_final2(find(asset_grid>=asset_grid(end)/10,1)))


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