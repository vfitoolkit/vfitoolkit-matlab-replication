% Gourio & Miao (2010) - Firm Heterogeneity and the Long-run Effects of Dividend Tax Reform
%
% Two small notational changes: I call 'hweight' the parameter in the
% utility function (which GM2010 call 'h'). I call 'capadjconstant' the
% parameter in the capital adjustment costs function that GM2010 call 'phi'.
%
% GM2010 uses pure discretization with n_a=300 (they call it nK, line 64 of their divGEsetparametersBASELINE.m)
% As far as I can tell my results differ mainly just because they are substantially more accurate (no big differences, but everything is a bit different)

n_d=901; % Dividend
n_a=701; % Capital
n_z=10; % Productivity (GM2010 use 10 points for productivity)

vfoptions.lowmemory=1;

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Parameters

% Representative HH Preferences
Params.beta=0.971;
Params.hweight=6.616;

% Production
Params.alpha_k=0.311; % diminishing returns to capital input
Params.alpha_l=0.650; % diminishing returns to labor input
Params.delta=0.095; % Deprecitation of physical capital

% Capital adjustment costs
Params.capadjconstant=1.080; % term in the capital adjustment cost (GM2010 call this phi)

% Tax
Params.tau_corp=0.34; % Tax rate on corporate earnings
Params.phi=0; % Fraction of capital adjustment costs that can be deducted from corporate earnings (GM2010 don't allow this, but in US law it is possible, I have left it here from some codes I adapted)
Params.tau_d=0.25; % Tax rate on dividends
Params.tau_cg=0.2; % Tax rate on capital gains
Params.tau_i=0.25;

% Idiosyncatic productivity shocks
Params.rho_z=0.767;
Params.sigma_z_e=0.211;

% Market prices
% beta(r(1-tau_i)+1)=1 (eqn 20 of GM2010).We make sure our choices of beta, r and tau_i satisfy this equation, so we know that capital markets will clear
Params.r=(1/Params.beta-1)/(1-Params.tau_i); % rate of return on shareholdings
Params.w=1.2; % Actual value of wage will be determined in general eqm (when I first ran it was 1, made it 1.2 to be a bit closer to general eqm value so code runs faster)

Params.firmbeta=1/(1+Params.r*(1-Params.tau_i)/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax


% Table 3
FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table3.tex', 'w');
fprintf(FID, 'Baseline Parametrization \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \n');
fprintf(FID, '  & Parameter & Value  \\\\ \\hline \n');
fprintf(FID, 'Corporate income tax & $\\tau_{corp}$ & %8.3f \\\\ \n', Params.tau_corp);
fprintf(FID, 'Personal income tax & $\\tau_i$ & %8.3f \\\\ \n', Params.tau_i);
fprintf(FID, 'Dividend tax & $\\tau_d$ & %8.3f \\\\ \n', Params.tau_d);
fprintf(FID, 'Capital gains tax & $\\tau_{cg}$ & %8.3f \\\\ \n', Params.tau_cg);
fprintf(FID, 'Exponent on capital & $\\alpha_k$ & %8.3f \\\\ \n', Params.alpha_k);
fprintf(FID, 'Exponent on labor & $\\alpha_l$ & %8.3f \\\\ \n', Params.alpha_l);
fprintf(FID, 'Shock persistence & $\\rho$ & %8.3f \\\\ \n', Params.rho_z);
fprintf(FID, 'Shock standard deviation & $\\sigma$ & %8.3f \\\\ \n', Params.sigma_z_e);
fprintf(FID, 'Depreciation rate & $\\delta$ & %8.3f \\\\ \n', Params.delta);
fprintf(FID, 'Discount factor & $\\beta$ & %8.3f \\\\ \n', Params.beta);
fprintf(FID, 'Weight on leisure & $h$ & %8.3f \\\\ \n', Params.hweight);
fprintf(FID, 'Adjustment cost & $\\phi$ & %8.3f \\\\ \n', Params.capadjconstant);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);


%% Grids
d_grid=8*linspace(0,1,n_d)';
dgridspacing=d_grid(2)-d_grid(1);

% a_grid=5*(linspace(0.001^(1/3),1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
% GM2010 say they use minimum capital of 0.001
a_grid=8*linspace(0.001,1,n_a)'; % Max of 5 (I originally tried 20, and that turned out to be way to big; GM2010 codes shows they use 8)
agridspacing=a_grid(2)-a_grid(1); % Use an evenly spaced grid so that this can be used later

% GM2010 use the following grid with n_a=300
% [lines 71-75, their divGEsetparametersBASELINE.m]
% kstar=1.08; % 0.3^(nu/(1-alpha))*kovery^(1/(1-alpha)); % kovery=2.31; % this is alpha*Anorm/(r+delta)
% kstarhigh=8*kstar;
% a_grid=[linspace(1e-3,1.9,200),linspace(1.905,kstarhigh,n_a-200)]';
% Note: cannot really use it with current codes as a_grid is not evenly spaced

% Because of how the codes work I use agridspacing to determine the regime,
% and thus it must be that dgridspacing<agridspacing
if dgridspacing>=agridspacing
    error('You have to have smaller grid spacing for dividends than capital (so increase n_d or decrease n_a)')
end

% GM2010 say they use Tauchen-Hussey. Nowadays this is a bad idea and you should use Farmer-Toda instead (it is much better).
tauchenhusseyoptions.baseSigma=Params.sigma_z_e; % By default Tauchen-Hussey method uses the Floden improvement, this forces just the basic/original Tauchen-Hussey
[z_grid,pi_z] = discretizeAR1_TauchenHussey(0,Params.rho_z,Params.sigma_z_e,n_z,tauchenhusseyoptions); 
z_grid=exp(z_grid);
% I double check their transition matrix against mine, it is identical [line 82, their divGEsetparametersBASELINE.m]

%% Now, create the return function 
DiscountFactorParamNames={'firmbeta'};

% Notice we use 'GourioMiao2010_ReturnFn'
ReturnFn=@(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg) GourioMiao2010_ReturnFn(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg);

%% Test value function calculation
% vfoptions=struct(); % just use defaults
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% V is value function
% Policy is policy function (but as an index of k_grid, not the actual values)

%% Test stationary distribution calculation
simoptions=struct(); % just use defaults
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

%% Define aggregates and general eqm conditions
% Create functions to be evaluated
FnsToEvaluate.L = @(d,kprime,k,z,w,alpha_k,alpha_l) (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % labor (demanded by firms)
FnsToEvaluate.Y = @(d,kprime,k,z,w,alpha_k,alpha_l) z*(k^alpha_k)*( ((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1))) ^alpha_l); % Output (have to substitute for l)
FnsToEvaluate.I = @(d,kprime,k,z,delta) kprime-(1-delta)*k; % Investment
FnsToEvaluate.CapitalAdjCosts = @(d,kprime,k,z,delta,capadjconstant) (capadjconstant/2)*((kprime-(1-delta)*k)^2)/k; % Capital Adjustment Costs


GEPriceParamNames={'w'}; % We don't need r because of the representative household means we can just clear the capital market analytically due to the (tax-adjusted version of) the usual requirement that r=1/beta-1. Namely, beta(r(1-tau_i)+1)=1 (eqn 20 of GM2010). (That is, we made sure our choices of beta, r and tau_i satisfy this equation, so we know that capital markets will clear)
% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
% Labor market clearance is evaluated using household FOC:
% -U_1(C,L)/U_2(C,L)=(1-tau_i)w
% With U(C,L)=log(c)-h/2 L^2 this becomes
% hCL=(1-tau_i)w
% Note that we have to calculate C via the aggregate resource constraint:
% Y=C+I+CapitalAdjCosts [note: this implicitly says that taxes are just being lump-sumed to households for consumption]
GeneralEqmEqns.LaborMarket = @(L,Y,I,CapitalAdjCosts,hweight,tau_i,w) hweight*(Y-I-CapitalAdjCosts)*L-(1-tau_i)*w; % 
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate
% [Looking at GM2010 codes I see on line 37 of divGEfindeq.m that this is essentially exact same condition that they used.]


%% Find the general equilibrium
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

Params.w=p_eqm.w; % GM2010 find that general eqm wage is 1.26

[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);


% % Code of GM2010 contains a plot of the first of the following graphs of agent distribution, mine and theirs look broadly similar
% figure(2)
% subplot(3,1,1); plot(sum(StationaryDist,2))
% subplot(3,1,2); plot(a_grid,sum(StationaryDist,2))
% subplot(3,1,3); plot(a_grid,cumsum(sum(StationaryDist,2)))
% title('Agent Distribution')


fprintf('Check for trying to leave top of grid')
[max(max(Policy(1,:,:))),n_d]
[max(max(Policy(2,:,:))),n_a]
sum(sum(StationaryDist(end-10:end,:)))

% % Look more closely at hitting top of d_grid
% figure(1)
% plot(1:1:n_a,Policy(1,:,1),1:1:n_a,Policy(1,:,n_z),1:1:n_a,n_d*ones(1,n_a))
% legend('z1','zend','maxvalueofd')


%% Create Table 4
% Definition of 'earnings' is unclear from GM2010, guessing it is just profits?
FnsToEvaluate.investmentrate = @(d,kprime,k,z,w,alpha_k,alpha_l,delta) (kprime-(1-delta)*k) / k; % investment/capital (this is what GM2010 define the investment rate as, I had originally thought it would be investment/output)
FnsToEvaluate.investmentoutputratio = @(d,kprime,k,z,w,alpha_k,alpha_l,delta) (kprime-(1-delta)*k) / (z*(k^alpha_k)*( ((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1))) ^alpha_l));
FnsToEvaluate.D = @(d,kprime,k,z) d; % dividends
FnsToEvaluate.S = @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % new equity: from firm budget constraint
FnsToEvaluate.earnings = @(dividend,kprime,k,z,w,alpha_k,alpha_l) GourioMiao2010_EarningsFn(dividend,kprime,k,z,w,alpha_k,alpha_l); % earnings
FnsToEvaluate.dividendearningsratio = @(dividend,kprime,k,z,w,alpha_k,alpha_l) dividend/GourioMiao2010_EarningsFn(dividend,kprime,k,z,w,alpha_k,alpha_l); % dividend/earnings
FnsToEvaluate.earningscapitalratio = @(dividend,kprime,k,z,w,alpha_k,alpha_l) GourioMiao2010_EarningsFn(dividend,kprime,k,z,w,alpha_k,alpha_l)/k; % earnings/capital
FnsToEvaluate.newequityinvestmentratio = @(dividend,kprime,k,z,w,alpha_k,alpha_l,capadjconstant,tau_corp,phi,delta) GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)/(kprime-(1-delta)*k); % new equity/investment
FnsToEvaluate.K = @(d,kprime,k,z) k; % capital
FnsToEvaluate.z = @(d,kprime,k,z) z; % needed for table 8



AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

% Note: autocorrelations have to done manually from simulated panel data
simoptions.simperiods=500; % Default is to do 1000 simulations, each with this many periods
SimPanelValues=SimPanelValues_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions);
tempIrate=SimPanelValues.investmentrate(2:end,:);
tempIratelag=SimPanelValues.investmentrate(1:end-1,:);
autocorr_Irate=corr(tempIrate(:),tempIratelag(:)); % This is going to include some where it is last period of one firm correlated with first period of next firm, but given the length of the time series is 500 periods this should be negligible
tempearningscapital=SimPanelValues.earningscapitalratio(2:end,:);
tempearningscapitallag=SimPanelValues.earningscapitalratio(1:end-1,:);
autocorr_earningscapital=corr(tempearningscapital(:),tempearningscapitallag(:)); % This is going to include some where it is last period of one firm correlated with first period of next firm, but given the length of the time series is 500 periods this should be negligible

% Table 4
FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table4.tex', 'w');
fprintf(FID, 'Aggregate and Cross-Sectional Moments in the Baseline Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \n');
fprintf(FID, 'Variable & Data & Model  \\\\ \\hline \n');
fprintf(FID, 'I/Y & - & %8.3f \\\\ \n', AllStats.I.Mean/AllStats.Y.Mean); % Note: we know in this will just equal Params.delta
fprintf(FID, 'Investment rate (I/K) & - & %8.3f \\\\ \n', AllStats.I.Mean/AllStats.K.Mean); % Note: this is what GM2010 actually report as the investment rate (clear from their codes)
fprintf(FID, 'Aggregate dividends/earnings & - & %8.3f \\\\ \n', AllStats.D.Mean/AllStats.earnings.Mean); % My impression is that GM2010 report the ratio of aggregates (rather than mean of ratios)
fprintf(FID, 'Aggregate new equity/investment & - & %8.3f \\\\ \n', AllStats.S.Mean/AllStats.I.Mean);
fprintf(FID, 'Volatility of investment rate & - & %8.3f \\\\ \n', AllStats.investmentrate.StdDev);
fprintf(FID, 'Autocorrelation of investment rate & - & %8.3f \\\\ \n', autocorr_Irate);
fprintf(FID, 'Volatility of earnings/capital & - & %8.3f \\\\ \n', AllStats.earningscapitalratio.StdDev); % Guessing 'volatility' means std dev based on top of page 154 says "The model also underpredicts the standard deviation of the ratio of earnings to capital."
fprintf(FID, 'Autocorrelation of earnings/capital & - & %8.3f \\\\ \n', autocorr_earningscapital);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: The autocorrelations are calculated at firm level (based on time series) and then averaged across firms. Volatility is calculated as standard deviation based on pooled observations. From codes of GM2010 it is clear that Investment Rate is actually reporting I/K, not I/Y. I include I/Y as reading the paper this is what I had initially assumed it to be. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);
% I looked through GM2010 codes, they define "Volatility of investment
% rate" to be the sqrt(varIoverK), and "volatility of earnings/capital" as
% sqrt(varEoverK). varEoverK is the cross sectional variance of
% earnings/capital-. varIoverK is the cross sectional variance of the
% investment/capital.


C=AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean; % Aggregate consumption, not used for anything just calculated to check it
Tobins_q=sum(nansum(V.*StationaryDist))/AllStats.K.Mean; % Tobin's q, not used for anything just calculated to check it [note: nansum is needed because of V being -Inf times 0 mass gives nan]
fprintf('Check some things (just post Table 4 in codes) \n')
[AllStats.investmentrate.Mean, AllStats.I.Mean/AllStats.K.Mean] % investment/capital
[AllStats.newequityinvestmentratio.Mean, AllStats.S.Mean/AllStats.I.Mean] % new equity/investment
[AllStats.dividendearningsratio.Mean, AllStats.D.Mean/AllStats.earnings.Mean] % dividends/earnings

% Just save the workspace
save ./SavedOutput/GM2010.mat

%% In baseline tau_d and tau_cg are not equal, so firms should not both issue new equity and pay dividends
% Look at values of grid for s and dividends, make sure they seem to be mutally exclusive.

%% Create Table 5
% Things won't be precisely zero due to numerical approximation error, so we define a tolerance for what we consider zero
% zerotol is space between grid points, following the original code of GM2010 (which also requires using an evenly spaced grid).
Params.zerotol=agridspacing; % Needs to be in Params so can pass it as a function input
FnsToEvaluate.dividendregime2=@(dividend,kprime,k,z,zerotol) (dividend>zerotol); % 0.001 rather than zero to allow numerical error in exact zeros
FnsToEvaluate.equityregime2= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)>zerotol); % new equity: from firm budget constraint
FnsToEvaluate.neitherregime2= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (abs(dividend)<zerotol)*(abs(GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi))<zerotol);
FnsToEvaluate.neitherregime_test1= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (abs(dividend)<zerotol);
FnsToEvaluate.neitherregime_test2= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (dividend<zerotol);
FnsToEvaluate.neitherregime_test3= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (abs(GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi))<zerotol);
% Comment: I looked through GM2010 codes to see how they deal with this
% issue (that dividends and new equity issuance will never exactly equal
% zero due to numerical approximation error). They do something like what I
% do here, but it is not exactly the same. You can find it in their "divGEstatdis" on lines 155-174.
% They calculate what the dividend would be if there is no equity issuance
% (if s=0), and call this divifszero. When divifszero<0 they call it
% 'is_constrained', which is a firm not issuing dividends ("cannot finance
% with internal funds an investment greater than the optimal one, i.e.
% d=0"). When s policy is greater than divifszero-divifszero(k_policy one
% grid point below), then they call this 'is_using_external' ("is really
% using external fund; i.e. s>0"). They then define
% constrained_notext=is_constraied-is_using_external (%d=0, s=0), and
% is_using_external=1-is_constrained (%d>0).

% Figure out the three financing regimes
simoptions.conditionalrestrictions.dividendregime=@(dividend,kprime,k,z,zerotol) (dividend>zerotol);
simoptions.conditionalrestrictions.equityregime= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)>zerotol); % new equity: from firm budget constraint
simoptions.conditionalrestrictions.neitherregime= @(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,zerotol) (abs(dividend)<zerotol)*(abs(GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi))<zerotol);

% Recalculate all the statistics, but this time also doing so for the conditional restrictions
AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

% Following two should give the same answer (this is just a double-check
% that conditional restrictions are being correctly calculated)
[AllStats.equityregime2.Mean,AllStats.neitherregime2.Mean,AllStats.dividendregime2.Mean]
[AllStats.equityregime.RestrictedSampleMass,AllStats.neitherregime.RestrictedSampleMass,AllStats.dividendregime.RestrictedSampleMass]

FnsToEvaluate2.firmvalue=@(dividend,kprime,k,z,V) V; % V has to refer to value function
simoptions.eval_valuefn=V; % If you want to use the value function as part of a function to evaluate you must put the value function into this option...
simoptions.eval_valuefnname={'V'}; % ...and the name you will use in the function to evaluate in this option
% To use the value function in a function to evaluate, it must be the first
% entry after z. And you must input it into EvalFnOnAgentDist_AggVars_Case1() in the position after simoptions.
ValueOfFirmStats=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate2,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,[],simoptions);


% Table 5
FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table5.tex', 'w');
fprintf(FID, 'Distribution of Firms across Finance Regimes in the Baseline Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & Equity & Liquidity & Dividend  \\\\ \n');
fprintf(FID, ' & issuance regine & constrained regime & distribution regime  \\\\ \\hline \n');
fprintf(FID, 'Share of firms & %8.3f & %8.3f & %8.3f \\\\ \n', AllStats.equityregime.RestrictedSampleMass, AllStats.neitherregime.RestrictedSampleMass, AllStats.dividendregime.RestrictedSampleMass); 
fprintf(FID, 'Share of capital & %8.3f & %8.3f & %8.3f \\\\ \n', AllStats.equityregime.K.Total/AllStats.K.Mean, AllStats.neitherregime.K.Total/AllStats.K.Mean, AllStats.dividendregime.K.Total/AllStats.K.Mean); 
fprintf(FID, 'Share of investment & %8.3f & %8.3f & %8.3f \\\\ \n', AllStats.equityregime.I.Total/AllStats.I.Mean, AllStats.neitherregime.I.Total/AllStats.I.Mean, AllStats.dividendregime.I.Total/AllStats.I.Mean); 
fprintf(FID, 'Earnings-Capital Ratio & %8.3f & %8.3f & %8.3f \\\\ \n', AllStats.equityregime.earnings.Total/AllStats.equityregime.K.Total, AllStats.neitherregime.earnings.Total/AllStats.neitherregime.K.Total, AllStats.dividendregime.earnings.Total/AllStats.dividendregime.K.Total);
fprintf(FID, 'Investment-Capital Ratio & %8.3f & %8.3f & %8.3f \\\\ \n', AllStats.equityregime.I.Total/AllStats.equityregime.K.Total, AllStats.neitherregime.I.Total/AllStats.neitherregime.K.Total, AllStats.dividendregime.I.Total/AllStats.dividendregime.K.Total);
fprintf(FID, 'Average Tobins q & %8.3f & %8.3f & %8.3f \\\\ \n', ValueOfFirmStats.equityregime.firmvalue.Mean/AllStats.equityregime.K.Mean, ValueOfFirmStats.neitherregime.firmvalue.Mean/AllStats.neitherregime.K.Mean, ValueOfFirmStats.dividendregime.firmvalue.Mean/AllStats.dividendregime.K.Mean); % 'market value of a company divided by its assets replacement cost'
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: The share of firms for a regime is equal to the total number of firms in that regime divided by the total number of firms in all regimes. The share of capital (investment) for a regime is equal to the total capital stock (resp. investment) of the firms in that regime divided by aggregate capital stock (investment) of all firms. The earnings-capital ratio for a regime is equal to the earnings of the firms in that regime divided by their total capital stock, i.e., it is capital-weighted. The investment-capital ratio and Tobins q are computed in a similar way. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Store some things for Table 8
Table8=zeros(4,3);
Table8(1,1)=AllStats.Y.Mean/((AllStats.K.Mean^Params.alpha_k)*(AllStats.L.Mean^Params.alpha_l)); % TFP
Table8(2,1)=AllStats.Y.Mean/AllStats.L.Mean; % TFP
temp=corrcoef(log(SimPanelValues.K(:)),log(SimPanelValues.z(:)));
Table8(3,1)=temp(1,2);
temp=regress(log(SimPanelValues.K(:)),[ones(numel(SimPanelValues.z),1),log(SimPanelValues.z(:))]);
Table8(4,1)=temp(2); % temp(1) is the constant coeff

% Store some things for Table 10
Table10prereform=zeros(1,4);
Table10prereform(1,1)=AllStats.K.Mean;
Table10prereform(1,2)=AllStats.Y.Mean;
Table10prereform(1,3)=AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean;
Table10prereform(1,4)=Params.w;

% For the welfare calculations of last row of Table 6 we need to calculate the utility of the representative household
C=AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean;
Ubaseline=log(C)-(Params.hweight/2)*AllStats.L.Mean^2;

Params_baseline=Params; % For use later
AllStats_baseline=AllStats; % For use later
C_baseline=C; % For use later
ValueOfFirmStats_baseline=ValueOfFirmStats; % For use later


%% Distribution of firms
% The codes of GM2010 produce a version of the following graph


%% Create Table 6
OtherParametrizeationsCheck=struct(); % just used to check everything seems to be working as intended

Table6=zeros(8,4);
for Reform=1:4
    if Reform==1
        Params.tau_d=0.22;
        Params.tau_cg=0.2;
    elseif Reform==2
        Params.tau_d=0.2;
        Params.tau_cg=0.2;
    elseif Reform==3
        Params.tau_d=0.15;
        Params.tau_cg=0.15;
    elseif Reform==4
        Params.tau_d=0;
        Params.tau_cg=0;
    end
    % Changing tau_cg means we have to recalculate firmbeta
    Params.firmbeta=1/(1+Params.r*(1-Params.tau_i)/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax

    % Solve the tax reform
    [p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    Params.w=p_eqm.w;

    [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

    % Create what we need for Table 6
    AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

    % For the welfare calculations of last row of Table 6 we need to calculate the utility of the representative household
    % GM2010, pg 155: "The welfare benefit can be measured by the equivalent increase in consumption holding leisure constant."
    C=AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean;
    Ucurrent=log(C)-(Params.hweight/2)*AllStats.L.Mean^2;
    CEV=exp(Ucurrent-Ubaseline)-1; % Rearranged version of: ln(1+CEV)=Ucurrent-Ubaseline [note: U(C*(1+CEV),L)=ln(C*(1+CEV))-h/2L^2=ln(C)+ln(1+CEV)-h/2L^2=U(C,L)+ln(1+CEV)]

    % Use value function to calculate the value of firm (value of firm is just the value function)
    simoptions.eval_valuefn=V; % If you want to use the value function as part of a function to evaluate you must put the value function into this option...
    ValueOfFirmStats=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate2,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,[],simoptions);

    
    Table6(1,Reform)=(AllStats.K.Mean-AllStats_baseline.K.Mean)/AllStats_baseline.K.Mean;
    Table6(2,Reform)=(AllStats.Y.Mean-AllStats_baseline.Y.Mean)/AllStats_baseline.Y.Mean;
    Table6(3,Reform)=((AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean)-C_baseline)/C_baseline; % Resource constraint: C+I+CapitalAdjCosts=Y
    Table6(4,Reform)=(AllStats.D.Mean-AllStats_baseline.D.Mean)/AllStats_baseline.D.Mean;
    Table6(5,Reform)=(AllStats.S.Mean-AllStats_baseline.S.Mean)/AllStats_baseline.S.Mean;
    Table6(6,Reform)=(Params.w-Params_baseline.w)/Params_baseline.w;
    Table6(7,Reform)=(ValueOfFirmStats.firmvalue.Mean-ValueOfFirmStats_baseline.firmvalue.Mean)/ValueOfFirmStats_baseline.firmvalue.Mean; % Value of firm
    Table6(8,Reform)=CEV; % Welfare: "The welfare benefit can be measured by the equivalent increase in consumption holding leisure constant, which is only 0.31 percent."
    
    % Create what we need for Table 7
    if Reform==1
        % Just a repeat of what we did for Table 5 (note AllStats already includes the conditional restrictions)
        Table7=zeros(6,3);
        Table7(1,:)=[AllStats.equityregime.RestrictedSampleMass, AllStats.neitherregime.RestrictedSampleMass, AllStats.dividendregime.RestrictedSampleMass];
        Table7(2,:)=[AllStats.equityregime.K.Total/AllStats.K.Mean, AllStats.neitherregime.K.Total/AllStats.K.Mean, AllStats.dividendregime.K.Total/AllStats.K.Mean];
        Table7(3,:)=[AllStats.equityregime.I.Total/AllStats.I.Mean, AllStats.neitherregime.I.Total/AllStats.I.Mean, AllStats.dividendregime.I.Total/AllStats.I.Mean];
        Table7(4,:)=[AllStats.equityregime.earnings.Total/AllStats.equityregime.K.Total, AllStats.neitherregime.earnings.Total/AllStats.neitherregime.K.Total, AllStats.dividendregime.earnings.Total/AllStats.dividendregime.K.Total];
        Table7(5,:)=[AllStats.equityregime.I.Total/AllStats.equityregime.K.Total, AllStats.neitherregime.I.Total/AllStats.neitherregime.K.Total, AllStats.dividendregime.I.Total/AllStats.dividendregime.K.Total];
        Table7(6,:)=[ValueOfFirmStats.equityregime.firmvalue.Mean/AllStats.equityregime.K.Mean, ValueOfFirmStats.neitherregime.firmvalue.Mean/AllStats.neitherregime.K.Mean, ValueOfFirmStats.dividendregime.firmvalue.Mean/AllStats.dividendregime.K.Mean];
    end
    
    % Create what we need for Table 8
    if Reform==1 || Reform==2
        SimPanelValues=SimPanelValues_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions);

        Table8(1,Reform+1)=AllStats.Y.Mean/((AllStats.K.Mean^Params.alpha_k)*(AllStats.L.Mean^Params.alpha_l)); % TFP
        Table8(2,Reform+1)=AllStats.Y.Mean/AllStats.L.Mean; % TFP
        temp=corrcoef(log(SimPanelValues.K(:)),log(SimPanelValues.z(:)));
        Table8(3,Reform+1)=temp(1,2);
        temp=regress(log(SimPanelValues.K(:)),[ones(numel(SimPanelValues.z),1),log(SimPanelValues.z(:))]);
        Table8(4,Reform+1)=temp(2); % temp(1) is the constant coeff

    end    
    
    % Some things that I just check to make sure everything seems to be doing what it should
    OtherParametrizeationsCheck(Reform).Params=Params;
    OtherParametrizeationsCheck(Reform).Reform=Reform;
    OtherParametrizeationsCheck(Reform).GeneralEqmCondn=GeneralEqmCondn;
    OtherParametrizeationsCheck(Reform).p_eqm=p_eqm;
    OtherParametrizeationsCheck(Reform).AllStats=AllStats;
end
% Table 8 reports percentage changes in first two rows
Table8(1,:)=(Table8(1,:)-Table8(1,1))/Table8(1,1);
Table8(2,:)=(Table8(2,:)-Table8(2,1))/Table8(2,1);

% Just save the workspace again
save ./SavedOutput/GM2010_v2.mat


FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table6.tex', 'w');
fprintf(FID, 'Aggregate Effects of the Dividend Tax Reform in the Baseline Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, ' & $\\tau_d=0.22$ & $\\tau_d=0.20$ & $\\tau_d=0.15$ & $\\tau_d=0$ \\\\ \n');
fprintf(FID, ' & $\\tau_{cg}=0.20$ & $\\tau_{cg}=0.20$ & $\\tau_{cg}=0.15$ & $\\tau_{cg}=0$ \\\\ \\hline \n');
fprintf(FID, 'Capital & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(1,:)); 
fprintf(FID, 'Output & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(2,:)); 
fprintf(FID, 'Consumption & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(3,:)); 
fprintf(FID, 'Dividends & %8.2f & N/A & N/A & N/A \\\\ \n', 100*Table6(4,1)); 
fprintf(FID, 'Equity Issuance & %8.2f & N/A & N/A & N/A \\\\ \n', 100*Table6(5,1)); 
fprintf(FID, 'Wage & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(6,:)); 
fprintf(FID, 'Firm Value & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(7,:)); 
fprintf(FID, 'Welfare & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table6(8,:)); % GM2010, pg 155, "The welfare benefit can be measured by the equivalent increase in consumption holding leisure constant, which is only 0.31 percent."
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: All results are measured in percentage change from the initial steady state before the reform. When dividend tax and capital gains tax are equal, firms are indifferent between dividends and equity issuance, hence the N/A. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table7.tex', 'w');
fprintf(FID, 'Distribution of Firms across Finance Regimes for $\\tau_d=0.22$ and $\\tau_{cg}=0.20$ \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & Equity & Liquidity & Dividend  \\\\ \n');
fprintf(FID, ' & issuance regine & constrained regime & distribution regime  \\\\ \\hline \n');
fprintf(FID, 'Share of firms & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(1,:)); 
fprintf(FID, 'Share of capital & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(2,:)); 
fprintf(FID, 'Share of investment & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(3,:)); 
fprintf(FID, 'Earnings-Capital Ratio & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(4,:)); 
fprintf(FID, 'Investment-Capital Ratio & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(5,:)); 
fprintf(FID, 'Average Tobins q & %8.3f & %8.3f & %8.3f \\\\ \n', Table7(6,:)); % 'market value of a company divided by its assets replacement cost'
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: The share of firms for a regime is equal to the total number of firms in that regime divided by the total number of firms in all regimes. The share of capital (investment) for a regime is equal to the total capital stock (resp. investment) of the firms in that regime divided by aggregate capital stock (investment) of all firms. The earnings-capital ratio for a regime is equal to the earnings of the firms in that regime divided by their total capital stock, i.e., it is capital-weighted. The investment-capital ratio and Tobins q are computed in a similar way. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table8.tex', 'w');
fprintf(FID, 'Productivty Gains from the Dividend Tax Cut \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & $\\tau_d=0.25$ & $\\tau_d=0.22$ & $\\tau_d=0.20$  \\\\ \\hline \n');
fprintf(FID, 'Percentage change in TFP & %8.3f & %8.3f & %8.3f \\\\ \n', 100*Table8(1,:)); 
fprintf(FID, 'Percentage change in Y/L & %8.3f & %8.3f & %8.3f \\\\ \n', 100*Table8(2,:)); 
fprintf(FID, 'Correlation between ln(k) and ln(z) & %8.3f & %8.3f & %8.3f \\\\ \n', Table8(3,:));
fprintf(FID, 'Regression coefficient of ln(k) on ln(z) & %8.3f & %8.3f & %8.3f \\\\ \n', Table8(4,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);


%% Redoing some of the above with different parameter values for Tables 9 and 10
Table9=zeros(6,7); % ignore data column
Table10=zeros(7,4);

OtherParametrizeationsCheck=struct();

z_grid_baseline=z_grid;
pi_z_baseline=pi_z;

for otherparametrizations=1:7
    Params=Params_baseline;
    if otherparametrizations==1
        % Just the baseline parametrization
    elseif otherparametrizations==2
        Params.rho_z=0.65;
    elseif otherparametrizations==3
        Params.rho_z=0.85;
    elseif otherparametrizations==4
        Params.sigma_z_e=0.1;
    elseif otherparametrizations==5
        Params.sigma_z_e=0.3;
    elseif otherparametrizations==6
        Params.capadjconstant=0.5;
    elseif otherparametrizations==7
        Params.capadjconstant=1.5;
    end
    
    % If changing z, need to recreate the z_grid and pi_z
    if otherparametrizations==2 || otherparametrizations==3 || otherparametrizations==4 || otherparametrizations==5
        % GM2010 say they use Tauchen-Hussey. Nowadays this is a bad idea and you should use Farmer-Toda instead (it is much better).
        tauchenhusseyoptions.baseSigma=Params.sigma_z_e; % By default Tauchen-Hussey method uses the Floden improvement, this forces just the basic/original Tauchen-Hussey
        [z_grid,pi_z] = discretizeAR1_TauchenHussey(0,Params.rho_z,Params.sigma_z_e,n_z,tauchenhusseyoptions);
        z_grid=exp(z_grid);
    else
        z_grid=z_grid_baseline;
        pi_z=pi_z_baseline;
    end

    % Original economy for table 9
    [p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

    Params.w=p_eqm.w;

    [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

    AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

    % Note: autocorrelations have to done manually from simulated panel data
    SimPanelValues=SimPanelValues_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions);
    tempIrate=SimPanelValues.investmentrate(2:end,:);
    tempIratelag=SimPanelValues.investmentrate(1:end-1,:);
    autocorr_Irate=corr(tempIrate(:),tempIratelag(:));  % This is going to include some where it is last period of one firm correlated with first period of next firm, but given the length of the time series is 500 periods this should be negligible
    tempearningscapital=SimPanelValues.earningscapitalratio(2:end,:);
    tempearningscapitallag=SimPanelValues.earningscapitalratio(1:end-1,:);
    autocorr_earningscapital=corr(tempearningscapital(:),tempearningscapitallag(:)); % This is going to include some where it is last period of one firm correlated with first period of next firm, but given the length of the time series is 500 periods this should be negligible

    Table9(1,otherparametrizations)=AllStats.D.Mean/AllStats.earnings.Mean;
    Table9(2,otherparametrizations)=AllStats.S.Mean/AllStats.I.Mean;
    Table9(3,otherparametrizations)=AllStats.investmentrate.StdDev;
    Table9(4,otherparametrizations)=autocorr_Irate;
    Table9(5,otherparametrizations)=AllStats.earningscapitalratio.StdDev;
    Table9(6,otherparametrizations)=autocorr_earningscapital;

    % The main tax reform for Table 10
    Params.tau_d=0.15; % These rates are specified at top of pg 160 of GM2010
    Params.tau_cg=0.15;
    % Changing tau_cg means we have to recalculate firmbeta
    Params.firmbeta=1/(1+Params.r*(1-Params.tau_i)/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax
    
    [p_eqm2,~,GeneralEqmCondn2]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    
    Params.w=p_eqm2.w;
    
    [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);
    
    AllStats2=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    
    Table10(otherparametrizations,1)=(AllStats2.K.Mean-AllStats.K.Mean)/AllStats.K.Mean;
    Table10(otherparametrizations,2)=(AllStats2.Y.Mean-AllStats.Y.Mean)/AllStats.Y.Mean;
    C=AllStats.Y.Mean-AllStats.I.Mean-AllStats.CapitalAdjCosts.Mean;
    C2=AllStats2.Y.Mean-AllStats2.I.Mean-AllStats2.CapitalAdjCosts.Mean;
    Table10(otherparametrizations,3)=(C2-C)/C;
    Table10(otherparametrizations,4)=(p_eqm2.w-p_eqm.w)/p_eqm.w;

    % Some things that I just check to make sure everything seems to be doing what it should
    OtherParametrizeationsCheck(otherparametrizations).Params=Params;
    OtherParametrizeationsCheck(otherparametrizations).otherparametrizations=otherparametrizations;
    OtherParametrizeationsCheck(otherparametrizations).GeneralEqmCondn=GeneralEqmCondn;
    OtherParametrizeationsCheck(otherparametrizations).GeneralEqmCondn2=GeneralEqmCondn2;
end

% Just save the workspace again
save ./SavedOutput/GM2010_v3.mat

% Table 9
FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table9.tex', 'w');
fprintf(FID, 'Moments for Different Parameter Values \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccc} \n \\hline \n');
fprintf(FID, ' & Data & $\\rho=0.65$ & $\\rho=0.85$ & $\\sigma=0.1$ & $\\sigma=0.3$ & $\\phi=0.5$ & $\\phi=1.5$  \\\\ \\hline \n');
fprintf(FID, 'Aggregate dividends/earnings & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(1,2:end)); % My impression is that GM2010 report the ratio of aggregates (rather than mean of ratios)
fprintf(FID, 'Aggregate new equity/investment & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(2,2:end));
fprintf(FID, 'Volatility of investment rate & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(3,2:end));
fprintf(FID, 'Autocorrelation of investment rate & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(4,2:end));
fprintf(FID, 'Volatility of earnings/capital & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(5,2:end)); 
fprintf(FID, 'Autocorrelation of earnings/capital & - & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table9(6,2:end));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);

% Table 10
FID = fopen('./SavedOutput/LatexInputs/GourioMiao2010_Table10.tex', 'w');
fprintf(FID, 'Moments for Different Parameter Values \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \n');
fprintf(FID, ' & Capital & Output & Consumption & Wage  \\\\ \\hline \n');
fprintf(FID, 'Baseline & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(1,:));
fprintf(FID, '$\\rho=0.65$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(2,:));
fprintf(FID, '$\\rho=0.85$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(3,:));
fprintf(FID, '$\\sigma=0.1$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(4,:));
fprintf(FID, '$\\sigma=0.3$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(5,:));
fprintf(FID, '$\\phi=0.5$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(6,:));
fprintf(FID, '$\\phi=1.5$ & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', 100*Table10(7,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: All results are measured in percentage change of model with tax-reform from model before tax-reform. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);






