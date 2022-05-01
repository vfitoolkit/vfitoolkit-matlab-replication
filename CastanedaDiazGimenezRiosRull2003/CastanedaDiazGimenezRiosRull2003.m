% Replicates the results (and more) of Castaneda, Diaz-Gimenez, and Rios-Rull (2003)
% Takes quite a while to run, simply because two of the moments (out of 40 odd) take a long time to simulate. Actual solution of the model itself is
% fast enough.
% Figures use plotly
%
% I make a slight notational change, namely using h for hours worked and l
% for effeciency-unit-hours worked.
%
% The presence of the estate tax (together with stochastic probability of
% death) means next periods endogenous state (assets) cannot be chosen with
% certainty and problem is therefore a 'Case 2' value function problem.
%
% Note: The algorithm used for calculating the general equilibrium is not what you would want to use normally for solving this model. It finds the general equilibrium using a discrete grid on interest rates, rather
% than just solving the fixed-point problem on interest rates directly by using optimization. This is done for robustness reasons; see my paper on BHA models.

% A line needed to run on the server
addpath(genpath('./MatlabToolkits/'))

SkipGE=0 % Just a placeholder I am using to work on codes without rerunning the GE step.

%% Set some basic variables

n_l=51; % 101
n_k=1501; % 2001
n_r=551; % Two General Eqm variables: interest rate r, tax rate a3.
% Note: d1 is l (labour supply), d2 is a (assets)

% Some Toolkit options
vfoptions=struct(); % Just use defaults
simoptions.parallel=3; % Sparse matrix for agent distribution
simoptions.ncores=feature('numcores'); % Number of CPU cores

%% Parameters

Params.J=4;

% From Table 3
% Parameters
Params.beta=0.942; % Time discount factor
Params.sigma1=1.5; % Curvature of consumption
Params.sigma2=1.016; % Curvature of leisure
Params.chi=1.138; % Relative share of consumption and leisure
Params.elle=3.2; % Productive time
% Age and employment process
Params.p_eg=0.022; % Common probability of retiring
Params.p_gg=1-0.066; % 1 - p_gg = Common probability of dying
Params.phi1=0.969; % Earnings life cycle controller
Params.phi2=0.525; % Intergenerational earnings persistence controller
% Technology
Params.theta=0.376; % Capital share
Params.delta=0.059; % Capital depreciation rate
% Government Policy
Params.G=0.296; % Government expenditures % Note: G is not really a 'parameter', it is determined residually so as to balance the government budget balance which is a 'market clearance' condition
Params.omega=0.696; % Normalized transfers to retirees
Params.a0=0.258; % Income tax function parameters 1 of 4
Params.a1=0.768; % Income tax function parameters 2 of 4
Params.a2=0.491; % Income tax function parameters 3 of 4
Params.a3=0.144; % Income tax function parameters 4 of 4
Params.zlowerbar=14.101; % Estate tax function parameter: tax exempt level
Params.tauE=0.16; % Estate tax function parameter: marginal tax rate

% Process on exogenous shocks
% From Table 5
Params.e1=1; Params.e2=3.15; Params.e3=9.78; Params.e4=1061.00; % Params.e1=1 is a normalization.
% From Javiers codes they actually do the nomalization on e(2). So e1=1/e2; e3=e3/e2; e4=e4/e2; e2=1; %No longer need to do this.
% From Table 4 (the diagonal elements of Gamma_ee can be considered as residually determined to make each row of Gamma (Gamma_ee) add up to 1 (to 1-p_eg); 
%    in principle it doesn't matter which element in each row is residual, but from perspective of calibration it is easiest to 
%    use the diagonals to avoid contraint that they as transition probabilities they must be >=0 from binding)
Params.Gamma_ee_12=0.0114; Params.Gamma_ee_13=0.0039; Params.Gamma_ee_14=0.0001; % Params.Gamma_ee_11=0.9624;
Params.Gamma_ee_21=0.0307; Params.Gamma_ee_23=0.0037; Params.Gamma_ee_24=0;       % Params.Gamma_ee_22=0.9433;
Params.Gamma_ee_31=0.015;  Params.Gamma_ee_32=0.0043; Params.Gamma_ee_34=0.0002;  % Params.Gamma_ee_33=0.9582;
Params.Gamma_ee_41=0.1066; Params.Gamma_ee_42=0.0049; Params.Gamma_ee_43=0.0611;  % Params.Gamma_ee_44=0.8051;


%% Create the grids

[e_grid,Gamma,gammastar,gammastarfull]=CastanedaDiazGimenezRiosRull2003_Create_Exog_Shock(Params,vfoptions);

l_grid=linspace(0,Params.elle,n_l)';
k_grid=2500*(linspace(0,1,n_k).^3)'; % Note that the upper limit on my asset grid is larger than the 1500 in CDGRR2003, but a look at the stationary distribution shows noone gets this high (although the policy functions would get them there the exogenous shock for s==4 is not persistent enough)
% The next line gives the values used by CDGRR2003
% k_grid=[0:0.02:1,1.05:0.05:2,2.1:0.1:10,10.5:0.5:100,104:4:1500]';
% n_k=length(k_grid);

% Bring model into the notational conventions used by the toolkit
d_grid=[l_grid; k_grid]; % Is a 'Case 2' value function problem
a_grid=k_grid;
z_grid=linspace(1,2*Params.J,2*Params.J)'; %(age (& determines retirement))
pi_z=Gamma;

% r_grid=linspace(0,1/Params.beta-1,n_r)';
% a3_grid=linspace(0.9*Params.a3,1.1*Params.a3,n_a3)';
% p_grid=[r_grid; a3_grid];

n_d=[n_l,n_k];
n_a=n_k;
n_z=length(z_grid);

disp('sizes')
n_d
n_a
n_z

%% Set up the model itself
GEPriceParamNames={'r','a3'};

DiscountFactorParamNames={'beta'};

ReturnFn=@(l, kprime, k, z_val,r,sigma1,sigma2,chi,elle,theta,delta,e1,e2,e3,e4,omega,a0,a1,a2,a3) CastanedaDiazGimenezRiosRull2003_ReturnFn(l, kprime, k, z_val,r,sigma1,sigma2,chi,elle,theta,delta,e1,e2,e3,e4,omega,a0,a1,a2,a3);

% Case 2 requires 'phiaprime' which determines next periods assets from
% this periods decisions.
Case2_Type=2;
vfoptions.phiaprimematrix=1;
PhiaprimeParamNames={};
Phi_aprimeMatrix=CastanedaDiazGimenezRiosRull2003_PhiaprimeMatrix(n_d,n_z,k_grid,Params.J,Params.zlowerbar,Params.tauE,vfoptions);


FnsToEvaluate.K = @(h,kprime,k,s) k; %K
FnsToEvaluate.L = @(h,kprime,k,s,e1,e2,e3,e4) h*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); % Efficiency hours worked: L
FnsToEvaluate.IncomeTaxRevenue = @(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluate.Pensions = @(h,kprime,k,s,J,omega) omega*(s>J); % If you are retired you earn pension omega (otherwise it is zero).
FnsToEvaluate.EstateTaxRevenue  = @(h,kprime,k,s,J,p_gg,zlowerbar,tauE) (s>J)*(1-p_gg)*tauE*max(kprime-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay

% Now define the functions for the General Equilibrium conditions
%     Should be written so that closer the value given by the function is to zero, the closer the general eqm condition is to being met.
%The requirement that the interest rate equals the marginal product of capital
GeneralEqmEqns.CapitalMarket = @(r,K,L,theta,delta) r-(theta*(K^(theta-1))*(L^(1-theta))-delta); 
% Government budget balance
GeneralEqmEqns.GovBudget = @(G,Pensions,IncomeTaxRevenue,EstateTaxRevenue) G+Pensions-IncomeTaxRevenue-EstateTaxRevenue; % The roles of 'a3' is captured in the total revenue of income taxes

%% Test a few commands out before getting into the main part of General equilibrium
Params.r=0.045; %Params.a3
[V, Policy]=ValueFnIter_Case2(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Phi_aprimeMatrix, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, [], PhiaprimeParamNames, vfoptions);
StationaryDist=StationaryDist_Case2(Policy,Phi_aprimeMatrix,Case2_Type,n_d,n_a,n_z,pi_z,simoptions);
AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, 2);

%% Solve the baseline model

if SkipGE==0
    % Find the competitive equilibrium
%     heteroagentoptions.pgrid=p_grid;
    Params.r=0.045; %Params.a3
    heteroagentoptions.verbose=1;
    [p_eqm,p_eqm_index,MarketClearance]=HeteroAgentStationaryEqm_Case2(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid,Phi_aprimeMatrix, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], PhiaprimeParamNames, [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    
    % Evaluate a few objects at the equilibrium
    Params.r=p_eqm.r;
    Params.a3=p_eqm.a3;
    Params.w=(1-Params.theta)*(((Params.r+Params.delta)/(Params.theta))^(Params.theta/(Params.theta-1)));
        
    [V, Policy]=ValueFnIter_Case2(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Phi_aprimeMatrix, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, [], PhiaprimeParamNames, vfoptions);
    
    StationaryDist=StationaryDist_Case2(Policy,Phi_aprimeMatrix,Case2_Type,n_d,n_a,n_z,pi_z,simoptions);
        
    save ./SavedOutput/CastanedaDiazGimenezRiosRull2003.mat p_eqm Params MarketClearance a_grid V Policy StationaryDist

else
    load ./SavedOutput/CastanedaDiazGimenezRiosRull2003.mat p_eqm Params MarketClearance a_grid V Policy StationaryDist
end

%% Reproduce Tables
% Tables 3, 4, & 5 simply report the calibrated parameters. 
% While not really part of replication I reproduce these anyway (combining Tables 4 & 5 into a single table)
CastanedaDiazGimenezRiosRull2003_Tables345

% First, calculate all of model statistics that appear in Table 6
FnsToEvaluate=struct(); % clear out what was previously
FnsToEvaluate.K = @(h,kprime,k,s) k; %K
FnsToEvaluate.I = @(h,kprime,k,s,delta) kprime-k*(1-delta); %I
FnsToEvaluate.L = @(h,kprime,k,s,e1,e2,e3,e4) h*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); % Efficiency hours worked: L
FnsToEvaluate.H = @(h,kprime,k,s) h; %H
FnsToEvaluate.IncomeTaxRevenue = @(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluate.Pensions = @(h,kprime,k,s,J,omega) omega*(s>J); % If you are retired you earn pension omega (otherwise it is zero).
FnsToEvaluate.EstateTaxRevenue  = @(h,kprime,k,s,J,p_gg,zlowerbar,tauE) (s>J)*(1-p_gg)*tauE*max(kprime-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay
FnsToEvaluate.Consumption = @(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_ConsumptionFn(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);


Y=(AggVars.K.Mean^Params.theta)*(AggVars.L.Mean^(1-Params.theta));

Table6variables(1)=AggVars.K.Mean/Y; % K/Y
Table6variables(2)=100*AggVars.I.Mean/Y; % I/Y as %age
Table6variables(3)=100*(AggVars.IncomeTaxRevenue.Mean+AggVars.EstateTaxRevenue.Mean-AggVars.Pensions.Mean)/Y; % G/Y  as %age: G=T-Tr=(Income Tax Revenue+ Estate Tax Revenue) - Pensions
Table6variables(4)=100*AggVars.Pensions.Mean/Y; % Tr/Y  as %age
Table6variables(5)=100*AggVars.EstateTaxRevenue.Mean/Y; % T_E/Y  as %age
Table6variables(6)=100*AggVars.H.Mean/Params.elle; % h  as %age

MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

Table6variables(7)=gather((MeanMedianStdDev.Consumption.StdDev/MeanMedianStdDev.Consumption.Mean)/(MeanMedianStdDev.H.StdDev/MeanMedianStdDev.H.Mean)); % Coefficient of Variation=std deviation divided by mean. 

NSimulations=10^6;
e=[Params.e1,Params.e2,Params.e3,Params.e4,0,0,0,0];
tic;
% Ratio of earnings of 40 year olds to 20 year olds. This is quite complicated to calculate and so required a dedicated script.
Table6variables(8)=CDGRR2003_RatioEarningsOldYoung(NSimulations, StationaryDist, Policy, Phi_aprimeMatrix, n_d,n_a,n_z, d_grid, Gamma, e,Params.w,Params.J);
toc
tic;
% Intergenerational correlation coefficient. This is quite complicated to calculate and so required a dedicated script.
Table6variables(9)=CDGRR2003_IntergenerationalEarnings(NSimulations,StationaryDist, Policy, Phi_aprimeMatrix, n_d,n_a,n_z,d_grid, Gamma, e,Params.w,Params.J);
toc

%Table 6
FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table6.tex', 'w');
fprintf(FID, 'Values of the Targeted Ratios and Aggregates in the United States and in the Benchmark Model Economies \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, ' & $K/Y$ & $I/Y$ & $G/Y$ & $Tr/Y$ & $T_E/Y$ & $mean(h)$ & $CV_C/CV_H$ & $e_{40/20}$ & $\rho(f,s)$ \\\\ \n \\hline \n');
fprintf(FID, ' Target (USA) & 3.13  & 18.6\\%%  & 20.2\\%%  & 4.9\\%%   & 0.20\\%%  & 30.0\\%%  & 3.00  &  1.30 & 0.40  \\\\ \n');
fprintf(FID, ' Benchmark              & %8.2f & %8.1f\\%% & %8.1f\\%% & %8.1f\\%% & %8.2f\\%% & %8.1f\\%% & %8.2f & %8.2f & %8.2f \\\\ \n', Table6variables);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Variable $mean(h)$ (column 6) denotes the average share of disposable time allocated to the market. The statistic $CV_c/CV_h$ (column 7) is the ratio of the coefficients of variation of consumption and of hours worked. \\\\ \n');
fprintf(FID, '$e_{40/20}$ is the ratio of average earnings of 40 year old to 20 year old. $\rho(f,s)$ the intergenerational correlation coefficient between lifetime earnings of father and so. Note that model actually has households, while data is individuals.');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


% Lorenz Curves needed for Tables 7 and 8
StationaryDist_LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);


% Calculate Distributions of Earnings and Wealth for Table 7
Table7variables=nan(2,9,'gpuArray');
%  Gini for Earnings
Table7variables(1,1)=Gini_from_LorenzCurve(StationaryDist_LorenzCurves.L);
%  Earnings Lorenz Curve: Quintiles (%) 
Table7variables(1,2:6)=100*(StationaryDist_LorenzCurves.L([20,40,60,80,100])-StationaryDist_LorenzCurves.L([1,21,41,61,81]));
%  Earnings Lorenz Curve: 90-95, 95-99, and 99-100 (%)
Table7variables(1,7:9)=100*(StationaryDist_LorenzCurves.L([95,99,100])-StationaryDist_LorenzCurves.L([90,95,99]));
%  Gini for Wealth
Table7variables(2,1)=Gini_from_LorenzCurve(StationaryDist_LorenzCurves.K);
%  Wealth Lorenz Curve: Quintiles (%)
Table7variables(2,2:6)=100*(StationaryDist_LorenzCurves.K([20,40,60,80,100])-StationaryDist_LorenzCurves.K([1,21,41,61,81]));
%  Wealth Lorenz Curve: 90-95, 95-99, and 99-100 (%)
Table7variables(2,7:9)=100*(StationaryDist_LorenzCurves.K([95,99,100])-StationaryDist_LorenzCurves.K([90,95,99]));

%Table 7
FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table7.tex', 'w');
fprintf(FID, 'Distributions of Earnings and of Wealth in the United States and in the Benchmark Model Economies (\\%%) \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{(TOP GROUPS} \\\\ \n');
fprintf(FID, '& & \\multicolumn{5}{c}{QUINTILE} & \\multicolumn{3}{c}{(Percentile)} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, 'ECONOMY        & GINI  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '& \\multicolumn{9}{c}{A. Distribution of Earnings} \\\\ \n \\cline{2-10} \n');
fprintf(FID, ' United States & 0.63  & -0.40  & 3.19  & 12.49 & 23.33 & 61.39 & 12.38 & 16.37 & 14.76 \\\\ \n');
fprintf(FID, ' Benchmark     & %8.2f & %8.2f  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n \\cline{2-10} \n', Table7variables(1,:));
fprintf(FID, '& \\multicolumn{9}{c}{B. Distribution of Wealth} \\\\ \n \\cline{2-10} \n');
fprintf(FID, ' United States & 0.78  & -0.39  & 1.74  & 5.72 & 13.43  & 79.49 & 12.62 & 23.95 & 29.55 \\\\ \n');
fprintf(FID, ' Benchmark     & %8.2f & %8.2f  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table7variables(2,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note:  \\\\ \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Calculate Distributions of Consumption for Table 8
Table8variables=nan(2,9,'gpuArray');
% Calculate cutoff for wealthiest 1%
temp=cumsum(sum(StationaryDist,2)); % cdf on wealth dimension alone
[thisshouldequal99,cutoff_wealth1percent_index]=max(temp>0.99); % index
cutoff_wealth1percent=k_grid(cutoff_wealth1percent_index); % value
Params.cutoff_wealth1percent=cutoff_wealth1percent;
FnsToEvaluate2.Consumption_ExWealthiest1percent = @(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3,cutoff_wealth1percent) CDGRR2003_ConsumptionFn_ExWealthiest1percent(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3,cutoff_wealth1percent);
evalfnoptions.exclude.condn='equal';
evalfnoptions.exclude.values=-Inf; % 'FnsToEvaluateFn_Consumption_ExWealthiest1percent' returns a value of -Inf for the wealthiest 1%, so we want to exclude these.
LorenzCurves_ExWealthiest1percent=EvalFnOnAgentDist_LorenzCurve_Case2(StationaryDist, Policy, FnsToEvaluate2, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,evalfnoptions);
%  Gini for Consumption_ExWealthiest1percent
Table8variables(1,1)=Gini_from_LorenzCurve(LorenzCurves_ExWealthiest1percent.Consumption_ExWealthiest1percent);
%  Consumption_ExWealthiest1percent Lorenz Curve: Quintiles (%) 
Table8variables(1,2:6)=100*(LorenzCurves_ExWealthiest1percent.Consumption_ExWealthiest1percent([20,40,60,80,100])-LorenzCurves_ExWealthiest1percent.Consumption_ExWealthiest1percent([1,21,41,61,81]));
%  Consumption_ExWealthiest1percent Lorenz Curve: 90-95, 95-99, and 99-100 (%)
Table8variables(1,7:9)=100*(LorenzCurves_ExWealthiest1percent.Consumption_ExWealthiest1percent([95,99,100])-LorenzCurves_ExWealthiest1percent.Consumption_ExWealthiest1percent([90,95,99]));
%  Gini for Consumption
Table8variables(2,1)=Gini_from_LorenzCurve(StationaryDist_LorenzCurves.Consumption);
%  Consumption Lorenz Curve: Quintiles (%) 
Table8variables(2,2:6)=100*(StationaryDist_LorenzCurves.Consumption([20,40,60,80,100])-StationaryDist_LorenzCurves.Consumption([1,21,41,61,81]));
%  Consumption Lorenz Curve: 90-95, 95-99, and 99-100 (%)
Table8variables(2,7:9)=100*(StationaryDist_LorenzCurves.Consumption([95,99,100])-StationaryDist_LorenzCurves.Consumption([90,95,99]));


%Table 8
FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table8.tex', 'w');
fprintf(FID, 'Distribution of Consumption in the United States and in the Benchmark Model Economies (\\%%) \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{(TOP GROUPS} \\\\ \n');
fprintf(FID, '& & \\multicolumn{5}{c}{QUINTILE} & \\multicolumn{3}{c}{(Percentile)} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, 'ECONOMY        & GINI  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '\\multicolumn{10}{l}{United States:} \\\\ \n \\cline{2-10} \n');
fprintf(FID, '\\quad Nondurables   & 0.32  & 6.87 & 12.27 & 17.27 & 23.33 & 40.27 & 9.71 & 10.30 & 4.83 \\\\ \n');
fprintf(FID, '\\quad Nondurables+* & 0.30  & 7.19 & 12.96 & 17.80 & 23.77 & 38.28 & 9.43 &  9.69 & 3.77 \\\\ \n');
fprintf(FID, '\\multicolumn{10}{l}{Benchmark:} \\\\ \n');
fprintf(FID, '\\multicolumn{10}{l}{\\quad Wealthiest} \\\\ \n');
fprintf(FID, '\\quad 1\\%% Excluded & %8.2f & %8.2f  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table8variables(1,:));
fprintf(FID, '\\quad Entire Sample & %8.2f & %8.2f  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table8variables(2,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, '*: Includes imputed services of consumer durables. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);



% Calculate Mobility Statistics for Table 9
Table9variables=nan(2,5,'gpuArray');
% Seems like ideal method for mobility would be based on cupolas, but for
% now just use simulation methods.

% Transition Probabilities needed for Table 9
FnsToEvaluate3.L=FnsToEvaluate.L;
FnsToEvaluate3.K=FnsToEvaluate.K;
% 5 period transtion probabilities
t=5;
% Quintiles, so
npoints=5;
% Number of simulations on which to base results
NSims=10^7;
TransitionProbabilities=EvalFnOnAgentDist_RankTransitionProbabilities_Case2(t,NSims,StationaryDist, Policy,Phi_aprimeMatrix, Case2_Type, FnsToEvaluate3, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, 2, npoints);

Table9variables=zeros(2,5,'gpuArray');
for ii=1:5
    Table9variables(1,ii)=TransitionProbabilities.L(ii,ii); % Probability of staying in same quintile, hence ii-ii entry.
    Table9variables(2,ii)=TransitionProbabilities.K(ii,ii);
end

%Table 9
FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table9.tex', 'w');
fprintf(FID, 'Earnings and Wealth Persistence in the United States and in the Benchmark Model Economies: Fraction of Households That Remain In The Same Quintile After Five Years \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{5}{c}{QUINTILE} \\\\ \n');
fprintf(FID, 'ECONOMY & First  & Second & Third & Fourth & Fifth \\\\ \n \\hline \n');
fprintf(FID, '& \\multicolumn{5}{c}{A. Earnings Persistence} \\\\ \n \\cline{2-6} \n');
fprintf(FID, 'United States  & 0.86  & 0.41  & 0.47 & 0.46 & 0.66 \\\\ \n');
fprintf(FID, 'Benchmark      & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table9variables(1,:));
fprintf(FID, '& \\multicolumn{5}{c}{A. Wealth Persistence} \\\\ \n \\cline{2-6} \n');
fprintf(FID, 'United States  & 0.67  & 0.47  &  0.45 & 0.50  & 0.71  \\\\ \n');
fprintf(FID, 'Benchmark      & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table9variables(2,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on %d simulations. \\\\ \n', NSims);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Take a look at a few things

% % Take a look at the policy functions
% plot(l_grid(shiftdim(Policy(1,:,1:4),1)))
% plot(a_grid(shiftdim(Policy(2,:,1:4),1)))
%
% % Where is everyone on the asset distribution (making sure no-one hits the top)
% plot(cumsum(sum(StationaryDist,2)))

%% Finished Actual Replication, now do some extra things

%% Comparison of Calibrations
% The following shows how to use the VFI Toolkit to implement a calibration of this kind. However because the original weights assigned to each
% moment in CDGRR2003 have been lost to the sands of time this will not actually return the parameter values in CDGRR2003, nor should it be expected to.

% Ordering of following is unimportant. (25 params)
ParamNamesToEstimate={'beta','sigma2','chi','G','phi1','phi2','omega','a2','zlowerbar','tauE','e2','e3','e4',...
    'Gamma_ee_12','Gamma_ee_13','Gamma_ee_14','Gamma_ee_21','Gamma_ee_23','Gamma_ee_24','Gamma_ee_31','Gamma_ee_32','Gamma_ee_34','Gamma_ee_41','Gamma_ee_42','Gamma_ee_43'};
% Additionally 'r' and 'a3' are determined by the general eqm conditions, rather than the calibration.

% Ordering of following is unimportant. (27 targets)
EstimationTargetNames={'CapitalOutputRatio','GovExpenditureToOutputRatio','TransfersToOutputRatio',...
    'ShareOfDisposableTimeAllocatedToMarket','EffectiveTaxRateOnAverageHHIncome', 'zlowerbarMinus10timesAverageIncome', 'EstateTaxRevenueAsFractionOfGDP',...
    'RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked','RatioOfEarningsOldtoYoung','CrossSectionalCorrelationOfIncomeBetweenFathersAndSons',...
    'EarningsGini', 'WealthGini','EarningsQuintileSharesAsFraction', 'WealthQuintileSharesAsFraction','EarningsTopSharesAsFraction','WealthTopSharesAsFraction'};

% B.2 Macroeconomic Aggregates
EstimationTargets.CapitalOutputRatio=3.13;
% EstimationTargets.CapitalIncomeShare=0.376;
% Params.theta=0.376; % Follows immediately from CapitalIncomeShare
% EstimationTargets.InvestmentToOutputRatio=0.186;
% Params.delta=0.0594; % Follows immediately from delta=I/K in stationary general eqm; hence delta=(I/Y)/(K/Y)
EstimationTargets.GovExpenditureToOutputRatio=0.202;
EstimationTargets.TransfersToOutputRatio=0.049;

% B.3 Allocation of Time and Consumption
% Params.elle=3.2;
EstimationTargets.ShareOfDisposableTimeAllocatedToMarket=0.3;
EstimationTargets.RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked=3.0;
% Params.sigma1=1.5; % Based on literature on risk aversion

% B.4 The Age Structure of the Population
% EstimationTargets.ExpectedDurationOfWorkingLife=45;
% EstimationTargets.ExpectedDurationOfRetirement=18;
% These lead us directly to
% Params.p_eg=0.022; % Note: 1/p_eg=45
% Params.p_gg=0.934; % Note: 1/(1-p_gg)=15 % Not the 18 that it should be.
% [Follows from theoretical results on 'survival analysis': the expected duration of process with constant-hazard-rate lambda is 1/lambda. 
% Here p_eg and (1-p_gg) are the hazard rates. See, e.g., example on middle of pg 3 of http://data.princeton.edu/wws509/notes/c7.pdf ]

% B.5 Life-Cycle Profile of Earnings
% RatioOfEarningsOldtoYoung: ratio of average earnings for households
% between ages of 41 & 60 to average earnings of households between ages of 21 & 40.
EstimationTargets.RatioOfEarningsOldtoYoung=1.303;

% B.6 The Intergenerational Transmission of Earnings Ability
EstimationTargets.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=0.4;

% B.7 Income Taxation
% Params.a0=0.258;
% Params.a1=0.768;
EstimationTargets.EffectiveTaxRateOnAverageHHIncome=0.0762;
% The 'EffectiveTaxRateOnAverageHHIncome' is not reported in Castañeda, Diaz-Gimenez, & Rios-Rull (2003). 
% The number used here is 
% According to the 1998 Economic Report of the President, Table B80, revenue from 'Individual Income Taxes' in 1992 was $476 billion.
% According to the 1998 Economic Report of the President, Table B1, GDP in 1992 was $6244.4 billion
% So use 0.0762=476/6224 as target.
% Alternatively,
% According to the 1998 Economic Report of the President, Table B80, total federal revenue in 1992 was $1091.3 billion.
% So use 0.1748=1091/6224 as target.
% Note that this is ratio of aggregate totals, rather than strictly being the mean 
% effective rate on average income HH which would be a cross-sectional concept.
% According to the 1992 Survey of Consumer Finances the average HH income was $58916 (pg 837 of CDGRR2003)
% EstimationTargets: government budget balance
% EstimationTargets.GovernmentBudgetBalance=0; % This is a General Eqm
% Condition, so no need to repeat it here.

% B.8 Estate Taxation
% Params.zlowerbar=10*AverageIncome;
EstimationTargets.zlowerbarMinus10timesAverageIncome=0;
EstimationTargets.EstateTaxRevenueAsFractionOfGDP=0.002;

% B.9 Normalization
% Params.e1=1; % Based on my own experience with variants of this model you are actually
% better of normalizing Params.e2=1 than Params.e1=1, but as this is a
% replication I follow them exactly.
% Normalize the diagonal elements of Gamma_ee (ie., Gamma_ee_11,
% Gamma_ee_22, Gamma_ee_33, Gamma_ee_44). Choose these as simply setting,
% e.g., Gamma_ee_11=1-Gamma_ee_12-Gamma_ee_13-Gamma_ee_14 is unlikely to
% lead us to a negative value of Gamma_ee_11. Thus in terms of maintaining
% the constraints on Gamma_ee (all elements between zero and one, rows sum
% to one) required for it to be a transition matrix, we are less likely to
% be throwing out parameter vectors when estimating because they failed to
% meet these constraints. This is just a numerical trick that works well in
% practice as we 'know' that most of the weight of the transition matrix is
% on the diagonal.
% Note that this normalization of the diagonal elements of Gamma_ee is
% actually hard-coded into the how we have written the codes that create
% the transition matrix.
% Note also that paper does not actually specify which elements of Gamma_ee were normalized.

% B.10 The Distributions of Earnings and Wealth
EstimationTargets.EarningsGini=0.63;
EstimationTargets.WealthGini=0.78;
EstimationTargets.EarningsQuintileSharesAsFraction=[-0.004,0.0319, 0.1249, 0.2333, 0.6139]; % Quintiles: Bottom to Top
EstimationTargets.WealthQuintileSharesAsFraction=[-0.0039, 0.0174, 0.0572, 0.1343, 0.7949];
EstimationTargets.EarningsTopSharesAsFraction=[0.1238,0.1637,0.1476]; % 90-95, 95-99, 99-100.
EstimationTargets.WealthTopSharesAsFraction=[0.1262,0.2395,0.2955]; % 90-95, 95-99, 99-100.

% The Pension Function
% Castañeda, Diaz-Gimenez, & Rio-Rull (2003) do not describe the
% calibration of omega(s). From Table 3 we have that 
% Params.omega=0.8;
% suggesting the idea was to target the replacement rate.
% Actually this is covered by 'Transfers to Output Ratio' as being a target.

% By default the VFI Toolkit estimation commands set bounds on
% parameter values (lower bound of 1/10th of initial value, upper bound of
% 10 times initial value). You can set these bounds manually where you wish to do
% so in the following manner. [First number is lower bound, Second number
% is upper bound].
estimationoptions.ParamBounds.beta=[0.8,0.99]; % Reasonable range for discount rate.
estimationoptions.ParamBounds.r=[0,0.15]; % Seems reasonable range for interest rate.
estimationoptions.ParamBounds.Gamma_ee_12=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_13=[0,0.2]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_14=[0,0.2]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_21=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_23=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_24=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_31=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_32=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_34=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_41=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_42=[0,0.3]; % Must be between 0 & 1 as is a probability.
estimationoptions.ParamBounds.Gamma_ee_43=[0,0.3]; % Must be between 0 & 1 as is a probability.

% By default the VFI Toolkit estimation commands assume that you want the
% distance for each of the targets to be measured as the square difference as a percentage of the
% target value. You can overrule these as follows.
estimationoptions.TargetDistanceFns.EarningsQuintileSharesAsFraction='absolute-difference';
estimationoptions.TargetDistanceFns.WealthQuintileSharesAsFraction='absolute-difference';

% By default the VFI Toolkit weights each of the targets equally (with a
% value of 1). You can manually increase or decrease these weights as follows.
estimationoptions.TargetWeights.CapitalIncomeRatio=20;
% EstimationTargets.TargetWeights.GovernmentBudgetBalance=100; % This is one of the general eqm conditions, so by 
      % default it gets a weight of 100 when we are using the (default) 'joint-fixed-pt' estimation algorithm.
% Targets include an excess of inequality stats, so decrease slightly the weights given to these.
estimationoptions.TargetWeights.EarningsQuintileSharesAsFraction=0.8;
% estimationoptions.TargetWeights.EarningsTopSharesAsFraction=1; % 90-95, 95-99, 99-100.
estimationoptions.TargetWeights.WealthQuintileSharesAsFraction=0.8;
estimationoptions.TargetWeights.WealthTopSharesAsFraction=1.5; % Increased these as they are important part of the purpose of model, and were otherwise being ignored during the calibration (in earlier runs)
% The data and link to model are not strongest for the following two, so I give them lower weights.
estimationoptions.TargetWeights.RatioOfEarningsOldtoYoung=0.7;
estimationoptions.TargetWeights.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=0.7;
% An early estimation attempt ended up going off-track and making almost nobody work. Following makes the fraction of time worked estimation target important.
estimationoptions.TargetWeights.ShareOfDisposableTimeAllocatedToMarket=10;


% VFI Toolkit uses CMA-ES algorithm to perform the calibration. You can
% manually set some of its options if you want.
estimationoptions.CMAES.MaxIter=1000;

%% Before estimation we need to set some things back to what they were for underlying model
clear FnsToEvaluate
FnsToEvaluate.K = @(h,kprime,k,s) k; %K
FnsToEvaluate.L = @(h,kprime,k,s,e1,e2,e3,e4) h*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); % Efficiency hours worked: L
FnsToEvaluate.IncomeTaxRevenue = @(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluate.Pensions = @(h,kprime,k,s,J,omega) omega*(s>J); % If you are retired you earn pension omega (otherwise it is zero).
FnsToEvaluate.EstateTaxRevenue  = @(h,kprime,k,s,J,p_gg,zlowerbar,tauE) (s>J)*(1-p_gg)*tauE*max(kprime-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay

%% Now we just need to create the 'ModelTargetsFn'. This will be a Matlab function
% that takes Params as an input and creates ModelTargets as an output.
% ModelTargets must be a structure containing the model values for the EstimationTargets.
ModelTargetsFn=@(Params) CDGRR2003_ModelTargetsFn(Params, n_d,n_a,n_z,a_grid,ReturnFn, DiscountFactorParamNames,Case2_Type,PhiaprimeParamNames,FnsToEvaluate,GEPriceParamNames,GeneralEqmEqns, vfoptions,simoptions)
% ModelTargets must also contain the model values for any General Equilibrium conditions.
GeneralEqmTargetNames={'GE_InterestRate','GE_GovBudgetBalance'};

%% Do the actual calibration

[Params1,fval,counteval,exitflag]=CalibrateFromModelTargetsFn(Params, ParamNamesToEstimate, EstimationTargets, ModelTargetsFn, estimationoptions, GEPriceParamNames, GeneralEqmTargetNames);

save ./SavedOutput/Calib/CDGRR2003_Calib1.mat Params1 fval counteval exitflag
% load ./SavedOutput/Calib/CDGRR2003_Calib1.mat Params1 fval counteval exitflag

%% Get the model estimation target values based on the estimated parameters.
ModelTargets=ModelTargetsFn(Params1);

save ./SavedOutput/Calib/CDGRR2003_Calib2.mat ModelTargets






