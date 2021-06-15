function OutputVector=Aiyagari1994_Fn(n_k,n_z,n_p,Params,vfoptions,simoptions)
%This code replicates the results of Aiyagari (1994) - Uninsured Idiosyncratic Risk and Aggregate Saving
disp('Running Aiyagari1994_Fn')
%% Set up


% %Parameters
% CreateIndividualParams(Params)
% whos

% beta=0.96; %Model period is one-sixth of a year
% alpha=0.36;
% delta=0.08;
% 
% mu=1; % {1,3,5}
% sigma=0.2; % {0.2,0.4}
% rho=0; % {0,0.3,0.6,0.9}

%Create markov process for the exogenous labour productivity, l.
% q=3; %Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3
[z_grid, pi_z]=TauchenMethod(0,(Params.sigma^2)*(1-Params.rho^2),Params.rho,n_z,Params.q); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)


%[s_grid, pi_z]=TauchenMethod_Param(0,(sigma^2)*(1-rho^2),rho,7,3); %This is the process for ln(l). Aiyagari uses 7 states, 
[z_mean,z_variance,z_corr,~]=MarkovChainMoments(z_grid,pi_z);
z_grid=exp(z_grid);
%Get some info on the markov process
[Expectation_l,~,~,~]=MarkovChainMoments(z_grid,pi_z); %Since l is exogenous, this will be it's eqm value 
%Note: Aiyagari (1994) actually then normalizes l by dividing it by
%Expectation_l (so that the resulting process has expectaion equal to 1
%(see Aiyagari (1993WP), footnote 33 pg 25-26).
%The following three lines do just this.
z_grid=z_grid./Expectation_l;
[Expectation_l,~,~,~]=MarkovChainMoments(z_grid,pi_z);


%In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

%% Grids
% Set grid for asset holdings
%Aiyagari uses 25 points, but with a piecewise-linear approx. see Aiyagari (1993WP, pg 28).
%His grid is not directly on k, but implicitly his k grid runs from zero up
%to k_max, where k_max is given by f(k_max,1)=delta*k_max
%k_max=delta^(1/(alpha-1));
%Doing this k_max is slightly less than 10*K_ss. But if I use this as the upper limit on
%the grid the results are wrong (in that increasing this to 15 or 20*K_ss
%gives different results). It might be that Aiyagari gets away with this
%due to his use of piecewise-linear approximation (ie. that policy fn is
%almost linear in this region anyway).
nk1=floor(n_k/3); nk2=floor(n_k/3); nk3=n_k-nk1-nk2;
k_grid=sort([linspace(0,K_ss,nk1),linspace(K_ss+0.0001,3*K_ss,nk2),linspace(3*K_ss+0.0001,15*K_ss,nk3)]');

%Set grid for interest rates (Aiyagari proves that with idiosyncratic
%uncertainty, the eqm interest rate is limited above by it's steady state value
%without idiosyncratic uncertainty, that is that r<r_ss).
%r_grid=linspace(0,r_ss,n_p)';
r_grid=sort([linspace(-Params.delta+0.0001,-0.0001,floor(1/3*n_p)),linspace(0,r_ss,ceil(2/3*n_p))])';

%Bring model into the notational conventions used by the toolkit
d_grid=0; %There is no d variable
a_grid=k_grid;
%pi_z;
%z_grid
p_grid=r_grid;


n_d=0;
n_a=n_k;
%n_z

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_z (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val) a_val; %We just want the aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
%Note: length(AggVars) is number_d_vars+number_a_vars and length(p) is number_p_vars
GeneralEqmEqnParamNames(1).Names={'alpha','delta'};
GeneralEqmEqn_1 = @(AggVars,GEprices,alpha,delta) GEprices-(alpha*(AggVars^(alpha-1))*(Expectation_l^(1-alpha))-delta); %The requirement that the interest rate corresponds to the agg capital level
GeneralEqmEqns={GeneralEqmEqn_1};

disp('sizes')
n_a
n_z
n_p


%% 
% Test value for prices, just to check codes are working
Params.r=0.04; %0.0485

DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, z_val,alpha,delta,mu,r) Aiyagari1994_ReturnFn(aprime_val, a_val, z_val,alpha,delta,mu,r);
ReturnFnParamNames={'alpha','delta','mu','r'}; %It is important that these are in same order as they appear in 'Aiyagari1994_ReturnFn'

% Test is currently commented out as not needed.
% tic;
% [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
% toc
% 
% tic;
% StationaryDist=SteadyState_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
% toc
% tic;
% AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid)
% toc
% MC=p_test - (alpha*(SSvalues_AggVars^(alpha-1))*(Expectation_l^(1-alpha))-delta);
% surf(StationaryDist)


%% Solve

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
heteroagentoptions.pgrid=p_grid;
[p_eqm,p_eqm_index, GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions,vfoptions);
% findeqmtime=toc
Params.r=p_eqm.r;

save ./SavedOutput/Aiyagari1994Market.mat p_eqm p_eqm_index GeneralEqmCondn

% Equilibrium wage
Params.w=(1-Params.alpha)*((Params.r+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

% %Now that we know what the equilibrium price is, lets calculate a bunch of
% %other things associated with the equilibrium
p_eqm_index
disp('Calculating various equilibrium objects')
[~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid, Parallel);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid)

eqm_MC=real(GeneralEqmConditions_Case1(AggVars,Params.r, GeneralEqmEqns, Params, GeneralEqmEqnParamNames)); % Note: is just giving GeneralEqmCondn(p_eqm_index)

save ./SavedOutput/Aiyagari1994SSObjects.mat p_eqm Policy StationaryDist

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}, and that L=1
% (exogeneous). Thus Y=K^{\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which
% equals delta*K. The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha})=delta*K^{1-\alpha}
aggsavingsrate=Params.delta*AggVars^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluate_Earnings = @(aprime_val,a_val,z_val,w) w*z_val;
FnsToEvaluateParamNames(2).Names={'r','w'};
FnsToEvaluate_Income = @(aprime_val,a_val,z_val,r,w) w*z_val+(1+r)*a_val;
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluate_Wealth = @(aprime_val,a_val,z_val) a_val;
FnsToEvaluateIneq={FnsToEvaluate_Earnings, FnsToEvaluate_Income, FnsToEvaluate_Wealth};
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluateIneq, Params,FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid);

% 3.5 The Distributions of Earnings and Wealth
%  Gini for Earnings
EarningsGini=Gini_from_LorenzCurve(LorenzCurves(1,:));
IncomeGini=Gini_from_LorenzCurve(LorenzCurves(2,:));
WealthGini=Gini_from_LorenzCurve(LorenzCurves(3,:));

% Calculate inverted Pareto coeff, b, from the top income shares as b=1/[log(S1%/S0.1%)/log(10)] (formula taken from Excel download of WTID database)
% No longer used: Calculate Pareto coeff from Gini as alpha=(1+1/G)/2; ( http://en.wikipedia.org/wiki/Pareto_distribution#Lorenz_curve_and_Gini_coefficient)
% Recalculte Lorenz curves, now with 1000 points
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluateIneq, Params,FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,[],1000);
EarningsParetoCoeff=1/((log(LorenzCurves(1,990))/log(LorenzCurves(1,999)))/log(10)); %(1+1/EarningsGini)/2;
IncomeParetoCoeff=1/((log(LorenzCurves(2,990))/log(LorenzCurves(2,999)))/log(10)); %(1+1/IncomeGini)/2;
WealthParetoCoeff=1/((log(LorenzCurves(3,990))/log(LorenzCurves(3,999)))/log(10)); %(1+1/WealthGini)/2;


%% Display some output about the solution
GraphString=['Aiyagari1994_MarketClearance_mu', num2str(Params.mu), 'rho', num2str(Params.rho), 'sigma', num2str(Params.sigma)];
GraphString=strrep(GraphString, '.', '_');
% Create a graph displaying the market clearance
fig1=figure(1);
plot(p_grid,GeneralEqmCondn, 'x',p_grid, zeros(n_p,1),p_eqm.r, 0, 'r o')
title(['Market clearance: mu=', num2str(Params.mu), ' rho=', num2str(Params.rho), ' sigma=', num2str(Params.sigma)],'FontSize',18);
xlabel('p','FontSize',16); ylabel('tilde(p)-p','FontSize',16);
set(fig1, 'Color', 'white');     % white bckgr
set(fig1, 'Unit', 'centimeters');  % set unit to cm
set(fig1,'position',[0 0 20 10]);  % set size
saveas(fig1,['./SavedOutput/Graphs/',GraphString,'.png'])

%plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values sigma=%.2f, mu=%.2f, rho=%.2f \n', [Params.sigma,Params.mu,Params.rho])
fprintf('The table 1 elements are sigma=%.4f, rho=%.4f \n',[sqrt(z_variance), z_corr])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm.r*100)
fprintf('The equilibrium value of the aggregate savings rate is r=%.4f \n', aggsavingsrate)
%fprintf('Time required to find the eqm was %.4f seconds \n',findeqmtime)

%% Outputs of the function
% Table 1: sqrt(z_variance), z_corr 
% Table 2: p_eqm*100, aggsavingsrate
% Inequality: EarningsGini, IncomeGini, WealthGini, EarningsParetoCoeff,
%  IncomeParetoCoeff, WealthParetoCoeff
% Share of top ten percent for all three?

disp('here1')
whos z_variance z_corr p_eqm aggsavingsrate EarningsGini IncomeGini WealthGini EarningsParetoCoeff IncomeParetoCoeff WealthParetoCoeff
disp('here2')
OutputVector=[sqrt(z_variance), z_corr, p_eqm.r*100, aggsavingsrate*100, EarningsGini, IncomeGini, WealthGini, EarningsParetoCoeff,IncomeParetoCoeff, WealthParetoCoeff];
disp('here3')
% Move OutputVector from GPU to CPU before returning it
OutputVector=gather(OutputVector);

end




