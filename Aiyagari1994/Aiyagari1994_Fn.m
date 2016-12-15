function OutputVector=Aiyagari1994_Fn(n_k,n_s,n_p,Params,Parallel, tauchenoptions, mcmomentsoptions, vfoptions, simoptions)
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
[s_grid, pi_s]=TauchenMethod(0,(Params.sigma^2)*(1-Params.rho^2),Params.rho,n_s,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)


%[s_grid, pi_s]=TauchenMethod_Param(0,(sigma^2)*(1-rho^2),rho,7,3); %This is the process for ln(l). Aiyagari uses 7 states, 
[s_mean,s_variance,s_corr,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions);
s_grid=exp(s_grid);
%Get some info on the markov process
[Expectation_l,~,~,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions); %Since l is exogenous, this will be it's eqm value 
%Note: Aiyagari (1994) actually then normalizes l by dividing it by
%Expectation_l (so that the resulting process has expectaion equal to 1
%(see Aiyagari (1993WP), footnote 33 pg 25-26).
%The following three lines do just this.
s_grid=s_grid./Expectation_l;
[Expectation_l,~,~,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions);


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
%pi_s;alpha, delta, gamma
%s_grid
p_grid=r_grid;


n_d=0;
n_a=n_k;
%n_s

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
SSvalueParamNames={};
SSvaluesFn_1 = @(aprime_val,a_val,s_val,p_val) a_val; %We just want the aggregate assets (which is this periods state)
SSvaluesFn={SSvaluesFn_1};

%Now define the functions for the Market Clearance conditions
    %Should be written as LHS of market clearance eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the market is to clearing.
%Note: length(AggVars) is number_d_vars+number_a_vars and length(p) is number_p_vars
MarketPriceParamNames={'alpha','delta'};
MarketPriceEqn_1 = @(AggVars,p,params) p-(params(1)*(AggVars^(params(1)-1))*(Expectation_l^(1-params(1)))-params(2)); %The requirement that the interest rate corresponds to the agg capital level
MarketPriceEqns={MarketPriceEqn_1};

disp('sizes')
n_a
n_s
n_p


%% 
% V0=ones([n_a,n_s],'gpuArray'); %(a,s)

% Test value for prices, just to check codes are working
Params.r=0.04; %0.0485
% n_z=n_s; z_grid=s_grid; pi_z=pi_s;

DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val,alpha,delta,mu,r) Aiyagari1994_ReturnFn(aprime_val, a_val, s_val,alpha,delta,mu,r);
ReturnFnParamNames={'alpha','delta','mu','r'}; %It is important that these are in same order as they appear in 'Aiyagari1994_ReturnFn'

% Test is currently commented out as not needed.
% tic;
% [V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, beta, ReturnFn,vfoptions,ReturnFnParams);
% toc
% 
% tic;
% StationaryDist=SteadyState_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
% toc
% tic;
% SSvalues_AggVars=SSvalues_AggVars_Case1(StationaryDist, Policy, SSvaluesFn, n_d, n_a, n_s, d_grid, a_grid,s_grid,pi_s,p_test, Parallel);
% toc
% MC=p_test - (alpha*(SSvalues_AggVars^(alpha-1))*(Expectation_l^(1-alpha))-delta);
% surf(StationaryDist)


%% Solve

V0=ones(n_a,n_s,'gpuArray'); %(a,s)
%Use the toolkit to find the equilibrium price index
PriceParamNames={'r'};

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
heteroagentoptions.pgrid=p_grid;
[p_eqm,p_eqm_index, MarketClearance]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_s, n_p, pi_s, d_grid, a_grid, s_grid, ReturnFn, SSvaluesFn, MarketPriceEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, MarketPriceParamNames, PriceParamNames,heteroagentoptions, simoptions, vfoptions);
% findeqmtime=toc
save ./SavedOutput/Aiyagari1994Market.mat p_eqm p_eqm_index MarketClearance

% Equilibrium wage
Params.w=(1-Params.alpha)*((p_eqm+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

% %Now that we know what the equilibrium price is, lets calculate a bunch of
% %other things associated with the equilibrium
p_eqm_index
disp('Calculating various equilibrium objects')
Params.r=p_eqm;
[~,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);


% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_s,pi_s, simoptions);

SSvalues_AggVars=SSvalues_AggVars_Case1(StationaryDist, Policy, SSvaluesFn,Params, SSvalueParamNames,n_d, n_a, n_s, d_grid, a_grid,s_grid,pi_s,p_eqm, Parallel)

eqm_MC=real(MarketClearance_Case1(SSvalues_AggVars,p_eqm,MarketPriceEqns, Params, MarketPriceParamNames));
save ./SavedOutput/Aiyagari1994SSObjects.mat p_eqm Policy StationaryDist

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}, and that L=1
% (exogeneous). Thus Y=K^{\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which
% equals delta*K. The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha})=delta*K^{1-\alpha}
aggsavingsrate=Params.delta*SSvalues_AggVars^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
SSvalueParamNames={'w'};
%  @(d_val,aprime_val,a_val,s_val,pi_s,p_val,param)
SSvalue_Earnings = @(aprime_val,a_val,s_val,p_val,param) param*s_val;
SSvalue_Income = @(aprime_val,a_val,s_val,p_val,param) param*s_val+(1+p_val)*a_val;
SSvalue_Wealth = @(aprime_val,a_val,s_val,p_val,param) a_val;
SSvaluesFnIneq={SSvalue_Earnings, SSvalue_Income, SSvalue_Wealth};
SSvalues_LorenzCurves=SSvalues_LorenzCurve_Case1(StationaryDist, Policy, SSvaluesFnIneq, Params,SSvalueParamNames, n_d, n_a, n_s, d_grid, a_grid, s_grid, pi_s,p_eqm,Parallel);

% 3.5 The Distributions of Earnings and Wealth
%  Gini for Earnings
EarningsGini=Gini_from_LorenzCurve(SSvalues_LorenzCurves(1,:));
IncomeGini=Gini_from_LorenzCurve(SSvalues_LorenzCurves(2,:));
WealthGini=Gini_from_LorenzCurve(SSvalues_LorenzCurves(3,:));

% Calculate inverted Pareto coeff, b, from the top income shares as b=1/[log(S1%/S0.1%)/log(10)] (formula taken from Excel download of WTID database)
% No longer used: Calculate Pareto coeff from Gini as alpha=(1+1/G)/2; ( http://en.wikipedia.org/wiki/Pareto_distribution#Lorenz_curve_and_Gini_coefficient)
% Recalculte Lorenz curves, now with 1000 points
SSvalues_LorenzCurves=SSvalues_LorenzCurve_Case1(StationaryDist, Policy, SSvaluesFnIneq, Params,SSvalueParamNames, n_d, n_a, n_s, d_grid, a_grid, s_grid, pi_s,p_eqm,Parallel,1000);
EarningsParetoCoeff=1/((log(SSvalues_LorenzCurves(1,990))/log(SSvalues_LorenzCurves(1,999)))/log(10)); %(1+1/EarningsGini)/2;
IncomeParetoCoeff=1/((log(SSvalues_LorenzCurves(2,990))/log(SSvalues_LorenzCurves(2,999)))/log(10)); %(1+1/IncomeGini)/2;
WealthParetoCoeff=1/((log(SSvalues_LorenzCurves(3,990))/log(SSvalues_LorenzCurves(3,999)))/log(10)); %(1+1/WealthGini)/2;


%% Display some output about the solution
GraphString=['Aiyagari1994_MarketClearance_mu', num2str(Params.mu), 'rho', num2str(Params.rho), 'sigma', num2str(Params.sigma)];
GraphString=strrep(GraphString, '.', '_');
% Create a graph displaying the market clearance
fig1=figure(1);
plot(p_grid,MarketClearance, 'x',p_grid, zeros(n_p,1),p_eqm, 0, 'r o')
title(['Market clearance: mu=', num2str(Params.mu), ' rho=', num2str(Params.rho), ' sigma=', num2str(Params.sigma)],'FontSize',18);
xlabel('p','FontSize',16); ylabel('tilde(p)-p','FontSize',16);
set(fig1, 'Color', 'white');     % white bckgr
set(fig1, 'Unit', 'centimeters');  % set unit to cm
set(fig1,'position',[0 0 20 10]);  % set size
export_fig(fig1, ...            % figure handle
    ['./SavedOutput/Graphs/',GraphString,'.pdf'],... % name of output file without extension   %  '-painters', ...            % renderer
    '-pdf', ...                 % file format
    '-r300' );                  % resolution

%plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values sigma=%.2f, mu=%.2f, rho=%.2f \n', [Params.sigma,Params.mu,Params.rho])
fprintf('The table 1 elements are sigma=%.4f, rho=%.4f \n',[sqrt(s_variance), s_corr])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm*100)
fprintf('The equilibrium value of the aggregate savings rate is r=%.4f \n', aggsavingsrate)
%fprintf('Time required to find the eqm was %.4f seconds \n',findeqmtime)

%% Outputs of the function
% Table 1: sqrt(s_variance), s_corr 
% Table 2: p_eqm*100, aggsavingsrate
% Inequality: EarningsGini, IncomeGini, WealthGini, EarningsParetoCoeff,
%  IncomeParetoCoeff, WealthParetoCoeff
% Share of top ten percent for all three?

disp('here1')
whos s_variance s_corr p_eqm aggsavingsrate EarningsGini IncomeGini WealthGini EarningsParetoCoeff IncomeParetoCoeff WealthParetoCoeff
disp('here2')
OutputVector=[sqrt(s_variance), s_corr, p_eqm*100, aggsavingsrate*100, EarningsGini, IncomeGini, WealthGini, EarningsParetoCoeff,IncomeParetoCoeff, WealthParetoCoeff];
disp('here3')
% Move OutputVector from GPU to CPU before returning it
OutputVector=gather(OutputVector);

end




