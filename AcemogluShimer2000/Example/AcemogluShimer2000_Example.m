%% Acemoglu & Shimer (2000) - Productivity gains from unemployment insurance
%
% This example solves the baseline model of Acemoglu & Shimer (2000) and computes a few model statistics.
% It is intended for anyone who just wants to see how the model is set up in VFI Toolkit.
% (A full replication is at: )
%
% The model: workers are infinitely lived, a period is a week, and they save in a riskless bond (a>=0) at an
% exogenous interest rate. A worker is in one of four employment states: a good job (high wage), a bad job
% (low wage), unemployed and receiving unemployment insurance (UI), or unemployed without UI. Unemployed workers
% choose how hard to search and whether to look for a good job (hard to find) or a bad job (easy to find).
% Jobs end with exogenous separation probability, after which the worker is insured unemployed; insured workers 
% randomly lose eligibility. UI is paid for by a tax on labour income, and the tax rate is set to balance the 
% government budget.
% The point of the paper: UI lets workers take the risk of searching for the good jobs.
%
% Everything is normalized by productivity growth Gamma, so next period assets are aprime=((1+r)/Gamma)*(a+y-c),
% and the discount factor is beta*Gamma^(1-theta).

%% How does VFI Toolkit think about this?
%
% Two decision variables: h and s, search hours, and a binary choice
%   (unemployed: s=1 look for a good job, s=0 look for a bad job; employed: s=1 quit, s=0 stay)
% One endogenous state variable: a, assets
% One semi-exogenous state: the employment state (good job, bad job, insured unemployed, uninsured unemployed)
% No (standard) exogenous states
%
% The employment state is a "semi-exogenous state": it is a Markov process, but its transition probabilities
% depend on the decision variables. Here both decision variables (h and s) influence the transitions:
% how likely you are to find a job depends on how many hours you search, h, and what kind of job you look for, s.
%
% Employed workers also choose hours of work, but these do not affect any transitions, so rather than make them
% a decision variable we solve for them inside the return function using the first-order condition (see the
% return function). This means work hours are exact, and not restricted to a grid.

%% Record everything printed to the screen in AcemogluShimer2000_Example_Diary.txt (fresh each run)
diary off
if exist('./AcemogluShimer2000_Example_Diary.txt','file'), delete('./AcemogluShimer2000_Example_Diary.txt'); end
diary ./AcemogluShimer2000_Example_Diary.txt

%% Begin setting up to use VFI Toolkit to solve

% Grid sizes to use
n_d=[101,2]; % search hours h, binary choice s
n_a=396; % assets
n_semiz=4; % employment state
n_z=0; % no (standard) exogenous state

%% Parameters (AS2000, Table 1)
% Interest rate and growth
Params.r=0.0006; % weekly interest rate (so the gross interest rate is R=1+r; AS2000 use R, this example uses r)
Params.Gamma=1.0003; % gross weekly growth rate of productivity
% Jobs
Params.wb=1; % bad job wage
Params.deltab=0.005; % bad job layoff rate
Params.mb=0.5; % bad job matching rate (prob of finding one is mb*h)
Params.wg=1.3; % good job wage
Params.deltag=0.005; % good job layoff rate
Params.mg=0.05; % good job matching rate (prob of finding one is mg*h)
% Preferences
Params.beta=0.998; % discount factor
Params.theta=4; % relative risk aversion
Params.hbar=1; % time endowment
Params.eta=1/6; % relative value of leisure
% Unemployment insurance
Params.uib=0.4; % unemployment benefit
Params.phi=0.04; % prob of losing eligibility each week (so benefits last 25 weeks on average)
Params.tau=0.022; % tax rate on labour income, initial guess as it will be determined in general eqm

% Growth changes the discount factor: it is beta*Gamma^(1-theta)
Params.Gammadisc=Params.Gamma^(1-Params.theta);

% A parameter that is not in the paper:
% A tiny income for the unemployed. Without it an uninsured worker with zero assets cannot consume anything,
% which is -Inf utility. Because you otherwise have 0 income in the uninsured&unemployed state, and because you
% can be there for any finite length of time, it gets numerically messy. 1e-7 of a wage has no effect on the results
% and resolves this numerical problem.
Params.uincfloor=10^(-7);

%% Grids
% Assets: the paper uses the grid {1.01^n-1}, which puts most points near zero
a_grid=1.01.^(0:n_a-1)'-1; % goes from 0 to about 50

% Decision variables
h_grid=linspace(0.5,1,n_d(1))'; % search hours
s_grid=[0;1]; % binary choice
d_grid=[h_grid; s_grid]; % stacked column vector, as always in VFI Toolkit

% There are no standard exogenous states
z_grid=[];
pi_z=[];

%% Semi-exogenous state: the employment state
% 1=good job, 2=bad job, 3=insured unemployed, 4=uninsured unemployed
vfoptions.n_semiz=n_semiz;
vfoptions.semiz_grid=(1:1:4)';
% The transition probabilities of the semi-exogenous state, as a function of (semiz,semizprime,d) and parameters.
% See AcemogluShimer2000_SemiExoStateFn.m for the probabilities themselves.
vfoptions.SemiExoStateFn=@(semiz,semizprime,h,s,deltag,deltab,mg,mb,phi) AcemogluShimer2000_SemiExoStateFn(semiz,semizprime,h,s,deltag,deltab,mg,mb,phi);
% By default only the last decision variable influences the semi-exogenous state, here it is the last two
vfoptions.l_dsemiz=2;

% We also need to tell simoptions about the semi-exogenous state
simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions.l_dsemiz=vfoptions.l_dsemiz;
simoptions.d_grid=d_grid; % the transitions depend on d, so simoptions needs d_grid

% Use the grid interpolation layer: aprime is chosen from 10 evenly-spaced points between each pair of
% a_grid points. This matters here: workers save only a little each week, less than the gap between grid points.
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=vfoptions.gridinterplayer; % grid interpolation layer must also be set in simoptions
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Return function
DiscountFactorParamNames={'beta','Gammadisc'}; % multiplied together to give the discount factor

% The inputs to the ReturnFn are first the decision variables (h,s), then next period endogenous state (aprime),
% then this period endogenous state (a), then the semi-exogenous state (semiz), and then the parameters.
% See AcemogluShimer2000_ReturnFn.m
ReturnFn=@(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,r,Gamma,hbar,eta,theta) ...
    AcemogluShimer2000_ReturnFn(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,r,Gamma,hbar,eta,theta);


%% We don't need to solve the value function, agent distribution, nor model stats before solving the stationary
% general eqm, but it is typically good to do so as that way you see they all work before putting it all together

%% Value fn
disp('Solve the value function')
tic;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
toc

%% Stationary distribution of agents
disp('Solve the stationary distribution')
tic;
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
toc

%% Model statistics
% FnsToEvaluate take the same first inputs as the ReturnFn (h,s,aprime,a,semiz), followed by any parameters.
% Hours are work hours when employed and search hours when unemployed; work hours of the employed are computed
% by AcemogluShimer2000_Hours.m in the same way as inside the ReturnFn.
FnsToEvaluate.TaxRevenue=@(h,s,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta) tau*AcemogluShimer2000_Hours(h,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta)*(wg*(semiz==1)+wb*(semiz==2));
FnsToEvaluate.UIspending=@(h,s,aprime,a,semiz,uib) uib*(semiz==3); % only the insured unemployed get benefits

AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% At our initial guess for the tax rate the budget does not balance
fprintf('At tau=%1.3f, tax revenue is %1.4f and UI spending is %1.4f \n',Params.tau,AllStats.TaxRevenue.Mean,AllStats.UIspending.Mean)

%% Solve the stationary general equilibrium
% The 'general eqm' here is just that the tax rate balances the government budget
heteroagentoptions.verbose=1; % give feedback on progress

GEPriceParamNames={'tau'}; % we want to find the tax rate

% A general eqm eqn is written so that it equals zero in general eqm
% Inputs to the GeneralEqmEqns can be any parameters or any AggVars (the Mean in AllStats)
GeneralEqmEqns.govbudget=@(TaxRevenue,UIspending) TaxRevenue-UIspending;

% The budget balances at more than one tax rate (a 'Laffer curve'); the paper uses the lowest one (fn 12),
% so our initial guess is a low tax rate
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

% The GEcondns are essentially zero
GEcondns

% p_eqm contains the equilibrium tax rate, we need to put it into the Params structure
Params.tau=p_eqm.tau;

%% Look at the equilibrium
% Add some other model stats we want to see
FnsToEvaluate.assets=@(h,s,aprime,a,semiz) a;
FnsToEvaluate.goodjob=@(h,s,aprime,a,semiz) (semiz==1);
FnsToEvaluate.badjob=@(h,s,aprime,a,semiz) (semiz==2);
FnsToEvaluate.unemployed=@(h,s,aprime,a,semiz) (semiz>=3);
FnsToEvaluate.searchhoursinsured=@(h,s,aprime,a,semiz) h*(semiz==3);
FnsToEvaluate.insured=@(h,s,aprime,a,semiz) (semiz==3);
FnsToEvaluate.wage=@(h,s,aprime,a,semiz,wg,wb) wg*(semiz==1)+wb*(semiz==2);
FnsToEvaluate.output=@(h,s,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta) AcemogluShimer2000_Hours(h,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta)*(wg*(semiz==1)+wb*(semiz==2));
FnsToEvaluate.consumption=@(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,r,Gamma,hbar,eta) a+(1-tau)*AcemogluShimer2000_Hours(h,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta)*(wg*(semiz==1)+wb*(semiz==2))+max(uib,uincfloor)*(semiz==3)+uincfloor*(semiz==4)-aprime*Gamma/(1+r);
FnsToEvaluate.hours=@(h,s,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta) AcemogluShimer2000_Hours(h,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta);

% We first solve for the value fn and agent dist at the equilibrium tax rate
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

% The FnsToEvaluate evaluated at the policy, on the grid, so we can plot them
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
ConsumptionPolicy=reshape(ValuesOnGrid.consumption,[n_a,n_semiz]); % assets-by-employment state
HoursPolicy=reshape(ValuesOnGrid.hours,[n_a,n_semiz]);
StationaryDist=reshape(StationaryDist,[n_a,n_semiz]);

%% Create some output based on the model results
fprintf(' \n')
fprintf('In the stationary general equilibrium (paper''s Table 2, column I, in brackets),... \n')
fprintf('Tax rate=%1.3f (0.022) \n', Params.tau)
fprintf('Mean assets=%1.1f (16.5) \n', AllStats.assets.Mean)
fprintf('Fraction of good jobs=%1.1f%% (51.5%%) \n', 100*AllStats.goodjob.Mean/(AllStats.goodjob.Mean+AllStats.badjob.Mean))
fprintf('Unemployment rate=%1.1f%% (6.3%%) \n', 100*AllStats.unemployed.Mean)
fprintf('Search hours by insured=%1.3f (0.816) \n', AllStats.searchhoursinsured.Mean/AllStats.insured.Mean)
fprintf('Average wage=%1.3f (1.155) \n', AllStats.wage.Mean/(1-AllStats.unemployed.Mean))
fprintf('Output=%1.3f (0.931) \n', AllStats.output.Mean)
fprintf('Gini of assets=%1.2f (not reported in paper) \n', AllStats.assets.Gini)

% Consumption policy in each employment state (AS2000 Figure 3)
plotind=(a_grid<=40);
StateNames={'Good job','Bad job','Unemployed, insured','Unemployed, uninsured'};
figure(1)
for ii=1:4
    subplot(2,2,ii); plot(a_grid(plotind),ConsumptionPolicy(plotind,ii))
    title(['Consumption: ',StateNames{ii}])
    xlabel('assets')
end
% Note: bad-job workers quit when their assets exceed about 19, so the bad-job consumption policy above 19 is never used

% Hours of work or search in each employment state (AS2000 Figure 4)
figure(2)
for ii=1:4
    subplot(2,2,ii); plot(a_grid(plotind),HoursPolicy(plotind,ii))
    title(['Hours: ',StateNames{ii}])
    xlabel('assets')
end

% Distribution of assets (AS2000 Figure 5 plots the pdf, we plot the cdf)
figure(3)
plot(a_grid(plotind),cumsum(sum(StationaryDist(plotind,:),2)))
title('cdf of workers over assets') % cdf=cumulative distribution function
xlabel('assets')

diary off
