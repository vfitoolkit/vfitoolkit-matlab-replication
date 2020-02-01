% Replication of Hopenhayn & Rogerson (1993) - Job Turnover and Policy Evaluation: A General Equilibrium Analysis

% QUESTIONS:
% What is the distribution of new agents entering (upsilon): answered (roughly, Marten Floden reverse engineered it, see below)
% Value of Tauchen q? (Paper does not specifically say tauchens, but does say discretize AR(1) with twenty grid points.)

vfoptions.parallel=0
simoptions.parallel=0

ImposeFootnote5=1 % If equal to one, then impose footnote 5 from Hopenhayn & Rogerson (1993):
    % Awkwardly, Hopenhayn & Rogerson (1993) state in footnote 5 on page 922
    % that "Note that we are assuming that a new entrant bears only the fixed
    % cost of entry and does not pay the cost cf". This is in direct conflict
    % with the statement on pg 919 at "Once this [entry] cost has been paid,
    % the entrant is in the same position as an incumbent that has chosen to
    % remain in the productive sector and had zero employees in the previous
    % period". It appears from the replication results that the footnote is the
    % final word on this, the following few lines implement footnote 5 under
    % the assumption that only a new entrant firm will have zero employees 'last period'
    % I implement this by adding an extra point to a_grid which takes value
    % 10^6. I then use ReturnFn to make it so that existing firms cannot
    % choose to move to this point (cost of doing so is -Inf; note this is also true of 
    % 'staying' at this point) and adjust the calculation of all model
    % statistics. This 'trick' to impose footnote 5 only works because we are only interested in
    % stationary equilibrium, it would not work for, e.g., a transition path.
Original_agrid=1; % When imposing footnote 5, can also impose a log-space 250 point grid on employment, as used by HR1993. (Is ignored when not imposing footnote 5)

ChrisEdmondCalib=0 % If equal to one, use alternative calibration/model from Chris Edmonds lecture notes, slide 15: http://www.chrisedmond.net/phd2014/90065_lecture4.pdf
% Note that Chris Edmonds calibration uses a different n_z
% Note that Hopenhayn & Rogerson (1993) normalize p to 1 and then choose ce in eqm, while Chris Edmonds instead sets ce and determines p in eqm.
if ChrisEdmondCalib==1
    ImposeFootnote5=0; % Chris Edmond (sensibly) does not follow Footnote 5.
end

% Comment: Hopenhayn & Rogerson (1993) use quite a small grid on 'employees' together with 
% pure discretization as a solution method. This has the effect that while
% the model looks like there is just a 'linear cost of adjustment' it acts,
% in their implementation, more like '(proportional to last period employment) fixed cost of adjustment plus linear cost of
% adjustment'. This is not true in this replication where many more grid
% points are used. Hence the quantitative differences in some results
% between this replication and the original.

%%

% Preferences: log(c)-AN (HR1993, pg 927. Pg 923 is typo (gets case of A incorrect)
Params.beta=0.8;
Params.A=0.6;

% Production fn
Params.alpha=0.64; % HR1993 call this theta
Params.cf=12; % Hopenhayn & Rogerson (1993) do not report, but Martin Flodén figures out the following (pg 5): http://martinfloden.net/files/macrolab.pdf
Params.ce=40; % Initial guess, this will be determined by general equilibrium condition.

% Firing cost (policy)
Params.tau=0; % Baseline

Params.p=1; % "value of ce is chosen so that [Free entry condition] is satisfied with p=1, pg 930, Hopenhayn & Rogerson (1993)
Params.w=1; % Normalization

% Exogenous AR(1) process on (log) productivity: logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
Params.rho=0.93; % Hopenhayn & Rogerson (1993)
Params.sigma_logz=sqrt(0.53); % Hopenhayn & Rogerson (1993)
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.078; % Hopenhayn & Rogerson (1993) do not report, but Martin Flodén figures out the following (pg 5): http://martinfloden.net/files/macrolab.pdf

tauchenoptions.parallel=vfoptions.parallel;
n_z=20; % I here call z, what Hopenhayn & Rogerson (1993) call s. The choice of n_z=20 follows them.
Params.q=4; % Hopenhayn & Rogerson (1993) do not report (based on Table 4 is seems something around q=4 is used, otherwise don't get values of z anywhere near as high as 27.3. (HR1993 have typo and call the column 'log(s)' when it should be 's') 
[z_grid, pi_z]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_z,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
z_grid=exp(z_grid);

if ChrisEdmondCalib==1
    % Parameters not listed in this if-statement just take same value as in Hopenhayn & Rogerson (1993).
    Params.beta=0.8;
    Params.alpha=2/3;
    Params.cf=20;
    
    Params.ce=39.659; % Exact value from CE codes (slides say 40)
    
    Params.rho=0.9;
    Params.sigma_epsilon=0.2026; % Exact value from CE codes (slides say 0.2)
    Params.logzbar=1.3855; % Exact value from CE codes (slides say 1.39)
    
    n_z=33; % The choice of n_z=33 follows Chris Edmonds lecture notes.
    Params.q=6; % Reverse engineered from Chris Edmonds code. He uses an implementation based on Tauchen-Hussey method that uses a quadrature routine from Miranda & Fackler which chooses q to make sure var(epsilon) takes the correct value.
    [z_grid, pi_z]=TauchenMethod((1-Params.rho)*Params.logzbar,Params.sigma_epsilon^2,Params.rho,n_z,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
    z_grid=exp(z_grid);
    
    % Edmonds calls 'theta' what Hopenhayn & Rogerson (1993) call A
    Params.A=100;
end

%%
n_a=401; % Hopenhayn & Rogerson (1993) use 250, this is intended as an approximation of continuous variable, so I use more.
a_grid=[linspace(0,100,101),((logspace(0,pi,n_a-101)-1)/(pi-1))*(5000-101)+101]'; % Put points on each of 0 to 100, then logarithmically spaced from 101 to 5000 (5000 is as reported by Hopenhayn & Rogerson (1993)).
if ImposeFootnote5==1
    a_grid=[linspace(0,100,101),((logspace(0,pi,n_a-101-1)-1)/(pi-1))*(5000-101)+101,10^6]'; % One less point in standard grid, instead add the 10^6 point to keep track of the new entrants.
    if Original_agrid==1
        n_a=250;
        a_grid=[exp(linspace(0,log(5001),n_a-1))-1,10^6]';
    end
end
n_d=0; % None.
d_grid=[];

% plot(1:1:n_a-1,a_grid(1:end-1))

%% Value Function Problem (with endogenous exit)
% Entry is irrelevant to the value function problem. Exit is relevant, and treated differently based on whether it is endogenous or exogenous.

% For exogenous exit, you would simply need to include the 'conditional surivival probability' as another 'DiscountFactorParamNames'
DiscountFactorParamNames={'beta'};
% For endogenous exit, you need to set
vfoptions.endogenousexit=1;
% We also need to create 'vfoptions.ReturnToExitFn' (and 'vfoptions.ReturnToExitFnParamNames'), as below.

ReturnFn=@(aprime_val, a_val, z_val, p, alpha, tau, cf) HopenhaynRogerson1993_ReturnFn(aprime_val, a_val, z_val, p, alpha, tau, cf);
ReturnFnParamNames={'p', 'alpha', 'tau', 'cf'}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnFn'
if ChrisEdmondCalib==1
    % ReturnFn is slightly different as the fixed cost of production is
    % denominated in dollars (so just -cf, not -p*cf inside return fn).
    ReturnFn=@(aprime_val, a_val, z_val, p, alpha, tau, cf) HopenhaynRogerson1993_ReturnFn_ChrisEdmonds(aprime_val, a_val, z_val, p, alpha, tau, cf);
end


% For endogenous exit, also need to define the 'return to exit'.
vfoptions.ReturnToExitFn=@(a_val, z_val,tau) -tau*a_val; % the exit cost is the cost of firing all remaining employees (pg 919 of Hopenhayn & Rogerson, 1993)
% Remark: if a more complex 'return to exit' was desired it could be
% created in much the same way as the ReturnFn is created, but depends only
% on (a,z) variables (and any parameters passed using ReturnToExitFnParamNames).
% The following commented out line provides an 'example'
% ReturnToExitFn=@(a_val, s_val) HopenhaynRogerson1993_ReturnToExitFn(a_val, s_val);
vfoptions.ReturnToExitFnParamNames={'tau'}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnToExitFn'

if ImposeFootnote5==1
    vfoptions.ReturnToExitFn=@(a_val, z_val,tau) -tau*a_val*(a_val~=10^6); % the exit cost is the cost of firing all remaining employees (pg 919 of Hopenhayn & Rogerson, 1993)
end

% Check that everything is working so far by solving the value function
if vfoptions.parallel==2
    V0=zeros(n_a,n_z,'gpuArray');
else
    V0=zeros(n_a,n_z);
end
[V,Policy,ExitPolicy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);

% When tau=0 there is a cut-off value of z below which all firms exit, and which is independent of n.
% This can be easily seen graphing the whole of the ExitPolicy, which takes
% a value of 1 for exit, and 0 when not choosing to exit.
figure(1)
surf(ExitPolicy)
title('Exit Decision (value of 1 indicates exit)')

% Can be informative to, e.g., change value of tau and look at how the
% following are changed as a result
% surf(shiftdim(Policy,1))
% surf(V)

%% Stationary Distribution of Agents with entry and exit
% Both entry and exit matter for stationary distribution of agents. But
% whether they are endogenous or exogenous is irrelevant.
simoptions.agententryandexit=1;
% In principle whether exit or entry is exogenous is not needed info, but
% for algorithm turns out to be useful info
simoptions.endogenousexit=1;


% Distribution of new agents:
EntryExitParamNames.DistOfNewAgents={'upsilon'};

Params.upsilon=[1;zeros(n_a-1,1)]*[ones(1,floor(0.65*n_z)),zeros(1,n_z-floor(0.65*n_z))];
Params.upsilon=Params.upsilon/sum(sum(Params.upsilon)); % Normalize to be pdf
% "We found that a uniform distribution on the lower part of the interval in which realizations of s 
% lie produced a reasonable fit.", pg 930 of Hopenhayn & Rogerson (1993), but no mention of what constitutes 'lower part'.
% Martin Flodén concludes that roughly the bottom 2/3 (0.65 to be precise) is a good definiting of 'lower' (pg 5): http://martinfloden.net/files/macrolab.pdf
% Also, "Once this cost has been paid the entrant is in the same position as an incumbent that has chosen 
% to remain in the productive sector and had zero employees in the previous period.", pg 919 of Hopenhayn 
% & Rogerson (1993) tells us that all entrants have zero 'lagged employees' (which is the lowest grid point).
if ChrisEdmondCalib==1
    % Uses the stationary distribution of the AR(1) process on z, together
    % with zero assets, as the initial distribution of agents.
    
    % Compute the stationary distribution of this Markov (this could be done more easily/directly/analytically)
    pistar_z=ones(size(z_grid))/n_z; % Initial guess
    dist=1;
    while dist>10^(-9)
        pistar_z_old=pistar_z;
        pistar_z=(pi_z')*pistar_z;
        dist=max(abs(pistar_z-pistar_z_old));
    end
    
    Params.upsilon=[1;zeros(n_a-1,1)]*pistar_z';
end

if ImposeFootnote5==1
    % Need to put the new entrants into the 'special value' used to
    % indicate new entrants so they can receive special treatment as per footnote.
    Params.upsilon=[zeros(n_a-1,1);1]*[ones(1,floor(0.65*n_z)),zeros(1,n_z-floor(0.65*n_z))]; % 'special value' is last point in a_grid
    Params.upsilon=Params.upsilon/sum(sum(Params.upsilon)); % Normalize to be pdf
end


% Note: VFI Toolkit requires the DistOfNewAgents to be a pdf (so unit mass), and then uses
% the 'MassOfNewAgents' to understand how many there will be entering
% relative to existing agents. (MassOfExistingAgents is kept track of.)
Params.Ne=0.5; % The initial guess for the mass of existing agents is always 1, this is implicitly hardcoded into how VFI Toolkit commands are implemented.
EntryExitParamNames.MassOfNewAgents={'Ne'};
% Note: this is just an initial guess, as it will anyway need to be determined in general equilibrium.

% Exit can be given in one of two forms. As a function, or as a matrix.
% With endogenous entry, it is always best to give as following matrix.
% simoptions.CondlProbOfSurvival=1-ExitPolicy;
EntryExitParamNames.CondlProbOfSurvival={'zeta'};
Params.zeta=1-ExitPolicy;
% This conditional probability can be a state-dependent parameter, in which case this input would be a vector or matrix, etc.

% Check that everything is working so far by solving the simulation of agent distribution to get the stationary distribution.
simoptions % Show which options are being set
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

% Note: When using models, such as entry and exit, where the mass of agents is not equal to 1
% the toolkit will automatically keep track of distributions as StationaryDist.pdf and StationaryDist.mass

% Note: when tau=0, conditional on z there is a value of aprime that all
% firms will choose. This means that next period all firms will appear in
% just one of n_z values on a_grid (the pdf over on a_grid will appear as
% just n_z point masses; actually n_z+1 because of new entrants).
figure(2)
if ImposeFootnote5==0
    plot(a_grid, cumsum(sum(StationaryDist.pdf,2))) % Graph of marginal cdf over assets (will look stupid when using ImposeFootnote5=1)
else % ImposeFootnote5==1
    temp=sum(StationaryDist.pdf,2);
    temp2=temp(1:end-1); temp2(1)=temp(1)+temp(end);
    plot(a_grid(1:end-1), cumsum(temp2))
end
title('Stationary Distribution over lagged employment (sum/integral over z)')

% StationaryDist.pdf(1:20,16:19)

% plot([a_grid(1:end-1)], cumsum(sum(StationaryDist.pdf(1:end-1,:),2)),[a_grid(1:end-1)], cumsum(sum(StationaryDist1.pdf(1:end-1,:),2)),[a_grid(1:end-1)], cumsum(sum(StationaryDist2.pdf(1:end-1,:),2)))
% legend('0','1','2')

%% General Equilibrium with endogenous entry and endogenous exit.
% The entry decision imposes a general equilibrium condition. In practice it determines the price level p.
% Notice that this depends on the value function of existing firms and on the 'DistOfNewAgents'.
% If entry is exogenous it would not impose any general equilibrium condition.
% Notice that exit (whether endogenous or exogenous) does not involve any
% general equilibrium conditions.
% Since the stationary equilibrium commands already need to know all the
% vfoptions and simoptions, we don't need to add any further info to highlight that 
% this problem includes entry and exit as these
% already contain everything the stationary equilibrium command needs to know.

%Use the toolkit to find the equilibrium prices
GEPriceParamNames={'ce','Ne'};

FnsToEvaluateParamNames(1).Names={'alpha'};%,'p'};
% Note: With entry-exit the mass of the distribution of agents often
% matters. So it becomes an extra input arguement in all functions to be evaluated.
% FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val,agentmass,alpha,p) p*z_val*(aprime_val^alpha); % Total output
FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val,agentmass,alpha,p) z_val*(aprime_val^alpha); % Real output
FnsToEvaluate={FnsToEvaluateFn_1};

% Just to test: (note, is same command as usual, just need to include the optional extra inputs 'simoptions' and 'EntryExitParamNames' which contains all the needed info about entry/exit)
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

% The general equilibrium condition is that the EV^e-ce=0.
% This does not fit standard format for general equilibrium conditions.
heteroagentoptions.specialgeneqmcondn={0,'entry'};
% Certain kinds of general equilibrium conditions that are non-standard can
% be used via heteroagentoptions.specialgeneqmcondn
GeneralEqmEqnParamNames(1).Names={'p','A'};
% SHOULD THE FOLLOWING BE MODIFIED TO C=Y-ce*Ne? (currently AggVars is Y, while condition comes from Rep HH and is on C; so currently just C=Y)
GeneralEqmEqn_1 = @(AggVars,GEprices,p,A) A/AggVars-p; %The requirement that the price is determined by the demand eqn (or equivalently, can think of this as goods market clearance). You can derive it from FOCs of standard consumption-leisure problem [it is the -U_c/U_N=p/w condition you often see in household problems; remember normalize w=1]: max_{c,N} log(c)-AN s.t. pC=wN+T
GeneralEqmEqnParamNames(2).Names={'p','beta'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,p,beta) beta*EValueFn-p*GEprices(1); % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
if ImposeFootnote5==1
    GeneralEqmEqnParamNames(2).Names={'p'};
    GeneralEqmEqn_Entry = @(EValueFn,GEprices,p) EValueFn-p*GEprices(1); % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
end
% The entry condition looks slightly different to more standard @(EValueFn,p,params)
% This is because 'p' is the name of a parameter, and so have used 'GEprices'
% instead of my usual 'p' to refer to the general equilibrium prices (here 'ce' and 'Ne')
GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_Entry};
% Note that GeneralEqmEqn_Entry needed to be pointed out as special because
% it depends on the distribution of entrants and not the distribution of
% existing agents (all standard general eqm conditions involve the later).

if ChrisEdmondCalib==1
    % Calibrates ce, and instead treats p as a general eqm price to be determined.
    GEPriceParamNames={'p','Ne'};
    GeneralEqmEqnParamNames(1).Names={'A'};
    GeneralEqmEqn_1 = @(AggVars,GEprices,A) AggVars/A-GEprices(1); %The requirement that the price is determined by the demand eqn (or equivalently, can think of this as goods market clearance)
    GeneralEqmEqnParamNames(2).Names={'beta','ce'};
    GeneralEqmEqn_Entry = @(ValueFn,GEprices,beta,ce) beta*ValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
    % The general eqm conditions look slightly different to more standard @(ValueFn,p,params)
    % This is because 'p' is the name of a parameter, and so have used 'GEprices'
    % instead of my usual 'p' to refer to the general equilibrium prices (parameter 'p' is actually GEprices(1))
    GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_Entry};
    
    % Just to test: (note, is same command as usual, just need to include the optional extra inputs 'simoptions' and 'EntryExitParamNames' which contains all the needed info about entry/exit)
    AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
end


heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
% findeqmtime=toc
if ChrisEdmondCalib==0
    Params.ce=p_eqm.ce;
    Params.Ne=p_eqm.Ne;
else
    Params.p=p_eqm.p;
    Params.Ne=p_eqm.Ne;
end

%% Calculate some relevant things in eqm
[V,Policy,ExitPolicy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.zeta=1-ExitPolicy;
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

if ChrisEdmondCalib==0
    save ./SavedOutput/HopenhaynRogerson1993.mat Params V Policy ExitPolicy StationaryDist
else
    save ./SavedOutput/HopenhaynRogerson1993_ChrisEdmonds.mat Params V Policy ExitPolicy StationaryDist
end

% load ./SavedOutput/HopenhaynRogerson1993.mat Params V Policy ExitPolicy StationaryDist
% load ./SavedOutput/HopenhaynRogerson1993_ChrisEdmonds.mat Params V Policy ExitPolicy StationaryDist

%% Now that the stationary equilibrium has been found, replicate Table 2

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_Emp = @(aprime_val,a_val,z_val,AgentDistMass) aprime_val; % Employment
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_Hiring = @(aprime_val,a_val,z_val,AgentDistMass) (aprime_val-a_val)*(aprime_val>a_val); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluateFn_Firing = @(aprime_val,a_val,z_val,AgentDistMass) -(aprime_val-a_val)*(aprime_val<a_val); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
FnsToEvaluate={FnsToEvaluateFn_Emp, FnsToEvaluateFn_Hiring, FnsToEvaluateFn_Firing};

if ImposeFootnote5==1 % Need to deal with the 'special value' in a_grid for new entrants
    % Is simply a matter of replacing a_val with a_val*(a_val~=10^6)
    FnsToEvaluateFn_Hiring = @(aprime_val,a_val,z_val,AgentDistMass) (aprime_val-a_val*(a_val~=10^6))*(aprime_val>a_val*(a_val~=10^6)); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
    FnsToEvaluateFn_Firing = @(aprime_val,a_val,z_val,AgentDistMass) -(aprime_val-a_val*(a_val~=10^6))*(aprime_val<a_val*(a_val~=10^6)); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
    FnsToEvaluate={FnsToEvaluateFn_Emp, FnsToEvaluateFn_Hiring, FnsToEvaluateFn_Firing};
end

% We will want the aggregate values of these. 
AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
% For much of Panel B we just need the pdf of the relevant measure (employment, hiring, firing)
ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
% We need a simulated panel based on whole distributions (for calculating
% variance of growth rates and serial correlation in log(n); for survivors).
% Note that because of these two moments we want to calculate it makes more
% sense to have a very large number of two period simulations, and since we
% just want survivors, we won't want entrants.
simoptions.noentryinpanel=1; % Don't want entry in this panel data simulation (we are just interested in 'survivors')
simoptions.simperiods=2;
simoptions.numbersims=10^4;
SimPanel=SimPanelValues_Case1(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions, EntryExitParamNames);
Survive_indicator=~isnan(shiftdim((SimPanel(2,2,:)),1));
SimPanel_Survivors=SimPanel(:,:,Survive_indicator);
GrowthRateEmploy=(SimPanel_Survivors(1,2,:)-SimPanel_Survivors(1,1,:))./SimPanel_Survivors(1,1,:);
VarianceOfGrowthRate_survivors=var(GrowthRateEmploy);
SerialCorrelationLogn_survivors=corr(log(shiftdim(SimPanel_Survivors(1,2,:),2)),log(shiftdim(SimPanel_Survivors(1,1,:),2)));
% We need a simulated panel based on new entrants for some (e.g., for stats by cohort)
simoptions.noentryinpanel=1; % Don't want further entry in this panel data simulation
simoptions.simperiods=20; % We anyway only need 10 for the stats being reported
simoptions.numbersims=10^4; % Default is 1000, this was not enough to get stable/smooth estimate of 'hazard rates by cohort'
EntrantDist.pdf=Params.upsilon;
EntrantDist.mass=Params.Ne;
SimPanel_Entrants=SimPanelValues_Case1(EntrantDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions, EntryExitParamNames);

% plot(sort(GrowthRateEmploy(:)))

EmploymentDecision_index=shiftdim(Policy,1);
EmploymentDecision=zeros(n_a,n_z);
a_grid_emp=a_grid; 
if ImposeFootnote5==1
    a_grid_emp(end)=0;
end
for a_c=1:n_a
    for z_c=1:n_z
%         if ~isnan(EmploymentDecision_index)
            if EmploymentDecision_index(a_c,z_c)>0
                EmploymentDecision(a_c,z_c)=a_grid_emp(EmploymentDecision_index(a_c,z_c));
            end
%         end
    end
end


% Start with Panel A.
% Average Firm Size (i.e., Average number of employees)
AvgFirmSize=AggValues(1)/StationaryDist.mass;
% Exit rate of Firms
MassOfExitingFirms=sum(sum(StationaryDist.pdf(logical(ExitPolicy))))*StationaryDist.mass;
ExitRateOfFirms=MassOfExitingFirms/StationaryDist.mass;
% In stationary eqm, firing must equal hiring, so can use either for
% turnover. (might need to adjust for entry???)
TurnoverRateOfJobs=AggValues(2)/AggValues(1); % the "/StationaryDist.mass" cancels top and bottom
% Fraction of hiring by new firms
TotalHiringByNewFirms=Params.Ne*sum(sum(Params.(EntryExitParamNames.DistOfNewAgents{1}).*EmploymentDecision)); %Note that all employment by new firms represents hiring (as they enter with zero employees).
FractionHiringByNewFirms=TotalHiringByNewFirms/AggValues(1);
% Average size of new firm
% This is awkward as Hopenhayn & Rogerson (1993) explictly state that new
% firms have zero employees in their first period of existence. Presumably
% this statistic means 'new firm after one period, conditional on survival'.
AvgSizeEnteringFirm=sum(sum(Params.(EntryExitParamNames.DistOfNewAgents{1}).*EmploymentDecision));
% Average size of exiting firm
DistOfExitingFirms=StationaryDist.pdf.*ExitPolicy/sum(sum(StationaryDist.pdf.*ExitPolicy));
AvgSizeExitingFirm=sum(sum(DistOfExitingFirms.*(a_grid_emp*ones(1,n_z)))); % Note, really this is 'lagged size', not current size.

% For Panel B of Table 1, first need to figure out the relevant partition of stationary distribution.
% The calculation is somewhat more complicated than at first glance, since 'current
% number of employees' is not the current state, it is the choice of next period state.
% % Following commented out lines give an equivalent way of calculating the partitions.
% [~,Firstcutoff]=min(abs(a_grid-20));
% [~,Secondcutoff]=min(abs(a_grid-100));
% [~,Thirdcutoff]=min(abs(a_grid-500));
% EmploymentDecision_index=shiftdim(Policy,1);
% FirstPartition=logical((EmploymentDecision_index>0).*(EmploymentDecision_index<Firstcutoff));
% SecondPartition=logical((EmploymentDecision_index>=Firstcutoff).*(EmploymentDecision_index<Secondcutoff));
% ThirdPartition=logical((EmploymentDecision_index>=Secondcutoff).*(EmploymentDecision_index<Thirdcutoff));
% FourthPartition=logical((EmploymentDecision_index>=Thirdcutoff));
FirstPartition=logical((EmploymentDecision>0).*(EmploymentDecision<20));
SecondPartition=logical((EmploymentDecision>=20).*(EmploymentDecision<100));
ThirdPartition=logical((EmploymentDecision>=100).*(EmploymentDecision<500));
FourthPartition=logical((EmploymentDecision>=500));

% Fraction of firm in each partition (conditional on not exiting)
FractionOfFirmsPerPartition=zeros(4,1);
pdfoffirms=StationaryDist.pdf;
FractionOfFirmsPerPartition(1)=sum(pdfoffirms(FirstPartition));
FractionOfFirmsPerPartition(2)=sum(pdfoffirms(SecondPartition));
FractionOfFirmsPerPartition(3)=sum(pdfoffirms(ThirdPartition));
FractionOfFirmsPerPartition(4)=sum(pdfoffirms(FourthPartition));

% Fraction of employment in each partition
FractionOfEmploymentPerPartition=zeros(4,1);
pdfofemploy=shiftdim(ProbDensityFns(1,:,:),1);
pdfofemploy2=(StationaryDist.pdf.*EmploymentDecision)/sum(sum(StationaryDist.pdf.*EmploymentDecision));
FractionOfEmploymentPerPartition(1)=sum(pdfofemploy(FirstPartition));
FractionOfEmploymentPerPartition(2)=sum(pdfofemploy(SecondPartition));
FractionOfEmploymentPerPartition(3)=sum(pdfofemploy(ThirdPartition));
FractionOfEmploymentPerPartition(4)=sum(pdfofemploy(FourthPartition));
% Fraction of hiring in each partition
FractionOfHiringPerPartition=zeros(4,1);
pdfofhiring=shiftdim(ProbDensityFns(2,:,:),1);
FractionOfHiringPerPartition(1)=sum(pdfofhiring(FirstPartition));
FractionOfHiringPerPartition(2)=sum(pdfofhiring(SecondPartition));
FractionOfHiringPerPartition(3)=sum(pdfofhiring(ThirdPartition));
FractionOfHiringPerPartition(4)=sum(pdfofhiring(FourthPartition));
% Fraction of firing in each partition
FractionOfFiringPerPartition=zeros(4,1);
pdfoffiring=shiftdim(ProbDensityFns(3,:,:),1);
FractionOfFiringPerPartition(1)=sum(pdfoffiring(FirstPartition));
FractionOfFiringPerPartition(2)=sum(pdfoffiring(SecondPartition));
FractionOfFiringPerPartition(3)=sum(pdfoffiring(ThirdPartition));
FractionOfFiringPerPartition(4)=sum(pdfoffiring(FourthPartition));

% Firms size by cohort: the current employment size of a firm is 'aprime'.
FractionOfFirmsPerPartition_CohortPeriod=zeros(4,size(SimPanel_Entrants,2));
for t=1:size(SimPanel_Entrants,2)
    temp=SimPanel_Entrants(1,t,:); % First variable in panel is 'aprime', which is the number of employees this period.
    temp=temp(~isnan(temp)); % For obvious reasons, calculation has to be done is conditional on survival (else won't sum to one)
    temp=temp(temp~=0); % Hopenhayn & Rogerson do this conditional on not choosing to exit (clear from fact the columns add to one in their numbers)
    NumberOfFirmsStillAlive_CurrentCohort=numel(temp);
    FractionOfFirmsPerPartition_CohortPeriod(1,t)=sum((temp>0).*(temp<20))/NumberOfFirmsStillAlive_CurrentCohort;
    FractionOfFirmsPerPartition_CohortPeriod(2,t)=sum((temp>=20).*(temp<100))/NumberOfFirmsStillAlive_CurrentCohort;
    FractionOfFirmsPerPartition_CohortPeriod(3,t)=sum((temp>=100).*(temp<500))/NumberOfFirmsStillAlive_CurrentCohort;
    FractionOfFirmsPerPartition_CohortPeriod(4,t)=sum((temp>=500))/NumberOfFirmsStillAlive_CurrentCohort;
end
% Hazard rates by cohort: Looking at SimPanel_Entrants a firm that exits
% will be represented by a number followed by a 'nan'.
NumberFirmsStillAlive=[size(SimPanel_Entrants,3),sum(~isnan(shiftdim(SimPanel_Entrants(1,:,:),1)),2)'];
PeriodHazardRate=1-NumberFirmsStillAlive(2:end)./NumberFirmsStillAlive(1:end-1);

%Table 1
FID = fopen('./SavedOutput/LatexInputs/HopenhaynRogerson1993_Table2.tex', 'w');
fprintf(FID, 'A: Summary Statistics for Benchmark Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lr} \n \\hline \\hline \n');
fprintf(FID, 'Average firm size & %8.2f \\\\ \n', AvgFirmSize);
fprintf(FID, 'Co-worker mean    & - \\\\ \n');
fprintf(FID, 'Variance of growth rates (survivors) & %8.2f \\\\ \n', VarianceOfGrowthRate_survivors);
fprintf(FID, 'Serial correlation in log(n) (survivors) & %8.2f \\\\ \n', SerialCorrelationLogn_survivors);
fprintf(FID, 'Exit rate of firms & %8.2f \\\\ \n', ExitRateOfFirms);
fprintf(FID, 'Turnover rate of jobs & %8.2f \\\\ \n', TurnoverRateOfJobs);
fprintf(FID, 'Fraction of hiring by new firms & %8.2f \\\\ \n', FractionHiringByNewFirms);
fprintf(FID, 'Average size of new firms & %8.2f \\\\ \n', AvgSizeEnteringFirm);
fprintf(FID, 'Average size of exiting firms & %8.2f \\\\ \n', AvgSizeExitingFirm);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, 'B: Size Distribution \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, ' & 1-19 & 20-99 & 100-499 & 500+  \\\\ \\hline \n');
fprintf(FID, 'Firms      & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFirmsPerPartition);
fprintf(FID, 'Employment & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfEmploymentPerPartition);
fprintf(FID, 'Hiring     & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfHiringPerPartition);
fprintf(FID, 'Firing     & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFiringPerPartition);
fprintf(FID, 'By cohort: & & & & \\\\ \n');
fprintf(FID, '\\, 1 period & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFirmsPerPartition_CohortPeriod(:,1));
fprintf(FID, '\\, 2 periods & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFirmsPerPartition_CohortPeriod(:,2));
fprintf(FID, '\\, 5 periods & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFirmsPerPartition_CohortPeriod(:,5));
fprintf(FID, '\\, 10 periods & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n',FractionOfFirmsPerPartition_CohortPeriod(:,10));
fprintf(FID, 'Hazard rates by cohort: & & & &  \\\\ \n');
fprintf(FID, '\\, 1 period & %8.2f & & &  \\\\ \n',PeriodHazardRate(1));
fprintf(FID, '\\, 2 periods & %8.2f & & &  \\\\ \n',PeriodHazardRate(2));
fprintf(FID, '\\, 5 periods & %8.2f & & &  \\\\ \n',PeriodHazardRate(5));
fprintf(FID, '\\, 10 periods & %8.2f & & &  \\\\ \n',PeriodHazardRate(10));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);



%% For comparison with Chris Edmonds slides, use following alternative cutoffs. 
% These results should be compared to those on slide 16. I have done this
% for tau=0, and get
% [~,CE_Firstcutoff]=min(abs(a_grid-20));
% [~,CE_Secondcutoff]=min(abs(a_grid-50));
% [~,CE_Thirdcutoff]=min(abs(a_grid-100));
% [~,CE_Fourthcutoff]=min(abs(a_grid-500));
CE_FirstPartition=logical((EmploymentDecision>0).*(EmploymentDecision<20));
CE_SecondPartition=logical((EmploymentDecision>=20).*(EmploymentDecision<50));
CE_ThirdPartition=logical((EmploymentDecision>=50).*(EmploymentDecision<100));
CE_FourthPartition=logical((EmploymentDecision>=100).*(EmploymentDecision<500));
CE_FifthPartition=logical((EmploymentDecision>=500));

% Fraction of firm in each partition
CE_FractionOfFirmsPerPartition=zeros(5,1);
CE_FractionOfFirmsPerPartition(1)=sum(pdfoffirms(CE_FirstPartition));
CE_FractionOfFirmsPerPartition(2)=sum(pdfoffirms(CE_SecondPartition));
CE_FractionOfFirmsPerPartition(3)=sum(pdfoffirms(CE_ThirdPartition));
CE_FractionOfFirmsPerPartition(4)=sum(pdfoffirms(CE_FourthPartition));
CE_FractionOfFirmsPerPartition(5)=sum(pdfoffirms(CE_FifthPartition));
% Fraction of employment in each partition
CE_FractionOfEmploymentPerPartition=zeros(5,1);
CE_FractionOfEmploymentPerPartition(1)=sum(pdfofemploy(CE_FirstPartition));
CE_FractionOfEmploymentPerPartition(2)=sum(pdfofemploy(CE_SecondPartition));
CE_FractionOfEmploymentPerPartition(3)=sum(pdfofemploy(CE_ThirdPartition));
CE_FractionOfEmploymentPerPartition(4)=sum(pdfofemploy(CE_FourthPartition));
CE_FractionOfEmploymentPerPartition(5)=sum(pdfofemploy(CE_FifthPartition));


%% Solving for various values of tau and creating Tables 3, 4 and 5.

fprintf('Current progress: Solving for various values of tau \n')
tau_vec=[0,0.1,0.2];
for tau_c=1:length(tau_vec)
    Params.tau=tau_vec(tau_c);
    Output.(['tau',num2str(tau_c)])=HopenhaynRogerson1993_Fn(ImposeFootnote5, Params, n_a, n_z,a_grid,z_grid,pi_z, ReturnFn, ReturnFnParamNames, DiscountFactorParamNames, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions);
end

save ./SavedOutput/HopenhaynRogerson1993_Output.mat Output
%% Table 3

% Most of the numbers in Table 3 are reported as relative to benchmark. So
% need to calculate these from the levels.
Table3_Cons=100*[Output.tau1.consumption, Output.tau2.consumption, Output.tau3.consumption]/Output.tau1.consumption;
Table3_AvgProductivity=100*[Output.tau1.averageproductivity, Output.tau2.averageproductivity, Output.tau3.averageproductivity]/Output.tau1.averageproductivity;
Table3_TotalEmp=100*[Output.tau1.totalemployment, Output.tau2.totalemployment, Output.tau3.totalemployment]/Output.tau1.totalemployment;
Table3_ConsEquiv=[100*exp(Output.tau1.utility_period-Output.tau1.utility_period),...
    100*exp(Output.tau2.utility_period-Output.tau1.utility_period),...
    100*exp(Output.tau3.utility_period-Output.tau1.utility_period)];

FID = fopen('./SavedOutput/LatexInputs/HopenhaynRogerson1993_Table3.tex', 'w');
fprintf(FID, 'Effect of Changes in $\\tau$ (Benchmark Model) \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & $\\tau=0$ & $\\tau=0.1$ & $\\tau=0.2$ \\\\ \\hline \n');
fprintf(FID, 'Price & %8.3f & %8.3f & %8.3f \\\\ \n', Output.tau1.price, Output.tau2.price, Output.tau3.price);
fprintf(FID, 'Consumption (output) & %8.1f & %8.1f & %8.1f \\\\ \n', Table3_Cons);
fprintf(FID, 'Average Productivity & %8.1f & %8.1f & %8.1f \\\\ \n', Table3_AvgProductivity);
fprintf(FID, 'Total Employment     & %8.1f & %8.1f & %8.1f \\\\ \n', Table3_TotalEmp);
fprintf(FID, 'Utility-adjusted consumption & %8.1f & %8.1f & %8.1f \\\\ \n', Table3_ConsEquiv);
fprintf(FID, 'Average firm size    & %8.1f & %8.1f & %8.1f \\\\ \n', Output.tau1.AvgFirmSize, Output.tau2.AvgFirmSize, Output.tau3.AvgFirmSize);
fprintf(FID, 'Layoff costs/wage bill & %8.1f & %8.3f & %8.3f \\\\ \n', Output.tau1.layoffcostsdivwagebill, Output.tau2.layoffcostsdivwagebill, Output.tau3.layoffcostsdivwagebill);
fprintf(FID, 'Job turnover rate    & %8.2f & %8.2f & %8.2f \\\\ \n', Output.tau1.TurnoverRateOfJobs, Output.tau2.TurnoverRateOfJobs, Output.tau3.TurnoverRateOfJobs);
fprintf(FID, 'Serial correlation in log(n) & %8.2f & %8.2f & %8.2f \\\\ \n', Output.tau1.SerialCorrelationLogn_survivors, Output.tau2.SerialCorrelationLogn_survivors, Output.tau3.SerialCorrelationLogn_survivors);
fprintf(FID, 'Variance in growth rate       & %8.2f & %8.2f & %8.2f \\\\ \n', Output.tau1.VarianceOfGrowthRate_survivors, Output.tau2.VarianceOfGrowthRate_survivors, Output.tau3.VarianceOfGrowthRate_survivors);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Consumption (output), Average Productivity, Total Employment and Utility-adjusted consumption are all reported relative to $\\tau=0$ (which is set equal to 100). Effectively also true of price which is normalized to one in the $\\tau=0$ case as part of calibration. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);



%% Table 4
zindex_Table4=[8,12,16,19,20];
FID = fopen('./SavedOutput/LatexInputs/HopenhaynRogerson1993_Table4.tex', 'w');
fprintf(FID, 'Effect of $\\tau$ on Decision Rules \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{2}{c}{$\\tau=0.1$} & \\multicolumn{2}{c}{$\\tau=0.2$} \\\\ \\cline{2-3} \\cline{4-5} \n');
fprintf(FID, 'z & $n_l$ & $n_u$ & $n_l$ & $n_u$ \\\\ \\hline \n');
fprintf(FID, '%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', z_grid(zindex_Table4(1)), Output.tau2.n_l(zindex_Table4(1)), Output.tau2.n_u(zindex_Table4(1)), Output.tau3.n_l(zindex_Table4(1)), Output.tau3.n_u(zindex_Table4(1)));
fprintf(FID, '%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', z_grid(zindex_Table4(2)), Output.tau2.n_l(zindex_Table4(2)), Output.tau2.n_u(zindex_Table4(2)), Output.tau3.n_l(zindex_Table4(2)), Output.tau3.n_u(zindex_Table4(2)));
fprintf(FID, '%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', z_grid(zindex_Table4(3)), Output.tau2.n_l(zindex_Table4(3)), Output.tau2.n_u(zindex_Table4(3)), Output.tau3.n_l(zindex_Table4(3)), Output.tau3.n_u(zindex_Table4(3)));
fprintf(FID, '%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', z_grid(zindex_Table4(4)), Output.tau2.n_l(zindex_Table4(4)), Output.tau2.n_u(zindex_Table4(4)), Output.tau3.n_l(zindex_Table4(4)), Output.tau3.n_u(zindex_Table4(4)));
fprintf(FID, '%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', z_grid(zindex_Table4(5)), Output.tau2.n_l(zindex_Table4(5)), Output.tau2.n_u(zindex_Table4(5)), Output.tau3.n_l(zindex_Table4(5)), Output.tau3.n_u(zindex_Table4(5)));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Hopenhayn \\& Rogerson (1993) call the first column log(s), but it is clear from the values it should be s, which I call z (the idiosyncratic productivity). Difference in the first column, z, between original and replication represent differences in the grids. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Table 5

FID = fopen('./SavedOutput/LatexInputs/HopenhaynRogerson1993_Table5.tex', 'w');
fprintf(FID, 'Absolute Deviations from MPL=1/p \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, '& \\multicolumn{2}{c}{Fraction of Firms} \\\\ \n');
fprintf(FID, '& \\multicolumn{2}{c}{within Interval} \\\\  \\cline{2-3} \n');
fprintf(FID, 'Size of Deviation (\\%%) & $\\tau=0.1$ & $\\tau=0.2$ \\\\ \\hline \n');
fprintf(FID, '0-3 & %8.2f & %8.2f \\\\ \n', Output.tau2.AbsDevsForTable5(1), Output.tau3.AbsDevsForTable5(1));
fprintf(FID, '3-5 & %8.2f & %8.2f \\\\ \n', Output.tau2.AbsDevsForTable5(2), Output.tau3.AbsDevsForTable5(2));
fprintf(FID, '5-10 & %8.2f & %8.2f \\\\ \n', Output.tau2.AbsDevsForTable5(3), Output.tau3.AbsDevsForTable5(3));
fprintf(FID, '10-15 & %8.2f & %8.2f \\\\ \n', Output.tau2.AbsDevsForTable5(4), Output.tau3.AbsDevsForTable5(4));
fprintf(FID, '$>$15 & %8.2f & %8.2f \\\\ \n', Output.tau2.AbsDevsForTable5(5), Output.tau3.AbsDevsForTable5(5));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: blah blah blah \\\\ \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);


