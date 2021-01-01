% Replication of Restuccia & Rogerson (2008) - Policy Distortions and Aggregate Productivity with Heterogeneous Establishments

% Essentially, everything replicates except Table 7. Good odds this is
% simply that I have been unable to understand exactly what the denominator
% is for the 'relative' calculations (the bottom row which is not
% 'relative' replicates cleanly).

% Note: rho is defined differently here to in the paper. Here it does not
% include the (exogenous) conditional probability of firm exit. In RR2008
% it does. This is done to emphasize how to include exit/death in models
% more generally when using VFI Toolkit.

% For explanation of model, the calibration choices, and what the policy
% experiments are aiming to understand, see paper. (This replication solves
% the model hundreds of times more than necessary, as solves many
% combinations of parameters/policies for which no results are reported in paper.)

% We will consider there to be 'two types of agents', namely firms and potential entrants.
% The importance of (existing) firms is obvious.
% The importance of the potential entrants to general equilibrium condition of the model is made clear in the paper.

% The code provided by RR2008 (namely RR_model_RED.m) contains an typo/error. It 
% sets variable i as interest rate on line 23. It then overwrites this as part 
% of using i in a for-loop on lines 43-46 so that value of i ends up equal to 
% 'length(datazupper2)'. On line 68 i is used to calculate rho, with the 
% intention that i is still the interest rate, but this is no longer true.
% Thanks to Denise Manfredini who first brought this to my attention. You
% can easily correct by simply editing for-loop on lines 43-46.

% Note: In principle the firm productivity and whether it faces
% tax/zero/subsidy are both 'agent permanent types' and could be modelled
% in this fashion by VFI Toolkit. Here we will instead model them as
% 'exogenous states', or 'z variables' in the nomenclature of the toolkit (they will be "made
% permanent" by setting the markov transition matrix for these states to be a diagonal matrix of zeros).
% Three reasons underlie this choice, the first is that setting them up as
% 'agent permanent types' would involve a lot more lines of code (because it is
% intended to allow for different agent permanent types to face genuinely
% different problems, rather than just different values of a single
% parameter as here), the second is the sheer number of different permanent types (300), 
% and the third is that because the problem to be solved
% is so simple we can easily still solve the problem. (The way they will be
% set up each permanent type will be calculating future expectations taking
% into account the (zero) probability that they become a different
% permanent type. This is wasteful, but not much because the small size of 
% the problem to be solved means we can just solve them all at the same time anyway.)

parpool
% PoolDetails=gcp;
% NCores=PoolDetails.NumWorkers;
Parallel=1; % 1 for (parallel) CPUs, 2 for GPU

%% Set some basic variables

n_s=100; % Firm-specific Productivity level
n_tau=3; % Negative (subsidy), zero, and positive (tax).

%Parameters
Params.beta=0.96; % Discount rate
% Firms
Params.alpha=0.283;
Params.gamma=0.567;
Params.delta=0.08; % Depreciation rate of physical capital
Params.cf=0; % Fixed cost of production
% Firm entry and exit
Params.ce=1; % Fixed cost of entry (this is a normalization)
Params.lambda=0.1; % Probability of firm exit

Params.Ne=1; % RR2008 call this E, the mass of new potential new entrant distribution.

% The actual 'distortionary policy rate'
Params.taurate=0; % This is the rate for the tax.
Params.subsidyrate=0; % This is the rate for the subsidy.

Params.w=1; % This is just an initial guess since w (wage) is determined in general eqm.

%% Grid on the productivity level, s, as well as the distribution over these h(s).
% The grid is here created by the next few lines of text.
% The original numbers which can be found in the codes provided by Restuccia & Rogerson (2008) for their model. 
% The original source is Rossi-Hansberg and Wright (AER, 2007) All industries for 2000,
% which in turn is created from data from the US Census of Businesses. The
% following commented out line does the import that creates what is stored here as text.
% load ./OriginalCodes/files/establishment_dist.txt
establishment_dist=[4, 0.482359949; 9, 0.21665686; 14, 0.092591623;...
    19, 0.049473226; 24, 0.031397291; 29, 0.021542855; 34, 0.015801151; 39, 0.011847888; 44, 0.009307797;...
    49, 0.007457313; 59, 0.011325843; 69, 0.007940958; 79, 0.005988763; 89, 0.004677303; 99, 0.003753013;...
    124, 0.00683308; 149, 0.004509741; 174, 0.003126242; 199, 0.002215598; 224, 0.001685936; 249, 0.001293212;...
    299, 0.001867145; 349, 0.001285596; 399, 0.000873831; 449, 0.000649621; 499, 0.000511097; 599, 0.000705634;...
    699, 0.000463018; 799, 0.000331475; 899, 0.0002469; 999, 0.000184065; 1249, 0.000311482; 1499, 0.000191522;...
    1749, 0.000130432; 1999, 9.10802E-05; 2249, 7.14044E-05; 2499, 5.42673E-05; 2999, 6.96589E-05; ...
    3499, 4.42707E-05; 3999, 3.09419E-05; 4499, 1.96759E-05; 4999, 1.60263E-05; 9999, 5.06178E-05; ...
    10000, 1.45982E-05];
% establishment_dist(:,1); %(upper range of number of employees)
% establishment_dist(:,2); %(fraction of establishment in each group)

% Note: Table 1 of RR2008 reports an 's range', this is actually implicit
% from the numbers in establishment_dist and the parameters gamma and alpha.
s_grid=exp(linspace(log(1),log(establishment_dist(end,1)^(1-Params.gamma-Params.alpha)),n_s));
% The stationary distribution, called h(s) by RR2008, will here be called pistar_s
cumsum_pistar_s=interp1(establishment_dist(:,1),cumsum(establishment_dist(:,2)),s_grid.^(1/(1-Params.gamma-Params.alpha)));
% The first few points of s_grid have to be extrapolated (as outside the
% range of the actual data). Following line implements this, and similar
% for the max value. I just divide equally across these first few (which is same as RR2008 do).
temp=~isnan(cumsum_pistar_s);
tempind=find(temp,1,'first');
cumsum_pistar_s(1:tempind)=cumsum((1/tempind)*cumsum_pistar_s(tempind)*ones(1,tempind));
cumsum_pistar_s(end)=1;
pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)])';
% Note: This stationary distribution is created by interpolating establishment_dist(:,2)
% (interpolation of the cumulative distribution tends to perform better in
% practice than interpolation of the pdf, as cdf tends to be smoother)
% Here this is done using 1D interpolation (the matlab default spline interpolation), RR2008 do this manually by evenly
% dividing the data across the s_grid points that are 'in-between' establishment_dist points.
% (You can see the difference if you use original RR2008 codes to create 'ns' and 'hs', and then run following line
% plot(establishment_dist(:,1),cumsum(establishment_dist(:,2)),s_grid.^(1/(1-Params.gamma-Params.alpha)), cumsum(pistar_s), z, cumsum(hs))
% The difference is not large, but the use of interpolation clearly provide better approximation.

% Note: if you want the actual establishment sizes implicit in s_grid to
% compare to data just run: s_grid.^(1/(1-Params.gamma-Params.alpha))

% Comment: The calibration of the model is awkward, in the sense that it is
% the distribution of potential entrants that is calibrated to match the US
% Census Bureau data on existing firms, while conceptually it is the
% stationary distribution of (exisiting) firms that should match the US
% Census Bureau data on existing firms. Although this later would require a
% much more complicated calibration procedure (moment matching based on simulation).

%% Grids
% Grids for d variables: there are none (because all the decisions being
% made by the firms have closed form solutions, so there are none that the
% toolkit actually needs to solve for)
n_d=0;
d_grid=[];
% Grids for a variables: there are none
n_a=1; % Still need to declare an a variable, the toolkit is not really designed for situation of none. Since it only has a single value and won't really be used anywhere it won't actually do anything anyway.
a_grid=1;
% Grids for z variables: there are two
% s_grid
tau_grid=[-1; 0; 1]; % Negative will be subsidy, zero, positive will be tax
n_z=[n_s,length(tau_grid)];
z_grid=[s_grid'; tau_grid];
pi_z=eye(prod(n_z),prod(n_z)); % Ones on the diagonal, zeros elsewhere. These are essentially permanent states, and there is zero probability of transition.

%% Figure 1

figure(1)
semilogx(s_grid.^(1/(1-Params.gamma-Params.alpha)),cumsum(pistar_s),'-',establishment_dist(:,1),cumsum(establishment_dist(:,2)),'o','LineWidth',2)
title('Distribution of Plants by Employment')
set(gca,'XTick',[1 10 100 1000 10000])
set(gca,'XTickLabel',{'1','10','100','1,000','10,000'})
xlabel('Number of Employees (log scale)')
ylabel('Cummulative Distribution of Establishments')
legend('Model','Data','Location','NorthWest')
axis([0 14000 0 1.1])
saveas(gcf,['./SavedOutput/Graphs/RestucciaRogerson2008_Figure1.pdf'])
% Remark: the 'sagging curves' of model between data points is a result of
% the graph being done with log scale, while the interpolation is done
% based on (linear scale) number of employees directly.

%% Odds and ends

% Because asset markets are complete the households problem leads to the
% standard result that the households can just be treated as a
% Representative Household, and equilibrium in the asset markets requires
% that the interest rate, i (referred to in RR2008 paper as R, but here
% follow their codes and refer to as i)
Params.i=1/Params.beta-1; % This is standard general eqm result in complete market models, comes from consumption euler eqn together with requirements of stationary eqm.
% The net return to capital in equilibrium will thus be
Params.r=Params.i+Params.delta; % That the gross return is just 1/beta-1 and equals i (that the gross return to capital equals the interest rate is a requirement of capital market clearance in model)

% This (i=1/beta-1) means HHs supply amount of K that clears market, and
% labour supply is perfectly inelastic (leisure is not valued by HHs) 
% means that any value of w clears labour market.

%% Solve the Value Function
% Entry is irrelevant to the value function problem. Exit is relevant, and treated differently based on whether it is endogenous or exogenous.

% For exogenous exit, you simply need to include the 'conditional surivival probability' as another 'DiscountFactorParamNames'
Params.oneminuslambda=1-Params.lambda; % This is now the conditional probability of survival.
Params.rho=1/(1+Params.i); % The factor at which firms discount the future depends on both the risk-free interest rate and the risk/probability of exit
DiscountFactorParamNames={'rho','oneminuslambda'};

ReturnFn=@(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf) RestucciaRogerson2008_ReturnFn(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf);
ReturnFnParamNames={'w','r','alpha','gamma','taurate','subsidyrate','cf'}; %It is important that these are in same order as they appear in 'RestucciaRogerson2008_ReturnFn'

% Check that everything is working so far by solving the value function
vfoptions.parallel=Parallel;
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);

% If you wanted to look at the value fn
% figure(2)
% surf(shiftdim(V,1))
% title('Value fn')

%% Simulate Agent Distribution, with Entry and Exit
% I use different notation to RR2008. The pdf of potential
% entrants I call upsilon (they used g). I call the 'conditional entry
% decision' ebar (they used xbar). I call the 'mass' by which this is
% multiplied Ne (they used E). Note that Ne is not the actual mass of
% entrants because of the role of ebar.
% So the actual distribution of entrants is then Ne*upsilon.*ebar (note that
% this is not a pdf). The actual mass of entrants can be calculated by
% repeatedly taking the sum of Ne*upsilon.*ebar.

% Need to 'activate' these options
simoptions.agententryandexit=1
% Note: it is possible to do just entry or just exit.
% To do just entry, simply set CondlProbOfExit=0;
% To do just exit, simply set DistOfNewAgents=0;
% (It is not envisaged that you are likely to want either of entry or exit without the other)
% Note that from the perspective of the simulation of agent distribution
% whether these are endogenous or exogenous entry/exit decisions is irrelevant.

pistar_tau=[0;1;0];

% Because they are not a default part of agent simulation, you need to pass the entry/exit aspects as part of simoptions.
EntryExitParamNames.DistOfNewAgents={'upsilon'};
%Params.upsilon=kron(pistar_tau,pistar_s); % Note: these should be in 'reverse order'
Params.upsilon=pistar_s.*(pistar_tau');
EntryExitParamNames.CondlEntryDecisions={'ebar'};
Params.ebar=ones([n_a,n_z]); % Takes value of one for enter, zero for not-enter. This is just an initial guess as the actual decisions are determined as part of general equilibrium.
% Note: VFI Toolkit requires the DistOfNewAgents to be a pdf (so unit mass), and then uses
% the 'MassOfNewAgents' to understand how many there will be entering relative to existing agents.
EntryExitParamNames.MassOfNewAgents={'Ne'}; % This is implied by the need for there to be a stationary distribution, so the mass of firm entry must equal the mass of firm exit, which is lambda*1

EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};
% This conditional probability can be a state-dependent parameter, in which case this input would be a vector or matrix, etc.

simoptions.parallel=Parallel;
% Check that everything is working so far by solving the simulation of
% agent distribution to get the stationary distribution.
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

% If you wanted to look at the pdf:
% surf(shiftdim(StationaryDist.pdf,1))

%% General equilibrium conditions (and aggregates needed to evaluation them)

% Endogenous entry adds a general equilibrium condition, specifically a
% 'free entry' condition.
% The only additional information needed for the general equilibrium in terms of entry
% and exit is to know which if any of the general equilibrium conditions are 'free entry'
% conditions and so need to be evaluated on potential entrants
% distribution.

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
% FnsToEvaluateFn_1 = @(aprime_val,a_val,z1_val,z2_val) a_val; %We just want the aggregate assets (which is this periods state)
% FnsToEvaluate={FnsToEvaluateFn_1};
FnsToEvaluate={};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
%Note: length(AggVars) is as for FnsToEvaluate and length(p) is length(n_p)
heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};

% A 'condlentry' general equilibrium condition will take values of greater
% than zero for firms that decide to enter, less than zero for first that
% decide not to enter (or more accurately, after entry decision they draw
% their state, and then decide to cancel/abort their entry).

GEPriceParamNames={'w'}; 
    % Note that this parameter does not directly appear in any of the general eqm conditions, only indirectly by it's effect on the value fn.
    % Note that 'ebar', the conditional entry decision, is also determined as part of general eqm, and so in some sense is a general eqm parameter. But since this is unavoidably the case for conditional entry there is no need to declare it.
GeneralEqmEqnParamNames(1).Names={'beta'};
GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0; % % Conditional entry condition
% (ValueFn>=0): should this be changed to (beta*ValueFn>=0)? Yes, because of timing.
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry};

% In principle, 'Ne', the mass of (potential) new entrants is also a
% parameter to be determined in general equilibrium, and hence would also
% appear in GEPriceParamNames. The relevant condition would be labour
% market clearance 1-Nbar=0. The labour supply in the model is equal to 1
% (households are endowed with a unit endowment of labour, and this is
% supplied perfectly inelastically because households do not value
% leisure). Nbar is the labour demand and comes from integral of nbar over
% the distribution of firms. To implement this we would need to add the
% following:
% GEPriceParamNames={'w','Ne'};
% FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate'};
% FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate) (((1-taurate*z2_val)*z1_val*gamma)/w)^(1/(1-gamma)) *((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-taurate*z2_val))^(1/(1-alpha-gamma)))^(alpha/(1-gamma)); % which evaluates to Nbar in the aggregate
% GeneralEqmEqnParamNames(3).Names={};
% GeneralEqmEqn_LabourMarket = @(AggVars,GEprices) 1-AggVars;
% The obvious changes to include these in 
% FnsToEvaluate={FnsToEvaluateFn_nbar};
% heteroagentoptions.specialgeneqmcondn={0,'condlentry','entry'};
% GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry, GeneralEqmEqn_LabourMarket}; 
% would also need to be made.
% Because RR2008 model is linear in Ne (not the case for most models of entry), we 
% can instead simply impose this afterwards, by setting Ne=1/Nbar. This is done
% here just after computing the general equilibrium.
% Obviously endogenizing labour supply would change this, and this 'short-cut' version would no longer work.

heteroagentoptions.verbose=0;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
% findeqmtime=toc
Params.w=p_eqm.w;
Params.ebar=p_eqm.ebar;

% RR2008 gets baseline solution of w=1.8955, xbar is all ones.
% I get w=1.9074, ebar is all ones.

% Calculate some things in the general eqm
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

% Impose the labour market clearance, which involves calculating Ne. See
% comments above (about 10-20 lines above).
FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*z1_val)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluate={FnsToEvaluateFn_nbar};
AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
InitialNe=Params.Ne;
Params.Ne=1/AggValues; % AggValues is presently equal to Nbar. This line is imposing/satisfying the labour market clearance condition.
StationaryDist.mass=StationaryDist.mass*(Params.Ne/InitialNe); % Take advantage of linearity of the stationary distribution in new entrants distribution.

%% Table 1
% Not exactly replicating anything as this is just the parameter values...

%Table 1
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table1.tex', 'w');
fprintf(FID, 'Benchmark calibration to US data \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lll} \n \\hline \n');
fprintf(FID, 'Parameters & Value & Target \\\\ \\hline \n');
fprintf(FID, '$\\alpha$ & %8.3f & Capital income share \\\\ \n', Params.alpha);
fprintf(FID, '$\\gamma$ & %8.3f & Labour income share \\\\ \n', Params.gamma);
fprintf(FID, '$\\beta$ & %8.2f & Real rate of return \\\\ \n', Params.beta);
fprintf(FID, '$\\delta$ & %8.2f & Investment to output rate \\\\ \n', Params.delta);
fprintf(FID, '$c_e$ & %8.1f & Normalization \\\\ \n', Params.ce);
fprintf(FID, '$c_f$ & %8.1f & Benchmark rate \\\\ \n', Params.cf);
fprintf(FID, '$\\lambda$ & %8.1f & Annual exit rate \\\\ \n', Params.lambda);
fprintf(FID, '$s$ range & [%d, %8.2f] & Relative establishment sizes \\\\ \n', s_grid(1), s_grid(end));
fprintf(FID, '$h(s)$ & see Fig. 1 & Size distribution of establishments \\\\ \n');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Table 2

% One could code this directly easily enough since there are closed form
% solutions for kbar and nbar in terms of (s,tau). But will use the VFI
% Toolkit commands simply to show how to apply them.

FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
FnsToEvaluateParamNames(2).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*z1_val)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluateParamNames(3).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_nbar, FnsToEvaluateFn_output};

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_Mass(StationaryDist.mass, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames,EntryExitParamNames, n_d, n_a, n_z, [], a_grid, z_grid, Parallel,simoptions);

ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, [], a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

% s_grid.^(1/(1-Params.gamma-Params.alpha))
nbarValues=shiftdim(ValuesOnGrid(2,:,:,:),1);
% THIS IS REALLY WEIRD. Model is setup using h(s) and s grid values to make
% it so that firms can choose realistic numbers of employees. They choose
% much less than in the data. RR2008 then just go and renormalize the
% choices so that it looks like the same magnitude as in the data again!!!
nbarValues=shiftdim(ValuesOnGrid(2,:,:,:),1);
normalize_employment=nbarValues(1,1,2); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./normalize_employment;

Partion1Indicator=logical(nbarValues<5);
Partion2Indicator=logical((nbarValues>=5).*(nbarValues<50));
Partion3Indicator=logical(nbarValues>=50);

% Check that the following is equal to prod(n_z), so 300
sum(sum(Partion1Indicator+Partion2Indicator+Partion3Indicator))

ShareOfEstablishments(1)=sum(sum(StationaryDist.pdf(Partion1Indicator)));
ShareOfEstablishments(2)=sum(sum(StationaryDist.pdf(Partion2Indicator)));
ShareOfEstablishments(3)=sum(sum(StationaryDist.pdf(Partion3Indicator)));

Output_pdf=shiftdim(ProbDensityFns(3,:,:,:),1);
ShareOfOutput(1)=sum(sum(sum(Output_pdf(Partion1Indicator))));
ShareOfOutput(2)=sum(sum(sum(Output_pdf(Partion2Indicator))));
ShareOfOutput(3)=sum(sum(sum(Output_pdf(Partion3Indicator))));

Labour_pdf=shiftdim(ProbDensityFns(2,:,:,:),1);
ShareOfLabour(1)=sum(sum(sum(Labour_pdf(Partion1Indicator))));
ShareOfLabour(2)=sum(sum(sum(Labour_pdf(Partion2Indicator))));
ShareOfLabour(3)=sum(sum(sum(Labour_pdf(Partion3Indicator))));

Capital_pdf=shiftdim(ProbDensityFns(1,:,:,:),1);
ShareOfCapital(1)=sum(sum(sum(Capital_pdf(Partion1Indicator))));
ShareOfCapital(2)=sum(sum(sum(Capital_pdf(Partion2Indicator))));
ShareOfCapital(3)=sum(sum(sum(Capital_pdf(Partion3Indicator))));

AverageEmployment(1)=sum(sum(nbarValues(Partion1Indicator).*StationaryDist.pdf(Partion1Indicator)))/sum(sum(StationaryDist.pdf(Partion1Indicator)));
AverageEmployment(2)=sum(sum(nbarValues(Partion2Indicator).*StationaryDist.pdf(Partion2Indicator)))/sum(sum(StationaryDist.pdf(Partion2Indicator)));
AverageEmployment(3)=sum(sum(nbarValues(Partion3Indicator).*StationaryDist.pdf(Partion3Indicator)))/sum(sum(StationaryDist.pdf(Partion3Indicator)));

%Table 2
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table2.tex', 'w');
fprintf(FID, 'Distribution statistics of benchmark economy \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llcr} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{3}{l}{Establishment size (number of employees)} \\\\ \\cline{2-4} \n');
fprintf(FID, ' & $<$5 & 5 to 49 & $\\geq$50 \\\\ \\hline \n');
fprintf(FID, 'Share of establishments & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfEstablishments);
fprintf(FID, 'Share of output         & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfOutput);
fprintf(FID, 'Share of labour         & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfLabour);
fprintf(FID, 'Share of capital        & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfCapital);
fprintf(FID, 'Share of employment     & %8.2f & %8.2f & %8.2f \\\\ \n', AverageEmployment);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Need to calculate the 'baseline K' (aggregate capital in the baseline case with no taxes/distortions)
% as this is used to determine the size of the subsidies (taurate_s) in the
% rest of the cases where the subsidy is set to ensure that the level of
% aggregate capital remains unchanged from the baseline.
FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','w','taurate','subsidyrate'};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
FnsToEvaluate={FnsToEvaluateFn_kbar};

Params.Kbaseline=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

%% We have now solved the baseline case. Will solve the rest of the cases using the following loop. 
% All the cases are essentially solved exactly the same way as the
% baseline, and contents of this function are largely just a copy-paste of
% the above.
% I solve the full 'cross product' of all possible combinations of cases.
% This is way more than actually required/reported in RR2008, but means I
% can just use a bunch of nested for-loops (obviously it takes much
% longer).

% TaxProductivityCorrelation: indexes economy to be computed:
% 1 - if iid plant-specific tax/subsidy
% 2 - if negative corr of tax/subsidy and productivity (tax high s)
% 3 - if positive corr of tax/subsidy and productivity  (tax low s)

% Five different tax rates are considered for output tax; the sixth is used for labour tax rates.
tauratevec=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5]; % taup as this is the tax/subsidy potentially correlated with productivity.
% Fraction taxed
fractiontaxedvec=0.1:0.1:0.9; % More than are actually needed (for Table 4).

taxbasevec={'outputtax','capitaltax','labourtax'};
% Following forumulae look complex, mainly because it has to contain two
% 'copies' to allow for taurate and subsidy rate to be different.
% (Note: the inputs remain unchanged in each case) 
clear FnsToEvaluateFn_kbar FnsToEvaluateFn_nbar FnsToEvaluateFn_output

FnsToEvaluateFn_kbar.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
FnsToEvaluateFn_kbar.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *z1_val^(1/(1-alpha-gamma));
FnsToEvaluateFn_kbar.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^(gamma/(1-gamma-alpha)) *z1_val^(1/(1-alpha-gamma));

FnsToEvaluateFn_nbar.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*z1_val)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluateFn_nbar.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) z1_val^(1/(1-alpha-gamma)) *(alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluateFn_nbar.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) z1_val^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate

FnsToEvaluateFn_output.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_output.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) z1_val^(1/(1-gamma-alpha)) *(alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_output.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^(gamma/(1-gamma-alpha));

FnsToEvaluateFn_subsidy.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*subsidyrate* ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)); % (z2_val<0)*subsidyrate)* output
FnsToEvaluateFn_subsidy.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*subsidyrate*r* (alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *z1_val^(1/(1-alpha-gamma)); % (z2_val<0)*subsidyrate)*r* capital
FnsToEvaluateFn_subsidy.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*subsidyrate*w* z1_val^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^((1-alpha)/(1-gamma-alpha)); % (z2_val<0)*subsidyrate)*w* nbar

% Following are just 'indicator for subsidised' times output, needed to calculate Ys
FnsToEvaluateFn_outputofsubsidised.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_outputofsubsidised.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*z1_val^(1/(1-gamma-alpha)) *(alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_outputofsubsidised.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<0)*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^(gamma/(1-gamma-alpha));

% For Table 7, RR2008 redefine Ys to be 'output of non-taxed', rather than 'output of subsidised'
% For dealing with this, need the following.
FnsToEvaluateFn_outputofnontaxed.outputtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<=0)*((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_outputofnontaxed.capitaltax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<=0)*z1_val^(1/(1-gamma-alpha)) *(alpha/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*r))^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
FnsToEvaluateFn_outputofnontaxed.labourtax = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,w,taurate,subsidyrate) (z2_val<=0)*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/((1+((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*w))^(gamma/(1-gamma-alpha));

% Use a slightly different ReturnFn for the different tax bases.
% (Note: the inputs remain unchanged in each case)
% You could alternatively have a 'tax base' input, and just use the same ReturnFn for all.
clear ReturnFn
ReturnFn.outputtax =@(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf) RestucciaRogerson2008_ReturnFn(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf);
ReturnFn.capitaltax =@(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf) RestucciaRogerson2008_ReturnFn_CapitalTax(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf);
ReturnFn.labourtax =@(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf) RestucciaRogerson2008_ReturnFn_LabourTax(aprime_val, a_val, z1_val,z2_val,w,r,alpha,gamma,taurate,subsidyrate,cf);

% Start with all the results for output tax (that is, everything except Tables 8 and 9)
taxbase_c=1; % Output tax
fixedsubsidy_c=0; % Setting to 1 allows to overrule standard behaviour and fix the subsidy rate rather than optimizing it as is the default.
for fractiontaxed_c=1:length(fractiontaxedvec) %[1,3,5,7,9]%    
    fractiontaxed=fractiontaxedvec(fractiontaxed_c);
    for index_taus=0:1 % Set to 1 if tax only experiments
        for TaxProductivityCorrelation=1:2 % Note: 3 is never actually reported in any results Tables.
            
            if TaxProductivityCorrelation==1 % Uncorrelated
                pistar_tau=[1-fractiontaxed;0;fractiontaxed]; % [Subsidy, nothing, tax]
                Params.upsilon=pistar_s.*(pistar_tau'); % upsilon is the 'potential entrants distribution' (actual entrants distribution is upsilon.*ebar)
            elseif TaxProductivityCorrelation==2
                [~,productivity_cutoff]=max((cumsum(pistar_s)>(1-fractiontaxed))); % '1-fractiontaxed', those below this will have tax set to zero, leaving fractiontaxed as those taxed.
                Params.upsilon=pistar_s.*[ones(n_s,1),zeros(n_s,1),ones(n_s,1)];
                Params.upsilon(1:(productivity_cutoff-1),3)=0; % No tax on low productivity
                Params.upsilon(productivity_cutoff:end,1)=0; % No subsidies on high productivity
            elseif TaxProductivityCorrelation==3
                [~,productivity_cutoff]=max((cumsum(pistar_s)>fractiontaxed));% 'fractiontaxed', those above this will have tax set to zero, leaving fractiontaxed as those taxed.
                Params.upsilon=pistar_s.*[ones(n_s,1),zeros(n_s,1),ones(n_s,1)];
                Params.upsilon(productivity_cutoff:end,3)=0; % No tax on high productivity
                Params.upsilon(1:(productivity_cutoff-1),1)=0; % No subsidy on low productivity
            end
            if index_taus==1 % No subsidies (and so not trying to keep capital equal to Kbaseline)
                Params.upsilon(:,2)=Params.upsilon(:,2)+Params.upsilon(:,1); % redistribute 'subsidy' firms to 'nothing'
                Params.upsilon(:,1)=0; % and eliminate subsidies
            end
            
            for tau_c=1:length(tauratevec)
                Params.taurate=tauratevec(tau_c);
                % Will use initial guess of
                Params.subsidyrate=0.5*Params.taurate; % Note: is zero when tax rate is zero
                if sum(Params.upsilon(:,1))==0 % If no firms are subsidied
                    Params.subsidyrate=0; % Irrelevant as no subsidized firms, but zero seems 'nice and intuitive' in this case.
                elseif index_taus==1 % 'Tax only experiments' (for Table 7)
                    Params.subsidyrate=0;
                end
                
                if fixedsubsidy_c==1
                    Params.subsidyrate=Params.taurate;
                end
                
                FullFnsToEvaluate{1}=FnsToEvaluateFn_kbar.(taxbasevec{taxbase_c});
                FullFnsToEvaluate{2}=FnsToEvaluateFn_nbar.(taxbasevec{taxbase_c});
                FullFnsToEvaluate{3}=FnsToEvaluateFn_output.(taxbasevec{taxbase_c});
                FullFnsToEvaluate{4}=FnsToEvaluateFn_subsidy.(taxbasevec{taxbase_c});
                FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofsubsidised.(taxbasevec{taxbase_c});
                if index_taus==1 % Table 7, need to redefine Ys to 'non-taxed' instead of 'subsidised' as everywhere else.
                    FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofnontaxed.(taxbasevec{taxbase_c});
                end
                
                fprintf('Now running RestucciaRogerson2008 for combination: \n')
                [fixedsubsidy_c+1, TaxProductivityCorrelation, fractiontaxed_c,index_taus+1, taxbase_c, tau_c]
                Params.subsidyrate
                sum(Params.upsilon,1)
                            
                FullResults(fixedsubsidy_c+1, TaxProductivityCorrelation,fractiontaxed_c,index_taus+1,taxbase_c,tau_c).Output=RestucciaRogerson2008_Fn(fixedsubsidy_c, normalize_employment, n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,ReturnFn.(taxbasevec{taxbase_c}),DiscountFactorParamNames, ReturnFnParamNames, FullFnsToEvaluate, GEPriceParamNames, GeneralEqmEqnParamNames, GeneralEqmEqns, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions);
            end
        end
    end
end

% Do the parts for the other tax bases seperately as these are quite different.
taxbase_c=2; % Capital tax
tauratevec_capitaltax=[0,0.5,1.0];
fixedsubsidy_c=0;
fractiontaxed_c=5; fractiontaxed=fractiontaxedvec(fractiontaxed_c);
index_taus=0;
for TaxProductivityCorrelation=1:2 % Note: 3 is never actually reported in any results Tables.
    if TaxProductivityCorrelation==1 % Uncorrelated
        pistar_tau=[1-fractiontaxed;0;fractiontaxed]; % [Subsidy, nothing, tax]
        Params.upsilon=pistar_s.*(pistar_tau'); % upsilon is the 'potential entrants distribution' (actual entrants distribution is upsilon.*ebar)
    elseif TaxProductivityCorrelation==2
        [~,productivity_cutoff]=max((cumsum(pistar_s)>(1-fractiontaxed))); % '1-fractiontaxed', those below this will have tax set to zero, leaving fractiontaxed as those taxed.
        Params.upsilon=pistar_s.*[ones(n_s,1),zeros(n_s,1),ones(n_s,1)];
        Params.upsilon(1:(productivity_cutoff-1),3)=0; % No tax on low productivity
        Params.upsilon(productivity_cutoff:end,1)=0; % No subsidies on high productivity
    end
    for tau_c=1:3
        Params.taurate=tauratevec_capitaltax(tau_c);
        % Will use initial guess of
        Params.subsidyrate=0.5*Params.taurate; % Note: if Params.subsidyrate=Params.taurate then codes go wrong when Params.taurate=1 because this leads to a division by zero.
        % Note: this leads to zero subsidy when tax rate is zero (as initial guess).
%         if sum(Params.upsilon(:,1))>0 % if some firms are actually receiving subsidies
%             Params.subsidyrate=Params.taurate; % From RR2008 we know this is way too large as an initial guess, but it will do.
%         end
        FullFnsToEvaluate{1}=FnsToEvaluateFn_kbar.(taxbasevec{taxbase_c});
        FullFnsToEvaluate{2}=FnsToEvaluateFn_nbar.(taxbasevec{taxbase_c});
        FullFnsToEvaluate{3}=FnsToEvaluateFn_output.(taxbasevec{taxbase_c});
        FullFnsToEvaluate{4}=FnsToEvaluateFn_subsidy.(taxbasevec{taxbase_c});
        FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofsubsidised.(taxbasevec{taxbase_c});
        if index_taus==1 % Table 7, need to redefine Ys to 'non-taxed' instead of 'subsidised' as everywhere else.
            FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofnontaxed.(taxbasevec{taxbase_c});
        end
        
        fprintf('Now running RestucciaRogerson2008 for comination: \n')
        [fixedsubsidy_c+1,TaxProductivityCorrelation, fractiontaxed_c,index_taus+1, taxbase_c, tau_c]
        FullResults(fixedsubsidy_c+1,TaxProductivityCorrelation,fractiontaxed_c,index_taus+1,taxbase_c,tau_c).Output=RestucciaRogerson2008_Fn(fixedsubsidy_c, normalize_employment, n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,ReturnFn.(taxbasevec{taxbase_c}),DiscountFactorParamNames, ReturnFnParamNames, FullFnsToEvaluate, GEPriceParamNames, GeneralEqmEqnParamNames, GeneralEqmEqns, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions);
    end
end

% Labour taxes
fixedsubsidy_c=1;
fractiontaxed_c=5; fractiontaxed=fractiontaxedvec(fractiontaxed_c);
index_taus=0;
for tau_c=[1,6] 
    Params.taurate=tauratevec(tau_c);
    for taxbase_c=[1,3]
        for TaxProductivityCorrelation=1:2 % Note: 3 is never actually reported in any results Tables.
            if TaxProductivityCorrelation==1 % Uncorrelated
                pistar_tau=[1-fractiontaxed;0;fractiontaxed]; % [Subsidy, nothing, tax]
                Params.upsilon=pistar_s.*(pistar_tau'); % upsilon is the 'potential entrants distribution' (actual entrants distribution is upsilon.*ebar)
            elseif TaxProductivityCorrelation==2
                [~,productivity_cutoff]=max((cumsum(pistar_s)>(1-fractiontaxed))); % '1-fractiontaxed', those below this will have tax set to zero, leaving fractiontaxed as those taxed.
                Params.upsilon=pistar_s.*[ones(n_s,1),zeros(n_s,1),ones(n_s,1)];
                Params.upsilon(1:(productivity_cutoff-1),3)=0; % No tax on low productivity
                Params.upsilon(productivity_cutoff:end,1)=0; % No subsidies on high productivity
            end
            Params.subsidyrate=Params.taurate;
            
            FullFnsToEvaluate{1}=FnsToEvaluateFn_kbar.(taxbasevec{taxbase_c});
            FullFnsToEvaluate{2}=FnsToEvaluateFn_nbar.(taxbasevec{taxbase_c});
            FullFnsToEvaluate{3}=FnsToEvaluateFn_output.(taxbasevec{taxbase_c});
            FullFnsToEvaluate{4}=FnsToEvaluateFn_subsidy.(taxbasevec{taxbase_c});
            FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofsubsidised.(taxbasevec{taxbase_c});
            if index_taus==1 % Table 7, need to redefine Ys to 'non-taxed' instead of 'subsidised' as everywhere else.
                FullFnsToEvaluate{5}=FnsToEvaluateFn_outputofnontaxed.(taxbasevec{taxbase_c});
            end
            
            fprintf('Now running RestucciaRogerson2008 for comination: \n')
            [fixedsubsidy_c+1,TaxProductivityCorrelation, fractiontaxed_c,index_taus+1, taxbase_c, tau_c]
            FullResults(fixedsubsidy_c+1,TaxProductivityCorrelation,fractiontaxed_c,index_taus+1,taxbase_c,tau_c).Output=RestucciaRogerson2008_Fn(fixedsubsidy_c, normalize_employment, n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,ReturnFn.(taxbasevec{taxbase_c}),DiscountFactorParamNames, ReturnFnParamNames, FullFnsToEvaluate, GEPriceParamNames, GeneralEqmEqnParamNames, GeneralEqmEqns, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions);
        end
    end
end

save ./SavedOutput/RR2008_FullResults.mat FullResults
% load ./SavedOutput/RR2008_FullResults.mat FullResults


%% Have all the results, now just create the Tables themselves

%Table 3
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table3.tex', 'w');
fprintf(FID, 'Effects of idiosyncratic distortions -- uncorrelated case \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, 'Variable & \\multicolumn{4}{l}{$\\tau_t$ (tax rate on output)} \\\\ \\cline{2-5} \n');
fprintf(FID, ' & 0.1 & 0.2 & 0.3 & 0.4 \\\\ \\hline \n');
fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.Y/FullResults(1,1,5,1,1,1).Output.Y, FullResults(1,1,5,1,1,3).Output.Y/FullResults(1,1,5,1,1,1).Output.Y, FullResults(1,1,5,1,1,4).Output.Y/FullResults(1,1,5,1,1,1).Output.Y, FullResults(1,1,5,1,1,5).Output.Y/FullResults(1,1,5,1,1,1).Output.Y);
fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,3).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,4).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,5).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP);
fprintf(FID, 'Relative $E$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.Params.Ne/FullResults(1,1,5,1,1,1).Output.Params.Ne, FullResults(1,1,5,1,1,3).Output.Params.Ne/FullResults(1,1,5,1,1,1).Output.Params.Ne, FullResults(1,1,5,1,1,4).Output.Params.Ne/FullResults(1,1,5,1,1,1).Output.Params.Ne, FullResults(1,1,5,1,1,5).Output.Params.Ne/FullResults(1,1,5,1,1,1).Output.Params.Ne);
fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.Ys_divY, FullResults(1,1,5,1,1,3).Output.Ys_divY, FullResults(1,1,5,1,1,4).Output.Ys_divY, FullResults(1,1,5,1,1,5).Output.Ys_divY);
fprintf(FID, '$S/Y$          & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.SdivY, FullResults(1,1,5,1,1,3).Output.SdivY, FullResults(1,1,5,1,1,4).Output.SdivY, FullResults(1,1,5,1,1,5).Output.SdivY);
fprintf(FID, '$\\tau_s$      & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.tau_s, FullResults(1,1,5,1,1,3).Output.tau_s, FullResults(1,1,5,1,1,4).Output.tau_s, FullResults(1,1,5,1,1,5).Output.tau_s);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 4
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table4.tex', 'w');
fprintf(FID, 'Relative TFP -- uncorrelated case \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, 'Fraction of  & \\multicolumn{4}{l}{$\\tau_t$ (tax rate on output)} \\\\ \\cline{2-5} \n');
fprintf(FID, 'estabilishments taxed (\\%%) & 0.1 & 0.2 & 0.3 & 0.4 \\\\ \\hline \n');
fprintf(FID, '90 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,1,1,2).Output.TFP/FullResults(1,1,9,1,1,1).Output.TFP, FullResults(1,1,9,1,1,3).Output.TFP/FullResults(1,1,9,1,1,1).Output.TFP, FullResults(1,1,9,1,1,4).Output.TFP/FullResults(1,1,9,1,1,1).Output.TFP, FullResults(1,1,9,1,1,5).Output.TFP/FullResults(1,1,9,1,1,1).Output.TFP);
fprintf(FID, '80 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,8,1,1,2).Output.TFP/FullResults(1,1,8,1,1,1).Output.TFP, FullResults(1,1,8,1,1,3).Output.TFP/FullResults(1,1,8,1,1,1).Output.TFP, FullResults(1,1,8,1,1,4).Output.TFP/FullResults(1,1,8,1,1,1).Output.TFP, FullResults(1,1,8,1,1,5).Output.TFP/FullResults(1,1,8,1,1,1).Output.TFP);
fprintf(FID, '60 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,6,1,1,2).Output.TFP/FullResults(1,1,6,1,1,1).Output.TFP, FullResults(1,1,6,1,1,3).Output.TFP/FullResults(1,1,6,1,1,1).Output.TFP, FullResults(1,1,6,1,1,4).Output.TFP/FullResults(1,1,6,1,1,1).Output.TFP, FullResults(1,1,6,1,1,5).Output.TFP/FullResults(1,1,6,1,1,1).Output.TFP);
fprintf(FID, '50 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,1,2).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,3).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,4).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP, FullResults(1,1,5,1,1,5).Output.TFP/FullResults(1,1,5,1,1,1).Output.TFP);
fprintf(FID, '40 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,4,1,1,2).Output.TFP/FullResults(1,1,4,1,1,1).Output.TFP, FullResults(1,1,4,1,1,3).Output.TFP/FullResults(1,1,4,1,1,1).Output.TFP, FullResults(1,1,4,1,1,4).Output.TFP/FullResults(1,1,4,1,1,1).Output.TFP, FullResults(1,1,4,1,1,5).Output.TFP/FullResults(1,1,4,1,1,1).Output.TFP);
fprintf(FID, '20 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,2,1,1,2).Output.TFP/FullResults(1,1,2,1,1,1).Output.TFP, FullResults(1,1,2,1,1,3).Output.TFP/FullResults(1,1,2,1,1,1).Output.TFP, FullResults(1,1,2,1,1,4).Output.TFP/FullResults(1,1,2,1,1,1).Output.TFP, FullResults(1,1,2,1,1,5).Output.TFP/FullResults(1,1,2,1,1,1).Output.TFP);
fprintf(FID, '10 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,1,1,1,2).Output.TFP/FullResults(1,1,1,1,1,1).Output.TFP, FullResults(1,1,1,1,1,3).Output.TFP/FullResults(1,1,1,1,1,1).Output.TFP, FullResults(1,1,1,1,1,4).Output.TFP/FullResults(1,1,1,1,1,1).Output.TFP, FullResults(1,1,1,1,1,5).Output.TFP/FullResults(1,1,1,1,1,1).Output.TFP);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table 5
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table5.tex', 'w');
fprintf(FID, 'Effects of idiosyncratic distortions -- correlated case \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{4}{l}{$\\tau_t$ (tax rate on output)} \\\\ \\cline{2-5} \n');
fprintf(FID, ' & 0.1 & 0.2 & 0.3 & 0.4 \\\\ \\hline \n');
fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.Y/FullResults(1,2,5,1,1,1).Output.Y, FullResults(1,2,5,1,1,3).Output.Y/FullResults(1,2,5,1,1,1).Output.Y, FullResults(1,2,5,1,1,4).Output.Y/FullResults(1,2,5,1,1,1).Output.Y, FullResults(1,2,5,1,1,5).Output.Y/FullResults(1,2,5,1,1,1).Output.Y);
fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,3).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,4).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,5).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP);
fprintf(FID, 'Relative $E$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.Params.Ne/FullResults(1,2,5,1,1,1).Output.Params.Ne, FullResults(1,2,5,1,1,3).Output.Params.Ne/FullResults(1,2,5,1,1,1).Output.Params.Ne, FullResults(1,2,5,1,1,4).Output.Params.Ne/FullResults(1,2,5,1,1,1).Output.Params.Ne, FullResults(1,2,5,1,1,5).Output.Params.Ne/FullResults(1,2,5,1,1,1).Output.Params.Ne);
fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.Ys_divY, FullResults(1,2,5,1,1,3).Output.Ys_divY, FullResults(1,2,5,1,1,4).Output.Ys_divY, FullResults(1,2,5,1,1,5).Output.Ys_divY);
fprintf(FID, '$S/Y$          & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.SdivY, FullResults(1,2,5,1,1,3).Output.SdivY, FullResults(1,2,5,1,1,4).Output.SdivY, FullResults(1,2,5,1,1,5).Output.SdivY);
fprintf(FID, '$\\tau_s$      & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.tau_s, FullResults(1,2,5,1,1,3).Output.tau_s, FullResults(1,2,5,1,1,4).Output.tau_s, FullResults(1,2,5,1,1,5).Output.tau_s);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 6
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table6.tex', 'w');
fprintf(FID, 'Relative TFP -- correlated case \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, 'Fraction of  & \\multicolumn{4}{l}{$\\tau_t$ (tax rate on output)} \\\\ \\cline{2-5} \n');
fprintf(FID, 'estabilishments taxed (\\%%) & 0.1 & 0.2 & 0.3 & 0.4 \\\\ \\hline \n');
fprintf(FID, '90 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,1,1,2).Output.TFP/FullResults(1,2,9,1,1,1).Output.TFP, FullResults(1,2,9,1,1,3).Output.TFP/FullResults(1,2,9,1,1,1).Output.TFP, FullResults(1,2,9,1,1,4).Output.TFP/FullResults(1,2,9,1,1,1).Output.TFP, FullResults(1,2,9,1,1,5).Output.TFP/FullResults(1,2,9,1,1,1).Output.TFP);
fprintf(FID, '80 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,8,1,1,2).Output.TFP/FullResults(1,2,8,1,1,1).Output.TFP, FullResults(1,2,8,1,1,3).Output.TFP/FullResults(1,2,8,1,1,1).Output.TFP, FullResults(1,2,8,1,1,4).Output.TFP/FullResults(1,2,8,1,1,1).Output.TFP, FullResults(1,2,8,1,1,5).Output.TFP/FullResults(1,2,8,1,1,1).Output.TFP);
fprintf(FID, '60 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,6,1,1,2).Output.TFP/FullResults(1,2,6,1,1,1).Output.TFP, FullResults(1,2,6,1,1,3).Output.TFP/FullResults(1,2,6,1,1,1).Output.TFP, FullResults(1,2,6,1,1,4).Output.TFP/FullResults(1,2,6,1,1,1).Output.TFP, FullResults(1,2,6,1,1,5).Output.TFP/FullResults(1,2,6,1,1,1).Output.TFP);
fprintf(FID, '50 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,5,1,1,2).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,3).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,4).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP, FullResults(1,2,5,1,1,5).Output.TFP/FullResults(1,2,5,1,1,1).Output.TFP);
fprintf(FID, '40 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,4,1,1,2).Output.TFP/FullResults(1,2,4,1,1,1).Output.TFP, FullResults(1,2,4,1,1,3).Output.TFP/FullResults(1,2,4,1,1,1).Output.TFP, FullResults(1,2,4,1,1,4).Output.TFP/FullResults(1,2,4,1,1,1).Output.TFP, FullResults(1,2,4,1,1,5).Output.TFP/FullResults(1,2,4,1,1,1).Output.TFP);
fprintf(FID, '20 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,2,1,1,2).Output.TFP/FullResults(1,2,2,1,1,1).Output.TFP, FullResults(1,2,2,1,1,3).Output.TFP/FullResults(1,2,2,1,1,1).Output.TFP, FullResults(1,2,2,1,1,4).Output.TFP/FullResults(1,2,2,1,1,1).Output.TFP, FullResults(1,2,2,1,1,5).Output.TFP/FullResults(1,2,2,1,1,1).Output.TFP);
fprintf(FID, '10 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,1,1,1,2).Output.TFP/FullResults(1,2,1,1,1,1).Output.TFP, FullResults(1,2,1,1,1,3).Output.TFP/FullResults(1,2,1,1,1,1).Output.TFP, FullResults(1,2,1,1,1,4).Output.TFP/FullResults(1,2,1,1,1,1).Output.TFP, FullResults(1,2,1,1,1,5).Output.TFP/FullResults(1,2,1,1,1,1).Output.TFP);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% More tables

% % % Uncorrelated version of Table 7
% % %Table 7 (Note, this has 'fraction exempt', as opposed to Table 6 which was 'fraction taxed')
% % FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table7.tex', 'w');
% % fprintf(FID, 'Taxing all but some exempt establishments ($\\tau_t = 0.40$) \\\\ \n');
% % fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllll} \n \\hline \\hline \n');
% % fprintf(FID, 'Variable & \\multicolumn{5}{l}{Establishments exempt (\\%%)} \\\\ \\cline{2-6} \n');
% % fprintf(FID, ' & 10 & 30 & 50 & 70 & 90 \\\\ \\hline \n');
% % fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.Y/FullResults(1,1,9,2,1,1).Output.Y,                 FullResults(1,1,7,2,1,5).Output.Y/FullResults(1,1,7,2,1,1).Output.Y,                 FullResults(1,1,5,2,1,5).Output.Y/FullResults(1,1,5,2,1,1).Output.Y,                 FullResults(1,1,3,2,1,5).Output.Y/FullResults(1,1,5,2,1,1).Output.Y,                 FullResults(1,1,1,2,1,5).Output.Y/FullResults(1,1,5,2,1,1).Output.Y);
% % fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.TFP/FullResults(1,1,9,2,1,1).Output.TFP,             FullResults(1,1,7,2,1,5).Output.TFP/FullResults(1,1,7,2,1,1).Output.TFP,             FullResults(1,1,5,2,1,5).Output.TFP/FullResults(1,1,5,2,1,1).Output.TFP,             FullResults(1,1,3,2,1,5).Output.TFP/FullResults(1,1,5,2,1,1).Output.TFP,             FullResults(1,1,1,2,1,5).Output.TFP/FullResults(1,1,5,2,1,1).Output.TFP);
% % fprintf(FID, 'Relative $E$   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.Params.Ne/FullResults(1,1,9,2,1,1).Output.Params.Ne, FullResults(1,1,7,2,1,5).Output.Params.Ne/FullResults(1,1,7,2,1,1).Output.Params.Ne, FullResults(1,1,5,2,1,5).Output.Params.Ne/FullResults(1,1,5,2,1,1).Output.Params.Ne, FullResults(1,1,3,2,1,5).Output.Params.Ne/FullResults(1,1,5,2,1,1).Output.Params.Ne, FullResults(1,1,1,2,1,5).Output.Params.Ne/FullResults(1,1,5,2,1,1).Output.Params.Ne);
% % fprintf(FID, 'Relative w     & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.Params.w/FullResults(1,1,9,2,1,1).Output.Params.w,   FullResults(1,1,7,2,1,5).Output.Params.w/FullResults(1,1,7,2,1,1).Output.Params.w,   FullResults(1,1,5,2,1,5).Output.Params.w/FullResults(1,1,5,2,1,1).Output.Params.w,   FullResults(1,1,3,2,1,5).Output.Params.w/FullResults(1,1,5,2,1,1).Output.Params.w,   FullResults(1,1,1,2,1,5).Output.Params.w/FullResults(1,1,5,2,1,1).Output.Params.w);
% % fprintf(FID, 'Relative K     & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.K/FullResults(1,1,9,2,1,1).Output.K,                 FullResults(1,1,7,2,1,5).Output.K/FullResults(1,1,7,2,1,1).Output.K,                 FullResults(1,1,5,2,1,5).Output.K/FullResults(1,1,5,2,1,1).Output.K,                 FullResults(1,1,3,2,1,5).Output.K/FullResults(1,1,5,2,1,1).Output.K,                 FullResults(1,1,1,2,1,5).Output.K/FullResults(1,1,5,2,1,1).Output.K);
% % fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,9,2,1,5).Output.Ys_divY, FullResults(1,1,7,2,1,5).Output.Ys_divY, FullResults(1,1,5,2,1,5).Output.Ys_divY, FullResults(1,1,3,2,1,5).Output.Ys_divY, FullResults(1,1,1,2,1,5).Output.Ys_divY);
% % fprintf(FID, '\\hline \n \\end{tabular*} \n');
% % % fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% % % fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% % % fprintf(FID, '}} \\end{minipage}');
% % fclose(FID);
% % % Note: Since no-one is receiving subsidies the subsidy rate is anyway zero
% % % and not optimized to ensure K equals Kbaseline. Hence there is no need/point to
% % % use the fixedsubsidy_c=1 case for creating this table.

% Correlated version of Table 7
%Table 7 (Note, this has 'fraction exempt', as opposed to Table 6 which was 'fraction taxed')
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table7.tex', 'w');
fprintf(FID, 'Taxing all but some exempt establishments ($\\tau_t = 0.40$) \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllll} \n \\hline \\hline \n');
fprintf(FID, 'Variable & \\multicolumn{5}{l}{Establishments exempt (\\%%)} \\\\ \\cline{2-6} \n');
fprintf(FID, ' & 10 & 30 & 50 & 70 & 90 \\\\ \\hline \n');
fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.Y/FullResults(1,2,5,2,1,1).Output.Y,                 FullResults(1,2,7,2,1,5).Output.Y/FullResults(1,2,7,2,1,1).Output.Y,                 FullResults(1,2,5,2,1,5).Output.Y/FullResults(1,2,5,2,1,1).Output.Y,                 FullResults(1,2,3,2,1,5).Output.Y/FullResults(1,2,5,2,1,1).Output.Y,                 FullResults(1,2,1,2,1,5).Output.Y/FullResults(1,2,5,2,1,1).Output.Y);
fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.TFP/FullResults(1,2,5,2,1,1).Output.TFP,             FullResults(1,2,7,2,1,5).Output.TFP/FullResults(1,2,7,2,1,1).Output.TFP,             FullResults(1,2,5,2,1,5).Output.TFP/FullResults(1,2,5,2,1,1).Output.TFP,             FullResults(1,2,3,2,1,5).Output.TFP/FullResults(1,2,5,2,1,1).Output.TFP,             FullResults(1,2,1,2,1,5).Output.TFP/FullResults(1,2,5,2,1,1).Output.TFP);
fprintf(FID, 'Relative $E$   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.Params.Ne/FullResults(1,2,5,2,1,1).Output.Params.Ne, FullResults(1,2,7,2,1,5).Output.Params.Ne/FullResults(1,2,7,2,1,1).Output.Params.Ne, FullResults(1,2,5,2,1,5).Output.Params.Ne/FullResults(1,2,5,2,1,1).Output.Params.Ne, FullResults(1,2,3,2,1,5).Output.Params.Ne/FullResults(1,2,5,2,1,1).Output.Params.Ne, FullResults(1,2,1,2,1,5).Output.Params.Ne/FullResults(1,2,5,2,1,1).Output.Params.Ne);
fprintf(FID, 'Relative w     & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.Params.w/FullResults(1,2,5,2,1,1).Output.Params.w,   FullResults(1,2,7,2,1,5).Output.Params.w/FullResults(1,2,7,2,1,1).Output.Params.w,   FullResults(1,2,5,2,1,5).Output.Params.w/FullResults(1,2,5,2,1,1).Output.Params.w,   FullResults(1,2,3,2,1,5).Output.Params.w/FullResults(1,2,5,2,1,1).Output.Params.w,   FullResults(1,2,1,2,1,5).Output.Params.w/FullResults(1,2,5,2,1,1).Output.Params.w);
fprintf(FID, 'Relative K     & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.K/FullResults(1,2,5,2,1,1).Output.K,                 FullResults(1,2,7,2,1,5).Output.K/FullResults(1,2,7,2,1,1).Output.K,                 FullResults(1,2,5,2,1,5).Output.K/FullResults(1,2,5,2,1,1).Output.K,                 FullResults(1,2,3,2,1,5).Output.K/FullResults(1,2,5,2,1,1).Output.K,                 FullResults(1,2,1,2,1,5).Output.K/FullResults(1,2,5,2,1,1).Output.K);
fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,2,9,2,1,5).Output.Ys_divY, FullResults(1,2,7,2,1,5).Output.Ys_divY, FullResults(1,2,5,2,1,5).Output.Ys_divY, FullResults(1,2,3,2,1,5).Output.Ys_divY, FullResults(1,2,1,2,1,5).Output.Ys_divY);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);
% Note: Since no-one is receiving subsidies the subsidy rate is anyway zero
% and not optimized to ensure K equals Kbaseline. Hence there is no need/point to
% use the fixedsubsidy_c=1 case for creating this table.

% Following lines (prior to Table 8) are just me playing around with Table
% 7 numbers trying to figure out if so other denominators seemed to improve
% things. And checking a few other outputs to see if there was anything
% obviously incorrect.
d1=[1,2,5,1,1,1]; % Allows you to easily change denominator of first column
% d2=d1; d3=d1; d4=d1; d5=d1;
d2=[1,2,5,1,1,1]; d3=[1,2,5,1,1,1]; d4=[1,2,5,1,1,1]; d5=[1,2,5,1,1,1]; % ...and the other columns
[FullResults(1,2,9,2,1,5).Output.Y/FullResults(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6)).Output.Y,                 FullResults(1,2,7,2,1,5).Output.Y/FullResults(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6)).Output.Y,                 FullResults(1,2,5,2,1,5).Output.Y/FullResults(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6)).Output.Y,                 FullResults(1,2,3,2,1,5).Output.Y/FullResults(d4(1),d4(2),d4(3),d4(4),d4(5),d4(6)).Output.Y,                 FullResults(1,2,1,2,1,5).Output.Y/FullResults(d5(1),d5(2),d5(3),d5(4),d5(5),d5(6)).Output.Y;...
 FullResults(1,2,9,2,1,5).Output.TFP/FullResults(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6)).Output.TFP,             FullResults(1,2,7,2,1,5).Output.TFP/FullResults(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6)).Output.TFP,             FullResults(1,2,5,2,1,5).Output.TFP/FullResults(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6)).Output.TFP,             FullResults(1,2,3,2,1,5).Output.TFP/FullResults(d4(1),d4(2),d4(3),d4(4),d4(5),d4(6)).Output.TFP,             FullResults(1,2,1,2,1,5).Output.TFP/FullResults(d5(1),d5(2),d5(3),d5(4),d5(5),d5(6)).Output.TFP;...
 FullResults(1,2,9,2,1,5).Output.Params.Ne/FullResults(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6)).Output.Params.Ne, FullResults(1,2,7,2,1,5).Output.Params.Ne/FullResults(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6)).Output.Params.Ne, FullResults(1,2,5,2,1,5).Output.Params.Ne/FullResults(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6)).Output.Params.Ne, FullResults(1,2,3,2,1,5).Output.Params.Ne/FullResults(d4(1),d4(2),d4(3),d4(4),d4(5),d4(6)).Output.Params.Ne, FullResults(1,2,1,2,1,5).Output.Params.Ne/FullResults(d5(1),d5(2),d5(3),d5(4),d5(5),d5(6)).Output.Params.Ne;...
 FullResults(1,2,9,2,1,5).Output.Params.w/FullResults(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6)).Output.Params.w,   FullResults(1,2,7,2,1,5).Output.Params.w/FullResults(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6)).Output.Params.w,   FullResults(1,2,5,2,1,5).Output.Params.w/FullResults(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6)).Output.Params.w,   FullResults(1,2,3,2,1,5).Output.Params.w/FullResults(d4(1),d4(2),d4(3),d4(4),d4(5),d4(6)).Output.Params.w,   FullResults(1,2,1,2,1,5).Output.Params.w/FullResults(d5(1),d5(2),d5(3),d5(4),d5(5),d5(6)).Output.Params.w;...
 FullResults(1,2,9,2,1,5).Output.K/FullResults(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6)).Output.K,                 FullResults(1,2,7,2,1,5).Output.K/FullResults(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6)).Output.K,                 FullResults(1,2,5,2,1,5).Output.K/FullResults(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6)).Output.K,                 FullResults(1,2,3,2,1,5).Output.K/FullResults(d4(1),d4(2),d4(3),d4(4),d4(5),d4(6)).Output.K,                 FullResults(1,2,1,2,1,5).Output.K/FullResults(d5(1),d5(2),d5(3),d5(4),d5(5),d5(6)).Output.K;...
 FullResults(1,2,9,2,1,5).Output.Ys_divY, FullResults(1,2,7,2,1,5).Output.Ys_divY, FullResults(1,2,5,2,1,5).Output.Ys_divY, FullResults(1,2,3,2,1,5).Output.Ys_divY, FullResults(1,2,1,2,1,5).Output.Ys_divY]

[FullResults(1,2,9,2,1,5).Output.SdivY, FullResults(1,2,7,2,1,5).Output.SdivY, FullResults(1,2,5,2,1,5).Output.SdivY, FullResults(1,2,3,2,1,5).Output.SdivY, FullResults(1,2,1,2,1,5).Output.SdivY]
[FullResults(1,2,9,2,1,5).Output.tau_s, FullResults(1,2,7,2,1,5).Output.tau_s, FullResults(1,2,5,2,1,5).Output.tau_s, FullResults(1,2,3,2,1,5).Output.tau_s, FullResults(1,2,1,2,1,5).Output.tau_s]
sum(sum(FullResults(1,2,9,2,1,5).Output.Params.ebar))
sum(sum(FullResults(1,2,5,2,1,5).Output.Params.ebar))
sum(sum(FullResults(1,2,1,2,1,5).Output.Params.ebar))

rate1_c=2;
rate2_c=3;
%Table 8
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table8.tex', 'w');
fprintf(FID, 'Idiosyncratic distortions to capital rental rates \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, 'Variable & \\multicolumn{2}{l}{Uncorrelated} & \\multicolumn{2}{l}{Correlated} \\\\ \\cline{2-3} \\cline{4-5} \n');
fprintf(FID, ' & $\\tau_t=%8.1f $ & $\\tau_t=%8.1f $ & $\\tau_t=%8.1f $ & $\\tau_t=%8.1f $ \\\\ \\hline \n', tauratevec_capitaltax(rate1_c), tauratevec_capitaltax(rate2_c), tauratevec_capitaltax(rate1_c), tauratevec_capitaltax(rate2_c));
fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,2,rate1_c).Output.Y/FullResults(1,1,5,1,2,1).Output.Y,                 FullResults(1,1,5,1,2,rate2_c).Output.Y/FullResults(1,1,5,1,2,1).Output.Y,                 FullResults(1,2,5,1,2,rate1_c).Output.Y/FullResults(1,2,5,1,2,1).Output.Y,                 FullResults(1,2,5,1,2,rate2_c).Output.Y/FullResults(1,2,5,1,2,1).Output.Y);
fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,2,rate1_c).Output.TFP/FullResults(1,1,5,1,2,1).Output.TFP,             FullResults(1,1,5,1,2,rate2_c).Output.TFP/FullResults(1,1,5,1,2,1).Output.TFP,             FullResults(1,2,5,1,2,rate1_c).Output.TFP/FullResults(1,2,5,1,2,1).Output.TFP,             FullResults(1,2,5,1,2,rate2_c).Output.TFP/FullResults(1,2,5,1,2,1).Output.TFP);
fprintf(FID, 'Relative $E$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,2,rate1_c).Output.Params.Ne/FullResults(1,1,5,1,2,1).Output.Params.Ne, FullResults(1,1,5,1,2,rate2_c).Output.Params.Ne/FullResults(1,1,5,1,2,1).Output.Params.Ne, FullResults(1,2,5,1,2,rate1_c).Output.Params.Ne/FullResults(1,2,5,1,2,1).Output.Params.Ne, FullResults(1,2,5,1,2,rate2_c).Output.Params.Ne/FullResults(1,2,5,1,2,1).Output.Params.Ne);
fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,2,rate1_c).Output.Ys_divY, FullResults(1,1,5,1,2,rate2_c).Output.Ys_divY, FullResults(1,2,5,1,2,rate1_c).Output.Ys_divY, FullResults(1,2,5,1,2,rate2_c).Output.Ys_divY);
fprintf(FID, '$S/Y$          & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(1,1,5,1,2,rate1_c).Output.SdivY,   FullResults(1,1,5,1,2,rate2_c).Output.SdivY,   FullResults(1,2,5,1,2,rate1_c).Output.SdivY, FullResults(1,2,5,1,2,rate2_c).Output.SdivY);
fprintf(FID, '$\\tau_s$     & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',  FullResults(1,1,5,1,2,rate1_c).Output.tau_s,   FullResults(1,1,5,1,2,rate2_c).Output.tau_s,   FullResults(1,2,5,1,2,rate1_c).Output.tau_s, FullResults(1,2,5,1,2,rate2_c).Output.tau_s);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table 9
FID = fopen('./SavedOutput/LatexInputs/RestucciaRogerson2008_Table9.tex', 'w');
fprintf(FID, 'Idiosyncratic distortions--outputs vs wages \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, 'Variable & \\multicolumn{2}{l}{Uncorrelated} & \\multicolumn{2}{l}{Correlated} \\\\ \\cline{2-3} \\cline{4-5} \n');
fprintf(FID, ' & Output & Wages & Output & Wages \\\\ \\hline \n');
fprintf(FID, 'Relative $Y$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',   FullResults(2,1,5,1,1,6).Output.Y/FullResults(2,1,5,1,1,1).Output.Y,                 FullResults(2,1,5,1,3,6).Output.Y/FullResults(2,1,5,1,3,1).Output.Y,                 FullResults(2,2,5,1,1,6).Output.Y/FullResults(2,2,5,1,1,1).Output.Y,                 FullResults(2,2,5,1,3,6).Output.Y/FullResults(2,2,5,1,3,1).Output.Y);
fprintf(FID, 'Relative TFP   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',   FullResults(2,1,5,1,1,6).Output.TFP/FullResults(2,1,5,1,1,1).Output.TFP,             FullResults(2,1,5,1,3,6).Output.TFP/FullResults(2,1,5,1,3,1).Output.TFP,             FullResults(2,2,5,1,1,6).Output.TFP/FullResults(2,2,5,1,1,1).Output.TFP,             FullResults(2,2,5,1,3,6).Output.TFP/FullResults(2,2,5,1,3,1).Output.TFP);
fprintf(FID, 'Relative $K$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',   FullResults(2,1,5,1,1,6).Output.K/FullResults(2,1,5,1,1,1).Output.K,                 FullResults(2,1,5,1,3,6).Output.K/FullResults(2,1,5,1,3,1).Output.K,                 FullResults(2,2,5,1,1,6).Output.K/FullResults(2,2,5,1,1,1).Output.K,                 FullResults(2,2,5,1,3,6).Output.K/FullResults(2,2,5,1,3,1).Output.K);
fprintf(FID, 'Relative $E$     & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(2,1,5,1,1,6).Output.Params.Ne/FullResults(2,1,5,1,1,1).Output.Params.Ne, FullResults(2,1,5,1,3,6).Output.Params.Ne/FullResults(2,1,5,1,3,1).Output.Params.Ne, FullResults(2,2,5,1,1,6).Output.Params.Ne/FullResults(2,2,5,1,1,1).Output.Params.Ne, FullResults(2,2,5,1,3,6).Output.Params.Ne/FullResults(2,2,5,1,3,1).Output.Params.Ne);
fprintf(FID, 'Relative $w$     & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(2,1,5,1,1,6).Output.Params.w/FullResults(2,1,5,1,1,1).Output.Params.w,   FullResults(2,1,5,1,3,6).Output.Params.w/FullResults(2,1,5,1,3,1).Output.Params.w,   FullResults(2,2,5,1,1,6).Output.Params.w/FullResults(2,2,5,1,1,1).Output.Params.w,   FullResults(2,2,5,1,3,6).Output.Params.w/FullResults(2,2,5,1,3,1).Output.Params.w);
fprintf(FID, '$Y_s /Y$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(2,1,5,1,1,6).Output.Ys_divY, FullResults(2,1,5,1,3,6).Output.Ys_divY, FullResults(2,2,5,1,1,6).Output.Ys_divY, FullResults(2,2,5,1,3,6).Output.Ys_divY);
fprintf(FID, '$S/Y$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', FullResults(2,1,5,1,1,6).Output.SdivY, FullResults(2,1,5,1,3,6).Output.SdivY, FullResults(2,2,5,1,1,6).Output.SdivY, FullResults(2,2,5,1,3,6).Output.SdivY);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);
















