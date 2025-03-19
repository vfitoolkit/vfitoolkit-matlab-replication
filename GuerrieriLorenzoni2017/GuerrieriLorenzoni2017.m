% Replicates the results of Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap

% QUESTIONS OUTSTANDING
% Do GL2017 change T for Fiscal Policy transitions? (I do, increase to T=50)
% Is my interpretation of the Fischer deflation experiment correct? (that B increases and then stays at this new higher level)
% Do GL2017 change tauchen q when reducing the number of grid points to 5?, in durable goods model? (probably, they appear to use Tauchen-Hussey, not Tauchen as reported in paper)
% Codes provided by GL2017 on website of Lorenzoni only cover baseline model for flexible and NK transitions. (so do not answer the above questions)

% To run on server I need the following line
addpath(genpath('./MatlabToolkits/'))

% Set all following to 0 to run everything, otherwise skip the parts set to 1
SkipInitialFinal=0
SkipBaselineTransPath=0
SkipPanelData=1
SkipNKTransPath=0
SkipAlternativeCalibrations=0
SkipFiscalTransPath=0
SkipFisherDebtDeflation=0
SkipDurables=1

% GL2017 take exogenous shock grid from a previous paper, and so it is just
% loaded from a file. Following lines enforce this (rather than try
% recreate it using quadrature methods).
OverwriteWithGL2017Grid=1
OverwriteWithGL2017CodeParams=1

CreateFigures=1 % Set to zero so can run on server

% Note: to run you need the following subfolders (to save results into)
% ./SavedOutput/Graphs
% ./SavedOutput/LatexInputs
% ./SavedOutput/GL2017Figs
% Note: the third of these is only needed if you set CreateFigures=0 (it
% will end up containing files that contain the objects you need to create
% the figures).

%% To translate Guerrieri & Lorenzoni (2017) into the standard setup of VFI Toolkit I use following:
% d variables: n_it
% aprime variables: b_it+1
% a variables: b_it
% z variables: theta_it, e_it


%% Set some basic variables

n_d=41;
n_a=901; % Guerrieri & Lorenzoni (2017) use 200 points for agent distribution; VFI is done the same but they use the EGM (endo grid method) and then use ConesaKrueger style probabilistic weights to nearest grid point for agent dist simulation
n_theta=13; % Guerrieri & Lorenzoni (2017), pg 1438, states that they use a "12 state markov chain, following the approach in Tauchen (1986)". This description can be misinterpreted as theta is in fact the combination of a two-state markov on employed/unemployed with a 12-state Tauchen approx of AR(1) applied to employment. (This is clear from their codes) [The precise wording of GL2017 is correct, just easily misread.]

n_r=1001; % for some graphs

%Parameters (mostly from G&L2017, Table 1)
Params.beta=0.9711; % Model period is one-quarter of a year.

Params.gamma=4; % Coefficient of relative risk aversion
Params.eta=1.5; % Curvature of utility from leisure
Params.psi=12.48; % Coefficient on leisure in utility

Params.pi_eu=0.0573; % Transition to unemployment (the last three comes from their code, not paper)
Params.pi_ue=0.882; % Transition to employment
% Note: you cannot change pi_eu or pi_ue without changing the return
% function because the implied uncondtional probability of being unemployed
% is hard coded there.
Params.rho=0.967; % Persistence of productivity shock (G&L2017 call this rho)
Params.sigmasq_epsilon=0.017; % Variance of the shock to the log-AR(1) process on labour productivity.
Params.tauchenq=2.1; % I have reverse engineered this value from the grid of GL2017 (copy of their grid is included in their codes). They themselves reverse engineered choice of roughly 2.1 so that the variance of the resulting process (variance of 12-state markov logz) is as close as possible to what it should be (variance of true AR(1) logz).
% This choice of tauchenq is crucial to the results/replication of GL2017. It means
% that the min and max productivity shocks in the model, while they have the correct variance, have a range which
% is roughly just +-standard deviation. Because this range of shocks is so (empirically unrealistically) small
% lots of agents in the model end up near the borrowing constraint (which is not true using, e.g., tauchenq=3).
% Without lots of agents near the borrowing constraint the credit crisis (change in phi) has very little effect on the model.

Params.v=0.1; % Unemployment benefit
Params.B=1.6; % Bond supply
Params.Bprime=Params.B; % Bond supply is unchanging (for most of the paper); this is needed as it is part of government budget constraint that determines lump-sum tax tau_t. (Obviously B=Bprime must hold in any stationary general eqm.)
% For creating graphs later it will be useful to store both the initial and final values of phi
Params.phi_initial=0.959;
Params.phi_final=0.525;
Params.phi=Params.phi_initial; % Borrowing limit

% Params.r % General equilibrium interest rate (I use this instead of q)
% Params.tau % Lump-sum taxes, determined in general equilibrium (can implement it directly inside the ReturnFn)

Params.omega=0; % This is not actually needed for anything until we get to the 'Sticky Wage' model (Section 4 of GL2017, pg 1450)
% In the New Keynesian Sticky Wages model this appears as a wedge that
% increases the utility of leisure (intended to capture a fall in labor
% demand as a result of wages not falling). See GL2017 for explanation.


%% GL2017 codes param values
% The codes 'compute_steady_states.m' provided on Lorenzoni's website
% actually use the following calibrated values (contents of par.mat that
% comes with the codes).
if OverwriteWithGL2017CodeParams==1
    Params.beta=0.9774;
    Params.v=0.1670;
    Params.B=2.6712;
    Params.Bprime=Params.B;
    Params.psi=15.8819;
    % For creating graphs later it will be useful to store both the initial and final values of phi
    Params.phi_initial=1.6005;
    Params.phi_final=0.8767;
    Params.phi=Params.phi_initial; % Borrowing limit
end
% Note that the values of psi, B, phi_initial, and phi_final in particular are
% substantially different from those reported in Table 1 of GL2017. The
% differences for B, phi_initial and phi_final are simply because GL2017 do
% not actually ever report these, they report these divided by annual GDP
% (roughly 1.64; e.g. in paper phi_initial=1.6/1.64=0.925). However this is
% never made clear in the paper (the targets are described as being in these 'as fraction of
% annual GDP' terms, but there is no indication the actual reported
% parameter values are).


%% Some Toolkit options

% Create markov process for the exogenous income (based on idea of employment and unemployment states, following Imrohoroglu, 1989).
[theta1_grid, pi_theta1]=discretizeAR1_Tauchen(0,Params.rho,sqrt(Params.sigmasq_epsilon),n_theta-1,Params.tauchenq); %[states, transmatrix], transmatix is (z,zprime)
z_grid=[0; exp(theta1_grid)];
% G&L2017, pg 1438 "when first employed, workers draw theta from its unconditional distribution"; so here compute the unconditional distribution
pistar_theta1=ones(n_theta-1,1)/(n_theta-1);
for ii=1:10^4 % There is a more efficient form to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_theta1=pi_theta1'*pistar_theta1; 
end

% Floden & Linde (2001) report annual values for the AR(1) on log(z) as
% rho_FL=0.9136, sigmasq_epsilon_FL=0.0426; can calculate sigmasq_z_FL=sigmasq_epsilon_FL/(1-rho_FL^2)=0.2577;
% GL2017, footnote 11 give the formulae for annual rho and sigmasqz in terms of the quarterly:
rho=Params.rho;
sigmasq_epsilon=Params.sigmasq_epsilon;
sigmasq_epsilon_annual=(1/(4^2))*(4+6*rho+4*rho^2+2*rho^3)*(sigmasq_epsilon/(1-rho^2)); % Gives 0.2572 which is very close to FL
autocovariance_annual=(1/(4^2))*(rho+2*rho^2+3*rho^3+4*rho^4+3*rho^5+2*rho^6+rho^7)*(sigmasq_epsilon/(1-rho^2));
rho_annual=autocovariance_annual/sigmasq_epsilon_annual; % Gives 0.9127, which is close to FL. Note that this depends only on quarterly rho (the quarterly sigmasq_epsilon disappears when doing the division of annual autocovariance by annual variance)
% These check out, so the values in GL2017 for rho and sigmasq_epsilon are correct.

pi_z=[(1-Params.pi_ue), Params.pi_ue*pistar_theta1'; Params.pi_eu*ones(n_theta-1,1),(1-Params.pi_eu)*pi_theta1];
% Rows did not sum to one due to rounding errors at order of 10^(-11), fix this
pi_z=pi_z./sum(pi_z,2);
pistar_z=ones(n_theta,1)/n_theta;
for ii=1:10^4 %  % There is a more efficient way to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_z=pi_z'*pistar_z; % Formula could be used to find stationary dist of the employment unemployment process, then just combine with stationary dist of theta1, which is already calculated
end
% "The average level of theta is chosen so that yearly output in the initial steady state is normalized to 1"
z_grid=z_grid/sum(z_grid.*pistar_z);
% Double-check that this is 1
% sum(z_grid.*pistar_z)

% That the "normalized to 1" refers to E[theta] and not E[n*theta] is clear from setting
% v=0.1 to satisfy "For the unemployment benefit, we also follow Shimer
% (2005) and set it to 40% of average labor income." (pg 1438)
% Note, they do not actually ever normalize to 1 in the codes, GL2017 has E[theta]=1.0718

% Either their reported parameter values are slightly incorrect or the
% Tauchen method code they used is not accurate, as they get slightly
% different grid, so I will just load theirs directly
if OverwriteWithGL2017Grid==1
    load inc_process.mat
    theta = [0; exp(x)];
    z_grid=theta;
    
    tol_dist = 1e-10; % tolerance lvl for distribution
    S     = length(theta);
    fin   = 0.8820;        % job-finding probability
    sep   = 0.0573;        % separation probability
    % new transition matrix
    Pr = [1-fin, fin*pr; sep*ones(S-1, 1), (1-sep)*Pr];
    % find new invariate distribution
    pr  = [0, pr];
    dif = 1;
    while dif > tol_dist
      pri = pr*Pr;
      dif = max(abs(pri-pr));
      pr  = pri;
    end
    
    pi_z=Pr;
    pistar_z=pr';
    
    % Note: GL2017 codes do not do the normalization of z_grid, they just have
    % E[theta]=1.07 rather than 1 as stated in the paper.
end

%% Grids
% Set grid for asset holdings
Params.alowerbar=-1.25*Params.phi; % This seems reasonable (No-one can go below -Params.phi in any case). Note that Fischer deflation experiment (second last part of paper) won't work unless this is below -1.1*phi
Params.aupperbar=20; % Not clear exactly what value is appropriate, have gone with this, and checked that increasing it makes no difference.
a_grid=(Params.aupperbar-Params.alowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a)')-1)+Params.alowerbar;
% GL2017 codes use alowerbar of -2 and aupperbar of 50, but since no-one
% goes anywhere near 50 this seems excessive (even for them 0.9995 of population ends up below 20.8690). 
% Of their 200 points, only 183 are actually above -phi=-0.925 (only 171 are above final -phi=-0.525)
% Because GL2017 use a weighted probablility to approach to put agents on
% grid (policy fn allows off grid choices), and do not enforce the budget constraint on this, they in fact end
% up with 0.19% of the initial stationary distribution violating the budget constraint! (0.9981 mass is above phi1)
% For final distributions this becomes 0.96% violating budget constraint (0.9904 mass is above phi1)
% Of their 200 points, only 115 are between -phi1 and 20.86.
% The actual magnitude of the numerical error induced by this is small enough not to matter in their final results.

% Set grid for interest rate, r
r_grid=linspace(-0.02,0.04,n_r)'; % This grid is substantially wider than the actual likely equilibrium values and so is somewhat overkill.
Params.r=0.02; % initial guess (for when not using grid)
% % Set grid for tax rate
% tau_grid=linspace(); % Can calculate this from the gov budget constraint

%Bring model into the notational conventions used by the toolkit
d_grid=linspace(0,1,n_d)'; % Labor supply
p_grid=r_grid;

n_z=n_theta;
n_p=n_r;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluate.A = @(d, aprime,a,z) a; % Aggregate assets (which is this periods state)

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqns.BondMarketClearence = @(A,B) A-B; %The requirement that the aggregate assets (lending and borrowing; by government and private) equal zero

%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d, aprime, a, z, r, gamma, psi, eta, phi, v,B,Bprime,omega) ...
    GuerrieriLorenzoni2017_ReturnFn(d, aprime, a, z,r, gamma, psi, eta, phi, v,B,Bprime,omega);

%%
vfoptions.divideandconquer=1; % for transition path, turn on divide-and-conquer
simoptions=struct(); % use defaults

%% Solve the initial stationary equilibrium

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};

heteroagentoptions.verbose=1;
heteroagentoptions.pgrid=p_grid;
if SkipInitialFinal==0
    disp('Calculating prices corresponding to the stationary eqm')
    p_eqm_initial=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params.r=p_eqm_initial.r;

    % To be able to draw Figure 2, going to do redo, but with a grid on r
    [~,~, MarketClearance_initial_grid]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

    disp('Calculating various equilibrium objects')
    [~,Policy_initial]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z, simoptions);
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

    save ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params p_eqm_initial MarketClearance_initial_grid Policy_initial StationaryDist_initial AggVars_initial n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat
end

% Some things you might want to take a look at just to see what is going on.
% p_eqm_initial
% AggVars_initial
% [Params.B, Params.Bprime]
% plot(p_grid, MarketClearance_initial)
% 
% Policy_initial(2,1000:end,11:13)
% plot(shiftdim(Policy_initial(2,:,:),1))
% plot(shiftdim(Policy_initial(2,:,:),1)-(1:1:n_a)'.*ones(n_a,n_z))
% plot(sum(StationaryDist_initial,2))
% plot(cumsum(sum(StationaryDist_initial,2)))

%% Table 1
FID = fopen('./SavedOutput/LatexInputs/GuerrieriLorenzoni2017_Table1.tex', 'w');
fprintf(FID, 'Parameters Values \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llll} \n \\hline \\hline \n');
fprintf(FID, 'Parameter & Explanation & Value & Target/Source \\\\ \n \\hline \n');
fprintf(FID, '$\\beta$ &  Discount Factor  & %8.4f  & Interest rate r=2.5\\%%  \\\\ \n', Params.beta);
fprintf(FID, '$\\gamma$ & Coefficient of relative risk aversion  & %i  &   \\\\ \n', Params.gamma);
fprintf(FID, '$\\eta$ & Curvature of utility of leisure  & %8.1f  & Average Frisch elasticity=1  \\\\ \n', Params.eta);
fprintf(FID, '$\\psi$ & Coefficient on leisure in utility & %8.2f  & Average hours worked 0.4 of endowment \\\\ \n', Params.psi);
fprintf(FID, '$\\rho$ & Persistence of productivity shock & %8.3f  & Persistence of wage process \\\\ \n', Params.rho);
fprintf(FID, '$\\sigma_{\\epsilon}$ & Variance of productivity shock & %8.3f  & Variance of wage process \\\\ \n', Params.sigmasq_epsilon);
fprintf(FID, '$\\pi_{e,u}$ & Transition to unemployment & %8.3f  & Shimer (2005)  \\\\ \n', Params.pi_eu);
fprintf(FID, '$\\pi_{u,e}$ & Transition to employment & %8.3f  & Shimer (2005)  \\\\ \n', Params.pi_ue);
fprintf(FID, '$\\upsilon$ & Unemployment benefit & %8.2f  & 40\\%% of average labor income \\\\ \n', Params.v);
fprintf(FID, '$B$ & Bond supply & %8.1f  & Liquid assets (flow of funds) \\\\ \n', Params.B);
fprintf(FID, '$phi$ & Borrowing limit & %8.3f  & Total gross debt (flow of funds) \\\\ \n', Params.phi);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: values of $B$ and $phi$ differ from those in Guerrieri \\& Lorenzoni (2017). The actual model parameters are reported here, while those in paper are the parameter divided by annual output (annual output equals 4 times quarterly output; model is quarterly). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Final stationary equilibrium
Params.phi=Params.phi_final;

if SkipInitialFinal==0
    disp('Calculating prices corresponding to the stationary eqm')
    p_eqm_final=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params.r=p_eqm_final.r;

    % To be able to draw Figure 2, going to do redo, but with a grid on r
    [~,~, MarketClearance_final_grid]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

    disp('Calculating various equilibrium objects')
    [V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
    StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z, simoptions);    
    AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    
    save ./SavedOutput/GuerrieriLorenzoni2017_final.mat Params p_eqm_final MarketClearance_final_grid V_final Policy_final StationaryDist_final AggVars_final n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_final.mat
end

%% Compute Annual GDP

% GL2017 describe Figure 1 as "Figure I shows the optimal values of consumption and labor supply as a function of the initial level of bond
% holdings". This is incorrect. The x-axis is not the level of bond holdings, it is bond-holdings as a fraction of annual output. I follow 
% GL2017 in what they plot, adding the footnote to bottom of figure to explain what is actually graphed.
% In fact this is true of many of the Figures of GL2017, which label the x-axis as either b or B, but are in fact reporting b (or B) divided by annual output.
% To do this I need to compute annual GDP: quarterly output is y=theta*n.
% Following few lines do this (together with multiplication by 4 to make it annual)
FnsToEvaluateExtra.output = @(d, aprime,a,z) d*z; % Output

QuarterlyOutput_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluateExtra,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
QuarterlyOutput_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluateExtra,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

AnnualOutput_initial=4*QuarterlyOutput_initial.output.Mean;
AnnualOutput_final=4*QuarterlyOutput_final.output.Mean;

%% Figure 1
l_a=length(n_a);
l_z=length(n_z);
Fig1FnsToEvaluate.consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

ConsumptionDecision=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy_initial, Fig1FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);
% ConsumptionDecision=EvalFnOnAgentDist_Grid_Case1(Fig1FnsToEvaluate,[Params.r,Params.v, Params.B, Params.Bprime],PolicyValuesPermute,n_d,n_a,n_z,a_grid,z_grid,2);
if CreateFigures==1
    figure(1)
    subplot(2,1,1); plot(a_grid,ConsumptionDecision.consumption(:,2),a_grid,ConsumptionDecision.consumption(:,8))
    % legend('Mean','10th','25th','Median')
    title({'Consumption'})
    xlabel('Bond holdings')
    % Labour supply is the d variable
    subplot(2,1,2); plot(a_grid,d_grid(Policy_initial(1,:,2)),a_grid,d_grid(Policy_initial(1,:,8)))
    legend('\theta^2','\theta^8');
    title({'Labor Supply'})
    xlabel('Bond holdings')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure1.pdf'])
    
    % Now a version that reproduces GL2017 exactly
    subplot(2,1,1); plot(a_grid/AnnualOutput_initial,ConsumptionDecision.consumption(:,2),a_grid/AnnualOutput_initial,ConsumptionDecision.consumption(:,8))
    % legend('Mean','10th','25th','Median')
    title({'Consumption'})
    xlabel('Bond holdings as fraction of annual output')
    ylabel('(Quarterly)')
    % Labour supply is the d variable
    subplot(2,1,2); plot(a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,2)),a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,8)))
    legend('\theta^2','\theta^8');
    title({'Labor Supply'})
    xlabel('Bond holdings as fraction of annual output')
    ylabel('(Quarterly)')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure1_GL2017.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig1.mat ConsumptionDecision Policy_initial a_grid d_grid AnnualOutput_initial
end

%% Figure 2
% This is a bit awkward to plot, but notice that MarketClearance is just
% measuring how large a change in B would be required to make that interest rate an equilibrium. 
% Need to change to annual interest rates (model period is quarterly)
annualinterestrate=100*((1+heteroagentoptions.pgrid).^4-1); % Note: GL2017 multiply by 4 rather than doing the explicit compounding.
if CreateFigures==1
    figure(2)
    plot(Params.B+MarketClearance_initial_grid, annualinterestrate, Params.B+MarketClearance_final_grid, annualinterestrate)
    legend('initial stationary eqm','final stationary eqm');
    % title({'Bond Market Equilibrum in Stationary Eqm'})
    xlabel('aggregate bond holdings')
    ylabel('(annual) interest rate')
    xlim([0.5,3])
    ylim([0,4])
    xline(Params.B,'-','Bond supply') % Requires at least Matlab R2018B
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure2.pdf'])
    
    % Now a version that reproduces GL2017 exactly
    plot((Params.B+MarketClearance_initial_grid)/AnnualOutput_initial, annualinterestrate, (Params.B+MarketClearance_final_grid)/AnnualOutput_final, annualinterestrate)
    legend('initial stationary eqm','final stationary eqm');
    % title({'Bond Market Equilibrum in Stationary Eqm'})
    xlabel('aggregate bond holdings')
    ylabel('(annual) interest rate')
    xlim([0.5,3])
    ylim([0,4])
    xline(Params.B/AnnualOutput_initial,'-','Bond supply') % Requires at least Matlab R2018B
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure2_GL2017.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig2.mat Params MarketClearance_initial annualinterestrate MarketClearance_final AnnualOutput_initial AnnualOutput_final
end
% I deliberately change the description to 'stationary equilibrium' instead of 'steady state'. I consider the later a slightly misleading
% description, albeit common in the literature, as the model contains plenty of shocks that are active (at agent level) and steady-state is 
% used in representative agent models to designate the equilibrium in the model with no  shocks (subtly, this is different to stationary 
% equilibrium, the model in which shocks exists, but there simply haven't been any for an infinite period of time). It is only the
% cross-section that is stationary, agents are moving around. But based on the existence of shocks in the model I consider stationary a 
% more accurate description/definition.


%% Figure 4
if CreateFigures==1
    figure(4)
    % Unclear from GL2017 paper what is plotted in the top panel.
    % From their codes it is clear that it is 'mean based on the unconditional
    % distribution of the exogenous shock' (note that this is not the same
    % as the 'mean based on the distribution of agents at that asset level').
    subplot(2,1,1); plot(a_grid,sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid,a_grid,sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)
    subplot(2,1,2); plot(a_grid,sum(StationaryDist_initial,2), a_grid,sum(StationaryDist_final,2))
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure4.pdf'])
    
    % Now a version that reproduces GL2017 exactly
    % Note: they are reporting x-axes that are 'as a fraction of initial annual
    % output' for the final equilibrium results, not the final annual output.
    % (Although in practice initial and final annual output are almost equal)
    % For the top panel in fact everything (both axes) is as a fraction of
    % initial annual output.
    subplot(2,1,1); plot(a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial,a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial)
    subplot(2,1,2); plot(a_grid/AnnualOutput_initial,sum(StationaryDist_initial,2), a_grid/AnnualOutput_initial,sum(StationaryDist_final,2))
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure4_GL2017.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig4.mat a_grid pistar_z Policy_initial Policy_final AnnualOutput_initial StationaryDist_initial StationaryDist_final
end

%%
% Free up space on gpu
clear Policy_initial
clear StationaryDist_final

% load ./SavedOutput/GuerrieriLorenzoni2017_final.mat V_final
% load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat StationaryDist_initial

%% Compute the transition path
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=50; % GL2017 use 25

% We want to look at a one off unanticipated path of phi. ParamPath & PathParamNames are thus given by
ParamPath.phi=Params.phi_final*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath.phi(1:6)=temp(2:7); % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
PricePath0.r=[linspace(-0.01, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % PricePath0 is matrix of size T-by-'number of prices'

% Rewrite the aggregate variable to be next period bonds rather than current bonds as this is the actual 
% timing of the decision which the interest rate (r) effects
TransPathFnsToEvaluate.Aprime = @(d, aprime,a,z) aprime; % Aggregate assets decisions
% Rewrite the General Eqm conditions as rules for updating the price
TransPathGeneralEqmEqns.BondMarket = @(Aprime,Bprime) Aprime-Bprime; % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)

transpathoptions.GEnewprice=3; % 3 allows you to say which general eqm eqn is used to update which general eqm price/param, and how
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'BondMarket','r',0,0.01};... % first GE condition will be positive if r is too big (as this will lead to too much Aprime), so subtract
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

% Now just run the TransitionPath_Case1 command (all of the other inputs are things we had already had to define to be able to solve for the initial and final equilibria)
transpathoptions.verbose=1;
if SkipBaselineTransPath==0
    tic;
    PricePath=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  TransPathFnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, vfoptions,simoptions);
    baselinetpathtime=toc
    save ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat PricePath n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat PricePath
end

[VPath,PolicyPath]=ValueFnOnTransPath_Case1(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPath=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath,n_d,n_a,n_z,pi_z,T);

% For later we will keep another copy
PricePath_Flex=PricePath;

%% Figure 3
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat
load ./SavedOutput/GuerrieriLorenzoni2017_final.mat
load ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat 

Fig3FnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig3FnsToEvaluate.debt = @(d, aprime,a,z) -a*(a<0); % debt is (minus of) negative assets

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig3FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1(Fig3FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);

Output_pch=([AggVars_initial.output.Mean; AggVarsPath.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;

if CreateFigures==1
    figure(3)
    % Borrowing limit
    subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath.phi])
    title('borrowing constraint')
    % household debt-to-GDP ratio
    subplot(2,2,2); plot(0:1:T,[AggVars_initial.debt.Mean./AggVars_initial.output.Mean; AggVarsPath.debt.Mean./AggVarsPath.output.Mean])
    title('household debt-to-annual-GDP ratio')
    % interest rate
    subplot(2,2,3); plot(0:1:T,100*(((1+[p_eqm_initial.r; PricePath.r]).^4)-1)) % converts to annual rate by compounding (not just multiplying by 4)
    title('annual interest rate')
    % output
    subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
    title('output')
    % ylabel('percent deviation from inital output in stationary eqm')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure3.pdf'])
    
    % Now a version that reproduces GL2017 exactly.
    % Borrowing limit
    subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath.phi]./AnnualOutput_initial)
    title('borrowing constraint as fraction-of-initial-annual-output')
    % household debt-to-GDP ratio
    subplot(2,2,2); plot(0:1:T,[AggVars_initial.debt.Mean; AggVarsPath.debt.Mean]./AggVars_initial.output.Mean)
    title('household debt-to-initial-annual-GDP ratio')
    % interest rate
    subplot(2,2,3); plot(0:1:T,100*4*[p_eqm_initial.r; PricePath.r]) % converts to annual rate by compounding (not just multiplying by 4)
    title('annual interest rate')
    % output
    subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
    % ylabel('percent deviation from inital output in stationary eqm')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure3_GL2017.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig4.mat T Params ParamPath AggVars_initial AggVarsPath p_eqm_initial PricePath Output_pch AnnualOutput_initial
end

%% Figure 5
Fig5FnsToEvaluate.consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig5FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
AggVarsPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
% Only debt limit reduction path
UnchangedPricePath.r=p_eqm_initial.r*ones(T,1);
[~,PolicyPath_partial_onlydebtlimit]=ValueFnOnTransPath_Case1(UnchangedPricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPath_partial_onlydebtlimit=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_partial_onlydebtlimit,n_d,n_a,n_z,pi_z,T);
AggVarsPath_partial_onlydebtlimit=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, AgentDistPath_partial_onlydebtlimit, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
% Only interest rate path
UnchangedParamPath.phi=Params.phi_initial*ones(T,1);
[~,PolicyPath_partial_onlyinterestrate]=ValueFnOnTransPath_Case1(PricePath, UnchangedParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPath_partial_onlyinterestrate=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_partial_onlyinterestrate,n_d,n_a,n_z,pi_z,T);
AggVarsPath_partial_onlyinterestrate=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, AgentDistPath_partial_onlyinterestrate, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);

AggVarsPath_GE_pch=100*(AggVarsPath_GE.consumption.Mean-AggVars_initial.consumption.Mean)./AggVars_initial.consumption.Mean;
AggVarsPath_partial_onlydebtlimit_pch=100*(AggVarsPath_partial_onlydebtlimit.consumption.Mean-AggVars_initial.consumption.Mean)./AggVars_initial.consumption.Mean;
AggVarsPath_partial_onlyinterestrate_pch=100*(AggVarsPath_partial_onlyinterestrate.consumption.Mean-AggVars_initial.consumption.Mean)./AggVars_initial.consumption.Mean;

if CreateFigures==1
    figure(5)
    plot(0:1:T,[0;AggVarsPath_GE_pch], 0:1:T,[0;AggVarsPath_partial_onlydebtlimit_pch], 0:1:T,[0;AggVarsPath_partial_onlyinterestrate_pch])
    title('Consumption Response Deviation')
    legend('General eqm response','Partial eqm response to debt limit reduction','Partial eqm response to interest rate changes')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure5.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig5.mat T AggVarsPath_GE_pch AggVarsPath_partial_onlydebtlimit_pch AggVarsPath_partial_onlyinterestrate_pch
end

%% Figure 6
if SkipPanelData==0
    % Unlike the other figures relating to the transition, which just require
    % the distribution of agents at each time period, this figure requires
    % following individuals along the transition path based on where they
    % started. Hence will use SimPanelValues_TransPath_Case1() rather than
    % EvalFnOnTransPath_AggVars_Case1(); actually the later could/should also be used
    % here, but want to show the different options available as part of VFI Toolkit.

    % For each of the four transition paths simulate 100 paths drawing from the relevant initial percentile, then take mean.
    % (Alternative would be a much bigger panel based on drawing from the actual inital stationary distribution, then
    % just restrict panel based on percentiles and take means, but this would involve much more computation time)
    simoptions.numbersims=100;

    % We are going to want the values for consumption
    Fig6FnsToEvaluate.consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

    % First, figure out the asset values that correspond to the percentiles
    assetdist=cumsum(sum(StationaryDist_initial,2));
    [~, prctileindexes]=min(abs(assetdist-(1:1:100)/100)); % prctilevalues_doublecheck should be approx equal to prctilevalues
    % prctilevalues_doublecheck=assetdist(prctileindexes); % should give [0.01,0.02, etc up to 1]

    % 1st percentile (note, the 1st percentile are going to be those closest to the borrowing constraint)
    % Set up the appropriate inital distribution and simulate a panel data set (over the transition) from this.
    InitialDist_1stpercentile=zeros(n_a,n_z,'gpuArray');
    InitialDist_1stpercentile(prctileindexes(1),:)=StationaryDist_initial(prctileindexes(1),:)./sum(StationaryDist_initial(prctileindexes(1),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
    % Everything else is just completely standard
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_1stpercentile, Policy_initial, Fig6FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    SimPanelValues=SimPanelValues_TransPath_Case1(PolicyPath, PricePath, ParamPath, T, InitialDist_1stpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, Fig6FnsToEvaluate, Params, transpathoptions,simoptions);
    % The line in the figure is just the mean for each time period of these (I am guessing), but expressed as
    % percent deviation from steady state. [Not obvious if I should take mean and then percent deviation, or
    % take percent deviation and then mean; have gone with the former.]
    % Fig6_1stPercentileTrace=(mean(SimPanelValues,2)-mean(SimPanelValues(1,:)))/mean(SimPanelValues(1,:));
    Fig6_1stPercentileTrace=(mean(SimPanelValues.consumption,2)-AggVars_initial.consumption.Mean)/AggVars_initial.consumption.Mean;

    % Now just repeat for 10th, 20th and 50th percentiles
    InitialDist_10thpercentile=zeros(n_a,n_z,'gpuArray');
    InitialDist_10thpercentile(prctileindexes(10),:)=StationaryDist_initial(prctileindexes(10),:)./sum(StationaryDist_initial(prctileindexes(10),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_10thpercentile, Policy_initial, Fig6FnsToEvaluate,Params,[],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    SimPanelValues=SimPanelValues_TransPath_Case1(PolicyPath, PricePath, ParamPath, T, InitialDist_10thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, Fig6FnsToEvaluate, Params, transpathoptions,simoptions);
    Fig6_10thPercentileTrace=(mean(SimPanelValues.consumption,2)-AggVars_initial.consumption.Mean)/AggVars_initial.consumption.Mean;
    InitialDist_20thpercentile=zeros(n_a,n_z,'gpuArray');
    InitialDist_20thpercentile(prctileindexes(20),:)=StationaryDist_initial(prctileindexes(20),:)./sum(StationaryDist_initial(prctileindexes(20),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_20thpercentile, Policy_initial, Fig6FnsToEvaluate,Params,[],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    SimPanelValues=SimPanelValues_TransPath_Case1(PolicyPath, PricePath, ParamPath, T, InitialDist_20thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, Fig6FnsToEvaluate, Params, transpathoptions,simoptions);
    Fig6_20thPercentileTrace=(mean(SimPanelValues.consumption,2)-AggVars_initial.consumption.Mean)/AggVars_initial.consumption.Mean;
    InitialDist_50thpercentile=zeros(n_a,n_z,'gpuArray');
    InitialDist_50thpercentile(prctileindexes(50),:)=StationaryDist_initial(prctileindexes(50),:)./sum(StationaryDist_initial(prctileindexes(50),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_50thpercentile, Policy_initial, Fig6FnsToEvaluate,Params,[],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    SimPanelValues=SimPanelValues_TransPath_Case1(PolicyPath, PricePath, ParamPath, T, InitialDist_50thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, Fig6FnsToEvaluate, Params, transpathoptions,simoptions);
    Fig6_50thPercentileTrace=(mean(SimPanelValues.consumption,2)-AggVars_initial.consumption.Mean)/AggVars_initial.consumption.Mean;

    if CreateFigures==1
        fig6=figure(6);
        plot(0:1:T, [0, Fig6_1stPercentileTrace'], 0:1:T, [0, Fig6_10thPercentileTrace'], 0:1:T, [0, Fig6_20thPercentileTrace'], 0:1:T, [0, Fig6_50thPercentileTrace'])
        title('Consumption Response by Percentile in Initial Wealth Distribution')
        legend('1st percentile','10th percentile','20th percentile','50th percentile')
        saveas(fig6,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure6.pdf'])
    else
        save ./SavedOutput/GL2017Figs/GL2017_Fig6.mat T Fig6_1stPercentileTrace Fig6_10thPercentileTrace Fig6_20thPercentileTrace Fig6_50thPercentileTrace
    end
end

%% Free up some memory
clear SimPanelValues

%% Figure 7
Fig7FnsToEvaluate.employment = @(d, aprime,a,z) d; %n_it in notation of GL2017

Employment_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig7FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
EmploymentPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig7FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
EmploymentPath_GE_pch=100*(EmploymentPath_GE.employment.Mean-Employment_initial.employment.Mean)./Employment_initial.employment.Mean;

if CreateFigures==1
    fig7=figure(7);
    plot(0:1:T,[0; EmploymentPath_GE_pch])
    title('Employment Response')
    saveas(fig7,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure7.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig7.mat T EmploymentPath_GE_pch
end

%% Sticky Wages
% What GL2017 refer to as sticky wages actually enters the model as an increase in
% the utility of leisure. (They justify this as decreasing hours worked
% without decreasing the actual real-age-per-unit-of-time-worked which is
% intended to capture a decrease in labour demand by (unmodelled) firms.)

% Only two changes: first is that the return function needs to be modified for this.
% Second is general equilibrium conditions, these now need to enforce r>=0
% and determine omega>0 as that which is needed to acheive that r=0 as a general
% eqm outcome, whenever omega=0 would otherwise mean r<0.
% Then compute transition path (omega=0 in initial and final states, so these remain unchanged).
% The key is all in writing the general equilibrium conditions in right way
% so that omega=0, except when this would lead to r<0, and in these cases
% need to pick omega so that r=0.

% We need to give an initial guess for the price path on interest rates and
% omega. Lets just start with the flexible prices general eqm path for
% interest rates, and with omega=0 for all periods as our initial guess.
PricePath0.r=PricePath.r;
PricePath0.omega=zeros(T,1);

% Note: I will 'subtract' BondMarket from r, but 'add' StickWages to omega
% Basic idea is do the same thing as before, but now when the interest rate
% is negative we set omega to be positive and increase omega until we reach
% that the interest rate reaches the zero lower bound.
NKTransPathGeneralEqmEqns.BondMarket = @(Aprime,Bprime,r,omega) (r>=0)*(omega<=0)*(Aprime-Bprime)-(r<0)*(omega>0)*(20*abs(r)+0.001);
% If r>=0: do the usual Aprime-Bprime
% If r<0 and omega=0: don't do anything (wait until omega>0 in a later iteration)
% If r<0 and omega>0: increase r by 20*abs(r)+0.001 (note, formula includes a minus as we subtract this from r)
NKTransPathGeneralEqmEqns.StickyWages = @(Aprime,Bprime,r,omega) (r<0)*(0.001+3*abs(r))+(r>=0)*(omega>0)*(-omega);
% If r>=0 (and omega>0) then send omega towards zero (in this case it returns -(0.0001+3*abs(omega)))
% If r<0 then increase omega (in this case it returns 0.001+3*abs(r))

% Remark: This approach to how to update omega is different from that used by GL2017. They instead force the ZLB on the interest rate r, and then use
% omega to ensure that goods markets clear (Y=C; remember this is an endowment economy). omega is updated based on how far the model is from
% goods market clearance.
% Remark: The above updating rules took a few attempts to come up with (the r>=0 constraint is a 'knife-edge' and so it was not trivial to find something
% that settles down rather than 'jump' back and forth between r<0 and r>0).

transpathoptions_NK.GEnewprice=3; % 3 allows you to say which general eqm eqn is used to update which general eqm price/param, and how
transpathoptions_NK.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'BondMarket','r',0,0.01;... % first GE condition will be positive if r is too big (as this will lead to too much Aprime) and omega is positive, so subtract
     'StickyWages','omega',1,0.03};
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

% Now just run the TransitionPath_Case1 command (all of the other inputs are things we had already had to define to be able to solve for the initial and final equilibria)
transpathoptions_NK.verbose=1;
transpathoptions_NK.tolerance=10^(-4); % will run until r and omega settle to four digits
if SkipNKTransPath==0
    PricePath_NK=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  TransPathFnsToEvaluate, NKTransPathGeneralEqmEqns, Params, DiscountFactorParamNames,transpathoptions_NK,vfoptions);
    save ./SavedOutput/GuerrieriLorenzoni2017_transpathNK.mat PricePath_NK n_d n_a n_z
else
    load ./SavedOutput/GuerrieriLorenzoni2017_transpathNK.mat    
end

[VPath_NK,PolicyPath_NK]=ValueFnOnTransPath_Case1(PricePath_NK, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPath_NK=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_NK,n_d,n_a,n_z,pi_z,T);

%% Figure 8 (I do an additional Figure 17 that shows the 'wedges' and employment)
Fig8FnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig8FnsToEvaluate.employment = @(d, aprime,a,z) d; %n_it in notation of GL2017

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig8FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
% Difference between following two lines is PricePath vs PricePath_NK
AggVarsPath_Flex=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
AggVarsPath_NK=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, AgentDistPath_NK, PolicyPath_NK, PricePath_NK, ParamPath, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);

OutputPath_pch_Flex=(AggVarsPath_Flex.output.Mean-AggVars_initial.output.Mean)./AggVars_initial.output.Mean;
OutputPath_pch_NK=(AggVarsPath_NK.output.Mean-AggVars_initial.output.Mean)./AggVars_initial.output.Mean;
% My extras
EmploymentPath_pch_Flex=(AggVarsPath_Flex.employment.Mean-AggVars_initial.employment.Mean)./AggVars_initial.employment.Mean; % This has already been calculated.
EmploymentPath_pch_NK=(AggVarsPath_NK.employment.Mean-AggVars_initial.employment.Mean)./AggVars_initial.employment.Mean;

if CreateFigures==1
    figure(8)
    % interest rate
    subplot(2,1,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_Flex.r], 0:1:T,4*100*[p_eqm_initial.r; PricePath_NK.r])
    title('annual interest rate')
    legend('flex price', 'NK fix price (ZLB)')
    % output
    subplot(2,1,2); plot(0:1:T,100*[0; OutputPath_pch_Flex],0:1:T,100*[0; OutputPath_pch_NK])
    title('output')
    % ylabel('percent deviation from inital output in stationary eqm')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure8.pdf'])
    
    % Not actually part of GL2017, but let's take a look at the path for omega (the wedge on labour supply that 'enforces' the zero lower bound on interest rates)
    figure(17)
    % wedges
    subplot(2,1,1); plot(0:1:T,[0; PricePath_NK.omega])
    title('omega (NK wedge on wages)')
    % employment
    subplot(2,1,2); plot(0:1:T,100*[0; EmploymentPath_pch_Flex],0:1:T,100*[0; EmploymentPath_pch_NK])
    title('employment')
    legend('flex price', 'NK fix price (ZLB)')
    % ylabel('percent deviation from inital output in stationary eqm')
    saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_ExtraFigure1.pdf'])
else
    save ./SavedOutput/GL2017Figs/GL2017_Fig8.mat T p_eqm_initial PricePath_Flex PricePath_NK OutputPath_pch_Flex OutputPath_pch_NK EmploymentPath_pch_Flex EmploymentPath_pch_NK
end

%%
ParamPath_baseline=ParamPath;


%% Do the alternative calibrations. There are four.
% The relevant parameter values are all given in the online appendix of GL2017 (pg 4 and 5).
% Note that some of the parameter values being set below are actually just the same as their baseline values.

% All of them graph subset of the interest rate, output, employment
% For each, need to calculate initial and final general eqm, then transition path for flex and NK
% The script 'GuerrieriLorezoni2017_altcalibration' is mostly just a copy
% paste of the above that generates the relevant outputs for whatever
% calibration is currently implicit in Ṕarams.

AnnualOutput_initial=gather(AnnualOutput_initial); % Note: model is quarterly, so this is four times the model period output

% Need to plot the baseline for some
Full_OutputPath_pch_Flex_baseline=100*[0; OutputPath_pch_Flex];
Full_PricePath_Flex_baseline=4*100*[p_eqm_initial.r; PricePath_Flex.r];

% Is not quite clear from GL2017 if the following are exactly what they do,
% in the sense that I just use the existing AnnualOutput_initial to figure
% out the parameter values for B and phi from the ones they report in the
% online appendix. They probably use the AnnualOutput_initial that would be
% relevant in each of these alternative calibrations. This could be backed
% out from use of their codes; if anyone wants to do so and send me the
% 'corrected' values I would be grateful, but I cannot be bothered right
% now.

% Figures 9 and 11 will 'fail' to replicate. This is 'intentional'. They
% calculate transition paths under a zero lower bound on interest rates in
% which the final equilibrium has a negative interest rate.
% Mathematically this can be accomadated using a wedge, as GL2017 do for
% the all the zero lower bound transitions, in the final equilibrium
% definition as well, but the economic justification underpinning this
% would become that fixed wages can never fall, even after an infinite
% amount of time passes. If you want to mathematically replicate these just
% add the wedge 'omega' and the 'second' general eqm condition (from
% transition paths) to the codes that determine the final eqm. I choose
% instead to 'fail' to replicate to highlight that these alternative calibrations
% actually involve modifications to the definition of the final eqm (or
% equivalently to the definition of an NK (fixed wage) transition path.

% Restore the original baseline calibration.
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

if SkipAlternativeCalibrations==0
    %% First alternative: Median Wealth Calibration
    fprintf('Doing first alternative calibration:  Median Wealth Calibration \n')
    Params.beta=0.9615;
    Params.eta=1.5;
    Params.psi=19.03;
    Params.v=0.1;
    Params.B=0.5241*AnnualOutput_initial;
    Params.Bprime=Params.B;
    Params.phi_initial=0.6031*AnnualOutput_initial;
    
    altcalib_figurenumber=9;
    
    GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.
    
    % Restore the original baseline calibration.
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
    
    %% Second alternative: Low phi Calibration
    fprintf('Doing second alternative calibration: Low phi Calibration \n')
    Params.beta=0.9755;
    Params.eta=0.25;
    Params.psi=1.389;
    Params.v=0.1;
    Params.B=1.6*AnnualOutput_initial;
    Params.Bprime=Params.B;
    Params.phi_initial=0.7077*AnnualOutput_initial;
    
    altcalib_figurenumber=10;
    
    GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.
    
    % Restore the original baseline calibration.
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
    
    %% Third alternative: Median wealth, low phi Calibration
    fprintf('Doing third alternative calibration: Median wealth, Low phi Calibration \n')
    Params.beta=0.9544;
    Params.eta=0.25;
    Params.psi=1.427;
    Params.v=0.1;
    Params.B=0.429*AnnualOutput_initial;
    Params.Bprime=Params.B;
    Params.phi_initial=0.7184*AnnualOutput_initial;
    
    altcalib_figurenumber=11;
    
    GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.
    
    % Restore the original baseline calibration.
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
    
    %% Fourth alternative: Calibration with gamma=6
    fprintf('Doing fourth alternative calibration:  Calibration with gamma=6 \n')
    % This one is kind of pointless. It turns up the risk aversion, but because
    % Tauchen q=2 (hyperparameter when doing Tauchen approx of the AR(1) log-income process)
    % the min and max values of the log-income process are only +/-1 standard
    % deviation anyway. So model economy has had all the actual risk removed
    % from it before starting.
    Params.beta=0.9756;
    Params.eta=1.5;
    Params.psi=101.22;
    Params.v=0.1;
    Params.B=1.6*AnnualOutput_initial;
    Params.Bprime=Params.B;
    Params.phi_initial=1.1*AnnualOutput_initial;
    % And presumably also (but not explicitly mentioned in GL2017 appendix)
    Params.gamma=6;
    
    altcalib_figurenumber=12;
    
    GuerrieriLorenzoni2017_altcalibration % Contents of this are largely just a copy paste of the above code, which is now rerun for each alternative calibration. Better practice would be to set it up as a seperate function with inputs and outputs.
    
    % Restore the original baseline calibration.
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
end


%%
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat StationaryDist_initial Policy_initial

%% Fiscal Policy
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

if SkipFiscalTransPath==0
    % Pg 1457 of GL2017: since taxes are lump-sum it is in principle
    % possible to fully offset the tightening of the borrowing limit by
    % transfering the 'same amount' to agents.
    % The experiment that will be performed here is to look at policies that
    % only partially offset the long run change in borrowing limit phi.
    
    % Note: After coding the PathForB below it becomes clear that T=25 will not
    % be enough to reach a new stationary recursive equilibrium (B_t has not yet reached Bprime). So have
    % increased to T=50.
    T=50;
    % Not perfect, but B_t now reaches very close to Bprime and is almost constant.
    % Based on Figures, it appears GL2017 kept T=25.
    
    % Two transfer policies are assessed, the first where government
    % temporarily reduces the lump-sum tax tau for all households, and the
    % second where the government temporarily raises the unemployment benefit.
    
    % Both transfer policies assume a fixed sequence of
    % government deficits implied by the bond supply following the path,
    %   B_t = rho_b^t Binitial + (1-rho_b^t)*Bfinal
    % where Bfinal is the new long-run level of B.
    % Set Bfinal=1.2*B, and rho_b=0.95.
    Params.rho_b=0.95;
    Params.Bfinal=1.2*Params.B;
    ParamPath_Fiscal1.B=(Params.rho_b.^(0:1:T-1)*Params.B+(1-Params.rho_b.^(0:1:T-1))*Params.Bfinal)';
    % Because of how the current lump-sum tax rate 'tau' is calculated we need
    % B_t over the whole path. We will actually also need B_tplus1 as well as
    % this is key to the general equilibrium condition. (See above discussion
    % when the general equilibrium condition was first described and introduced
    % when computing the initial general equilibrium.)
    ParamPath_Fiscal1.Bprime=[ParamPath_Fiscal1.B(2:end); Params.rho_b^(T)*Params.B+(1-Params.rho_b^(T))*Params.Bfinal];
    % Force that Bprime=B in the final period.
    ParamPath_Fiscal1.Bprime(end)=ParamPath_Fiscal1.B(end);
    
    % Before we get into the transition paths themselves we need to calculate the final general equilibrium.
    Params.phi=Params.phi_final;
    % Params.Binitial=Params.B; % So can restore after calculation. % Just load the inital params again instead; I do this roughly 10-15 lines below
    Params.B=ParamPath_Fiscal1.Bprime(end); % Use ParamPath_Fiscal1.Bprime(end) rather than Params.Bfinal
    Params.Bprime=Params.B;
    FnsToEvaluate.A = @(d, aprime,a,z) a; % Aggregate assets (which is this periods state)
    GeneralEqmEqn.BondMarket = @(A,p,B) A-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
    p_eqm_final_Fiscal=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params.r=p_eqm_final_Fiscal.r;
    [V_final_Fiscal,~]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, []);
    % Params.phi=Params.phi_initial;
    % Params.B=Params.Binitial;
    % Params.Bprime=Params.B;
    % Done, now back to the transition path.
    load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params
    
    % The two transfer policies are then different ways of spending the associated deficit implied by the bond supply path.
    
    % First, we consider a policy where the tax tau_t adjusts to balance government budget to
    % generate the relevant deficit (period-by-period). tau_t is determined inside the return function,
    % and is already determined in this way.
    
    % The second, transition considers a policy in which unemployment benefits are increased by 50%
    % for the first two years after the shock. Afterward the unemployment
    % benefit reverts to its initial value and throughout the transition path
    % tau_t adjusts to satisfy the government budget constraint.
    % This will simply involve adding the unemployment benefit path for v to
    % ParamPath, and will be done below after we finish the first case.
    
    % Set the new path for parameters to be modelled. In addition to the same
    % path on phi as before, we now have a path for B.
    ParamPath_Fiscal1.phi=Params.phi_final*ones(T,1);
    temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath_Fiscal1.phi(1:6)=temp(2:7)'; % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
    
    % The rest of the setup is just same thing as the baseline case with flexible prices. So most of the following lines are just copy-paste from there.
    
    % We need to give an initial guess for the price path on interest rates
    PricePath0_Fiscal.r=[linspace(p_eqm_initial.r, p_eqm_final_Fiscal.r, floor(T/2))'; p_eqm_final_Fiscal.r*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
    % Note: tau is also being determined, but this is done inside the ReturnFn
    
    % Rewrite the aggregate variable to be next period bonds rather than current bonds as this is the
    % actual timing of the decision which the interest rate (r) effects
    FnsToEvaluate.Aprime = @(d, aprime,a,z) aprime; % Aggregate assets decisions
    % Rewrite the General Eqm conditions as rules for updating the price
    FiscalTransPathGeneralEqmEqns.BondMarket = @(Aprime,Bprime) Aprime-Bprime; % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
    
    transpathoptions_Fiscal.GEnewprice=3; % 3 allows you to say which general eqm eqn is used to update which general eqm price/param, and how
    transpathoptions_Fiscal.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'BondMarket','r',0,0.01}; % first GE condition will be positive if r is too big (as this will lead to too much Aprime), so subtract
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

    % Now just run the TransitionPath_Case1 command (all of the other inputs are things we had already had to define to be able to solve for the initial and final equilibria)
    transpathoptions_Fiscal.weightscheme=1
    transpathoptions_Fiscal.verbose=1;
    PricePath_Fiscal1=TransitionPath_Case1(PricePath0_Fiscal, ParamPath_Fiscal1, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, FiscalTransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_Fiscal,vfoptions);
    
    save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fiscalpolicy.mat PricePath_Fiscal1

    [~,PolicyPath_Fiscal1]=ValueFnOnTransPath_Case1(PricePath_Fiscal1, ParamPath_Fiscal1, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPath_Fiscal1=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_Fiscal1,n_d,n_a,n_z,pi_z,T);

    % Now for the second transition.
    % The second, transition considers a policy in which unemployment benefits are increased by 50% for the first two years after the shock.
    % Afterward the unemployment benefit reverts to its initial value and throughout the transition path tau_t adjusts to satisfy the
    % government budget constraint. This will simply involve adding the unemployment benefit path for v to ParamPath, and will be done below after we finish the first case.
    ParamPath_Fiscal2=ParamPath_Fiscal1; % Use same path for phi,B,Bprime
    ParamPath_Fiscal2.v=Params.v*[1.5*ones(8,1); ones(T-8,1)];
    % Final equilibrium is same as for the other Fiscal Policy transition path.
    PricePath_Fiscal2=TransitionPath_Case1(PricePath0_Fiscal, ParamPath_Fiscal2, T, V_final_Fiscal, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, FiscalTransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_Fiscal,vfoptions);
    
    save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fiscalpolicy.mat PricePath_Fiscal1 PricePath_Fiscal2

    [~,PolicyPath_Fiscal2]=ValueFnOnTransPath_Case1(PricePath_Fiscal2, ParamPath_Fiscal2, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPath_Fiscal2=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_Fiscal2,n_d,n_a,n_z,pi_z,T);
    
    % Now we have done the Fiscal Policy, just need to create the relevant Figure.
    FiscPolFnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
    
    % Note that Params is currently all the initial values.
    
    % Get the
    AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FiscPolFnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    AggVarsPath_Fiscal1=EvalFnOnTransPath_AggVars_Case1(FiscPolFnsToEvaluate, AgentDistPath_Fiscal1, PolicyPath_Fiscal1, PricePath_Fiscal1, ParamPath_Fiscal1, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
    AggVarsPath_Fiscal2=EvalFnOnTransPath_AggVars_Case1(FiscPolFnsToEvaluate, AgentDistPath_Fiscal2, PolicyPath_Fiscal2, PricePath_Fiscal2, ParamPath_Fiscal2, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
    
    Output_pch_Fiscal1=([AggVars_initial.output.Mean; AggVarsPath_Fiscal1.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
    Output_pch_Fiscal2=([AggVars_initial.output.Mean; AggVarsPath_Fiscal2.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
    
    if CreateFigures==1
        % Figure 13
        figure(13)
        % interest rate
        subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_Flex.r; ones(T-length(PricePath_Flex.r),1)*PricePath_Flex.r(end)],0:1:T,4*100*[p_eqm_initial.r; PricePath_Fiscal1.r],0:1:T,4*100*[p_eqm_initial.r; PricePath_Fiscal2.r])
        title('annual interest rate')
        % output
        subplot(1,2,2); plot(0:1:T,100*[0; OutputPath_pch_Flex; ones(T-length(OutputPath_pch_Flex),1)*OutputPath_pch_Flex(end)], 0:1:T,100*Output_pch_Fiscal1, 0:1:T,100*Output_pch_Fiscal2)
        title('output')
        saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure13.pdf'])
    else
        save ./SavedOutput/GL2017Figs/GL2017_Fig13.mat T p_eqm_initial PricePath_Flex PricePath_Fiscal1 PricePath_Fiscal2 OutputPath_pch_Flex Output_pch_Fiscal1 Output_pch_Fiscal2
    end
end

%% Fischer Debt Deflation
load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

if SkipFisherDebtDeflation==0
    % Page 1459, "A simple experiment in which we assume that all debt (private and public)
    % is in nominal terms and there is an unexpected 10% reduction in the price
    % level at date t=0. After that, the debt level is constant."
    % Not entirely clear but this seems to mean they simply adjust the initial
    % stationary distribution of agents, and change 'B'. Note that it seems
    % like this means that B will remain fixed at this new level as Figure XIV does not show the
    % economy returning to its initial steady state, nor the baseline final steady state but instead
    % to a new higher level, but this is not clear from the text alone.
    
    % Inflation factor
    InflationFactor=0.9; % An unexpected 10% reduction in the price level (=1 means no change; 1.1 would be a 10% inflation, 0.9 would be a 10% deflation)
    
    % Start with the public debt.
    Params.B=(1/InflationFactor)*Params.B;
    Params.Bprime=(1/InflationFactor)*Params.Bprime;
    % Now the private debt, figure out where on a_grid people get shifted to.
    % (Note, the way I do this will only work if no-one is in the very ends of the grid (deflation would push them off the grid, but this is not possible so they would instead be 'distorted' back onto the grid).)
    newa_grid=(1/InflationFactor)*a_grid; % (absolute value of) debts decrease in real terms, savings decrease in real terms. [Opposite if it is a deflation, as in this example.]
    [~,indextoshiftto]=min(abs(a_grid-newa_grid'));
    checkthislookslikeagrid=a_grid(indextoshiftto); % It does
    newStationaryDist_initial=zeros(size(StationaryDist_initial),'gpuArray');
    for ii=1:length(a_grid)
        newStationaryDist_initial(indextoshiftto(ii),:)=newStationaryDist_initial(indextoshiftto(ii),:)+StationaryDist_initial(ii,:);
    end % Note: could probably have done this better by shifting fraction to each of the nearest two points based on the relative distance to each.
    % The rest is just to recalculate the final distribution (based on new B & Bprime)
    Params.phi=Params.phi_final;
    FnsToEvaluate.A = @(d, aprime,a,z) a; % Aggregate assets (which is this periods state)
    GEPriceParamNames={'r'};
    GeneralEqmEqns.AssetMarket = @(A,B) A-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
    p_eqm_final_Fisher=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params.r=p_eqm_final_Fisher.r;
    [V_final,~]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    Params.phi=Params.phi_initial;
    % Then the transition path. The param path is trivial, as it is just the
    % new parameters, together with the usual path on borrowing limit (phi).
    ParamPath_Fisher.phi=Params.phi_final*ones(T,1);
    temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath_Fisher.phi(1:6,1)=temp(2:7)'; % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
    ParamPath_Fisher.B=Params.B*ones(T,1);
    ParamPath_Fisher.Bprime=Params.Bprime*ones(T,1);
    
    PricePath0_Fisher.r=p_eqm_final_Fisher.r*ones(T,1);
    % Rewrite the aggregate variable to be next period bonds rather than
    % current bonds as this is the actual timing of the decision which the
    % interest rate (r) effects
    FnsToEvaluate.Aprime = @(d,aprime,a,z) aprime; % Aggregate assets decisions
    
    FisherTransPathGeneralEqmEqns.BondMarket = @(Aprime,r,Bprime) Aprime-Bprime; 
       % Note: if r is too big, then Aprime is too big
    
    transpathoptions_Fisher.GEnewprice=3; % 3 allows you to say which general eqm eqn is used to update which general eqm price/param, and how
    transpathoptions_Fisher.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'BondMarket','r',0,0.01}; % first GE condition will be positive if r is too big (as this will lead to too much Aprime), so subtract
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

    % Now just run the TransitionPath_Case1 command
    transpathoptions_Fisher.verbose=1;
    PricePath_Fisher=TransitionPath_Case1(PricePath0_Fisher, ParamPath_Fisher, T, V_final, newStationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, FisherTransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_Fisher,vfoptions);
    
    save ./SavedOutput/GuerrieriLorenzoni2017_transpath_fischerdeflation.mat PricePath_Fisher V_final p_eqm_final_Fisher
    
    [~,PolicyPath_Fisher]=ValueFnOnTransPath_Case1(PricePath_Fisher, ParamPath_Fisher, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPath_Fisher=AgentDistOnTransPath_Case1(StationaryDist_initial,PolicyPath_Fisher,n_d,n_a,n_z,pi_z,T);

    % Now we have done the Fischer Deflation, just need to create the relevant Figure.
    FisherDefFnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
    
    % Already have the AggVars_initial from previously, so just use those.
    AggVarsPath_Fisher=EvalFnOnTransPath_AggVars_Case1(FisherDefFnsToEvaluate, AgentDistPath_Fisher, PolicyPath_Fisher, PricePath_Fisher, ParamPath_Fisher, Params, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, simoptions);
    
    Output_pch_Fisher=([AggVars_initial.output.Mean; AggVarsPath_Fisher.output.Mean]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;
    
    if CreateFigures==1
        % Figure 14
        figure(14)
        % interest rate
        subplot(1,2,1); plot(0:1:T,4*100*[p_eqm_initial.r; PricePath_Flex,r;ones(T-length(PricePath_Flex.r),1)*PricePath_Flex.r(end)], 0:1:T,4*100*[p_eqm_initial.r; PricePath_Fisher.r])
        legend('Baseline','Fisher Deflation')
        title('annual interest rate')
        % output
        subplot(1,2,2); plot(0:1:T,100*[0; OutputPath_pch_Flex; ones(T-length(OutputPath_pch_Flex),1)*OutputPath_pch_Flex(end)], 0:1:T,100*Output_pch_Fisher)
        title('output')
        saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure14.pdf'])
    else
        save ./SavedOutput/GL2017Figs/GL2017_Fig14.mat T p_eqm_initial PricePath_Flex PricePath_Fisher OutputPath_pch_Flex Output_pch_Fisher
    end
end

%% Durable Goods
% Extend the model to introduce durable goods. These will be modelled as a
% second asset. They are less liquid (captured in model by an adjustment
% cost) and can be used as collateral. Consumers will want durable goods as 
% the durable goods are simply added directly into the utility function.
% The extended model all includes a spread between borrowing and lending
% interest rates.
% Notice that because there is no supply side to the durable goods (the
% supply of durable goods is perfectly elastic) the addition of durable
% goods does not involve adding any complication whatsoever to the general
% equilibrium aspects of the model.

load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params

% To be able to solve this we will need
vfoptions.lowmemory=1;
transpathoptions.lowmemory=1;
n_d=21;
% n_a=2^9;
% n_k=2^9; % Number of grid points for the durable asset.
n_a=2^5;
n_k=2^5; % Number of grid points for the durable asset.

% Need to update the grids for d and a based on new number of grid points
d_grid=linspace(0,1,n_d)'; % Labor supply
a_grid=(Params.aupperbar-Params.alowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a)')-1)+Params.alowerbar;


if SkipDurables==0
    % Based on graphs it appears GL2017 use T=50.
    T=50;
    
    % Because this involves solving a whole other model with different
    % dimensions, return function, etc., I have done it in a seperate function
    % that just returns the few things that are needed for the graphs.
    % As inputs it takes the parameters, and the grid sizes for those parts of
    % the model that are not changed.
    
    % The following lines are the options relating to the durable goods model
    % that can be set, and which are inputs to the function.
    Params.kupperbar=10; % (max value for grid, the min is hard-coded as zero)
    
    Params.tauchenq=1; % Reduces Tauchen approximation of theta to 5-states. This is done by GL2017 to ease computation (footnote 40 of online appendix).
    
    %Parameters that are changed or newly introduced for the durable goods model (G&L2017, Table 1)
    Params.beta=0.9713; %Model period is one-quarter of a year.
    Params.psi=2.54; % Coefficient on leisure in utility
    Params.alpha=0.7; % Coefficient on non-durables (Newly introduced)
    Params.delta=0.0129; % Durables depreciation rate (Newly introduced)
    Params.zeta=0.15; % Proportional loss on durable sales (Newly introduced)
    Params.chi=0.01; % Intermediation cost (Newly introduced)
    Params.v=0.16; % Unemployment benefit
    Params.phi_k=0.8; % Borrowing limit (Replaces phi: is now a collateral constraint on borrowing against durable goods)
    % NOTE: These come from Table A.4 of online appendix. It contains footnote
    % "The quantities v and B are expressed in terms of yearly aggregate
    % output." No such footnote appears on Table 1 which gives their
    % (identical) values in baseline calibration. It is unclear what to make of
    % this (codes appear to treat v & B values as those for model period,
    % namely quarterly).
    
    % Do two experiments.
    % The first it to change phi_k to 0.56 gradually over a linear path that lasts six quarters.
    ParamPath_Durables1.phi_k=0.56*ones(T,1);
    ParamPath_Durables1.phi_k(1:6)=linspace(Params.phi_k,0.56,6)';
    % The second is a transitory increase in chi of 0.06 (to 0.07), this is done
    % instantly and then the transitory increase decays back to zero at a rate of decay of 0.6.
    ParamPath_Durables2.chi=(Params.chi+0.06*(0.6.^(0:1:T-1)))';
        
    % We need to give an initial guess for the price path on interest rates.
    % This is done inside the durable goods function once initial and final
    % stationary general eqm interest rates are known. (Both paths have same
    % initial general eqm, but I am going to be lazy and just calculate it
    % twice, this allows the durable goods fn to be more generalized.)
    
    % The general equilibrium conditions are unchanged. These are declared inside the durable good function.
    figurenumber=15;
    Output_DurableGoods1=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a,n_k,d_grid,a_grid, T, ParamPath_Durables1, heteroagentoptions,transpathoptions, vfoptions,simoptions);
    figurenumber=16;
    Output_DurableGoods2=GuerrieriLorenzoni2017_DurableGoods(figurenumber,Params,n_d,n_a,n_k,d_grid,a_grid, T, ParamPath_Durables2, heteroagentoptions,transpathoptions, vfoptions,simoptions);
    
    
    if CreateFigures==1
        for figurenumber=15:16
            % This was created and saved by GuerrieriLorenzoni2017_DurableGoods()
            load (['./SavedOutput/GL2017_fig',num2str(figurenumber),'.mat'],'initialvalforParamPathParams','ParamPath','p_eqm_initial','PricePath','Output_pch','durablespurchases_pch','nondurablespurchases_pch', 'T','AggVars_initial','AggVarsPath')
            
            figure(figurenumber)
            % First parameter in ParamPath
            subplot(2,2,1); plot(0:1:T,[initialvalforParamPathParams(1); ParamPath.phi_k])
            % interest rate
            subplot(2,2,2); plot(0:1:T,[p_eqm_initial.r; PricePath.r])
            % output
            subplot(2,2,3); plot(0:1:T,4*100*Output_pch)
            % durable and nondurable purchases
            subplot(2,2,4); plot(0:1:T,durablespurchases_pch,0:1:T,nondurablespurchases_pch)
            saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni2017_Figure',num2str(figurenumber),'.pdf'])
        end
    end
end



