% Replication of Kitao (2008) - Entrepreneurship, taxation and capital investment
%
% assets: assets
% e: entrepreneur/worker
% eta: labor productivity
% theta: entrepreneurial ability

% Record the full run (incl. any error) to Kitao2008_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./Kitao2008_Diary.txt','file'), delete('./Kitao2008_Diary.txt'); end
diary ./Kitao2008_Diary.txt

% Subcodes (the no-entrepreneur model scripts, and the ReturnFn/ConsumptionFn/FnsToEvaluate family)
% live in ./Kitao2008subcodes/; MATLAB does not search subfolders of the working directory, so put
% it on the path. Note this does NOT change cwd, so the './SavedOutput/...' paths further down still
% resolve to here.
addpath('./Kitao2008subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
addpath(genpath('./VFIToolkit-matlab/'));

close all % close any figures, make sure they are all cleanly built from scratch
% Make sure the subfolders to save output exist
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end
% Note: HuggettVentura2000 keeps this as Params.headlessFigures, but here it has to be a plain
% variable. The doPart blocks below save() the whole workspace and load() it straight back, and a
% load() overwrites Params in its entirety, so a new field on Params would silently disappear the
% first time a doPart(i)=0 branch fires and the guards further down would then error on a
% non-existent field. A plain variable is not in the old .mat files, so it survives the load.
% (Nothing needs this inside the subcodes: set(0,'DefaultFigureVisible','off') is global.)

% Model has two endogenous states and two exogenous markov states (you can change n_asset, but the rest are hardcoded)
n_asset=1200; % assets (Kitao, 2008 uses 3000 points, does pure discretization) [also using vfoptions.ngridinterp]
n_e=2; % 1 is entrepreneur, 0 is worker
n_eta=5; % Kitao (2008) uses 5
n_theta=4; % Kitao (2004) uses 4

% Tax rates for which a transition path is computed (a subset of the tax rates for which the final
% stationary general eqm is computed; those vectors are tau_k_vec and tau_E1_vec, set in doPart(2)
% and doPart(6)). Set here rather than further down because doPart(2) and doPart(6) need to know
% which of their stationary equilibria will later be used as the endpoint of a transition path.
tau_k_transition_vec=[0,5,10,15,20,25,30,35,40]/100;
tau_E1_transition_vec=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4];

doPart=[1,1,1,1,1,1,1];
% Note: each doPart block save()s the whole workspace and the else branch load()s it straight back,
% so running with any doPart(i)=0 restores the Params, grids and p_eqm values from whenever that
% block was last run. After any change to the model codes you must run with doPart all ones (and
% doNoEntrepreneurs=1), otherwise you will silently be mixing old and new results.
% doPart(1): baseline (includes the No-entrepreneur economies)
doNoEntrepreneurs=1; % If zero, just loads results for the No-Entrepreneur economies. If one, compute and save them.
% doPart(2): tax capital income
% doPart(3): transition paths (x2)
% doPart(4): nine transition paths for Figure 7
% doPart(5): just draws a Table
% doPart(6): tax entrepreneurial business income
% doPart(7): transition paths for tax entrepreneurial business income

% Comment: Kitao (2008) uses tau_I to refer to the proportional tax on
% income. It remains a proportional tax on income in the two alternative
% tax regimes, but the definition of income (or more accurately taxable
% income) changes, so the interpretation of tau_I is subtlely different
% accross the different tax systems. Just something to be aware of.

% Table 5: I am not certain on how to compute the 'average leverage ratio (%)' in the last column. So have just gone with a best guess.

% What exactly is the earnings process for economy 'No Entrepreneurs 2'?
% Currently I just recalibrate the top point of eta grid targeting the
% wealth gini and use the rest of eta grid and transitions based on
% Castaneda, Diaz-Gimenez, & Rios-Rull (2003) but taking only the working age
% and dropping retirement. Seems likely to be correct, but not certain.

% For the final stationary equilibrium, Kitao (2008) considers tax rates
% from 0 to 40%, in intervals of 5%, so 9 of them. She does the same for transition
% paths. Here I do the final stationary equilibrium for 0 to 40% in
% intervals of 1%, so 41 of them. But I only do the same 9 transition
% paths.

%% Parameters

Params.sigma=2; % CES utility
Params.beta=0.9428; % time preference

% Corporate sector
Params.Z=1; % Technology level (normalized to 1) (this is never actually used, Kitao (2008) calls it A but I want to use A for aggregate assets)
Params.alpha=0.36; % Cobb-Douglas prodn fn, capital share
Params.delta=0.06; % depreciation rate

% Non-corporate sector
Params.upsilon=0.88;
Params.upsilon1=Params.alpha*Params.upsilon; % pg 50 of Kitao (2008)
Params.upsilon2=(1-Params.alpha)*Params.upsilon; % pg 50 of Kitao (2008)
% the way upsilon1 and upsilon2 two are done is so that capital has same
% 'importance' in non-corporate sector as it does in corporate sector

% Borrowing
Params.phi=0.05; % additional cost of borrowing for non-corporate sector
Params.d=0.5; % max borrowing leverage

% Government
Params.tau_a0=0.258; % Three parameters for non-linear income tax schedule
Params.tau_a1=0.768;
Params.tau_a2=0.438;
Params.tau_I=0.0316; % Proportional tax on income
Params.tau_c=0.0567; % consumption tax
Params.taxincome=1; % 1 is tax income, 2 is to tax capital income and labor income seperately (tau_k), 3 is tax on entrepreneurial business income (tau_E1)
% Following two are needed for the two alternative tax systems
Params.tau_k=0;
Params.tau_E1=0;
% Note: tau_k is only used when taxincome=2, tau_E1 is only used when taxincome=3
% Note: The interpretation of tau_I changes with taxincome. It is always a
% flat tax on taxable income, but the definition of taxable income changes.

%% Set up the exogenous shock processes

% mu: labor productivity (from Kitao (2008), Appendix B)
eta_grid=[0.646; 0.798; 0.966; 1.169; 1.444];
pi_eta=[0.731, 0.253, 0.016, 0.000, 0.000; 0.192, 0.555, 0.236, 0.017, 0.000; 0.011, 0.222, 0.533, 0.222, 0.011; 0.000, 0.017, 0.236, 0.555, 0.192; 0.000, 0.000, 0.016, 0.253, 0.731];
% Note: third row of pi_eta actually sums to 0.999, so need to normalize it to 1
pi_eta=pi_eta./sum(pi_eta,2);

% Kitao (2008) provides grid and transition matrix in Appendix B (in body of article
% it says use 5-state Tauchen-Hussey, which following commented out lines would implement)
% % Params.rho_eta=0.94;
% % Params.sigmasq_eta_epsilon=0.02;
% % tauchenhusseyoptions.baseSigma = sqrt(Params.sigmasq_eta_epsilon); % This is what original Tauchen-Hussey method uses
% % [eta_grid,pi_eta] = discretizeAR1_TauchenHussey(0,Params.rho_eta,sqrt(Params.sigmasq_eta_epsilon),n_eta,tauchenhusseyoptions);
% % eta_grid=exp(eta_grid);
% % % Normalize to unity
% % [E_eta,~,~,statdist_eta]=MarkovChainMoments(eta_grid,pi_eta);
% % eta_grid=eta_grid/E_eta; % Normalize so unconditional mean of eta is unity (Kitao, pg 50)

% theta: entrepreneurial ability (from Kitao (2008), Appendix B)
theta_grid=[0.000; 0.706; 1.470; 2.234];
pi_theta=[0.780, 0.220, 0.000, 0.000; 0.430, 0.420, 0.150, 0.000; 0.000, 0.430, 0.420, 0.150; 0.000, 0.000, 0.220, 0.780];

%% Grids

% In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); % The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
assetmaxfactor=40; % This times K_ss is the max assets
% Note: was 20, which put the top of the grid at 110.3 model units ($4.7m). That is below the
% assets at which a theta4 entrepreneur stops being collateral constrained, k*/(1+d)=149.9 model
% units ($6.4m) at the benchmark prices, so every theta4 entrepreneur on the grid sat at k=(1+d)a
% and the average leverage ratio in Table 5 came out at exactly 50.0% against the 32.6% in Kitao
% (2008). Her Figure 3 has theta4 leaving the max leverage line at about $6.3m, consistent with
% that. 40 puts the top of the grid at 220.7 model units ($9.4m), which clears the kink with room
% to spare as prices move across the tax experiments.
asset_grid=assetmaxfactor*K_ss*(linspace(0,1,n_asset).^3)'; % linspace ^3 puts more points near zero, where the curvature of value and policy functions is higher and where model spends more time
% Note: K_ss is about 5.5, so max assets is about 100

% Economy 'no entrepreneurs 2' needs its own, much larger, asset grid. The Castaneda,
% Diaz-Gimenez & Rios-Rull (2003) process used there has a top labor productivity state of around 300
% times the mean, so a household in that state earns about w*302=330 per year, which is three times
% the entire asset grid above. Their assets then have nowhere to go, which truncates precisely the
% part of the wealth distribution that this economy exists in order to generate. (With the grid above,
% the top 1% wealth share came out at 6.6% against the 35.4% that Kitao (2008) reports, and the Gini
% calibration could not reach its target, because it was fighting the grid rather than the economics.)
% The grid is in two pieces:
%  - below 20*K_ss: exactly the grid above, so that the three economies are directly comparable over
%    the range where essentially the whole population lives (61% of the population of that economy is
%    in the lowest eta state, near the borrowing constraint)
%  - above 20*K_ss: geometric spacing, which is cheap and gives constant proportional resolution
%    through the top tail
assetmax_noE2=20000; % in model units; times wealthscalingfactor/1000 this is about $850 million
% Note: that sounds enormous, but the model is internally consistent about it. The median worker in
% that economy earns about $13,000 per year and the top state earns about $14 million per year, a
% ratio of 1061 to 1, which is Kitao's "more than 1000 times the median".
n_asset_high=350;
a_grid_high=(assetmaxfactor*K_ss)*exp(linspace(0,log(assetmax_noE2/(assetmaxfactor*K_ss)),n_asset_high+1))';
asset_grid_noE2=[asset_grid; a_grid_high(2:end)]; % drop a_grid_high(1), it duplicates the top point of asset_grid
n_asset_noE2=length(asset_grid_noE2);
% Note: the tail here is fat, each extra five years worth of grid room contributes about 85% of what
% the previous five did, so the results are sensitive to assetmax_noE2. Re-run with assetmax_noE2
% doubled and check that the Gini and top wealth shares do not move. Note also that e4 is calibrated,
% and if it rises then the grid requirement rises with it, so check Fig14 after it has converged.

e_grid=[0;1]; % 1 is entrepreneur, 0 is worker

%% Get into form for VFI toolkit
n_d=0;
n_a=[n_asset,n_e];
n_z=[n_eta, n_theta];
d_grid=[];
a_grid=[asset_grid; e_grid];
z_grid=[eta_grid; theta_grid];
pi_z=kron(pi_theta, pi_eta); % in reverse order

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,eprime,a,e,eta,theta,sigma,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)...
    Kitao2008_ReturnFn(aprime,eprime,a,e,eta,theta,sigma,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1);
% The first inputs must be: decision variables, next period endogenous state, endogenous state, exogenous state. Followed by any parameters


%% Aggregates

% Create functions to be evaluated
FnsToEvaluate.K_noncorp = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_kFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % Assets used in non-corporate sector (=entrepreneurs)
FnsToEvaluate.A = @(aprime,eprime,a,e,eta,theta) a; % Total assets of households (workers and entrepreneurs)
FnsToEvaluate.N_noncorp = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_nFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % Labor used in non-corporate sector (=entrepreneurs)
FnsToEvaluate.L = @(aprime,eprime,a,e,eta,theta) eta; % Total labor supply
FnsToEvaluate.TaxRevenue = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)...
    Kitao2008_TaxFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1); % Tax Revenue


%% General equilbrium
GEPriceParamNames={'r','w','G'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but not in transition)

GeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,L,alpha,delta) r-(alpha*((A-K_noncorp)^(alpha-1))*((L-N_noncorp)^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,L,alpha) w-(1-alpha)*((A-K_noncorp)^(alpha))*((L-N_noncorp)^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; %Government runs balanced budget
% % When first solving this using just these GE conditions it was giving me a solution in which L<N_noncorp.
% % So I added the following which is essentially a penalty term on L<N_noncorp.
% GeneralEqmEqns.LaborPenalty = @(L,N_noncorp) 10*abs(L-N_noncorp)*(L<=N_noncorp+0.03); % Add a penalty whenever L<N_noncorp+0.03 (seems like any solution should have at least 0.03 labor in the corporate sector)

%% Some model output it reported in 1000s of dollars, we need a constant scaling factor to multiply model output by to get into dollars
wealthscalingfactor=42618; % Converts model units to dollars
% I cannot find the number by which to multiply model results to get dollars in the paper. 
% I emailed Sagiri Kitao to ask if she still has the constant anywhere. She replied that 
% digging through her old codes "It appears that I rescaled the model unit using the average 
% labor income in 2002, $42,618."

%%
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=25;

simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

% Use different copy of vfoptions for transition paths as there we want to use divide-and-conquer
vfoptionstpath=vfoptions;
vfoptionstpath.divideandconquer=1; % for transition path, also use divide-and-conquer

heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.fminalgo=[8,1]; % fast but not really high accuracy, then a higher accuracy
heteroagentoptions.toleranceGEcondns=[1e-4,1e-5]; % high accuracy on final solve
% Note on MATLAB versions: fminalgo=8 is lsqnonlin(), and the toolkit calls it with the ten-argument
% form lsqnonlin(fun,x0,lb,ub,A,b,Aeq,beq,nonlcon,options). That signature only exists in recent
% MATLAB. On MATLAB 2021a arguments six onwards are instead treated as extra parameters to pass
% through to the objective function, which is a one-input anonymous function, so it fails
% immediately with "Too many input arguments". If you ever have to run this on an old MATLAB, set
% heteroagentoptions.fminalgo=1 (fminsearch) and heteroagentoptions.toleranceGEcondns=1e-5 instead
% (the tolerance has to become a scalar alongside it, because with a scalar fminalgo the toolkit
% skips its multi-stage branch and passes toleranceGEcondns straight to optimset('TolFun',...)).

%%
if doPart(1)==1
    Params.taxincome=1;
    % Set initial value for general eqm
    Params.r=0.04;
    Params.w=1.5;
    Params.G=0.3;
    % Note: I had originally put r=0.04, w=1, G=0.3 but this actually let L-N_noncorp<0. Switched to these.
    
    %% Solve for the stationary general equilbirium
    [p_eqm_initial,GeneralEqmCondn]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    
    p_eqm_initial % The equilibrium values of the GE prices
    % We will want this eqm later for transition paths

    % For later, when doing the transtion paths, we will be using tau_I instead of G as the general eqm parameter.
    p_eqm_initial.tau_I=Params.tau_I;
    
    Params.r=p_eqm_initial.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_initial.w;
    Params.G=p_eqm_initial.G;
    
    %% Now that we have the GE, let's calculate a bunch of related objects   
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    
    % PolicyValues=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_s,d_grid,a_grid, vfoptions); % This will give you the policy in terms of values rather than index
    
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);

    %% Check that we don't hit the top of asset grid (Fig12 is not one from Kitao 2008, just something I want to see)
    Fig12=figure(12);
    assetdist=sum(sum(StationaryDist,4),3);
    assetdist=cumsum(assetdist,1);
    plot(asset_grid,assetdist(:,1),asset_grid,assetdist(:,2))
    title('cdf of HHs over assets')
    legend('worker','entrepreneur')
    xlabel('assets (model units)')
    saveas(Fig12,'./SavedOutput/Graphs/Kitao2008_Fig12.png')
    % The picture is suggestive, the number is the actual test. If any noticeable mass has piled up in
    % the top few gridpoints then the top of the wealth distribution is being clipped and assetmaxfactor
    % needs raising. (Note this is also what stops theta4 entrepreneurs reaching the level of investment
    % at which they would stop being collateral constrained.)
    massattopofassetgrid=sum(sum(sum(sum(StationaryDist(end-9:end,:,:,:)))));
    fprintf('Benchmark: mass in the top 10 asset gridpoints = %e   (top of grid is %.2f model units, $%.0f thousand) \n', massattopofassetgrid, asset_grid(end), asset_grid(end)*wealthscalingfactor/1000)
    if massattopofassetgrid>10^(-5)
        warning('Benchmark: %.4e mass in the top 10 asset gridpoints, raise assetmaxfactor',massattopofassetgrid)
    end

    %% Calculate various statistics related to eqm for Table 2
    FnsToEvaluate.Y_noncorp =  @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)...
        Kitao2008_yFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % output of non-corporate sector
    FnsToEvaluate.IncomeTaxRevenue = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1)...
        Kitao2008_IncomeTaxFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1); % Income Tax Revenue (all of it)
    FnsToEvaluate.IncomeTaxRevenue_NonLinearPartOnly = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1)...
        Kitao2008_IncomeTaxFn_NonLinearPartOnly(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1); % Income Tax Revenue raised by the non-linear tau_a0*() part alone
    FnsToEvaluate.FractionEntrepreneur =  @(aprime,eprime,a,e,eta,theta) (e==1); % entrepreneurs
    FnsToEvaluate.EntrepreneurExit =  @(aprime,eprime,a,e,eta,theta) (eprime==0)*(e==1); % exit of entrepreneurs
    FnsToEvaluate.EntrepreneurEntry =  @(aprime,eprime,a,e,eta,theta) (eprime==1)*(e==0); % exit of entrepreneurs
    FnsToEvaluate.A_entrepreneurs = @(aprime,eprime,a,e,eta,theta) a*(e==1); % Assets owned by entrepreneurs
    FnsToEvaluate.A_workers = @(aprime,eprime,a,e,eta,theta) a*(e==0); % Assets owned by workers
    FnsToEvaluate.Income_workers = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) (e==0)*Kitao2008_IncomeFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % Income of workers
    FnsToEvaluate.Income_entrepreneurs = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) (e==1)*Kitao2008_IncomeFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % Income of entrepreneurs
    
    simoptions.conditionalrestrictions.entrepreneur =  @(aprime,eprime,a,e,eta,theta) (e==1); % entrepreneurs
    simoptions.conditionalrestrictions.worker =  @(aprime,eprime,a,e,eta,theta) (e==0); % worker
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    simoptions=rmfield(simoptions,'conditionalrestrictions');

    % Just give some feedback so we can see nothing looks odd
    AllStats.FractionEntrepreneur.Mean
    [AllStats.A.Mean, AllStats.L.Mean]
    [AllStats.A_entrepreneurs.Mean, AllStats.A_workers.Mean]
    [AllStats.A.Mean, AllStats.K_noncorp.Mean, AllStats.A.Mean-AllStats.K_noncorp.Mean] % assets, capital to entrepreneurs, capital to corporate
    [AllStats.L.Mean, AllStats.N_noncorp.Mean,AllStats.L.Mean-AllStats.N_noncorp.Mean] % labor supply, labor to entrepreneurs, labor to corporate
    
    Output_corp=((AllStats.A.Mean-AllStats.K_noncorp.Mean)^Params.alpha)*((AllStats.L.Mean-AllStats.N_noncorp.Mean)^(1-Params.alpha));
    Y=Output_corp+AllStats.Y_noncorp.Mean;
    [Y,AllStats.Y_noncorp.Mean, Output_corp]

    %% Check the aggregate resource constraint. This is a cheap and very informative diagnostic:
    % it will pick up almost any inconsistency between the household budget constraint (ReturnFn),
    % the aggregation (FnsToEvaluate), and the general eqm conditions.
    % In the stationary eqm:  Y = C + delta*A + G + phi*(borrowing by entrepreneurs)
    % The last term is the intermediation cost, which Kitao (2008) assumes is a pure waste
    % ("thrown away into the ocean"), so it is a use of output but not part of anyones income.
    FnsToEvaluate_ResourceCheck.C = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)...
        Kitao2008_ConsumptionFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1); % consumption
    FnsToEvaluate_ResourceCheck.Borrowing = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)...
        max(Kitao2008_kFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)-a,0); % borrowing is max(k-a,0)
    % Note: this used to be a*Kitao2008_leverageFn(...), which was only max(k-a,0) while leverage was
    % defined as (k-a)/a. It is now (k-a)/k, so compute the borrowing directly and keep the two apart.
    FnsToEvaluate_ResourceCheck.LaborShortfall = @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)...
        (e==1)*min(Kitao2008_nFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)-eta,0); % min(n-eta,0) for entrepreneurs, which should always be exactly zero
    AggVars_ResourceCheck=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate_ResourceCheck,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    ResourceConstraint=Y-AggVars_ResourceCheck.C.Mean-Params.delta*AllStats.A.Mean-Params.G-Params.phi*AggVars_ResourceCheck.Borrowing.Mean;
    fprintf('Aggregate resource constraint, Y-C-delta*A-G-phi*Borrowing = %e  (as a fraction of Y: %e) \n', ResourceConstraint, ResourceConstraint/Y)
    % Entrepreneurs always use all of their own labor endowment in their own project, so the labor
    % they hire, n-eta, must be non-negative for every one of them. If it is not, then the labor market
    % clearing condition is wrong: it computes the corporate sector's labor as L-N_noncorp, which is only
    % equal to (labor supplied by workers)-(labor hired by entrepreneurs) when n>=eta for all entrepreneurs.
    fprintf('Entrepreneurs with n<eta (should be exactly zero), mass-weighted sum of min(n-eta,0) = %e \n', AggVars_ResourceCheck.LaborShortfall.Mean)

    %%
    Table2.CapitalOutputRatio=AllStats.A.Mean/Y; % Capital-Output ratio
    Table2.GdivY=Params.G/Y; % Goverment expenditures/GDP
    % Kitao (2008), Section 3.2, calibrates tau_a2 so that the share of government expenditures raised
    % by the NON-LINEAR PART of the income tax function is 65%. So that, and not the whole income tax,
    % is what should be compared against the 65% target in Table 2. (Her Table 2 row is labelled
    % "income tax/total tax revenue", but the text is explicit that it is the non-linear part.)
    % Note that total tax revenue equals G in equilibrium, since the government budget is balanced.
    Table2.NonLinearIncomeTaxAsShareOfTaxRevenue=AllStats.IncomeTaxRevenue_NonLinearPartOnly.Mean/AllStats.TaxRevenue.Mean;
    Table2.IncomeTaxAsShareOfTaxRevenue=AllStats.IncomeTaxRevenue.Mean/AllStats.TaxRevenue.Mean; % (reported as an extra row, is not the calibration target)
    % Cross-check on the split: with taxincome=1 the income tax is (non-linear part)+tau_I*I for
    % everyone, so the gap between the two must be exactly tau_I times aggregate taxable income
    fprintf('Check of the income tax split, (whole-nonlinear)-tau_I*(aggregate taxable income) = %e \n', (AllStats.IncomeTaxRevenue.Mean-AllStats.IncomeTaxRevenue_NonLinearPartOnly.Mean)-Params.tau_I*(AllStats.Income_workers.Mean+AllStats.Income_entrepreneurs.Mean))
    Table2.FractionEntrepreneur=AllStats.FractionEntrepreneur.Mean;
    Table2.EntrepreneursShareOfIncome=AllStats.Income_entrepreneurs.Mean/(AllStats.Income_workers.Mean+AllStats.Income_entrepreneurs.Mean);
    Table2.ExitRateOfEntrepreneurs=AllStats.EntrepreneurExit.Mean/AllStats.FractionEntrepreneur.Mean;
    % exit rate of new entrants (done seperately below)
    Table2.CapitalUsedByEntrepreneurs=AllStats.K_noncorp.Mean/AllStats.A.Mean;
    Table2.AssetsOwnedByEntrepreneurs=AllStats.A_entrepreneurs.Mean/AllStats.A.Mean;
    Table2.RatioMedianAssets=AllStats.entrepreneur.A_entrepreneurs.Median/AllStats.worker.A_workers.Median;
    
    % To calculate the exit rate of new entrants we need to simulate some panel data and calculate from that
    FnsToEvaluate_e.Entrepreneur = @(aprime,eprime,a,e,eta,theta) e; % Entrepreneurs
    simoptions.simperiods=500;% 500
    simoptions.numbersims=1000; % 1000
    SimPanelValues=SimPanelValues_InfHorz(StationaryDist,Policy,FnsToEvaluate_e,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions);
    % First, count 'new entrant who stays'
    NewEntrantStays=0;
    NewEntrantExits=0;
    for ii=1:simoptions.numbersims
        for tt=2:simoptions.simperiods-1
            if SimPanelValues.Entrepreneur(tt-1,ii)==0 && SimPanelValues.Entrepreneur(tt,ii)==1 % new entrant
                if SimPanelValues.Entrepreneur(tt+1,ii)==1 % stays
                    NewEntrantStays=NewEntrantStays+1;
                elseif  SimPanelValues.Entrepreneur(tt+1,ii)==0 % exits
                    NewEntrantExits=NewEntrantExits+1;
                end
            end
        end
    end
    Table2.ExitRateOfNewEntrants=NewEntrantExits/(NewEntrantExits+NewEntrantStays);
    
    % Table 2
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table2.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lll} \n');
    fprintf(FID, '\\multicolumn{3}{l}{Calibration targets and results} \\\\ \\hline  \n');
    fprintf(FID, '  & Target & Model  \\\\ \\hline \n');
    fprintf(FID, 'capital-output ratio  & 2.65 & %8.3f \\\\ \n', Table2.CapitalOutputRatio);
    fprintf(FID, 'government expenditures/GDP  & 18\\%% & %8.2f \\%% \\\\ \n', 100*Table2.GdivY);
    fprintf(FID, 'non-linear part of income tax/total tax revenue  & 65\\%% & %8.2f \\%% \\\\ \n', 100*Table2.NonLinearIncomeTaxAsShareOfTaxRevenue);
    fprintf(FID, '\\emph{(memo: whole income tax/total tax revenue)}  & - & %8.2f \\%% \\\\ \n', 100*Table2.IncomeTaxAsShareOfTaxRevenue);
    fprintf(FID, ' & &  \\\\ \n');
    fprintf(FID, 'fraction of entrepreneurs  & 12\\%% & %8.2f \\%% \\\\ \n', 100*Table2.FractionEntrepreneur);
    fprintf(FID, 'share of entrepreneurs income  & 27\\%% & %8.2f \\%% \\\\ \n', 100*Table2.EntrepreneursShareOfIncome);
    fprintf(FID, 'exit rate (overall) & 20\\%% & %8.2f \\%% \\\\ \n', 100*Table2.ExitRateOfEntrepreneurs);
    fprintf(FID, 'exit rate (new entrants)  & 40\\%% & %8.2f \\%% \\\\ \n', 100*Table2.ExitRateOfNewEntrants);
    fprintf(FID, 'capital used by entrepreneurs  & 35\\%% & %8.2f \\%% \\\\ \n', 100*Table2.CapitalUsedByEntrepreneurs);
    fprintf(FID, 'assets owned by entrepreneurs  & 40\\%% & %8.2f \\%% \\\\ \n', 100*Table2.AssetsOwnedByEntrepreneurs);
    fprintf(FID, 'ratio of median assets (entrepreneur to worker)  & 8 & %8.2f \\\\ \n', Table2.RatioMedianAssets);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: targets not checked as part of replication. Definition of new entrants assumed to be those that entered last period. \n');
    fprintf(FID, 'The 65\\%% target is for the non-linear part of the income tax function alone: Kitao (2008), Section 3.2, states that $a_2$ is pinned down in equilibrium so that the share of government expenditures raised by the non-linear part of the function equals 65\\%%. The whole income tax as a share of total tax revenue is reported as a memo line, and is not a target. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    % Keep some things for Table 7
    Table7.Entrepreneurs.initial.r=p_eqm_initial.r;
    Table7.Entrepreneurs.initial.w=p_eqm_initial.w;
    Table7.Entrepreneurs.initial.tau_I=p_eqm_initial.tau_I;
    Table7.Entrepreneurs.initial.A=AllStats.A.Mean;
    Table7.Entrepreneurs.initial.Y=Y;
    
    % Keep some things for Figure 8
    Figure8.initial.entry=AllStats.EntrepreneurEntry.Mean;
    Figure8.initial.exit=AllStats.EntrepreneurExit.Mean;
    Figure8.initial.entrepreneurs=AllStats.FractionEntrepreneur.Mean;

    % Keep some things from the benchmark model to use in Figure 4
    Fig4_benchmark.w=p_eqm_initial.w;
    Fig4_benchmark.r=p_eqm_initial.r;
    Fig4_benchmark.tau_I=p_eqm_initial.tau_I;
    Fig4_benchmark.A=AllStats.A.Mean;
    Fig4_benchmark.Y=Y;
    Fig4_benchmark.K_noncorp=AllStats.K_noncorp.Mean;
    Fig4_benchmark.N_noncorp=AllStats.N_noncorp.Mean;
    Fig4_benchmark.Y_noncorp=AllStats.Y_noncorp.Mean;
    Fig4_benchmark.WealthGini=AllStats.A.Gini;
    Fig4_benchmark.FractionEntrepreneur=AllStats.FractionEntrepreneur.Mean;
    
    %% Solve the two "no entrepreneur" economies that are used as comparisons
    % Done here as it requires the GE value of G from the main model (I have inferred this as tau_I is allowed to change)
    % Model without entrepreneurs uses same parameters, except the discount factor (beta) which is recalibrated to target the capital-output ratio
    % First model keeps the exact same process on worker productivity (eta)
    % Second model changes to process on worker productivity to be like that in Castaneda, Diaz-Gimenez & Rios-Rull (2003) with small change so we get same Gini of weath as in the model with entrepreneurs
    if doNoEntrepreneurs==1
        KdivYtarget=AllStats.A.Mean/Y; % K2008 reports 2.65;
        GdivYtarget=Params.G/Y; % K2008 sets G as a fixed fraction of GDP, at 18%
        % Note: it is this RATIO, not the level of G, that the two 'no entrepreneur' economies inherit
        % from the benchmark. They have no non-corporate sector and so a substantially smaller GDP,
        % and holding the level of G fixed would give them a much larger government than 18% of GDP.
        GiniWealthTarget=AllStats.A.Gini; % K2008 reports 0.801
        
        % First model just uses eta_grid and pi_eta
        Output_NoEntrepreneurs1=Kitao2008_NoEntrepreneurs1(Params,n_asset, n_eta, asset_grid, eta_grid, pi_eta,KdivYtarget,GdivYtarget,vfoptions,simoptions,vfoptionstpath);
        Table3.NoEntrepreneurs1=Output_NoEntrepreneurs1.Table3.NoEntrepreneurs1;
        Table7.NoEntrepreneurs1=Output_NoEntrepreneurs1.Table7.NoEntrepreneurs1;
        
        % Second model creates n_eta, eta_grid and pi_eta internally based on Castaneda, Diaz-Gimenez & Rios-Rull (2003)
        % Note: use n_eta=8, as they have 4 points for working age and 4 for retirement.
        % But Kitao (2008) does not have any kind of pension, so I am guessing that she only used the 4 points for working age, which is therefore what I do here.
        % Note this means this Kitao2008_NoEntrepreneurs2() has n_eta=4, rather than the n_eta=5 that Kitao (2008) uses for the main model.
        Output_NoEntrepreneurs2=Kitao2008_NoEntrepreneurs2(Params,n_asset_noE2, asset_grid_noE2,KdivYtarget,GdivYtarget,GiniWealthTarget,vfoptions,simoptions,vfoptionstpath);
        Table3.NoEntrepreneurs2=Output_NoEntrepreneurs2.Table3.NoEntrepreneurs2;
        Table7.NoEntrepreneurs2=Output_NoEntrepreneurs2.Table7.NoEntrepreneurs2;
        save ./SavedOutput/Kitao2008_NoE.mat Table3 Table7
    else
        load ./SavedOutput/Kitao2008_NoE.mat Table3 Table7
    end
    
    %% Inequality statistics related to eqm for Table 3
    FnsToEvaluateIneq.Wealth=FnsToEvaluate.A;
    simoptions.npoints=100; % Use 100 points for the lorenz curve
    AllStats_OnlyWealth=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluateIneq, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);
    
    Table3.benchmark.WealthGini=AllStats_OnlyWealth.Wealth.Gini;
    Table3.benchmark.WealthShareTop1=1-AllStats_OnlyWealth.Wealth.LorenzCurve(99);
    Table3.benchmark.WealthShareTop5=1-AllStats_OnlyWealth.Wealth.LorenzCurve(95);
    Table3.benchmark.WealthShareTop10=1-AllStats_OnlyWealth.Wealth.LorenzCurve(90);
    Table3.benchmark.WealthShareTop20=1-AllStats_OnlyWealth.Wealth.LorenzCurve(80);
    Table3.benchmark.WealthShareTop40=1-AllStats_OnlyWealth.Wealth.LorenzCurve(60);
    Table3.benchmark.WealthShareTop60=1-AllStats_OnlyWealth.Wealth.LorenzCurve(40);
    
    %Table 3
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table3.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllllll} \n');
    fprintf(FID, '\\multicolumn{7}{l}{Wealth distribution: data and models} \\\\  \\hline \n');
    fprintf(FID, '  & Wealth & \\multicolumn{6}{c}{Percentage wealth in the top}  \\\\ \n');
    fprintf(FID, '  & Gini   &  1\\%% & 5\\%% & 10\\%% & 20\\%% & 40\\%% & 60\\%%  \\\\ \\hline \n');
    fprintf(FID, 'US data  & 0.803 & 34.7 & 57.8 & 69.1 & 81.7 & 93.9 & 98.9 \\\\ \n');
    fprintf(FID, 'benchmark model  & %8.3f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',             Table3.benchmark.WealthGini, 100*Table3.benchmark.WealthShareTop1,        100*Table3.benchmark.WealthShareTop5,        100*Table3.benchmark.WealthShareTop10,        100*Table3.benchmark.WealthShareTop20,        100*Table3.benchmark.WealthShareTop40,        100*Table3.benchmark.WealthShareTop60);
    fprintf(FID, 'no entrepreneurs (1) & %8.3f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  Table3.NoEntrepreneurs1.WealthGini, 100*Table3.NoEntrepreneurs1.WealthShareTop1, 100*Table3.NoEntrepreneurs1.WealthShareTop5, 100*Table3.NoEntrepreneurs1.WealthShareTop10, 100*Table3.NoEntrepreneurs1.WealthShareTop20, 100*Table3.NoEntrepreneurs1.WealthShareTop40, 100*Table3.NoEntrepreneurs1.WealthShareTop60);
    fprintf(FID, 'no entrepreneurs (2) & %8.3f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  Table3.NoEntrepreneurs2.WealthGini, 100*Table3.NoEntrepreneurs2.WealthShareTop1, 100*Table3.NoEntrepreneurs2.WealthShareTop5, 100*Table3.NoEntrepreneurs2.WealthShareTop10, 100*Table3.NoEntrepreneurs2.WealthShareTop20, 100*Table3.NoEntrepreneurs2.WealthShareTop40, 100*Table3.NoEntrepreneurs2.WealthShareTop60);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: US data were not checked as part of replication (they come from another paper). "no entrepreneurs (1)" indicates the model with the labor income process as in the benchmark, and in "no entrepreneurs (2)", the process is calibrated to match the wealth Gini of the benchmark. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    %% Figure 1
    
    cdfofassets_workers=cumsum(sum(sum(StationaryDist(:,1,:,:),4),3));
    massofentrepreneurs=sum(sum(sum(StationaryDist(:,2,:,:))));
    cdfofassets_etheta2=cumsum(sum(StationaryDist(:,2,:,2),3));
    cdfofassets_etheta3=cumsum(sum(StationaryDist(:,2,:,3),3));
    cdfofassets_etheta4=cumsum(sum(StationaryDist(:,2,:,4),3));
    fig1=figure(1);
    subplot(1,2,1); plot(asset_grid*wealthscalingfactor/1000,100*cdfofassets_workers);
    title('workers')
    ylabel('percentage (%)')
    xlabel('assets in $1,000')
    xlim([0,500])
    ylim([0,100])
    subplot(1,2,2); plot(asset_grid*wealthscalingfactor/1000,100*(cdfofassets_etheta4+cdfofassets_etheta3+cdfofassets_etheta2)/massofentrepreneurs);
    hold on
    subplot(1,2,2); plot(asset_grid*wealthscalingfactor/1000,100*(cdfofassets_etheta4+cdfofassets_etheta3)/massofentrepreneurs);
    subplot(1,2,2); plot(asset_grid*wealthscalingfactor/1000,100*cdfofassets_etheta4/massofentrepreneurs);
    hold off
    title('entrepreneurs')
    ylabel('percentage (%)')
    xlabel('assets in $1,000')
    xlim([0,2000])
    ylim([0,100])
    legend({'theta 2,3 and 4','theta 3 and 4','theta 4'},'Location','northeast')
    saveas(fig1,'./SavedOutput/Graphs/Kitao2008_Fig1.png')

    %% Table 4
    cdfofassets=cumsum(sum(sum(sum(StationaryDist(:,:,:,:),4),3),2));
    % Top 1 percent that are entrepreneurs
    [~,index]=min(abs(cdfofassets-0.99));
    mass_eintop1=sum(sum(sum(StationaryDist(index:end,2,:,:),4),3),1);
    Table4.entrepreneurinwealthpercTop1=mass_eintop1/0.01;
    % Top 5 percent that are entrepreneurs
    [~,index]=min(abs(cdfofassets-0.95));
    mass_eintop5=sum(sum(sum(StationaryDist(index:end,2,:,:),4),3),1);
    Table4.entrepreneurinwealthpercTop5=mass_eintop5/0.05;
    % Top 10 percent that are entrepreneurs
    [~,index]=min(abs(cdfofassets-0.9));
    mass_eintop10=sum(sum(sum(StationaryDist(index:end,2,:,:),4),3),1);
    Table4.entrepreneurinwealthpercTop10=mass_eintop10/0.1;
    % Top 20 percent that are entrepreneurs
    [~,index]=min(abs(cdfofassets-0.8));
    mass_eintop20=sum(sum(sum(StationaryDist(index:end,2,:,:),4),3),1);
    Table4.entrepreneurinwealthpercTop20=mass_eintop20/0.2;
    
    %Table 4
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table4.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n');
    fprintf(FID, '\\multicolumn{5}{l}{Entrepreneurs and wealth distribution} \\\\  \\hline \n');
    fprintf(FID, '    & \\multicolumn{2}{c}{US data} & \\multicolumn{2}{c}{Benchmark model}  \\\\ \n');
    fprintf(FID, 'top & \\%% of wealth  & \\%% of entrep. & \\%% of wealth  & \\%% of entrep.  \\\\ \n');
    fprintf(FID, '    & held & in percentile & held & in percentile  \\\\ \\hline \n');
    fprintf(FID, '1\\%%  & 30 & 63 & %8.0f & %8.0f \\\\ \n', 100*Table3.benchmark.WealthShareTop1, 100*Table4.entrepreneurinwealthpercTop1);
    fprintf(FID, '5\\%%  & 54 & 49 & %8.0f & %8.0f \\\\ \n', 100*Table3.benchmark.WealthShareTop5, 100*Table4.entrepreneurinwealthpercTop5);
    fprintf(FID, '10\\%% & 67 & 39 & %8.0f & %8.0f \\\\ \n', 100*Table3.benchmark.WealthShareTop10, 100*Table4.entrepreneurinwealthpercTop10);
    fprintf(FID, '20\\%% & 81 & 28 & %8.0f & %8.0f \\\\ \n', 100*Table3.benchmark.WealthShareTop20, 100*Table4.entrepreneurinwealthpercTop20);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\  \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: US data were not checked as part of replication (they come from Cagetti and De Nardi (2006)). \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    %% Worker and Entrepreneur theta prevalence for Figure 2
    Params.theta1=theta_grid(1);
    Params.theta2=theta_grid(2);
    Params.theta3=theta_grid(3);
    Params.theta4=theta_grid(4);
    FnsToEvaluate_Fig2.etheta1=@(aprime,eprime,a,e,eta,theta,theta1) (theta==theta1)*e;
    FnsToEvaluate_Fig2.etheta2=@(aprime,eprime,a,e,eta,theta,theta2) (theta==theta2)*e;
    FnsToEvaluate_Fig2.etheta3=@(aprime,eprime,a,e,eta,theta,theta3) (theta==theta3)*e;
    FnsToEvaluate_Fig2.etheta4=@(aprime,eprime,a,e,eta,theta,theta4) (theta==theta4)*e;
    FnsToEvaluate_Fig2.theta1=@(aprime,eprime,a,e,eta,theta,theta1) (theta==theta1);
    FnsToEvaluate_Fig2.theta2=@(aprime,eprime,a,e,eta,theta,theta2) (theta==theta2);
    FnsToEvaluate_Fig2.theta3=@(aprime,eprime,a,e,eta,theta,theta3) (theta==theta3);
    FnsToEvaluate_Fig2.theta4=@(aprime,eprime,a,e,eta,theta,theta4) (theta==theta4);
    AggVars_Fig2=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate_Fig2,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    
    Figure2values=[AggVars_Fig2.etheta1.Mean, AggVars_Fig2.etheta2.Mean, AggVars_Fig2.etheta3.Mean, AggVars_Fig2.etheta4.Mean;...
        AggVars_Fig2.theta1.Mean, AggVars_Fig2.theta2.Mean, AggVars_Fig2.theta3.Mean, AggVars_Fig2.theta4.Mean];
    
    fig2=figure(2);
    subplot(1,1,1); plot(theta_grid,100*Figure2values(2,:),'b--o')
    hold on
    plot(theta_grid,100*Figure2values(1,:),'g-o')
    hold off
    legend({'workers and entrepreneurs','entrepreneurs'},'Location','northeast')
    ylabel('percentage (%)')
    xlabel('entrepreneurial ability \theta')
    xlim([theta_grid(1), theta_grid(end)])
    saveas(fig2,'./SavedOutput/Graphs/Kitao2008_Fig2.png')
    
    
    %% Table 5 info about entrepreneurs
    FnsToEvaluate_Table5.k=@(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_kFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
    FnsToEvaluate_Table5.a=@(aprime,eprime,a,e,eta,theta) a;
    FnsToEvaluate_Table5.leverage=@(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_leverageFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
    
    simoptions.conditionalrestrictions.etheta2=@(aprime,eprime,a,e,eta,theta,theta2) (theta==theta2)*e;
    simoptions.conditionalrestrictions.etheta3=@(aprime,eprime,a,e,eta,theta,theta3) (theta==theta3)*e;
    simoptions.conditionalrestrictions.etheta4=@(aprime,eprime,a,e,eta,theta,theta4) (theta==theta4)*e;
    simoptions.conditionalrestrictions.entrepreneur=@(aprime,eprime,a,e,eta,theta,theta2) e;
    AllStats_Table5=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate_Table5,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    simoptions=rmfield(simoptions,'conditionalrestrictions');

    Table5.thetavalue=theta_grid;
    Table5.thetainpopulation=Figure2values(1,:);
    Table5.thetainentrepreneurs=Figure2values(1,:)/AllStats.FractionEntrepreneur.Mean; % Note: AggVars.FractionEntrepreneur.Mean is just same as sum(Figure2values(1,:))
    % Table5.avginvest=[0, AggVars_Table5.etheta2_k.Mean/AggVars_Fig2.etheta2.Mean, AggVars_Table5.etheta3_k.Mean/AggVars_Fig2.etheta3.Mean, AggVars_Table5.etheta4_k.Mean/AggVars_Fig2.etheta4.Mean];
    % Table5.avgasset=[0, AggVars_Table5.etheta2_a.Mean/AggVars_Fig2.etheta2.Mean, AggVars_Table5.etheta3_a.Mean/AggVars_Fig2.etheta3.Mean, AggVars_Table5.etheta4_a.Mean/AggVars_Fig2.etheta4.Mean];
    % Table5.avgleverageratio=[1, AggVars_Table5.etheta2_lev.Mean/AggVars_Fig2.etheta2.Mean, AggVars_Table5.etheta3_lev.Mean/AggVars_Fig2.etheta3.Mean, AggVars_Table5.etheta4_lev.Mean/AggVars_Fig2.etheta4.Mean];
    Table5.avginvest=[0, AllStats_Table5.etheta2.k.Mean, AllStats_Table5.etheta3.k.Mean, AllStats_Table5.etheta4.k.Mean];
    Table5.avgasset=[0, AllStats_Table5.etheta2.a.Mean, AllStats_Table5.etheta3.a.Mean, AllStats_Table5.etheta4.a.Mean];
    Table5.avgleverageratio=[1, AllStats_Table5.etheta2.leverage.Mean, AllStats_Table5.etheta3.leverage.Mean, AllStats_Table5.etheta4.leverage.Mean];
    Table5.avgleverageratiopct=100*Table5.avgleverageratio; % turn into percentage
    % Bottom two rows of table 5
    Table5.total=[100*AllStats.FractionEntrepreneur.Mean,100];
    Table5.average=[AllStats_Table5.entrepreneur.k.Mean,AllStats_Table5.entrepreneur.a.Mean,AllStats_Table5.entrepreneur.leverage.Mean];
    
    % Table 5
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table5.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n');
    fprintf(FID, '\\multicolumn{7}{l}{Benchmark model: entrepreneurial activities by ability $\\theta$ } \\\\  \\hline \n');
    fprintf(FID, ' $\\theta$ & $\\theta$ & \\%%    & \\%%       & avg. invstm. & avg. assets  & avg. lev. ratio  \\\\ \n');
    fprintf(FID, ' grid      & value     & in pop. & in entrep. & in \\$1,000   & in \\$1,000 & \\%%    \\\\ \\hline \n');
    fprintf(FID, '$\\theta_1$  & %8.3f & %8.2f\\%% & %8.2f\\%% & - & - & - \\\\ \n',            Table5.thetavalue(1), 100*Table5.thetainpopulation(1), 100*Table5.thetainentrepreneurs(1));
    fprintf(FID, '$\\theta_2$  & %8.3f & %8.2f\\%% & %8.2f\\%% & %8.0f & %8.0f & %8.1f\\%% \\\\ \n',  Table5.thetavalue(2), 100*Table5.thetainpopulation(2), 100*Table5.thetainentrepreneurs(2),Table5.avginvest(2)*wealthscalingfactor/1000,Table5.avgasset(2)*wealthscalingfactor/1000,Table5.avgleverageratiopct(2));
    fprintf(FID, '$\\theta_3$  & %8.3f & %8.2f\\%% & %8.2f\\%% & %8.0f & %8.0f & %8.1f\\%% \\\\ \n',  Table5.thetavalue(3), 100*Table5.thetainpopulation(3), 100*Table5.thetainentrepreneurs(3),Table5.avginvest(3)*wealthscalingfactor/1000,Table5.avgasset(3)*wealthscalingfactor/1000,Table5.avgleverageratiopct(3));
    fprintf(FID, '$\\theta_4$  & %8.3f & %8.2f\\%% & %8.2f\\%% & %8.0f & %8.0f & %8.1f\\%% \\\\ \n',  Table5.thetavalue(4), 100*Table5.thetainpopulation(4), 100*Table5.thetainentrepreneurs(4),Table5.avginvest(4)*wealthscalingfactor/1000,Table5.avgasset(4)*wealthscalingfactor/1000,Table5.avgleverageratiopct(4));
    fprintf(FID, 'total        & -     & %8.2f\\%% & %8.2f\\%% & -  & -  & - \\\\ \n',          Table5.total(1), Table5.total(2));
    fprintf(FID, 'average      & -     & -         & -         & %8.1f & %8.1f & %8.1f\\%% \\\\ \n',  Table5.average(1)*wealthscalingfactor/1000, Table5.average(2)*wealthscalingfactor/1000, 100*Table5.average(3));
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: average refers to conditional on being an entrepreneur. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    
    % Keep some things for Table 8
    Table8.avginvest=[Table5.avginvest(2:4),AllStats_Table5.k.Mean/AllStats.FractionEntrepreneur.Mean];
    
    %% Figure 3 requires us to plot the capital k
    FnsToEvaluate_k.investment=@(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_kFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate_k,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    
    investmentpolicy_theta2=ValuesOnGrid.investment(:,2,1,2); % all a, second e is entrepreneur, first eta, second theta
    investmentpolicy_theta3=ValuesOnGrid.investment(:,2,1,3); % all a, second e is entrepreneur, first eta, third theta
    investmentpolicy_theta4=ValuesOnGrid.investment(:,2,1,4); % all a, second e is entrepreneur, first eta, fourth theta
    % Note: Kitao (2008) does not mention value of eta, presumably this policy is independent of eta? (it is looking at policy; thinking about it, it will be independent because is static decision and eta is irrelevant to that static decision)
    
    % plot the investment policy against asset_grid
    fig3=figure(3);
    plot(asset_grid*wealthscalingfactor/1000,(1+Params.d)*asset_grid*wealthscalingfactor/1000,'k-.') % max leverage line
    hold on
    plot(asset_grid*wealthscalingfactor/1000,asset_grid*wealthscalingfactor/1000,'k:') % 45-degree line
    plot(asset_grid*wealthscalingfactor/1000,investmentpolicy_theta2*wealthscalingfactor/1000,'b--')
    plot(asset_grid*wealthscalingfactor/1000,investmentpolicy_theta3*wealthscalingfactor/1000,'g--')
    plot(asset_grid*wealthscalingfactor/1000,investmentpolicy_theta4*wealthscalingfactor/1000,'r--')
    hold off
    legend({'max leverage','45 degree line','theta2','theta3','theta4'},'Location','northwest')
    ylabel('investment in $1,000')
    xlabel('assets in $1,000')
    xlim([0,8000])
    saveas(fig3,'./SavedOutput/Graphs/Kitao2008_Fig3.png')

    
    %% Save the initial general eqm for later transition paths
    save ./SavedOutput/Kitao2008_InitialGE.mat p_eqm_initial
    % Save only what the later doPart blocks need, rather than the whole workspace. Saving the whole
    % workspace made the .mat files cascade (the last one reached 580MB) and, worse, silently kept
    % stale variables in scope: Figure 9 in doPart(7) was plotting savingsrate_path1 and
    % savingsrate_path2 left over from doPart(3)'s tau_k transitions. Listing what is handed on turns
    % that class of mistake into an undefined-variable error. Anything a block computes that nothing
    % later needs (Vpath, PolicyPath, AgentDistPath, ...) simply stays out, so it is stored once in
    % its own block's .mat rather than being copied into every later one.
    save ./SavedOutput/Kitao2008_doPart1.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8
else
    load ./SavedOutput/Kitao2008_doPart1.mat
end
% Create a Params_init so I can use it more easily later
Params_init=Params;
Params_init.taxincome=1;
Params_init.r=p_eqm_initial.r;
Params_init.w=p_eqm_initial.w;
Params_init.G=p_eqm_initial.G;
Params_init.tau_I=p_eqm_initial.tau_I;



%% Done with the baseline setup. Now we turn to looking at a capital income tax (together with seperate labor income tax; rather than the income tax that was in baseline)

GEPriceParamNames={'r','w','tau_I'}; % tau_I instead of G
% Note: the general eqm conditions are unchanged, is just that the
% government budget balance is now about finding tau_I, where previously it
% was about finding G

if doPart(2)==1
    Params.taxincome=2; % 1 is tax income, 2 is to tax capital income and labor income seperately
    % Kitao (2008) does tau_k from 0 to 40 in intervals of 5. Here do intervals of 1
    tau_k_vec=0:0.01:0.4;
    % Work out which of these stationary equilibria doPart(4) will later want V and Policy for. Uses the
    % same lookup as doPart(4) does, so the two cannot disagree, and errors out if they do not line up.
    keepVPolicy_tauk=false(1,length(tau_k_vec));
    for kk=1:length(tau_k_transition_vec)
        [~,jj]=min(abs(tau_k_vec-tau_k_transition_vec(kk))); % nearest value, not exact equality (see note at the doPart(4) lookup)
        keepVPolicy_tauk(jj)=true;
    end
    if sum(keepVPolicy_tauk)~=length(tau_k_transition_vec)
        error('tau_k_transition_vec does not line up with tau_k_vec, so V and Policy would not be kept for every transition')
    end
    %% Diagnostic: reproduce the failing EvalFnOnAgentDist call and report its inputs.
    % Kitao2008_CompileTest compiles every one of these Fns on gpu01/R2025b, through the toolkit's
    % own EvalFnOnAgentDist_Grid, via the anonymous wrapper, with the real Params and the real grid
    % sizes, and they all pass. Yet the same call fails here at doPart(2) iteration 1. So rather
    % than keep guessing at the difference, print exactly what the toolkit derives and passes,
    % one FnsToEvaluate field at a time, and let the diary settle it.
    % Set doDiagnoseFns=0 once this is resolved.
    doDiagnoseFns=0;
    if doDiagnoseFns==1
        Params.tau_k=0; % the configuration of iteration 1, which is where it fails
        Params.r=p_eqm_initial.r; Params.w=p_eqm_initial.w; Params.tau_I=p_eqm_initial.tau_I;
        fprintf('\n=== diagnose: taxincome=%d tau_k=%g, solving V and dist at the initial eqm prices \n', Params.taxincome, Params.tau_k)
        [V_d,Policy_d]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
        StationaryDist_d=StationaryDist_InfHorz(Policy_d,n_d,n_a,n_z,pi_z, simoptions);
        fprintf('=== diagnose: V and dist solved (so the ReturnFn compiles); now each Fn on its own \n')
        l_d_d=0;
        if n_d(1)~=0
            l_d_d=length(n_d);
        end
        ndrop=l_d_d+length(n_a)+length(n_a)+length(n_z); % l_d+l_aprime+l_a+l_z, as the toolkit does it
        fnames_d=fieldnames(FnsToEvaluate);
        for ffd=1:length(fnames_d)
            nm=fnames_d{ffd};
            temp=getAnonymousFnInputNames(FnsToEvaluate.(nm));
            pn=cell(0);
            if length(temp)>ndrop
                pn={temp{ndrop+1:end}};
            end
            fprintf('  %-34s %2d inputs, %2d params after dropping %d: %s \n', nm, length(temp), length(pn), ndrop, strjoin(pn,','));
            for pp=1:length(pn)
                if isfield(Params,pn{pp})
                    vv=Params.(pn{pp});
                    fprintf('      %-12s %-10s %-8s %s \n', pn{pp}, class(vv), mat2str(size(vv)), mat2str(gather(vv)));
                else
                    fprintf('      %-12s NOT IN Params \n', pn{pp});
                end
            end
            oneFn=struct();
            oneFn.(nm)=FnsToEvaluate.(nm);
            try
                AV=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_d, Policy_d, oneFn, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);
                fprintf('    -> PASS   %s = %.6g \n', nm, AV.(nm).Mean);
            catch ME
                fprintf('    -> FAIL   %s \n', ME.message);
            end
        end
        fprintf('=== diagnose: done \n\n')
    end

    for ii=1:length(tau_k_vec)
        Params.tau_k=tau_k_vec(ii);
        % Report which experiment we are in BEFORE the solve. The bare ii further down only echoes
        % once HeteroAgentStationaryEqm_InfHorz has returned, so if a solve fails you cannot tell
        % which of the 41 it was without counting echoes.
        fprintf('doPart(2): iteration %d of %d, tau_k=%.2f \n', ii, length(tau_k_vec), Params.tau_k)
        
        [p_eqm_final,GeneralEqmCondn]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
        
        ii
        p_eqm_final % The equilibrium values of the GE prices
        
        % Now that we have the GE, let's calculate a bunch of related objects
        Params.r=p_eqm_final.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
        Params.w=p_eqm_final.w;
        Params.tau_I=p_eqm_final.tau_I;
        % Note: this will also get used as the initial guess next iteration
        
        [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
        
        StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);
        
        AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
        
        Output_corp=((AllStats.A.Mean-AllStats.K_noncorp.Mean)^Params.alpha)*((AllStats.L.Mean-AllStats.N_noncorp.Mean)^(1-Params.alpha));
        Y=Output_corp+AllStats.Y_noncorp.Mean;
        [Y,AllStats.Y_noncorp.Mean, Output_corp]
        
        Figure4data(ii).w=p_eqm_final.w;
        Figure4data(ii).r=p_eqm_final.r;
        Figure4data(ii).tau_I=p_eqm_final.tau_I;
        Figure4data(ii).aggcapital=AllStats.A.Mean;
        Figure4data(ii).aggoutput=Y;
        Figure4data(ii).Y_noncorp=AllStats.Y_noncorp.Mean;
        Figure4data(ii).K_noncorp=AllStats.K_noncorp.Mean;
        Figure4data(ii).N_noncorp=AllStats.N_noncorp.Mean;
        Figure4data(ii).WealthGini=AllStats.A.Gini;
        % To be able to do an alternative 'per entrepreneur' version I need
        Figure4data(ii).FractionEntrepreneur=AllStats.FractionEntrepreneur.Mean;


        % Only the tau_k values that get a transition path in doPart(4) ever have their V and Policy read
        % back, so only keep those (all 41 would be about 65MB of the .mat for nothing). StationaryDist is
        % not kept at all: nothing reads it back (Table 4 uses the benchmark distribution from doPart(1)).
        if keepVPolicy_tauk(ii)
            Figure4data(ii).V=V;
            Figure4data(ii).Policy=Policy;
        end
        
        % Mainly just in case we want to check how they went
        Figure4data(ii).p_eqm_final=p_eqm_final;
        Figure4data(ii).GeneralEqmCondn=GeneralEqmCondn;
        
        if Params.tau_k==0
            p_eqm_final1=p_eqm_final;
        elseif Params.tau_k==0.4
            p_eqm_final2=p_eqm_final;
        end
        
        Figure4data(ii).MassAtTopOfAssetGrid=sum(sum(sum(sum(StationaryDist(end-9:end,:,:,:)))));

        save ./SavedOutput/Kitao2008_Figure4data.mat Figure4data
    end

    % Across all 41 stationary equilibria, has anyone piled up at the top of the asset grid? Aggregate
    % capital is highest when tau_k is lowest, so the binding case is normally tau_k=0.
    [worstmass,worstii]=max(cell2mat({Figure4data(:).MassAtTopOfAssetGrid}));
    fprintf('tau_k experiments: worst mass in the top 10 asset gridpoints = %e, at tau_k=%.2f \n', worstmass, tau_k_vec(worstii))
    if worstmass>10^(-5)
        warning('tau_k experiments: %.4e mass in the top 10 asset gridpoints at tau_k=%.2f, raise assetmaxfactor',worstmass,tau_k_vec(worstii))
    end
    
    save ./SavedOutput/Kitao2008_FinalGE.mat p_eqm_final1  p_eqm_final2
    

    % Figure 4
    fig4=figure(4);
    subplot(2,3,1); plot(100*tau_k_vec,Fig4_benchmark.w*ones(1,length(tau_k_vec)),'b.',100*tau_k_vec,cell2mat({Figure4data(:).w}),'b-');
    title('wage')
    subplot(2,3,2); plot(100*tau_k_vec,100*Fig4_benchmark.r*ones(1,length(tau_k_vec)),'b.',100*tau_k_vec,100*cell2mat({Figure4data(:).r}),'b-');
    title('interest rate (%)')
    subplot(2,3,3); plot(100*tau_k_vec,100*Fig4_benchmark.tau_I*ones(1,length(tau_k_vec)),'b.',100*tau_k_vec,100*cell2mat({Figure4data(:).tau_I}),'b-');
    title('proportional tax \tau_I (%)')
    subplot(2,3,4); yyaxis left 
    plot(100*tau_k_vec,Fig4_benchmark.A*ones(1,length(tau_k_vec)),'b.',100*tau_k_vec,cell2mat({Figure4data(:).aggcapital}),'b-');
    subplot(2,3,4); yyaxis right 
    plot(100*tau_k_vec,Fig4_benchmark.Y*ones(1,length(tau_k_vec)),'r.',100*tau_k_vec, cell2mat({Figure4data(:).aggoutput}),'r-');
    yyaxis left 
    title('aggregate activities')
    legend('','agg. capital','','agg. output','northwest')
    subplot(2,3,5); plot(100*tau_k_vec,cell2mat({Figure4data(:).Y_noncorp})/Fig4_benchmark.Y_noncorp,'b-',100*tau_k_vec,cell2mat({Figure4data(:).K_noncorp})/Fig4_benchmark.K_noncorp,'g-',100*tau_k_vec, cell2mat({Figure4data(:).N_noncorp})/Fig4_benchmark.N_noncorp,'r-');
    title('entrepreneurial activities (normalized)')
    legend('output','capital','labor','northeast')
    subplot(2,3,6); plot(100*tau_k_vec,Fig4_benchmark.WealthGini*ones(1,length(tau_k_vec)),'b.',100*tau_k_vec,cell2mat({Figure4data(:).WealthGini}),'b-');
    title('Wealth Gini')
    saveas(fig4,'./SavedOutput/Graphs/Kitao2008_Fig4.png')
    
    % Do an alternative version of Fig 4, where only change is that in the fifth panel, I have 'per entrepreneur'
    Y_noncorp_perentrepreneur=(cell2mat({Figure4data(:).Y_noncorp})./cell2mat({Figure4data(:).FractionEntrepreneur}))/(Fig4_benchmark.Y_noncorp/Fig4_benchmark.FractionEntrepreneur);
    K_noncorp_perentrepreneur=(cell2mat({Figure4data(:).K_noncorp})./cell2mat({Figure4data(:).FractionEntrepreneur}))/(Fig4_benchmark.K_noncorp/Fig4_benchmark.FractionEntrepreneur);
    N_noncorp_perentrepreneur=(cell2mat({Figure4data(:).N_noncorp})./cell2mat({Figure4data(:).FractionEntrepreneur}))/(Fig4_benchmark.N_noncorp/Fig4_benchmark.FractionEntrepreneur);
    subplot(2,3,5); plot(100*tau_k_vec,Y_noncorp_perentrepreneur,'b-',100*tau_k_vec,K_noncorp_perentrepreneur,'g-',100*tau_k_vec, N_noncorp_perentrepreneur,'r-');
    title('entrepreneurial activities (normalized)')
    legend('output','capital','labor','northeast')
    saveas(fig4,'./SavedOutput/Graphs/Kitao2008_Fig4B.png')

    [~,ii_tau_k_0]=min(abs(tau_k_vec-0)); % find tau_k=0
    [~,ii_tau_k_0p4]=min(abs(tau_k_vec-0.4)); % find tau_k=0.4

    
    % Keep some things for Table 7
    Table7.Entrepreneurs.tau_k_final1.r=Figure4data(ii_tau_k_0).r;
    Table7.Entrepreneurs.tau_k_final1.w=Figure4data(ii_tau_k_0).w;
    Table7.Entrepreneurs.tau_k_final1.tau_I=Figure4data(ii_tau_k_0).tau_I;
    Table7.Entrepreneurs.tau_k_final1.A=Figure4data(ii_tau_k_0).aggcapital;
    Table7.Entrepreneurs.tau_k_final1.Y=Figure4data(ii_tau_k_0).aggoutput;
    
    Table7.Entrepreneurs.tau_k_final2.r=Figure4data(ii_tau_k_0p4).r;
    Table7.Entrepreneurs.tau_k_final2.w=Figure4data(ii_tau_k_0p4).w;
    Table7.Entrepreneurs.tau_k_final2.tau_I=Figure4data(ii_tau_k_0p4).tau_I;
    Table7.Entrepreneurs.tau_k_final2.A=Figure4data(ii_tau_k_0p4).aggcapital;
    Table7.Entrepreneurs.tau_k_final2.Y=Figure4data(ii_tau_k_0p4).aggoutput;
    
    % As for doPart(1) above, plus what this block adds: Params_init, p_eqm_final1, p_eqm_final2, tau_k_vec, keepVPolicy_tauk, Figure4data
    save ./SavedOutput/Kitao2008_doPart2.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8 Params_init p_eqm_final1 p_eqm_final2 tau_k_vec keepVPolicy_tauk Figure4data
else
    load ./SavedOutput/Kitao2008_doPart2.mat
end

%% Test a 'nothing happens' transition [works just fine, one iteration then stops, shows nothing really changing]
doTest=0
if doTest==1
    Params=Params_init;    
    
    % Initial eqm: Params_init
    [V_init,Policy_init]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params_init, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_init=StationaryDist_InfHorz(Policy_init,n_d,n_a,n_z,pi_z, simoptions);

    % Final eqm: same
    V_final=V_init;
    Policy_final=Policy_init;
    StationaryDist_final=StationaryDist_init;

    T=120 % Kitao (2008) graphs suggest she uses 50 periods, but this is not enough
    
    FnsToEvaluate_TransPath.K_noncorp=FnsToEvaluate.K_noncorp;
    FnsToEvaluate_TransPath.A=FnsToEvaluate.A;
    FnsToEvaluate_TransPath.N_noncorp=FnsToEvaluate.N_noncorp;
    FnsToEvaluate_TransPath.L=FnsToEvaluate.L;
    FnsToEvaluate_TransPath.TaxRevenue=FnsToEvaluate.TaxRevenue;
 
    % Add some things to FnsToEvaluate_TransPath that we are interested in but didn't need for solving the transition path
    FnsToEvaluate2_TransPath=FnsToEvaluate_TransPath;
    FnsToEvaluate2_TransPath.Y_noncorp =  @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_yFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % output of non-corporate sector
    FnsToEvaluate2_TransPath.C =  @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)...
        Kitao2008_ConsumptionFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1); % consumption
  
    
    TransPathGeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,L,alpha,delta) r-(alpha*((A-K_noncorp)^(alpha-1))*((L-N_noncorp)^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,L,alpha) w-(1-alpha)*((A-K_noncorp)^(alpha))*((L-N_noncorp)^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; % Government runs balanced budget
    % Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.
    
    transpathoptions.GEnewprice=2; % Anderson acceleration of the shooting update
    % Need to explain to transpathoptions how to use the GeneralEqmEqns to
    % update the general eqm transition prices (in PricePath).
    transpathoptions.GEnewprice2.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'CapitalMarket','r',0,0.1;... % CapitalMarket is positive is r is to large, so subtract
        'LaborMarket','w',0,0.1;... % LaborMarket is positive is r is to large, so subtract
        'GovBudget','tau_I',1,0.1}; % GovBudget is positive if tau_I is too small, so add
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.
    % Note: GEnewprice=2 is Anderson acceleration, which accelerates exactly this shooting update, so it
    % needs the same howtoupdate instructions (just under GEnewprice2 rather than GEnewprice3).
    
    % Now just run the TransitionPath_InfHorz command (all of the other inputs
    % are things we had already had to define to be able to solve for the initial and final equilibria)
    transpathoptions.verbose=1;

    transpathoptions.graphpricepath=1; % 1: creates a graph of the 'current' price path which updates each iteration.
    transpathoptions.graphaggvarspath=1; % 1: creates a graph of the 'current' aggregate variables which updates each iteration.
    if headlessFigures==1 % these two redraw every iteration, which is pointless offscreen and costs real time on a server
        transpathoptions.graphpricepath=0;
        transpathoptions.graphaggvarspath=0;
    end
    
    %% Compute the transition path to where nothing happens

    % We want to look at a one off unanticipated path of tau_k.
    ParamPath.tau_k=Params_init.tau_k*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
    % (the way ParamPath is set is designed to allow for a series of changes in the parameters)
    
    % We need to give an initial guess for the price path
    PricePath0.r=p_eqm_initial.r*ones(T,1); % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
    PricePath0.w=p_eqm_initial.w*ones(T,1);
    PricePath0.tau_I=p_eqm_initial.tau_I*ones(T,1);
    
    [PricePath,GEcondnPath1]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions,vfoptionstpath);
    
    % overrule
    PricePath=PricePath0;

    [Vpath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
    AgentDistPath=AgentDistOnTransPath_InfHorz(StationaryDist_init, PricePath, ParamPath, PolicyPath,n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
    AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2_TransPath,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    

    % Nothing changes, so every one of these should be (near) zero. If the first one is not, then the
    % CEV calculations further below are meaningless, as they compare Vpath(...,1) against V_init and
    % those two are solved by different codes (V_init by ValueFnIter_InfHorz with vfoptions, Vpath by
    % ValueFnOnTransPath_InfHorz with vfoptionstpath).
    % Which AgentDistOnTransPath_InfHorz did this run actually load, and does that copy carry the
    % L2flag override? Dropbox sync to the server can lag, so a run can silently use the old file
    % and produce output byte-identical to the pre-fix run.
    whichAD=which('AgentDistOnTransPath_InfHorz');
    fprintf('doTest: AgentDistOnTransPath_InfHorz loaded from %s \n', whichAD)
    fprintf('doTest: that copy carries the L2flag override: %d \n', contains(fileread(whichAD),'L2flag'))
    % Firing census: if the flag is 2 everywhere then the override is a no-op here by construction,
    % and identical output proves nothing about whether the fix is present.
    L2flag_init=gather(Policy_init(end,:));
    fprintf('doTest: size(Policy_init,1)=%d (last channel is the L2flag) \n', size(Policy_init,1))
    fprintf('doTest: L2flag fires at %d of %d states (%.4f%%): force-lower %d, force-upper %d \n', ...
        sum(L2flag_init~=2), numel(L2flag_init), 100*sum(L2flag_init~=2)/numel(L2flag_init), ...
        sum(L2flag_init==1), sum(L2flag_init==3))

    fprintf('doTest: max abs deviation of the value fn from the initial stationary eqm value fn \n')
    fprintf('  period 1:   %e \n', max(abs(Vpath(:,:,:,:,1)-V_init),[],'all'))
    fprintf('  period T-1: %e \n', max(abs(Vpath(:,:,:,:,T-1)-V_init),[],'all'))
    fprintf('  period T:   %e \n', max(abs(Vpath(:,:,:,:,T)-V_init),[],'all'))
    fprintf('doTest: max abs deviation of the value fn between period 1 and period 11: %e \n', max(abs(Vpath(:,:,:,:,1)-Vpath(:,:,:,:,11)),[],'all'))
    fprintf('doTest: number of policy indexes that differ from the initial stationary eqm in period 1: %d \n', sum(PolicyPath(:,:,:,:,:,1)~=Policy_init,'all'))
    fprintf('doTest: max abs deviation of the agent dist from the initial stationary dist \n')
    fprintf('  period 1: %e \n', max(abs(AgentDistPath(:,:,:,:,1)-StationaryDist_init),[],'all'))
    fprintf('  period T: %e \n', max(abs(AgentDistPath(:,:,:,:,T)-StationaryDist_init),[],'all'))
    % The elementwise max understates this: the distribution has ~48,000 cells summing to one, so a
    % typical cell holds ~2e-5 and an elementwise max of 1e-4 is already a large relative error. Report
    % the summed absolute deviation, and the aggregate that Figures 5 and 9 actually plot.
    fprintf('doTest: sum abs deviation of the agent dist from the initial stationary dist \n')
    fprintf('  period 1: %e \n', sum(abs(AgentDistPath(:,:,:,:,1)-StationaryDist_init),'all'))
    fprintf('  period T: %e \n', sum(abs(AgentDistPath(:,:,:,:,T)-StationaryDist_init),'all'))
    AggVars_doTest=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_init, Policy_init, FnsToEvaluate2_TransPath,Params_init, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    fprintf('doTest: aggregate assets, none of which should move \n')
    fprintf('  initial stationary dist: %.6f \n', AggVars_doTest.A.Mean)
    for tt=[1,2,3,5,10,20,50,T]
        fprintf('  period %3d: A=%.6f   sum|dist-init|=%.6e \n', tt, AggVarsPath.A.Mean(tt), ...
            sum(abs(AgentDistPath(:,:,:,:,tt)-StationaryDist_init),'all'))
    end
    % Total mass. If this is not 1 at every t the iteration is destroying or creating agents,
    % which is a different bug from merely putting them in the wrong place.
    fprintf('doTest: total mass at t=1,2,T: %.10f %.10f %.10f \n', ...
        sum(AgentDistPath(:,:,:,:,1),'all'), sum(AgentDistPath(:,:,:,:,2),'all'), sum(AgentDistPath(:,:,:,:,T),'all'))
    % StationaryDist_init is invariant under the true iteration, so whatever deviation appears at
    % period 2 IS the single-step error of AgentDistOnTransPath_InfHorz. Localise it.
    dev2=gather(AgentDistPath(:,:,:,:,2)-StationaryDist_init);
    [~,ord]=sort(abs(dev2(:)),'descend');
    fprintf('doTest: the 8 largest period-2 deviations \n')
    for kk=1:8
        [i1,i2,i3,i4]=ind2sub(size(dev2),ord(kk));
        fprintf('  a=%4d(of %d)  e=%d  eta=%d  theta=%d   dev=%+.6e \n', i1,n_a(1),i2,i3,i4, dev2(ord(kk)))
    end
    fprintf('doTest: period-2 deviation by asset region: bottom 10 gridpts %+.6e, top 10 %+.6e, middle %+.6e \n', ...
        sum(dev2(1:10,:,:,:),'all'), sum(dev2(end-9:end,:,:,:),'all'), sum(dev2(11:end-10,:,:,:),'all'))
    % a-marginal of the period-2 deviation. The signed version largely cancels (the error is a
    % redistribution within the asset direction), so report |dev| too, and how concentrated it is:
    % a handful of gridpoints means a specific index bug, hundreds means a weight formula.
    dev2_a=squeeze(sum(sum(sum(dev2,2),3),4));      % signed, over a
    adev  =squeeze(sum(sum(sum(abs(dev2),2),3),4)); % |dev|, over a
    fprintf('doTest: a-gridpoints carrying the period-2 error: %d above 1e-6, %d above 1e-5, %d above 1e-4 (of %d) \n', ...
        sum(adev>1e-6), sum(adev>1e-5), sum(adev>1e-4), n_a(1))
    fprintf('doTest: total |deviation| %.6e, of which the top 20 a-gridpoints hold %.2f%% \n', ...
        sum(adev), 100*sum(maxk(adev,20))/sum(adev))
    [~,aord]=sort(adev,'descend');
    fprintf('doTest: top 15 a-gridpoints by |dev| (index, a in $1,000, signed dev, |dev|) \n')
    for kk=1:15
        jj=aord(kk);
        fprintf('  a=%4d  $%9.1f  signed %+.6e  abs %.6e \n', jj, a_grid(jj)*wealthscalingfactor/1000, dev2_a(jj), adev(jj))
    end
    % The z-marginal is the stationary distribution of pi_z and must not move at all. If it does,
    % the exogenous transition applied along the path is not the pi_z the stationary solve used.
    zm=@(D) squeeze(sum(sum(gather(D),1),2)); % sum over a and e -> [n_z(1),n_z(2)]
    zm_init=zm(StationaryDist_init); zm_2=zm(AgentDistPath(:,:,:,:,2)); zm_T=zm(AgentDistPath(:,:,:,:,T));
    fprintf('doTest: z-marginal max|change|: t=2 %.6e, t=T %.6e \n', ...
        max(abs(zm_2-zm_init),[],'all'), max(abs(zm_T-zm_init),[],'all'))
    fprintf('doTest: eta marginal   init/t=2: %s / %s \n', mat2str(round(sum(zm_init,2)',6)), mat2str(round(sum(zm_2,2)',6)))
    fprintf('doTest: theta marginal init/t=2: %s / %s \n', mat2str(round(sum(zm_init,1),6)), mat2str(round(sum(zm_2,1),6)))
    % Mass must not move between the two occupations beyond what the policy itself sends.
    fprintf('doTest: mass by occupation: init %.6f/%.6f, t=2 %.6f/%.6f, t=T %.6f/%.6f \n', ...
        sum(StationaryDist_init(:,1,:,:),'all'), sum(StationaryDist_init(:,2,:,:),'all'), ...
        sum(AgentDistPath(:,1,:,:,2),'all'), sum(AgentDistPath(:,2,:,:,2),'all'), ...
        sum(AgentDistPath(:,1,:,:,T),'all'), sum(AgentDistPath(:,2,:,:,T),'all'))

    % Save the inputs to one distribution step so both operators can be rebuilt and diffed offline
    % on CPU (48,000 states, ~96,000 nonzeros). Small: Policy is 4x1200x2x5x4, the dist 48,000.
    AgentDistPath_t2=AgentDistPath(:,:,:,:,2);
    save ./SavedOutput/Kitao2008_doTestStep.mat Policy_init StationaryDist_init AgentDistPath_t2 pi_z a_grid n_a n_z simoptions

    % Same again, but with lowmemory, as a cross-check that the two give the same answer
    vfoptionstpath.lowmemory=1;
    [Vpath_lowmem,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
    vfoptionstpath=rmfield(vfoptionstpath,'lowmemory');
    fprintf('doTest: max abs deviation between the lowmemory=0 and lowmemory=1 value fn paths: %e \n', max(abs(Vpath(:)-Vpath_lowmem(:))))

end



%% Transition: two of them, to tau_k=0 and tau_k=0.4
if doPart(3)==1
    % Note: currently have Params.taxincome=2 from above    
    tau_k_final1=0;
    tau_k_final2=0.4;
    
    % Initial eqm: Params_init
    [V_init,Policy_init]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params_init, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_init=StationaryDist_InfHorz(Policy_init,n_d,n_a,n_z,pi_z, simoptions);
    
    % Final eqm 1 (tau_k=0)
    Params.taxincome=2;
    Params.r=p_eqm_final1.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_final1.w;
    Params.tau_I=p_eqm_final1.tau_I;
    Params.tau_k=tau_k_final1;
    
    [V_final1,Policy_final1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final1=StationaryDist_InfHorz(Policy_final1,n_d,n_a,n_z,pi_z, simoptions);

    % Final eqm 2 (tau_k=0.4)
    Params.taxincome=2;
    Params.r=p_eqm_final2.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
    Params.w=p_eqm_final2.w;
    Params.tau_I=p_eqm_final2.tau_I;
    Params.tau_k=tau_k_final2;
    
    [V_final2,Policy_final2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final2=StationaryDist_InfHorz(Policy_final2,n_d,n_a,n_z,pi_z, simoptions);


    %% The following comment out section is how I realised that using 50 periods for transitions was insufficient, even with Policy_final1 directly T=50 was not enough to converge to final distribution, with T=100 it worked better
    % figure()
    % plot(asset_grid,cumsum(sum(sum(sum(StationaryDist_final1,4),3),2)), asset_grid,cumsum(sum(sum(sum(StationaryDist_init,4),3),2)))
    % legend('final','initial')
    % 
    % AgentDistPathTest=AgentDistOnTransPath_InfHorz(StationaryDist_init, PricePath, ParamPath, Policy_final1.*ones(1,1,1,1,1,T),n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
    % 
    % figure()
    % plot(asset_grid,cumsum(sum(sum(sum(StationaryDist_final1,4),3),2)), asset_grid,cumsum(sum(sum(sum(StationaryDist_init,4),3),2)), asset_grid,cumsum(sum(sum(sum(AgentDistPathTest(:,:,:,:,end),4),3),2)) )
    % legend('final','initial')

    %% Setup for the transition path to tau_k (to p_eqm_final)
    % For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
    % (already calculated V_final & StationaryDist_init above)
    
    % Number of time periods to allow for the transition (if you set T too low
    % it will cause problems, too high just means run-time will be longer).
    T=120 % Kitao (2008) graphs suggest she uses 50 periods, but this is not enough
    Tgraph=50;
    
    FnsToEvaluate_TransPath.K_noncorp=FnsToEvaluate.K_noncorp;
    FnsToEvaluate_TransPath.A=FnsToEvaluate.A;
    FnsToEvaluate_TransPath.N_noncorp=FnsToEvaluate.N_noncorp;
    FnsToEvaluate_TransPath.L=FnsToEvaluate.L;
    FnsToEvaluate_TransPath.TaxRevenue=FnsToEvaluate.TaxRevenue;
 
    % Add some things to FnsToEvaluate_TransPath that we are interested in but didn't need for solving the transition path
    FnsToEvaluate2_TransPath=FnsToEvaluate_TransPath;
    FnsToEvaluate2_TransPath.Y_noncorp =  @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k) Kitao2008_yFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k); % output of non-corporate sector
    FnsToEvaluate2_TransPath.C =  @(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)...
        Kitao2008_ConsumptionFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1); % consumption
  
    
    TransPathGeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,L,alpha,delta) r-(alpha*((A-K_noncorp)^(alpha-1))*((L-N_noncorp)^(1-alpha))-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,L,alpha) w-(1-alpha)*((A-K_noncorp)^(alpha))*((L-N_noncorp)^(-alpha)); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
    TransPathGeneralEqmEqns.GovBudget = @(G,TaxRevenue) G-TaxRevenue; %Government runs balanced budget
    % Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.
    
    transpathoptions.GEnewprice=2; % Anderson acceleration of the shooting update
    % Need to explain to transpathoptions how to use the GeneralEqmEqns to
    % update the general eqm transition prices (in PricePath).
    transpathoptions.GEnewprice2.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'CapitalMarket','r',0,0.1;... % CapitalMarket is positive is r is to large, so subtract
        'LaborMarket','w',0,0.1;... % LaborMarket is positive is r is to large, so subtract
        'GovBudget','tau_I',1,0.1}; % GovBudget is positive if tau_I is too small, so add
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.
    % Note: GEnewprice=2 is Anderson acceleration, which accelerates exactly this shooting update, so it
    % needs the same howtoupdate instructions (just under GEnewprice2 rather than GEnewprice3).
    
    % Now just run the TransitionPath_InfHorz command (all of the other inputs
    % are things we had already had to define to be able to solve for the initial and final equilibria)
    transpathoptions.verbose=1;

    transpathoptions.graphpricepath=1; % 1: creates a graph of the 'current' price path which updates each iteration.
    transpathoptions.graphaggvarspath=1; % 1: creates a graph of the 'current' aggregate variables which updates each iteration.
    if headlessFigures==1 % these two redraw every iteration, which is pointless offscreen and costs real time on a server
        transpathoptions.graphpricepath=0;
        transpathoptions.graphaggvarspath=0;
    end
    
    %% Compute the transition path to tau_k=0 (to p_eqm_final1)
    Params.taxincome=2;

    % We want to look at a one off unanticipated path of tau_k.
    ParamPath1.tau_k=tau_k_final1*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
    % (the way ParamPath is set is designed to allow for a series of changes in the parameters)
    
    % We need to give an initial guess for the price path
    PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final1.r, floor(T/3))'; p_eqm_final1.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
    PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final1.w, floor(T/3))'; p_eqm_final1.w*ones(T-floor(T/3),1)];
    PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final1.tau_I, floor(T/3))'; p_eqm_final1.tau_I*ones(T-floor(T/3),1)];
    
    [PricePath1,GEcondnPath1]=TransitionPath_InfHorz(PricePath0, ParamPath1, T, V_final1, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions,vfoptionstpath);
    
    [Vpath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePath1, ParamPath1, T, V_final1, Policy_final1, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
    AgentDistPath1=AgentDistOnTransPath_InfHorz(StationaryDist_init, PricePath1, ParamPath1, PolicyPath1,n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
    AggVarsPath1=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2_TransPath,AgentDistPath1,PolicyPath1,PricePath1,ParamPath1, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    
    % CEV calculations (CES utility fn makes this easy)
    CEV1=(Vpath1(:,:,:,:,1)./V_init).^(1/(1-Params.sigma))-1;
    % CEV can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove
    % them (otherwise NaN*0 poisons the mass-weighted averages used for Figure 6)
    CEV1((StationaryDist_init==0))=0;
       
    % Calculate the fraction that gain CEV
    gainCEV=sum(StationaryDist_init(CEV1>0)); % Note: must use the initial distribution, as the CEV is evaluated at the initial states
    Table7.Entrepreneurs.tau_k_final1.gainCEV=gainCEV;
    
    %% Compute the transition path to tau_k=0.4 (to p_eqm_final2)
    Params.taxincome=2;

    % We want to look at a one off unanticipated path of tau_k.
    ParamPath2.tau_k=tau_k_final2*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
    % (the way ParamPath is set is designed to allow for a series of changes in the parameters)
    
    % We need to give an initial guess for the price path
    PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final2.r, floor(T/3))'; p_eqm_final2.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
    PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final2.w, floor(T/3))'; p_eqm_final2.w*ones(T-floor(T/3),1)];
    PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final2.tau_I, floor(T/3))'; p_eqm_final2.tau_I*ones(T-floor(T/3),1)];
    
    [PricePath2,GEcondnPath2]=TransitionPath_InfHorz(PricePath0, ParamPath2, T, V_final2, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions,vfoptionstpath);
    
    [Vpath2,PolicyPath2]=ValueFnOnTransPath_InfHorz(PricePath2, ParamPath2, T, V_final2, Policy_final2, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
    AgentDistPath2=AgentDistOnTransPath_InfHorz(StationaryDist_init, PricePath2, ParamPath2, PolicyPath2,n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
    AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2_TransPath,AgentDistPath2,PolicyPath2,PricePath2,ParamPath2, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    
    % CEV calculations (CES utility fn makes this easy)
    CEV2=(Vpath2(:,:,:,:,1)./V_init).^(1/(1-Params.sigma))-1;
    % CEV can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove
    % them (otherwise NaN*0 poisons the mass-weighted averages used for Figure 6)
    CEV2((StationaryDist_init==0))=0;
        
    % Calculate the fraction that gain CEV
    gainCEV=sum(StationaryDist_init(CEV2>0)); % Note: must use the initial distribution, as the CEV is evaluated at the initial states
    Table7.Entrepreneurs.tau_k_final2.gainCEV=gainCEV;
    
    %% Calculate these same FnsToEvaluate2_TransPath for the initial eqm    
    AggVars_init=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_init, Policy_init, FnsToEvaluate2_TransPath,Params_init, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
    
    % Need to calculate the savings rate for Figure 5.
    % Savings is the change in assets, plus the depreciation.
    asset_path1=[AggVars_init.A.Mean, AggVarsPath1.A.Mean];
    dasset_path1=[0, asset_path1(2:end)-asset_path1(1:end-1)];
    savings_path1=Params.delta*asset_path1+dasset_path1;
    
    asset_path2=[AggVars_init.A.Mean, AggVarsPath2.A.Mean];
    dasset_path2=[0, asset_path2(2:end)-asset_path2(1:end-1)];
    savings_path2=Params.delta*asset_path2+dasset_path2;
    
    % Now need output, so that we can do savings rate
    A_path1=[AggVars_init.A.Mean, AggVarsPath1.A.Mean];
    K_noncorp_path1=[AggVars_init.K_noncorp.Mean, AggVarsPath1.K_noncorp.Mean];
    L_path1=[AggVars_init.L.Mean, AggVarsPath1.L.Mean];
    N_noncorp_path1=[AggVars_init.N_noncorp.Mean, AggVarsPath1.N_noncorp.Mean];
    Output_corp_path1=((A_path1-K_noncorp_path1).^Params.alpha).*((L_path1-N_noncorp_path1).^(1-Params.alpha));
    Y_path1=Output_corp_path1+[AggVars_init.Y_noncorp.Mean, AggVarsPath1.Y_noncorp.Mean];
    
    A_path2=[AggVars_init.A.Mean, AggVarsPath2.A.Mean];
    K_noncorp_path2=[AggVars_init.K_noncorp.Mean, AggVarsPath2.K_noncorp.Mean];
    L_path2=[AggVars_init.L.Mean, AggVarsPath2.L.Mean];
    N_noncorp_path2=[AggVars_init.N_noncorp.Mean, AggVarsPath2.N_noncorp.Mean];
    Output_corp_path2=((A_path2-K_noncorp_path2).^Params.alpha).*((L_path2-N_noncorp_path2).^(1-Params.alpha));
    Y_path2=Output_corp_path2+[AggVars_init.Y_noncorp.Mean,AggVarsPath2.Y_noncorp.Mean];
    
    % Savings rate is savings as fraction of income
    savingsrate_path1=savings_path1./Y_path1;
    savingsrate_path2=savings_path2./Y_path2;
    
    %% Figure 5
    fig5=figure(5);
    subplot(2,4,1); plot(0:1:T,[p_eqm_initial.w, PricePath1.w],'b-',0:1:T,[p_eqm_initial.w, PricePath2.w],'r-.')
    title('wage')
    subplot(2,4,2); plot(0:1:T,100*[p_eqm_initial.r, PricePath1.r],'b-',0:1:T,100*[p_eqm_initial.r, PricePath2.r],'r-.')
    title('interest rate (%)')
    subplot(2,4,3); plot(0:1:T,[AggVars_init.A.Mean, AggVarsPath1.A.Mean],'b-',0:1:T,[AggVars_init.A.Mean, AggVarsPath2.A.Mean],'r-.')
    title('agg. capital')
    subplot(2,4,4); plot(0:1:T,[AggVars_init.C.Mean, AggVarsPath1.C.Mean],'b-',0:1:T,[AggVars_init.C.Mean, AggVarsPath2.C.Mean],'r-.')
    title('agg. consumption')
    subplot(2,4,5); plot(0:1:T,100*savingsrate_path1,'b-',0:1:T,100*savingsrate_path2,'r-.')
    title('savings rate')
    subplot(2,4,6); plot(0:1:T,100*[p_eqm_initial.tau_I, PricePath1.tau_I],'b-',0:1:T,100*[p_eqm_initial.tau_I, PricePath2.tau_I],'r-.')
    title('proportional tax \tau_I (%)')
    subplot(2,4,7); plot(0:1:T,[AggVars_init.K_noncorp.Mean, AggVarsPath1.K_noncorp.Mean],'b-',0:1:T,[AggVars_init.K_noncorp.Mean, AggVarsPath2.K_noncorp.Mean],'r-.')
    title('entrep. capital')
    subplot(2,4,8); plot(0:1:T,[AggVars_init.Y_noncorp.Mean, AggVarsPath1.Y_noncorp.Mean],'b-',0:1:T,[AggVars_init.Y_noncorp.Mean, AggVarsPath2.Y_noncorp.Mean],'r-.')
    title('entrep. output')
    saveas(fig5,'./SavedOutput/Graphs/Kitao2008_Fig5.png')
    
    
    %% Figure 6
    
    fig6=figure(6);
    % tau_k=0
    workermass=StationaryDist_init(:,1,:,:); % 1 is workers
    workerCEV1=sum(CEV1(:,1,:,:).*workermass,4)./sum(workermass,4); % average over theta, weighted by the initial distribution (Kitao (2008), footnote 15)
    % Note: this is NaN wherever there is no mass, which is what we want, as those points then simply do not get plotted
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,1),'b--')
    hold on
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,3),'g.')
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,5),'r-')
    hold off
    legend('eta 1','eta 3', 'eta 5')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('\tau_k=0: workers')
    
    entrepreneurmass=StationaryDist_init(:,2,:,:); % 2 is entrepreneurs
    entrepreneurCEV1=sum(CEV1(:,2,:,:).*entrepreneurmass,3)./sum(entrepreneurmass,3); % average over eta, weighted by the initial distribution (Kitao (2008), footnote 15)
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,2),'b--')
    hold on
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,3),'g.')
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,4),'r-')
    hold off
    legend('theta 2','theta 3', 'theta 4')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('\tau_k=0: entrepreneurs')
    
    % tau_k=0.4
    workermass=StationaryDist_init(:,1,:,:); % 1 is workers
    workerCEV2=sum(CEV2(:,1,:,:).*workermass,4)./sum(workermass,4); % average over theta, weighted by the initial distribution (Kitao (2008), footnote 15)
    % Note: this is NaN wherever there is no mass, which is what we want, as those points then simply do not get plotted
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,1),'b--')
    hold on
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,3),'g.')
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,5),'r-')
    hold off
    legend('eta 1','eta 3', 'eta 5')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('\tau_k=0.4: workers')
    
    entrepreneurmass=StationaryDist_init(:,2,:,:); % 2 is entrepreneurs
    entrepreneurCEV2=sum(CEV2(:,2,:,:).*entrepreneurmass,3)./sum(entrepreneurmass,3); % average over eta, weighted by the initial distribution (Kitao (2008), footnote 15)
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,2),'b--')
    hold on
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,3),'g.')
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,4),'r-')
    hold off
    legend('theta 2','theta 3', 'theta 4')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('\tau_k=0.4: entrepreneurs')
    
    saveas(fig6,'./SavedOutput/Graphs/Kitao2008_Fig6.png')
    
    
    %% Table 6
    % CEV again, but report by wealth levels
    wealthcutoffs=[10000,50000,100000,250000,500000,1000000,2000000]; % wealth cutoffs in dollars
    wealthcutoffs=wealthcutoffs/wealthscalingfactor; % wealth cutoffs in model units
    
    % First, deal with
    Table6data_mass=zeros(9,3); % Mass of agents in each wealth band
    % Stationary dist for the tau_k=0 eqm (CEV calculations all made from initial distribution)
    StationaryDist_temp=StationaryDist_init;
    % We don't want eta or theta for any of table 2
    StationaryDist_temp=sum(sum(StationaryDist_temp,4),3); % This is why I am calling it _temp, so I don't accidently use it elsewhere
    % Worker (1 in second dim of StationaryDist)
    Table6data_mass(1,1)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),1));
    Table6data_mass(2,1)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),1));
    Table6data_mass(3,1)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),1));
    Table6data_mass(4,1)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),1));
    Table6data_mass(5,1)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),1));
    Table6data_mass(6,1)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),1));
    Table6data_mass(7,1)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),1));
    Table6data_mass(8,1)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),1));
    Table6data_mass(9,1)=sum(StationaryDist_temp(:,1));
    % Entrepreneur (2 in second dim of StationaryDist)
    Table6data_mass(1,2)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),2));
    Table6data_mass(2,2)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),2));
    Table6data_mass(3,2)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),2));
    Table6data_mass(4,2)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),2));
    Table6data_mass(5,2)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),2));
    Table6data_mass(6,2)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),2));
    Table6data_mass(7,2)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),2));
    Table6data_mass(8,2)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),2));
    Table6data_mass(9,2)=sum(StationaryDist_temp(:,2));
    % All (sum over second dimension)
    Table6data_mass(1,3)=sum(sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),:),2));
    Table6data_mass(2,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),:),2));
    Table6data_mass(3,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),:),2));
    Table6data_mass(4,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),:),2));
    Table6data_mass(5,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),:),2));
    Table6data_mass(6,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),:),2));
    Table6data_mass(7,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),:),2));
    Table6data_mass(8,3)=sum(sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),:),2));
    Table6data_mass(9,3)=sum(sum(StationaryDist_temp(:,:),2));
    
    % Repeat, but this time only the mass for agents with positive CEV
    Table6data_tauk0_gainmass=zeros(9,3); % Mass of agents in each wealth band with positive CEV
    % Stationary dist for the tau_k=0 eqm
    StationaryDist_temp=StationaryDist_init;
    CEV_tauk01=CEV1;
    gainCEV_tauk1=(CEV_tauk01>0);
    % So the distribution mass of those who gained from the reform is
    StationaryDist_temp=StationaryDist_temp.*gainCEV_tauk1;
    % Now get the mass of those who gained
    StationaryDist_temp=sum(sum(StationaryDist_temp,4),3);
    % Worker (1 in second dim of StationaryDist)
    Table6data_tauk0_gainmass(1,1)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),1));
    Table6data_tauk0_gainmass(2,1)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),1));
    Table6data_tauk0_gainmass(3,1)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),1));
    Table6data_tauk0_gainmass(4,1)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),1));
    Table6data_tauk0_gainmass(5,1)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),1));
    Table6data_tauk0_gainmass(6,1)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),1));
    Table6data_tauk0_gainmass(7,1)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),1));
    Table6data_tauk0_gainmass(8,1)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),1));
    Table6data_tauk0_gainmass(9,1)=sum(StationaryDist_temp(:,1));
    % Entrepreneur (2 in second dim of StationaryDist)
    Table6data_tauk0_gainmass(1,2)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),2));
    Table6data_tauk0_gainmass(2,2)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),2));
    Table6data_tauk0_gainmass(3,2)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),2));
    Table6data_tauk0_gainmass(4,2)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),2));
    Table6data_tauk0_gainmass(5,2)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),2));
    Table6data_tauk0_gainmass(6,2)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),2));
    Table6data_tauk0_gainmass(7,2)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),2));
    Table6data_tauk0_gainmass(8,2)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),2));
    Table6data_tauk0_gainmass(9,2)=sum(StationaryDist_temp(:,2));
    % All (sum over second dimension)
    Table6data_tauk0_gainmass(1,3)=sum(sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),:),2));
    Table6data_tauk0_gainmass(2,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),:),2));
    Table6data_tauk0_gainmass(3,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),:),2));
    Table6data_tauk0_gainmass(4,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),:),2));
    Table6data_tauk0_gainmass(5,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),:),2));
    Table6data_tauk0_gainmass(6,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),:),2));
    Table6data_tauk0_gainmass(7,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),:),2));
    Table6data_tauk0_gainmass(8,3)=sum(sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),:),2));
    Table6data_tauk0_gainmass(9,3)=sum(sum(StationaryDist_temp(:,:),2));
    
    
    % Redo all the calculations for tau_k=0.4
    % Repeat, but this time only the mass for agents with positive CEV
    Table6data_tauk04_gainmass=zeros(9,3); % Mass of agents in each wealth band with positive CEV
    % Stationary dist for the tau_k=0.4 eqm
    StationaryDist_temp=StationaryDist_init;
    CEV_tauk04=CEV2;
    gainCEV_tauk04=(CEV_tauk04>0);
    % So the distribution mass of those who gained from the reform is
    StationaryDist_temp=StationaryDist_temp.*gainCEV_tauk04;
    % Now get the mass of those who gained
    StationaryDist_temp=sum(sum(StationaryDist_temp,4),3);
    % Worker (1 in second dim of StationaryDist)
    Table6data_tauk04_gainmass(1,1)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),1));
    Table6data_tauk04_gainmass(2,1)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),1));
    Table6data_tauk04_gainmass(3,1)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),1));
    Table6data_tauk04_gainmass(4,1)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),1));
    Table6data_tauk04_gainmass(5,1)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),1));
    Table6data_tauk04_gainmass(6,1)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),1));
    Table6data_tauk04_gainmass(7,1)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),1));
    Table6data_tauk04_gainmass(8,1)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),1));
    Table6data_tauk04_gainmass(9,1)=sum(StationaryDist_temp(:,1));
    % Entrepreneur (2 in second dim of StationaryDist)
    Table6data_tauk04_gainmass(1,2)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),2));
    Table6data_tauk04_gainmass(2,2)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),2));
    Table6data_tauk04_gainmass(3,2)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),2));
    Table6data_tauk04_gainmass(4,2)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),2));
    Table6data_tauk04_gainmass(5,2)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),2));
    Table6data_tauk04_gainmass(6,2)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),2));
    Table6data_tauk04_gainmass(7,2)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),2));
    Table6data_tauk04_gainmass(8,2)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),2));
    Table6data_tauk04_gainmass(9,2)=sum(StationaryDist_temp(:,2));
    % All (sum over second dimension)
    Table6data_tauk04_gainmass(1,3)=sum(sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),:),2));
    Table6data_tauk04_gainmass(2,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),:),2));
    Table6data_tauk04_gainmass(3,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),:),2));
    Table6data_tauk04_gainmass(4,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),:),2));
    Table6data_tauk04_gainmass(5,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),:),2));
    Table6data_tauk04_gainmass(6,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),:),2));
    Table6data_tauk04_gainmass(7,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),:),2));
    Table6data_tauk04_gainmass(8,3)=sum(sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),:),2));
    Table6data_tauk04_gainmass(9,3)=sum(sum(StationaryDist_temp(:,:),2));
    
    
    
    
    
    %Table 6
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table6.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n');
    fprintf(FID, '\\multicolumn{7}{l}{Fraction of households with welfare gains: capital income tax } \\\\  \\hline \n');
    fprintf(FID, ' Wealth (in \\$1,000) & Workers & & Entrepreneurs & & All &  \\\\ \\hline \n');
    fprintf(FID, '\\multicolumn{7}{c}{(a) Capital income tax $\\tau_k =0$\\%%}    \\\\ \n');
    fprintf(FID, '0-10       & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(1,1),100*Table6data_tauk0_gainmass(1,1)/Table6data_mass(1,1), 100*Table6data_tauk0_gainmass(1,2),100*Table6data_tauk0_gainmass(1,2)/Table6data_mass(1,2), 100*Table6data_tauk0_gainmass(1,3),100*Table6data_tauk0_gainmass(1,3)/Table6data_mass(1,3));
    fprintf(FID, '10-50      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(2,1),100*Table6data_tauk0_gainmass(2,1)/Table6data_mass(2,1), 100*Table6data_tauk0_gainmass(2,2),100*Table6data_tauk0_gainmass(2,2)/Table6data_mass(2,2), 100*Table6data_tauk0_gainmass(2,3),100*Table6data_tauk0_gainmass(2,3)/Table6data_mass(2,3));
    fprintf(FID, '50-100     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(3,1),100*Table6data_tauk0_gainmass(3,1)/Table6data_mass(3,1), 100*Table6data_tauk0_gainmass(3,2),100*Table6data_tauk0_gainmass(3,2)/Table6data_mass(3,2), 100*Table6data_tauk0_gainmass(3,3),100*Table6data_tauk0_gainmass(3,3)/Table6data_mass(3,3));
    fprintf(FID, '100-250    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(4,1),100*Table6data_tauk0_gainmass(4,1)/Table6data_mass(4,1), 100*Table6data_tauk0_gainmass(4,2),100*Table6data_tauk0_gainmass(4,2)/Table6data_mass(4,2), 100*Table6data_tauk0_gainmass(4,3),100*Table6data_tauk0_gainmass(4,3)/Table6data_mass(4,3));
    fprintf(FID, '250-500    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(5,1),100*Table6data_tauk0_gainmass(5,1)/Table6data_mass(5,1), 100*Table6data_tauk0_gainmass(5,2),100*Table6data_tauk0_gainmass(5,2)/Table6data_mass(5,2), 100*Table6data_tauk0_gainmass(5,3),100*Table6data_tauk0_gainmass(5,3)/Table6data_mass(5,3));
    fprintf(FID, '500-1000   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(6,1),100*Table6data_tauk0_gainmass(6,1)/Table6data_mass(6,1), 100*Table6data_tauk0_gainmass(6,2),100*Table6data_tauk0_gainmass(6,2)/Table6data_mass(6,2), 100*Table6data_tauk0_gainmass(6,3),100*Table6data_tauk0_gainmass(6,3)/Table6data_mass(6,3));
    fprintf(FID, '1000-2000  & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(7,1),100*Table6data_tauk0_gainmass(7,1)/Table6data_mass(7,1), 100*Table6data_tauk0_gainmass(7,2),100*Table6data_tauk0_gainmass(7,2)/Table6data_mass(7,2), 100*Table6data_tauk0_gainmass(7,3),100*Table6data_tauk0_gainmass(7,3)/Table6data_mass(7,3));
    fprintf(FID, '$>$2000      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(8,1),100*Table6data_tauk0_gainmass(8,1)/Table6data_mass(8,1), 100*Table6data_tauk0_gainmass(8,2),100*Table6data_tauk0_gainmass(8,2)/Table6data_mass(8,2), 100*Table6data_tauk0_gainmass(8,3),100*Table6data_tauk0_gainmass(8,3)/Table6data_mass(8,3));
    fprintf(FID, 'all        & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk0_gainmass(9,1),100*Table6data_tauk0_gainmass(9,1)/Table6data_mass(9,1), 100*Table6data_tauk0_gainmass(9,2),100*Table6data_tauk0_gainmass(9,2)/Table6data_mass(9,2), 100*Table6data_tauk0_gainmass(9,3),100*Table6data_tauk0_gainmass(9,3)/Table6data_mass(9,3));
    fprintf(FID, '\\multicolumn{7}{c}{ }    \\\\ \n');
    fprintf(FID, '\\multicolumn{7}{c}{(b) Capital income tax $\\tau_k =40$\\%%}    \\\\ \n');
    fprintf(FID, '0-10       & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(1,1),100*Table6data_tauk04_gainmass(1,1)/Table6data_mass(1,1), 100*Table6data_tauk04_gainmass(1,2),100*Table6data_tauk04_gainmass(1,2)/Table6data_mass(1,2), 100*Table6data_tauk04_gainmass(1,3),100*Table6data_tauk04_gainmass(1,3)/Table6data_mass(1,3));
    fprintf(FID, '10-50      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(2,1),100*Table6data_tauk04_gainmass(2,1)/Table6data_mass(2,1), 100*Table6data_tauk04_gainmass(2,2),100*Table6data_tauk04_gainmass(2,2)/Table6data_mass(2,2), 100*Table6data_tauk04_gainmass(2,3),100*Table6data_tauk04_gainmass(2,3)/Table6data_mass(2,3));
    fprintf(FID, '50-100     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(3,1),100*Table6data_tauk04_gainmass(3,1)/Table6data_mass(3,1), 100*Table6data_tauk04_gainmass(3,2),100*Table6data_tauk04_gainmass(3,2)/Table6data_mass(3,2), 100*Table6data_tauk04_gainmass(3,3),100*Table6data_tauk04_gainmass(3,3)/Table6data_mass(3,3));
    fprintf(FID, '100-250    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(4,1),100*Table6data_tauk04_gainmass(4,1)/Table6data_mass(4,1), 100*Table6data_tauk04_gainmass(4,2),100*Table6data_tauk04_gainmass(4,2)/Table6data_mass(4,2), 100*Table6data_tauk04_gainmass(4,3),100*Table6data_tauk04_gainmass(4,3)/Table6data_mass(4,3));
    fprintf(FID, '250-500    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(5,1),100*Table6data_tauk04_gainmass(5,1)/Table6data_mass(5,1), 100*Table6data_tauk04_gainmass(5,2),100*Table6data_tauk04_gainmass(5,2)/Table6data_mass(5,2), 100*Table6data_tauk04_gainmass(5,3),100*Table6data_tauk04_gainmass(5,3)/Table6data_mass(5,3));
    fprintf(FID, '500-1000   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(6,1),100*Table6data_tauk04_gainmass(6,1)/Table6data_mass(6,1), 100*Table6data_tauk04_gainmass(6,2),100*Table6data_tauk04_gainmass(6,2)/Table6data_mass(6,2), 100*Table6data_tauk04_gainmass(6,3),100*Table6data_tauk04_gainmass(6,3)/Table6data_mass(6,3));
    fprintf(FID, '1000-2000  & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(7,1),100*Table6data_tauk04_gainmass(7,1)/Table6data_mass(7,1), 100*Table6data_tauk04_gainmass(7,2),100*Table6data_tauk04_gainmass(7,2)/Table6data_mass(7,2), 100*Table6data_tauk04_gainmass(7,3),100*Table6data_tauk04_gainmass(7,3)/Table6data_mass(7,3));
    fprintf(FID, '$>$2000      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(8,1),100*Table6data_tauk04_gainmass(8,1)/Table6data_mass(8,1), 100*Table6data_tauk04_gainmass(8,2),100*Table6data_tauk04_gainmass(8,2)/Table6data_mass(8,2), 100*Table6data_tauk04_gainmass(8,3),100*Table6data_tauk04_gainmass(8,3)/Table6data_mass(8,3));
    fprintf(FID, 'all        & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table6data_tauk04_gainmass(9,1),100*Table6data_tauk04_gainmass(9,1)/Table6data_mass(9,1), 100*Table6data_tauk04_gainmass(9,2),100*Table6data_tauk04_gainmass(9,2)/Table6data_mass(9,2), 100*Table6data_tauk04_gainmass(9,3),100*Table6data_tauk04_gainmass(9,3)/Table6data_mass(9,3));
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: In parentheses are the fractions of households with welfare gains conditional on their occupation and their wealth category. \n');
    fprintf(FID, 'A NaN in parentheses indicates a wealth category in which the occupation in question has zero mass, so the conditional fraction is undefined. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);

    %% Check that we don't hit the top of asset grid (Fig15 is not one from Kitao 2008, just something I want to see)
    Fig15=figure(15);
    assetdist_init=sum(sum(StationaryDist_init,4),3);
    assetdist_init=cumsum(assetdist_init,1);
    assetdist_final1=sum(sum(StationaryDist_final1,4),3);
    assetdist_final1=cumsum(assetdist_final1,1);
    assetdist_final2=sum(sum(StationaryDist_final2,4),3);
    assetdist_final2=cumsum(assetdist_final2,1);
    plot(asset_grid,assetdist_init,asset_grid,assetdist_final1,asset_grid,assetdist_final2)
    title('cdf of HHs over assets')
    legend('init','final tauk=0','final tauk=0.4')
    % along with this, just calculate the mass in the top 10 and top 5 grid points of asset grid
    assetdist_initB=sum(assetdist_init,2);
    assetdist_final1B=sum(assetdist_final1,2);
    assetdist_final2B=sum(assetdist_final2,2);
    temp=[1-assetdist_initB(end-9), 1-assetdist_initB(end-4);...
        1-assetdist_final1B(end-9), 1-assetdist_final1B(end-4);...
        1-assetdist_final2B(end-9), 1-assetdist_final2B(end-4)];
    fprintf('For initial stationary general eqm, mass in top 10 grid points is %8.6f, and mass in top 5 grid points is %8.6f \n' , temp(1,1), temp(1,2) );
    fprintf('For tauk=0 stationary general eqm, mass in top 10 grid points is %8.6f, and mass in top 5 grid points is %8.6f \n' , temp(2,1), temp(2,2) );
    fprintf('For tauk=40 stationary general eqm, mass in top 10 grid points is %8.6f, and mass in top 5 grid points is %8.6f \n' , temp(3,1), temp(3,2) );
    % In all three cases it is of the order of 1/1000th of the mass in these top few points. Seems like asset_grid has a suitable max (high enough to not bind, but not so high as to waste heaps of points)



    % As for doPart(1) above, plus what this block adds: T, V_init, StationaryDist_init, AggVars_init, Table6data_mass, FnsToEvaluate_TransPath, FnsToEvaluate2_TransPath, TransPathGeneralEqmEqns, transpathoptions
    save ./SavedOutput/Kitao2008_doPart3.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8 Params_init p_eqm_final1 p_eqm_final2 tau_k_vec keepVPolicy_tauk Figure4data T V_init StationaryDist_init AggVars_init Table6data_mass FnsToEvaluate_TransPath FnsToEvaluate2_TransPath TransPathGeneralEqmEqns transpathoptions
else
    load ./SavedOutput/Kitao2008_doPart3.mat
end



%% Figure 7
% For figure 7 we need nine transition paths (we have already solved the first and ninth, but will redo them anyway)
% Note: tau_k_transition_vec is now set up at the top, alongside tau_E1_transition_vec
% Note: we already solved the final stationary general eqm for all of these (and many more)

if doPart(4)==1
    Params.taxincome=2;
    
    Figure7checks=struct();
    Figure7data=zeros(length(tau_k_transition_vec),2);
    for ii=1:length(tau_k_transition_vec)
        ii

        tau_k=tau_k_transition_vec(ii);
        % Note: matched on nearest value rather than on exact equality. tau_k_vec is built by the colon
    % operator, so its element for 35% is 35*0.01, while tau_k_transition_vec has 35/100, and those
    % two are different doubles (by 5.6e-17). With == the match failed and max() of an all-false
    % vector returns index 1, so the tau=0 equilibrium was silently used instead.
    [~,jj]=min(abs(tau_k_vec-tau_k)); % jj is the index from when the stationary general eqm was computed earlier
        V_final=Figure4data(jj).V;
        Policy_final=Figure4data(jj).Policy;
        p_eqm_final=Figure4data(jj).p_eqm_final;

        % We want to look at a one off unanticipated path of tau_k.
        ParamPath.tau_k=tau_k*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
        % (the way ParamPath is set is designed to allow for a series of changes in the parameters)
        
        % We need to give an initial guess for the price path
        PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
        PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final.w, floor(T/3))'; p_eqm_final.w*ones(T-floor(T/3),1)];
        PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final.tau_I, floor(T/3))'; p_eqm_final.tau_I*ones(T-floor(T/3),1)];

        [PricePath,GEcondnPath]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions,vfoptionstpath);
        
        [Vpath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
        
        % CEV calculations (CES utility fn makes this easy)
        CEVpath=(Vpath(:,:,:,:,1)./V_init).^(1/(1-Params.sigma))-1;
        % CEVpath can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove them
        CEVpath((StationaryDist_init==0))=0;
        CEVpath_behindtheveil=sum(sum(sum(sum(CEVpath.*StationaryDist_init))));
        
        % CEV calculations for steady-state comparison
        CEV=(V_final./V_init).^(1/(1-Params.sigma))-1;
        % CEV can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove them
        CEV((StationaryDist_init==0))=0;
        CEV_behindtheveil=sum(sum(sum(sum(CEV.*StationaryDist_init))));
        
        Figure7data(ii,1)=CEVpath_behindtheveil; % transition
        Figure7data(ii,2)=CEV_behindtheveil; % compare initial and final (ignoring transtion path)
        
        % Keep so stuff so can double check things
        Figure7checks(ii).PricePath=PricePath;
        Figure7checks(ii).GEcondnPath=GEcondnPath;
        % Note: Vpath used to be kept here too. It is nine full value function paths, about 350MB, and
        % nothing ever read it back. CEVpath below is what is actually wanted out of it.
    end
    
    % Figure 7
    fig7=figure(7);
    subplot(1,2,1); plot(100*tau_k_transition_vec,100*Figure7data(:,1),'b-o')
    title('transition')
    ylabel('percentage (%)')
    subplot(1,2,2); plot(100*tau_k_transition_vec,100*Figure7data(:,2),'b-o')
    title('steady state')
    ylabel('percentage (%)')
    saveas(fig7,'./SavedOutput/Graphs/Kitao2008_Fig7.png')

    % As for doPart(1) above, plus what this block adds: Figure7data
    save ./SavedOutput/Kitao2008_doPart4.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8 Params_init p_eqm_final1 p_eqm_final2 tau_k_vec keepVPolicy_tauk Figure4data T V_init StationaryDist_init AggVars_init Table6data_mass FnsToEvaluate_TransPath FnsToEvaluate2_TransPath TransPathGeneralEqmEqns transpathoptions Figure7data
else
    load ./SavedOutput/Kitao2008_doPart4.mat
end


%% Table 7

if doPart(5)==1
    % Table7 is computed as part of doPart(3)
    
    %Table 7
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table7.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n');
    fprintf(FID, '\\multicolumn{7}{l}{Flat capital income tax: economy without entrepreneurs} \\\\  \\hline \n');
    fprintf(FID, ' & $\\tau_k=0$\\%% & & & $\\tau_k=40$\\%% & & \\\\ \n');
    fprintf(FID, ' & bnch. & no E(1) & no E(2) & bnch. & no E(1) & no E(2) \\\\ \\hline \n');
    fprintf(FID, 'interest rate                                & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',  Table7.Entrepreneurs.tau_k_final1.r/Table7.Entrepreneurs.initial.r, Table7.NoEntrepreneurs1.tau_k_final1.r/Table7.NoEntrepreneurs1.initial.r, Table7.NoEntrepreneurs2.tau_k_final1.r/Table7.NoEntrepreneurs2.initial.r, Table7.Entrepreneurs.tau_k_final2.r/Table7.Entrepreneurs.initial.r, Table7.NoEntrepreneurs1.tau_k_final2.r/Table7.NoEntrepreneurs1.initial.r, Table7.NoEntrepreneurs2.tau_k_final2.r/Table7.NoEntrepreneurs2.initial.r);
    fprintf(FID, 'wage                                         & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',  Table7.Entrepreneurs.tau_k_final1.w/Table7.Entrepreneurs.initial.w, Table7.NoEntrepreneurs1.tau_k_final1.w/Table7.NoEntrepreneurs1.initial.w, Table7.NoEntrepreneurs2.tau_k_final1.w/Table7.NoEntrepreneurs2.initial.w, Table7.Entrepreneurs.tau_k_final2.w/Table7.Entrepreneurs.initial.w, Table7.NoEntrepreneurs1.tau_k_final2.w/Table7.NoEntrepreneurs1.initial.w, Table7.NoEntrepreneurs2.tau_k_final2.w/Table7.NoEntrepreneurs2.initial.w);
    fprintf(FID, 'agg. capital                                 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',  Table7.Entrepreneurs.tau_k_final1.A/Table7.Entrepreneurs.initial.A, Table7.NoEntrepreneurs1.tau_k_final1.A/Table7.NoEntrepreneurs1.initial.A, Table7.NoEntrepreneurs2.tau_k_final1.A/Table7.NoEntrepreneurs2.initial.A, Table7.Entrepreneurs.tau_k_final2.A/Table7.Entrepreneurs.initial.A, Table7.NoEntrepreneurs1.tau_k_final2.A/Table7.NoEntrepreneurs1.initial.A, Table7.NoEntrepreneurs2.tau_k_final2.A/Table7.NoEntrepreneurs2.initial.A);
    fprintf(FID, 'agg. output                                  & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',  Table7.Entrepreneurs.tau_k_final1.Y/Table7.Entrepreneurs.initial.Y, Table7.NoEntrepreneurs1.tau_k_final1.Y/Table7.NoEntrepreneurs1.initial.Y, Table7.NoEntrepreneurs2.tau_k_final1.Y/Table7.NoEntrepreneurs2.initial.Y, Table7.Entrepreneurs.tau_k_final2.Y/Table7.Entrepreneurs.initial.Y, Table7.NoEntrepreneurs1.tau_k_final2.Y/Table7.NoEntrepreneurs1.initial.Y, Table7.NoEntrepreneurs2.tau_k_final2.Y/Table7.NoEntrepreneurs2.initial.Y);
    fprintf(FID, 'proportional tax $\\tau_I$ (change in \\%%)  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',  100*(Table7.Entrepreneurs.tau_k_final1.tau_I-Table7.Entrepreneurs.initial.tau_I), 100*(Table7.NoEntrepreneurs1.tau_k_final1.tau_I-Table7.NoEntrepreneurs1.initial.tau_I), 100*(Table7.NoEntrepreneurs2.tau_k_final1.tau_I-Table7.NoEntrepreneurs2.initial.tau_I), 100*(Table7.Entrepreneurs.tau_k_final2.tau_I-Table7.Entrepreneurs.initial.tau_I), 100*(Table7.NoEntrepreneurs1.tau_k_final2.tau_I-Table7.NoEntrepreneurs1.initial.tau_I), 100*(Table7.NoEntrepreneurs2.tau_k_final2.tau_I-Table7.NoEntrepreneurs2.initial.tau_I));
    fprintf(FID, '\\%% with CEV>0                                 & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  100*Table7.Entrepreneurs.tau_k_final1.gainCEV, 100*Table7.NoEntrepreneurs1.tau_k_final1.gainCEV, 100*Table7.NoEntrepreneurs2.tau_k_final1.gainCEV, 100*Table7.Entrepreneurs.tau_k_final2.gainCEV, 100*Table7.NoEntrepreneurs1.tau_k_final2.gainCEV, 100*Table7.NoEntrepreneurs2.tau_k_final2.gainCEV);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: The variables are normalized by their benchmark level, except for the proportional tax rate $\\tau_I$ and the percentage of the population with a positive CEV. In the initial steady state, $\\tau_I=$%8.2f, %8.2f, and %8.2f\\%% in these economies. \n',100*Table7.Entrepreneurs.initial.tau_I,100*Table7.NoEntrepreneurs1.initial.tau_I,100*Table7.NoEntrepreneurs2.initial.tau_I);
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
end


%% Now to do the 'tax entrepreneurial business income'

Params.taxincome=3; % 1 is tax income, 2 is to tax capital income and labor income seperately, 3 is tax entrepreneurial business income
Params.tau_k=0; % Note: Params.tau_k is now anyway irrelevant (as taxincome=3)
Params.tau_E1=0; % Initially

GEPriceParamNames={'r','w','tau_I'}; % tau_I instead of G
% Note: the general eqm conditions are unchanged, is just that the
% government budget balance is now about finding tau_I, where previously it
% was about finding G

% For some reason the general eqm optimization insists on trying out negative r and w, so going to stop it doing that
% heteroagentoptions.constrainpositive={'r','w'};
% Turned this off as the codes kept getting stuck with r being zero and failing for find the general eqm, not sure why.

if doPart(6)==1
    Params.taxincome=3; % 1 is tax income, 2 is to tax capital income and labor income seperately, 3 is tax entrepreneurial business income

    % ---- keep the general eqm conditions real-valued in this experiment ----
    % tau_E1=0 leaves entrepreneurial business income untaxed, and with the corrected static problem
    % the entrepreneurs expand until they are using essentially all of the labour force: at the first
    % iterate N_noncorp=0.9396 against L=0.9998, so 94%% of all labour is already in the non-corporate
    % sector. One finite-difference probe on r (the toolkit uses FiniteDifferenceStepSize=1e-2, and
    % dN_noncorp/dr is about +6.2 here, so a step of 0.01 is almost exactly the remaining slack) then
    % takes N_noncorp past L, and (L-N_noncorp)^(-alpha) with a NEGATIVE base and a non-integer
    % exponent is COMPLEX. The diary hides this, because fprintf('%%8.4f',...) prints only the real
    % part: the reported LaborMarket residual of -4.735969 was really -4.8017+13.1634i.
    % lsqnonlin then gets a complex residual, takes a complex step, and hands complex r and w to the
    % ReturnFn, whereupon (upsilon2/w)^(upsilon2/onegg) is complex and gpuArray/arrayfun refuses to
    % compile with "Variable 'k' changed type".
    % Two changes, only needed for the tau_E1 experiments:
    %   the max(...) floors keep the conditions real everywhere, and
    %   LaborPenalty restores a gradient in the region where the floor has clamped them (without it
    %   the clamped conditions are flat in N_noncorp, so the Jacobian would see no way back out).
    % Note this cannot be done with parameter constraints instead: the prices at the failing probe
    % (r=0.0404, w=1.3925) are perfectly ordinary and positive. The forbidden region is a diagonal in
    % (r,w) whose position is only known after solving the model, and it moves with tau_E1.
    GeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,L,alpha,delta) r-(alpha*(max(A-K_noncorp,10^(-6))^(alpha-1))*(max(L-N_noncorp,10^(-2))^(1-alpha))-delta);
    GeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,L,alpha) w-(1-alpha)*(max(A-K_noncorp,10^(-6))^(alpha))*(max(L-N_noncorp,10^(-2))^(-alpha));
    GeneralEqmEqns.LaborPenalty = @(L,N_noncorp) 100*max(N_noncorp-L+0.03,0); % zero unless the corporate sector is down to under 3%% of labour
    % Note: a fourth general eqm eqn against three prices is fine here. fminalgo=8 is lsqnonlin,
    % which is a least-squares solve and so is happy over-identified, and fminalgo=1 (fminsearch)
    % works on the weighted sum of squares. Neither uses the howtoupdate mapping.

    % Kitao (2008) does tau_k from 0 to 40 in intervals of 5. Here do intervals of 1
    Figure8data=struct();
    tau_E1_vec=0:0.01:0.4;
    % As above: which of these will doPart(7) later want V and Policy for
    keepVPolicy_tauE1=false(1,length(tau_E1_vec));
    for kk=1:length(tau_E1_transition_vec)
        [~,jj]=min(abs(tau_E1_vec-tau_E1_transition_vec(kk))); % nearest value, not exact equality (see note at the doPart(4) lookup)
        keepVPolicy_tauE1(jj)=true;
    end
    if sum(keepVPolicy_tauE1)~=length(tau_E1_transition_vec)
        error('tau_E1_transition_vec does not line up with tau_E1_vec, so V and Policy would not be kept for every transition')
    end

    % Print the GE prices being tried BEFORE the value fn iteration rather than after it. With
    % heteroagentoptions.verbose=1 the toolkit reports prices at the end of the objective function
    % (HeteroAgentStationaryEqm_InfHorz_subfn line 109), so when the value fn iteration itself fails
    % the prices that caused it never reach the diary. With verbose=2 they are reported first
    % (subfn line 6). The signature to look for in the diary is a "Current GE prices" block with no
    % "Current aggregate variables" block after it: those are the prices that broke it.
    heteroagentoptions.verbose=2;

    % initial guess for first run
    Params.r=p_eqm_initial.r;
    Params.w=p_eqm_initial.w;
    Params.tau_I=p_eqm_initial.tau_I;
    

    for ii=1:length(tau_E1_vec)
        Params.tau_E1=tau_E1_vec(ii);
        % As in doPart(2): report which experiment we are in before the solve, not after it
        fprintf('doPart(6): iteration %d of %d, tau_E1=%.2f \n', ii, length(tau_E1_vec), Params.tau_E1)
        
        [p_eqm_final,GeneralEqmCondn]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
        
        ii
        p_eqm_final % The equilibrium values of the GE prices
        
        % Now that we have the GE, let's calculate a bunch of related objects
        Params.r=p_eqm_final.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters
        Params.w=p_eqm_final.w;
        Params.tau_I=p_eqm_final.tau_I;
        
        [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
        
        StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);
        
        AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
        
        Output_corp=((AllStats.A.Mean-AllStats.K_noncorp.Mean)^Params.alpha)*((AllStats.L.Mean-AllStats.N_noncorp.Mean)^(1-Params.alpha));
        Y=Output_corp+AllStats.Y_noncorp.Mean;
        [Y,AllStats.Y_noncorp.Mean, Output_corp]
        
        Figure8data(ii).w=p_eqm_final.w;
        Figure8data(ii).r=p_eqm_final.r;
        Figure8data(ii).tau_I=p_eqm_final.tau_I;
        Figure8data(ii).aggcapital=AllStats.A.Mean;
        Figure8data(ii).aggoutput=Y;
        Figure8data(ii).entrepreneur_output=AllStats.Y_noncorp.Mean;
        Figure8data(ii).entrepreneur_capital=AllStats.K_noncorp.Mean;
        Figure8data(ii).entrepreneur_labor=AllStats.N_noncorp.Mean;
        Figure8data(ii).WealthGini=AllStats.A.Gini;
        
        % As for Figure4data above: only the tau_E1 values that get a transition path in doPart(7) have
        % their V and Policy read back, and nothing ever reads StationaryDist back
        if keepVPolicy_tauE1(ii)
            Figure8data(ii).V=V;
            Figure8data(ii).Policy=Policy;
        end
        
        Figure8data(ii).p_eqm_final=p_eqm_final;
        Figure8data(ii).GeneralEqmCondn=GeneralEqmCondn;
        
        Figure8data(ii).MassAtTopOfAssetGrid=sum(sum(sum(sum(StationaryDist(end-9:end,:,:,:)))));
        Figure8data(ii).entry=AllStats.EntrepreneurEntry.Mean;
        Figure8data(ii).exit=AllStats.EntrepreneurExit.Mean;
        Figure8data(ii).entrepreneurs=AllStats.FractionEntrepreneur.Mean;

        AggVars2=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate_Fig2,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
        Figure8data(ii).AggVars2=AggVars2;

        simoptions.conditionalrestrictions.etheta2=@(aprime,eprime,a,e,eta,theta,theta2) (theta==theta2)*e;
        simoptions.conditionalrestrictions.etheta3=@(aprime,eprime,a,e,eta,theta,theta3) (theta==theta3)*e;
        simoptions.conditionalrestrictions.etheta4=@(aprime,eprime,a,e,eta,theta,theta4) (theta==theta4)*e;
        AllStats2=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate_Table5,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
        simoptions=rmfield(simoptions,'conditionalrestrictions');
        Figure8data(ii).AllStats2=AllStats2;
        
        if Params.tau_E1==0.1
            p_eqm_final1=p_eqm_final;
        elseif Params.tau_E1==0.4
            p_eqm_final2=p_eqm_final;
        end
        
        save ./SavedOutput/Kitao2008_Figure8data.mat Figure8data
    end

    % As for the tau_k experiments above
    [worstmass,worstii]=max(cell2mat({Figure8data(:).MassAtTopOfAssetGrid}));
    fprintf('tau_E1 experiments: worst mass in the top 10 asset gridpoints = %e, at tau_E1=%.2f \n', worstmass, tau_E1_vec(worstii))
    if worstmass>10^(-5)
        warning('tau_E1 experiments: %.4e mass in the top 10 asset gridpoints at tau_E1=%.2f, raise assetmaxfactor',worstmass,tau_E1_vec(worstii))
    end
    
    save ./SavedOutput/Kitao2008_FinalGE2.mat p_eqm_final1  p_eqm_final2

    % Figure 8
    fig8=figure(8);
    subplot(2,3,1); plot(100*tau_E1_vec,Fig4_benchmark.w*ones(1,length(tau_E1_vec)),'b.',100*tau_E1_vec,cell2mat({Figure8data(:).w}),'b-');
    title('wage')
    subplot(2,3,2); plot(100*tau_E1_vec,100*Fig4_benchmark.r*ones(1,length(tau_E1_vec)),'b.',100*tau_E1_vec,100*cell2mat({Figure8data(:).r}),'b-');
    title('interest rate (%)')
    subplot(2,3,3); plot(100*tau_E1_vec,100*Fig4_benchmark.tau_I*ones(1,length(tau_E1_vec)),'b.',100*tau_E1_vec,100*cell2mat({Figure8data(:).tau_I}),'b-');
    title('proportional tax \tau_I (%)')
    subplot(2,3,4); yyaxis left 
    plot(100*tau_E1_vec,Fig4_benchmark.A*ones(1,length(tau_E1_vec)),'b.',100*tau_E1_vec,cell2mat({Figure8data(:).aggcapital}),'b-');
    subplot(2,3,4); yyaxis right 
    plot(100*tau_E1_vec,Fig4_benchmark.Y*ones(1,length(tau_E1_vec)),'r.',100*tau_E1_vec, cell2mat({Figure8data(:).aggoutput}),'r-');
    yyaxis left 
    title('aggregate activities')
    legend('','agg. capital','','agg. output','northwest')
    subplot(2,3,5); plot(100*tau_E1_vec,cell2mat({Figure8data(:).entrepreneur_capital})/Fig4_benchmark.K_noncorp,'b-',100*tau_E1_vec,cell2mat({Figure8data(:).entrepreneur_output})/Fig4_benchmark.Y_noncorp,'g-',100*tau_k_vec, cell2mat({Figure8data(:).entrepreneur_labor})/Fig4_benchmark.N_noncorp,'r-');
    title('entrepreneurial activities (normalized)')
    legend('capital','output','labor','northeast') % Note: Order in Figure 4 was output,capital,labor, but here is capital,output,labor [following K2008]
    subplot(2,3,6); plot(100*tau_E1_vec,Fig4_benchmark.WealthGini*ones(1,length(tau_E1_vec)),'b.',100*tau_E1_vec,cell2mat({Figure8data(:).WealthGini}),'b-');
    title('Wealth Gini')
    saveas(fig8,'./SavedOutput/Graphs/Kitao2008_Fig8.png')

    
    %% Table 8
    [~,ii_tau_E1_0p1]=min(abs(tau_E1_vec-0.1)); % find tau_k=0.1
    [~,ii_tau_E1_0p4]=min(abs(tau_E1_vec-0.4)); % find tau_k=0.4
    
    Figure8_tau_E1_1=[Figure8data(ii_tau_E1_0p1).AggVars2.etheta1.Mean,Figure8data(ii_tau_E1_0p1).AggVars2.etheta2.Mean,Figure8data(ii_tau_E1_0p1).AggVars2.etheta3.Mean,Figure8data(ii_tau_E1_0p1).AggVars2.etheta4.Mean];
    Figure8_tau_E1_2=[Figure8data(ii_tau_E1_0p4).AggVars2.etheta1.Mean,Figure8data(ii_tau_E1_0p4).AggVars2.etheta2.Mean,Figure8data(ii_tau_E1_0p4).AggVars2.etheta3.Mean,Figure8data(ii_tau_E1_0p4).AggVars2.etheta4.Mean];
    
    Figure8_tau_E1_1a=[Figure8data(ii_tau_E1_0p1).AllStats2.etheta2.k.Mean, Figure8data(ii_tau_E1_0p1).AllStats2.etheta3.k.Mean,Figure8data(ii_tau_E1_0p1).AllStats2.etheta4.k.Mean,Figure8data(ii_tau_E1_0p1).AllStats2.k.Mean/Figure8data(ii_tau_E1_0p1).entrepreneurs];
    Figure8_tau_E1_2a=[Figure8data(ii_tau_E1_0p4).AllStats2.etheta2.k.Mean, Figure8data(ii_tau_E1_0p4).AllStats2.etheta3.k.Mean,Figure8data(ii_tau_E1_0p4).AllStats2.etheta4.k.Mean,Figure8data(ii_tau_E1_0p4).AllStats2.k.Mean/Figure8data(ii_tau_E1_0p4).entrepreneurs];

    % Table 8
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table8.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n');
    fprintf(FID, '\\multicolumn{7}{l}{Flat tax on entrepreneurs business income} \\\\  \\hline \n');
    fprintf(FID, ' \\multicolumn{7}{c}{(a) Entrepreneurs in the population} \\\\ \\hline \n');
    fprintf(FID, ' & \\multicolumn{4}{l}{\\%% of entrepreneurs} & exit rate & entry rate \\\\ \n');
    fprintf(FID, ' & $\\theta_2$ & $\\theta_3$ & $\\theta_4$ & all & (\\%%) & (\\%%) \\\\ \\hline \n');
    fprintf(FID, 'benchmark            & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  100*Figure2values(1,2),  100*Figure2values(1,3), 100*Figure2values(1,4), 100*sum(Figure2values(1,2:4)),  100*Figure8.initial.exit/Figure8.initial.entrepreneurs,   100*Figure8.initial.entry/(1-Figure8.initial.entrepreneurs));
    fprintf(FID, '$\\tau_{E1}=10\\%%$  & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  100*Figure8_tau_E1_1(2), 100*Figure8_tau_E1_1(3),100*Figure8_tau_E1_1(4),100*sum(Figure8_tau_E1_1(2:4)), 100*Figure8data(ii_tau_E1_0p1).exit/Figure8data(ii_tau_E1_0p1).entrepreneurs,  100*Figure8data(ii_tau_E1_0p1).entry/(1-Figure8data(ii_tau_E1_0p1).entrepreneurs));
    fprintf(FID, '$\\tau_{E1}=40\\%%$  & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n',  100*Figure8_tau_E1_2(2), 100*Figure8_tau_E1_2(3),100*Figure8_tau_E1_2(4),100*sum(Figure8_tau_E1_2(2:4)), 100*Figure8data(ii_tau_E1_0p4).exit/Figure8data(ii_tau_E1_0p4).entrepreneurs,  100*Figure8data(ii_tau_E1_0p4).entry/(1-Figure8data(ii_tau_E1_0p4).entrepreneurs));
    fprintf(FID, ' & \\multicolumn{5}{c}{(b) Investment by entrepreneurs} & \\\\ \\hline \n');
    fprintf(FID, ' & & \\multicolumn{4}{l}{avg. investment in \\$1,000} & \\\\ \\hline \n');
    fprintf(FID, ' & & $\\theta_2$ & $\\theta_3$ & $\\theta_4$ & all & \\\\ \\hline \n');
    fprintf(FID, ' & benchmark            & %8.0f & %8.0f & %8.0f & %8.0f &  \\\\ \n',  Table8.avginvest*wealthscalingfactor/1000);
    fprintf(FID, ' & $\\tau_{E1}=10\\%%$  & %8.0f & %8.0f & %8.0f & %8.0f & \\\\ \n',  Figure8_tau_E1_1a*wealthscalingfactor/1000);
    fprintf(FID, ' & $\\tau_{E1}=40\\%%$  & %8.0f & %8.0f & %8.0f & %8.0f & \\\\ \n',  Figure8_tau_E1_2a*wealthscalingfactor/1000);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: exit rate is the mass of entrepreneurs who exit as a percentage of entrepreneurs; entry rate is the mass of workers who become entrepreneurs as a percentage of workers (in a stationary eqm these two masses are equal, they just have different denominators). \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    % As for doPart(1) above, plus what this block adds: tau_E1_vec, keepVPolicy_tauE1, Figure8data
    save ./SavedOutput/Kitao2008_doPart6.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8 Params_init p_eqm_final1 p_eqm_final2 tau_k_vec keepVPolicy_tauk Figure4data T V_init StationaryDist_init AggVars_init Table6data_mass FnsToEvaluate_TransPath FnsToEvaluate2_TransPath TransPathGeneralEqmEqns transpathoptions Figure7data tau_E1_vec keepVPolicy_tauE1 Figure8data
else
    load ./SavedOutput/Kitao2008_doPart6.mat
end



%% Compute 9 transitions for tau_e, including tau_e=0.1 and 0.4 and do some additional analysis of this last two
if doPart(7)==1
    Params.taxincome=3; % 1 is tax income, 2 is to tax capital income and labor income seperately, 3 is tax entrepreneurial business income

    % Same floors as in doPart(6) above, for the same reason. NO LaborPenalty here though: the
    % transition path update needs exactly one row of GEnewprice2.howtoupdate per general eqm eqn
    % (setupGEnewprice3_shooting errors otherwise) and there is no fourth price to give it. The
    % floors alone keep the conditions real, and the shooting update that Anderson is accelerating
    % already pushes w up hard when the LaborMarket residual goes large and negative, which is the
    % direction that reduces entrepreneurial labour demand.
    TransPathGeneralEqmEqns.CapitalMarket = @(r,K_noncorp,A,N_noncorp,L,alpha,delta) r-(alpha*(max(A-K_noncorp,10^(-6))^(alpha-1))*(max(L-N_noncorp,10^(-2))^(1-alpha))-delta);
    TransPathGeneralEqmEqns.LaborMarket = @(w,K_noncorp,A,N_noncorp,L,alpha) w-(1-alpha)*(max(A-K_noncorp,10^(-6))^(alpha))*(max(L-N_noncorp,10^(-2))^(-alpha));
    
    Figure10and11=struct();
    
    for ii=1:length(tau_E1_transition_vec) % tau_E1_transition_vec is set up at the top
        Params.tau_E1=tau_E1_transition_vec(ii);
        
        [~,jj]=min(abs(tau_E1_vec-Params.tau_E1)); % nearest value, not exact equality (see note at the doPart(4) lookup); jj is the index from when the stationary general eqm was computed earlier
        V_final=Figure8data(jj).V;
        Policy_final=Figure8data(jj).Policy;
        p_eqm_final=Figure8data(jj).p_eqm_final;
        
        % We want to look at a one off unanticipated path of tau_k.
        ParamPath=struct(); % (clear out any tau_k path left over from the earlier experiments)
        ParamPath.tau_E1=Params.tau_E1*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
        % (the way ParamPath is set is designed to allow for a series of changes in the parameters)

        % We need to give an initial guess for the price path
        PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
        PricePath0.w=[linspace(p_eqm_initial.w, p_eqm_final.w, floor(T/3))'; p_eqm_final.w*ones(T-floor(T/3),1)];
        PricePath0.tau_I=[linspace(p_eqm_initial.tau_I, p_eqm_final.tau_I, floor(T/3))'; p_eqm_final.tau_I*ones(T-floor(T/3),1)];

        
        [PricePath,GEcondnPath]=TransitionPath_InfHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_TransPath, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions,simoptions,vfoptionstpath);
        
        [Vpath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
        AgentDistPath=AgentDistOnTransPath_InfHorz(StationaryDist_init, PricePath, ParamPath, PolicyPath,n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
        AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2_TransPath,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

        % CEV calculations (CES utility fn makes this easy)
        CEVpath=(Vpath(:,:,:,:,1)./V_init).^(1/(1-Params.sigma))-1;
        % CEVpath can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove them
        CEVpath((StationaryDist_init==0))=0;
        CEVpath_behindtheveil=sum(sum(sum(sum(CEVpath.*StationaryDist_init))));
        
        % CEV calculations for steady-state comparison
        CEV=(V_final./V_init).^(1/(1-Params.sigma))-1;
        % CEV can have some NaN from -Inf/-Inf, but these correspond to points where no-one is, so remove them
        CEV((StationaryDist_init==0))=0;
        CEV_behindtheveil=sum(sum(sum(sum(CEV.*StationaryDist_init))));
        
        Figure10and11.(['ii',num2str(ii)]).CEVpath=CEVpath;
        Figure10and11.CEVpath_behindtheveil(ii)=CEVpath_behindtheveil;
        Figure10and11.CEV_behindtheveil(ii)=CEV_behindtheveil;
       
        % So we can double check them
        Figure10and11.(['ii',num2str(ii)]).PricePath=PricePath;
        Figure10and11.(['ii',num2str(ii)]).GEcondnPath=GEcondnPath;

        if Params.tau_E1==0.1
            PricePath1=PricePath;
            AggVarsPath1=AggVarsPath;
            Output_corp_path1=((AggVarsPath1.A.Mean-AggVarsPath1.K_noncorp.Mean).^Params.alpha).*((AggVarsPath1.L.Mean-AggVarsPath1.N_noncorp.Mean).^(1-Params.alpha));
            CEV1=CEVpath;
        elseif Params.tau_E1==0.4
            PricePath2=PricePath;
            AggVarsPath2=AggVarsPath;
            Output_corp_path2=((AggVarsPath2.A.Mean-AggVarsPath2.K_noncorp.Mean).^Params.alpha).*((AggVarsPath2.L.Mean-AggVarsPath2.N_noncorp.Mean).^(1-Params.alpha));
            CEV2=CEVpath;
        end
        
    end
    
    %% Figure 9
    
    % Savings is the change in assets, plus the depreciation; savings rate is that over income.
    % These have to be computed here from this experiment's own paths. They used to be inherited from
    % doPart(3), which meant this figure was plotting the savings rate of the tau_k transitions rather
    % than of the tau_E1 transitions that every other panel is drawing.
    % Note: all of the paths below are padded with the initial stationary eqm value, so they are length
    % T+1 and line up with the 0:1:T used on the horizontal axis. (Output_corp_path1 and
    % Output_corp_path2 as set inside the loop above are length T, not T+1, so they are not reused here.)
    A_path1=[AggVars_init.A.Mean, AggVarsPath1.A.Mean];
    K_noncorp_path1=[AggVars_init.K_noncorp.Mean, AggVarsPath1.K_noncorp.Mean];
    L_path1=[AggVars_init.L.Mean, AggVarsPath1.L.Mean];
    N_noncorp_path1=[AggVars_init.N_noncorp.Mean, AggVarsPath1.N_noncorp.Mean];
    Y_path1=((A_path1-K_noncorp_path1).^Params.alpha).*((L_path1-N_noncorp_path1).^(1-Params.alpha))+[AggVars_init.Y_noncorp.Mean, AggVarsPath1.Y_noncorp.Mean];
    savings_path1=Params.delta*A_path1+[0, A_path1(2:end)-A_path1(1:end-1)];
    savingsrate_path1=savings_path1./Y_path1;

    A_path2=[AggVars_init.A.Mean, AggVarsPath2.A.Mean];
    K_noncorp_path2=[AggVars_init.K_noncorp.Mean, AggVarsPath2.K_noncorp.Mean];
    L_path2=[AggVars_init.L.Mean, AggVarsPath2.L.Mean];
    N_noncorp_path2=[AggVars_init.N_noncorp.Mean, AggVarsPath2.N_noncorp.Mean];
    Y_path2=((A_path2-K_noncorp_path2).^Params.alpha).*((L_path2-N_noncorp_path2).^(1-Params.alpha))+[AggVars_init.Y_noncorp.Mean, AggVarsPath2.Y_noncorp.Mean];
    savings_path2=Params.delta*A_path2+[0, A_path2(2:end)-A_path2(1:end-1)];
    savingsrate_path2=savings_path2./Y_path2;

    fig5=figure(9);
    subplot(2,4,1); plot(0:1:T,[p_eqm_initial.w, PricePath1.w],'b-',0:1:T,[p_eqm_initial.w, PricePath2.w],'r-.')
    title('wage')
    subplot(2,4,2); plot(0:1:T,100*[p_eqm_initial.r, PricePath1.r],'b-',0:1:T,100*[p_eqm_initial.r, PricePath2.r],'r-.')
    title('interest rate (%)')
    subplot(2,4,3); plot(0:1:T,[AggVars_init.A.Mean, AggVarsPath1.A.Mean],'b-',0:1:T,[AggVars_init.A.Mean, AggVarsPath2.A.Mean],'r-.')
    title('agg. capital')
    subplot(2,4,4); plot(0:1:T,[AggVars_init.C.Mean, AggVarsPath1.C.Mean],'b-',0:1:T,[AggVars_init.C.Mean, AggVarsPath2.C.Mean],'r-.')
    title('agg. consumption')
    subplot(2,4,5); plot(0:1:T,100*savingsrate_path1,'b-',0:1:T,100*savingsrate_path2,'r-.')
    title('savings rate')
    subplot(2,4,6); plot(0:1:T,100*[p_eqm_initial.tau_I, PricePath1.tau_I],'b-',0:1:T,100*[p_eqm_initial.tau_I, PricePath2.tau_I],'r-.')
    title('proportional tax \tau_I (%)')
    subplot(2,4,7); plot(0:1:T,[AggVars_init.K_noncorp.Mean, AggVarsPath1.K_noncorp.Mean],'b-',0:1:T,[AggVars_init.K_noncorp.Mean, AggVarsPath2.K_noncorp.Mean],'r-.')
    title('entrep. capital')
    subplot(2,4,8); plot(0:1:T,[AggVars_init.Y_noncorp.Mean, AggVarsPath1.Y_noncorp.Mean],'b-',0:1:T,[AggVars_init.Y_noncorp.Mean, AggVarsPath2.Y_noncorp.Mean],'r-.')
    title('entrep. output')
    saveas(fig5,'./SavedOutput/Graphs/Kitao2008_Fig9.png')
    
    %% Figure 10
    
    fig10=figure(10);
    % tau_E1=0.1
    workermass=StationaryDist_init(:,1,:,:); % 1 is workers
    workerCEV1=sum(CEV1(:,1,:,:).*workermass,4)./sum(workermass,4); % average over theta, weighted by the initial distribution (Kitao (2008), footnote 15)
    % Note: this is NaN wherever there is no mass, which is what we want, as those points then simply do not get plotted
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,1),'b--')
    hold on
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,3),'g.')
    subplot(2,2,1); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV1(:,1,5),'r-')
    hold off
    legend('eta 1','eta 3', 'eta 5')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('workers')
    
    entrepreneurmass=StationaryDist_init(:,2,:,:); % 2 is entrepreneurs
    entrepreneurCEV1=sum(CEV1(:,2,:,:).*entrepreneurmass,3)./sum(entrepreneurmass,3); % average over eta, weighted by the initial distribution (Kitao (2008), footnote 15)
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,2),'b--')
    hold on
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,3),'g.')
    subplot(2,2,2); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV1(:,1,1,4),'r-')
    hold off
    legend('theta 2','theta 3', 'theta 4')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('entrepreneurs')
    
    % tau_E1=0.4
    workermass=StationaryDist_init(:,1,:,:); % 1 is workers
    workerCEV2=sum(CEV2(:,1,:,:).*workermass,4)./sum(workermass,4); % average over theta, weighted by the initial distribution (Kitao (2008), footnote 15)
    % Note: this is NaN wherever there is no mass, which is what we want, as those points then simply do not get plotted
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,1),'b--')
    hold on
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,3),'g.')
    subplot(2,2,3); plot(asset_grid*wealthscalingfactor/1000,100*workerCEV2(:,1,5),'r-')
    hold off
    legend('eta 1','eta 3', 'eta 5')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('workers')
    
    entrepreneurmass=StationaryDist_init(:,2,:,:); % 2 is entrepreneurs
    entrepreneurCEV2=sum(CEV2(:,2,:,:).*entrepreneurmass,3)./sum(entrepreneurmass,3); % average over eta, weighted by the initial distribution (Kitao (2008), footnote 15)
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,2),'b--')
    hold on
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,3),'g.')
    subplot(2,2,4); plot(asset_grid*wealthscalingfactor/1000,100*entrepreneurCEV2(:,1,1,4),'r-')
    hold off
    legend('theta 2','theta 3', 'theta 4')
    xlabel('assets in $1,000')
    ylabel('percentage (%)')
    title('entrepreneurs')
    
    saveas(fig10,'./SavedOutput/Graphs/Kitao2008_Fig10.png')
    
    
    %% Figure 11
    
    fig11=figure(11);
    subplot(1,2,1); plot(100*tau_E1_transition_vec,100*Figure10and11.CEVpath_behindtheveil,'b-o')
    title('transition')
    ylabel('percentage (%)')
    subplot(1,2,2); plot(100*tau_E1_transition_vec,100*Figure10and11.CEV_behindtheveil,'b-o')
    title('steady state')
    ylabel('percentage (%)')
    saveas(fig11,'./SavedOutput/Graphs/Kitao2008_Fig11.png')
    
    %% Table 9
    
    
    % CEV again, but report by wealth levels
    wealthcutoffs=[10000,50000,100000,250000,500000,1000000,2000000]; % wealth cutoffs in dollars
    wealthcutoffs=wealthcutoffs/wealthscalingfactor; % wealth cutoffs in model units
    
    % First, deal with
    Table9data_mass=Table6data_mass; % Mass of agents in each wealth band
    % Is just the same initial distribution, so no change
    
    % Repeat, but this time only the mass for agents with positive CEV
    Table9data_tauE101_gainmass=zeros(9,3); % Mass of agents in each wealth band with positive CEV
    % Stationary dist for the tau_k=0 eqm
    StationaryDist_temp=StationaryDist_init;
    CEV_tauE101=CEV1;
    gainCEV_tauE101=(CEV_tauE101>0);
    % So the distribution mass of those who gained from the reform is
    StationaryDist_temp=StationaryDist_temp.*gainCEV_tauE101;
    % Now get the mass of those who gained
    StationaryDist_temp=sum(sum(StationaryDist_temp,4),3);
    % Worker (1 in second dim of StationaryDist)
    Table9data_tauE101_gainmass(1,1)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),1));
    Table9data_tauE101_gainmass(2,1)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),1));
    Table9data_tauE101_gainmass(3,1)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),1));
    Table9data_tauE101_gainmass(4,1)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),1));
    Table9data_tauE101_gainmass(5,1)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),1));
    Table9data_tauE101_gainmass(6,1)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),1));
    Table9data_tauE101_gainmass(7,1)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),1));
    Table9data_tauE101_gainmass(8,1)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),1));
    Table9data_tauE101_gainmass(9,1)=sum(StationaryDist_temp(:,1));
    % Entrepreneur (2 in second dim of StationaryDist)
    Table9data_tauE101_gainmass(1,2)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),2));
    Table9data_tauE101_gainmass(2,2)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),2));
    Table9data_tauE101_gainmass(3,2)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),2));
    Table9data_tauE101_gainmass(4,2)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),2));
    Table9data_tauE101_gainmass(5,2)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),2));
    Table9data_tauE101_gainmass(6,2)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),2));
    Table9data_tauE101_gainmass(7,2)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),2));
    Table9data_tauE101_gainmass(8,2)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),2));
    Table9data_tauE101_gainmass(9,2)=sum(StationaryDist_temp(:,2));
    % All (sum over second dimension)
    Table9data_tauE101_gainmass(1,3)=sum(sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),:),2));
    Table9data_tauE101_gainmass(2,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),:),2));
    Table9data_tauE101_gainmass(3,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),:),2));
    Table9data_tauE101_gainmass(4,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),:),2));
    Table9data_tauE101_gainmass(5,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),:),2));
    Table9data_tauE101_gainmass(6,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),:),2));
    Table9data_tauE101_gainmass(7,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),:),2));
    Table9data_tauE101_gainmass(8,3)=sum(sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),:),2));
    Table9data_tauE101_gainmass(9,3)=sum(sum(StationaryDist_temp(:,:),2));
    
    
    % Redo all the calculations for tau_k=0.4
    % Repeat, but this time only the mass for agents with positive CEV
    Table9data_tauE104_gainmass=zeros(9,3); % Mass of agents in each wealth band with positive CEV
    % Stationary dist for the tau_k=0.4 eqm
    StationaryDist_temp=StationaryDist_init;
    CEV_tauE104=CEV2;
    gainCEV_tauE104=(CEV_tauE104>0);
    % So the distribution mass of those who gained from the reform is
    StationaryDist_temp=StationaryDist_temp.*gainCEV_tauE104;
    % Now get the mass of those who gained
    StationaryDist_temp=sum(sum(StationaryDist_temp,4),3);
    % Worker (1 in second dim of StationaryDist)
    Table9data_tauE104_gainmass(1,1)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),1));
    Table9data_tauE104_gainmass(2,1)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),1));
    Table9data_tauE104_gainmass(3,1)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),1));
    Table9data_tauE104_gainmass(4,1)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),1));
    Table9data_tauE104_gainmass(5,1)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),1));
    Table9data_tauE104_gainmass(6,1)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),1));
    Table9data_tauE104_gainmass(7,1)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),1));
    Table9data_tauE104_gainmass(8,1)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),1));
    Table9data_tauE104_gainmass(9,1)=sum(StationaryDist_temp(:,1));
    % Entrepreneur (2 in second dim of StationaryDist)
    Table9data_tauE104_gainmass(1,2)=sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),2));
    Table9data_tauE104_gainmass(2,2)=sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),2));
    Table9data_tauE104_gainmass(3,2)=sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),2));
    Table9data_tauE104_gainmass(4,2)=sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),2));
    Table9data_tauE104_gainmass(5,2)=sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),2));
    Table9data_tauE104_gainmass(6,2)=sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),2));
    Table9data_tauE104_gainmass(7,2)=sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),2));
    Table9data_tauE104_gainmass(8,2)=sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),2));
    Table9data_tauE104_gainmass(9,2)=sum(StationaryDist_temp(:,2));
    % All (sum over second dimension)
    Table9data_tauE104_gainmass(1,3)=sum(sum(StationaryDist_temp(asset_grid<wealthcutoffs(1),:),2));
    Table9data_tauE104_gainmass(2,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(1)<=asset_grid).*(asset_grid<wealthcutoffs(2))),:),2));
    Table9data_tauE104_gainmass(3,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(2)<=asset_grid).*(asset_grid<wealthcutoffs(3))),:),2));
    Table9data_tauE104_gainmass(4,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(3)<=asset_grid).*(asset_grid<wealthcutoffs(4))),:),2));
    Table9data_tauE104_gainmass(5,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(4)<=asset_grid).*(asset_grid<wealthcutoffs(5))),:),2));
    Table9data_tauE104_gainmass(6,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(5)<=asset_grid).*(asset_grid<wealthcutoffs(6))),:),2));
    Table9data_tauE104_gainmass(7,3)=sum(sum(StationaryDist_temp(logical((wealthcutoffs(6)<=asset_grid).*(asset_grid<wealthcutoffs(7))),:),2));
    Table9data_tauE104_gainmass(8,3)=sum(sum(StationaryDist_temp((wealthcutoffs(7)<=asset_grid),:),2));
    Table9data_tauE104_gainmass(9,3)=sum(sum(StationaryDist_temp(:,:),2));
    
    
    
    
    %Table 9
    FID = fopen('./SavedOutput/LatexInputs/Kitao2008_Table9.tex', 'w');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n \\hline\n');
    fprintf(FID, '\\multicolumn{7}{l}{Fraction of households with welfare gains: flat tax on entrepreneurs business income} \\\\  \\hline \n');
    fprintf(FID, ' Wealth (in \\$1,000) & Workers & & Entrepreneurs & & All &  \\\\ \\hline \n');
    fprintf(FID, '\\multicolumn{7}{c}{(a) Entrepreneurial business tax $\\tau_{E1} =10$\\%%}    \\\\ \n');
    fprintf(FID, '0-10      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(1,1),100*Table9data_tauE101_gainmass(1,1)/Table9data_mass(1,1), 100*Table9data_tauE101_gainmass(1,2),100*Table9data_tauE101_gainmass(1,2)/Table9data_mass(1,2), 100*Table9data_tauE101_gainmass(1,3),100*Table9data_tauE101_gainmass(1,3)/Table9data_mass(1,3));
    fprintf(FID, '10-50     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(2,1),100*Table9data_tauE101_gainmass(2,1)/Table9data_mass(2,1), 100*Table9data_tauE101_gainmass(2,2),100*Table9data_tauE101_gainmass(2,2)/Table9data_mass(2,2), 100*Table9data_tauE101_gainmass(2,3),100*Table9data_tauE101_gainmass(2,3)/Table9data_mass(2,3));
    fprintf(FID, '50-100    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(3,1),100*Table9data_tauE101_gainmass(3,1)/Table9data_mass(3,1), 100*Table9data_tauE101_gainmass(3,2),100*Table9data_tauE101_gainmass(3,2)/Table9data_mass(3,2), 100*Table9data_tauE101_gainmass(3,3),100*Table9data_tauE101_gainmass(3,3)/Table9data_mass(3,3));
    fprintf(FID, '100-250   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(4,1),100*Table9data_tauE101_gainmass(4,1)/Table9data_mass(4,1), 100*Table9data_tauE101_gainmass(4,2),100*Table9data_tauE101_gainmass(4,2)/Table9data_mass(4,2), 100*Table9data_tauE101_gainmass(4,3),100*Table9data_tauE101_gainmass(4,3)/Table9data_mass(4,3));
    fprintf(FID, '250-500   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(5,1),100*Table9data_tauE101_gainmass(5,1)/Table9data_mass(5,1), 100*Table9data_tauE101_gainmass(5,2),100*Table9data_tauE101_gainmass(5,2)/Table9data_mass(5,2), 100*Table9data_tauE101_gainmass(5,3),100*Table9data_tauE101_gainmass(5,3)/Table9data_mass(5,3));
    fprintf(FID, '500-1000  & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(6,1),100*Table9data_tauE101_gainmass(6,1)/Table9data_mass(6,1), 100*Table9data_tauE101_gainmass(6,2),100*Table9data_tauE101_gainmass(6,2)/Table9data_mass(6,2), 100*Table9data_tauE101_gainmass(6,3),100*Table9data_tauE101_gainmass(6,3)/Table9data_mass(6,3));
    fprintf(FID, '1000-2000 & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(7,1),100*Table9data_tauE101_gainmass(7,1)/Table9data_mass(7,1), 100*Table9data_tauE101_gainmass(7,2),100*Table9data_tauE101_gainmass(7,2)/Table9data_mass(7,2), 100*Table9data_tauE101_gainmass(7,3),100*Table9data_tauE101_gainmass(7,3)/Table9data_mass(7,3));
    fprintf(FID, '$>$2000     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(8,1),100*Table9data_tauE101_gainmass(8,1)/Table9data_mass(8,1), 100*Table9data_tauE101_gainmass(8,2),100*Table9data_tauE101_gainmass(8,2)/Table9data_mass(8,2), 100*Table9data_tauE101_gainmass(8,3),100*Table9data_tauE101_gainmass(8,3)/Table9data_mass(8,3));
    fprintf(FID, 'all       & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE101_gainmass(9,1),100*Table9data_tauE101_gainmass(9,1)/Table9data_mass(9,1), 100*Table9data_tauE101_gainmass(9,2),100*Table9data_tauE101_gainmass(9,2)/Table9data_mass(9,2), 100*Table9data_tauE101_gainmass(9,3),100*Table9data_tauE101_gainmass(9,3)/Table9data_mass(9,3));
    fprintf(FID, '\\multicolumn{7}{c}{ }    \\\\ \n');
    fprintf(FID, '\\multicolumn{7}{c}{(b) Entrepreneurial business tax $\\tau_{E1} =40$\\%%}    \\\\ \n');
    fprintf(FID, '0-10      & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(1,1),100*Table9data_tauE104_gainmass(1,1)/Table9data_mass(1,1), 100*Table9data_tauE104_gainmass(1,2),100*Table9data_tauE104_gainmass(1,2)/Table9data_mass(1,2), 100*Table9data_tauE104_gainmass(1,3),100*Table9data_tauE104_gainmass(1,3)/Table9data_mass(1,3));
    fprintf(FID, '10-50     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(2,1),100*Table9data_tauE104_gainmass(2,1)/Table9data_mass(2,1), 100*Table9data_tauE104_gainmass(2,2),100*Table9data_tauE104_gainmass(2,2)/Table9data_mass(2,2), 100*Table9data_tauE104_gainmass(2,3),100*Table9data_tauE104_gainmass(2,3)/Table9data_mass(2,3));
    fprintf(FID, '50-100    & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(3,1),100*Table9data_tauE104_gainmass(3,1)/Table9data_mass(3,1), 100*Table9data_tauE104_gainmass(3,2),100*Table9data_tauE104_gainmass(3,2)/Table9data_mass(3,2), 100*Table9data_tauE104_gainmass(3,3),100*Table9data_tauE104_gainmass(3,3)/Table9data_mass(3,3));
    fprintf(FID, '100-250   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(4,1),100*Table9data_tauE104_gainmass(4,1)/Table9data_mass(4,1), 100*Table9data_tauE104_gainmass(4,2),100*Table9data_tauE104_gainmass(4,2)/Table9data_mass(4,2), 100*Table9data_tauE104_gainmass(4,3),100*Table9data_tauE104_gainmass(4,3)/Table9data_mass(4,3));
    fprintf(FID, '250-500   & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(5,1),100*Table9data_tauE104_gainmass(5,1)/Table9data_mass(5,1), 100*Table9data_tauE104_gainmass(5,2),100*Table9data_tauE104_gainmass(5,2)/Table9data_mass(5,2), 100*Table9data_tauE104_gainmass(5,3),100*Table9data_tauE104_gainmass(5,3)/Table9data_mass(5,3));
    fprintf(FID, '500-1000  & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(6,1),100*Table9data_tauE104_gainmass(6,1)/Table9data_mass(6,1), 100*Table9data_tauE104_gainmass(6,2),100*Table9data_tauE104_gainmass(6,2)/Table9data_mass(6,2), 100*Table9data_tauE104_gainmass(6,3),100*Table9data_tauE104_gainmass(6,3)/Table9data_mass(6,3));
    fprintf(FID, '1000-2000 & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(7,1),100*Table9data_tauE104_gainmass(7,1)/Table9data_mass(7,1), 100*Table9data_tauE104_gainmass(7,2),100*Table9data_tauE104_gainmass(7,2)/Table9data_mass(7,2), 100*Table9data_tauE104_gainmass(7,3),100*Table9data_tauE104_gainmass(7,3)/Table9data_mass(7,3));
    fprintf(FID, '$>$2000     & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(8,1),100*Table9data_tauE104_gainmass(8,1)/Table9data_mass(8,1), 100*Table9data_tauE104_gainmass(8,2),100*Table9data_tauE104_gainmass(8,2)/Table9data_mass(8,2), 100*Table9data_tauE104_gainmass(8,3),100*Table9data_tauE104_gainmass(8,3)/Table9data_mass(8,3));
    fprintf(FID, 'all       & %8.1f & (%8.1f) & %8.1f & (%8.1f) & %8.1f & (%8.1f) \\\\ \n', 100*Table9data_tauE104_gainmass(9,1),100*Table9data_tauE104_gainmass(9,1)/Table9data_mass(9,1), 100*Table9data_tauE104_gainmass(9,2),100*Table9data_tauE104_gainmass(9,2)/Table9data_mass(9,2), 100*Table9data_tauE104_gainmass(9,3),100*Table9data_tauE104_gainmass(9,3)/Table9data_mass(9,3));
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\\\ \\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Note: In parentheses are the fractions of households with welfare gains conditional on their occupation and their wealth category. \n');
    fprintf(FID, 'A NaN in parentheses indicates a wealth category in which the occupation in question has zero mass, so the conditional fraction is undefined. \n');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
    
    
    % As for doPart(1) above, plus what this block adds: Figure10and11
    save ./SavedOutput/Kitao2008_doPart7.mat Params simoptions FnsToEvaluate p_eqm_initial Fig4_benchmark Figure2values Figure8 FnsToEvaluate_Fig2 FnsToEvaluate_Table5 Table2 Table3 Table4 Table5 Table7 Table8 Params_init p_eqm_final1 p_eqm_final2 tau_k_vec keepVPolicy_tauk Figure4data T V_init StationaryDist_init AggVars_init Table6data_mass FnsToEvaluate_TransPath FnsToEvaluate2_TransPath TransPathGeneralEqmEqns transpathoptions Figure7data tau_E1_vec keepVPolicy_tauE1 Figure8data Figure10and11
end

%% End of run
diary off
