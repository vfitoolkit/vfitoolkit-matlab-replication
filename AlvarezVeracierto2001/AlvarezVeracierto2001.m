% Alvarez & Veracierto (2001) - Severance payments in an economy with frictions
% Journal of Monetary Economics 47, 477-498. https://doi.org/10.1016/S0304-3932(01)00054-8
%
% Establishment dynamics (Hopenhayn-Rogerson style, with entry and exit),
% households who search for jobs, save in a riskless asset, retire and die,
% rigid wage contracts and no insurance markets. The government sets a
% firing tax epsilon (per worker, in units of the before-tax wage wstar) and a
% one-off severance payment lambda (per fired worker, in units of the after-tax
% wage w), balancing its budget with a payroll tax/subsidy wstar-w.
%
% Model period is half a quarter. Stationary equilibria only (no transitions).
% The equilibrium is one equation in the interest rate i (the paper's
% appendix); everything else nests inside: wstar from free entry, the
% establishment distribution up to scale, households at w=1, the scale nu
% from labour market clearing, w from the government budget, and goods market
% clearing as the residual. See AV2001subcodes/AV2001_StationaryEqm.m; the
% root-find on i is AV2001_SolveGE.m.
%
% Household exogenous state is (establishment state (e,s), employment status),
% restricted to the establishment states reachable from entry under the
% establishment policy; see AV2001_StationaryEqm.m for the layout and
% AV2001_HouseholdReturnFn.m for the encoding of status+income into one value.
%
% Two versions of the model:
%  - ENDOGENOUS job-finding rate (Params.endogsearch=1): search intensity eta
%    is a decision on a grid, and the unemployed-to-employed transition
%    probability eta^xi depends on it: a semi-exogenous shock
%    (vfoptions.n_semiz, n_d=n_eta, n_z=0). Tables 1, 3, 5.
%  - EXOGENOUS job-finding rate (Params.endogsearch=0): eta is fixed so that
%    eta^xi equals the average job-finding rate of the endogenous laissez-faire
%    economy (the paper's 61.0%; here the model's own value from doPart(1)).
%    Households have no decision variable and their shocks are a plain Markov
%    chain. Tables 2 and 4.
%
% Two typos in the published paper, both resolved from the 1998 WP (WP-98-2):
%  - eq (3.5) prints the unemployed budget as (1+i)a+w; eq (13) of the WP has
%    no wage income while unemployed, only interest plus severance/UI.
%  - p=0.037 should be p=0.0037 (WP Table 1); see the comment at Params.p.

% Record the full run (incl. any error) to AV2001_Diary.txt, fresh each run
diary off
if exist('./AV2001_Diary.txt','file'), delete('./AV2001_Diary.txt'); end
diary ./AV2001_Diary.txt

addpath('./AV2001subcodes');
addpath(genpath('./VFIToolkit-matlab/')); % for the server

if ~exist('./SavedOutput','dir'), mkdir('./SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'), mkdir('./SavedOutput/Graphs'); end
if ~exist('./SavedOutput/LatexInputs','dir'), mkdir('./SavedOutput/LatexInputs'); end

doPart=[1,1,1,1,1,1,1,1,1];
% doPart(1): the two laissez-faire economies: endogenous search first (its average job-finding rate sets eta for the exogenous economies), then exogenous
% doPart(2): Table 2, firing taxes rebated as employment subsidies (exogenous job-finding rate)
% doPart(3): Table 4, one-period unemployment insurance financed by payroll tax (exogenous job-finding rate)
% doPart(4): Table 1, severance payments (endogenous job-finding rate); column 1 is the endogenous laissez-faire from doPart(1)
% doPart(5): Table 3, firing taxes rebated as employment subsidies (endogenous job-finding rate)
% doPart(6): Table 5, one-period unemployment insurance financed by payroll tax (endogenous job-finding rate)
% doPart(7): calibrated U.S. economy with unemployment insurance, first column of Table 1 (rewrites Table 1 with it; run after doPart(4))
% doPart(8): DIAGNOSTICS for the capital/interest-rate side (Table 2: the paper's interest rate falls 4.0->3.8% with
%            firing taxes and capital rises; here it barely moves and capital falls). Needs doPart(1),(2),(4). The record
%            behind the write-up's section on this gap; about an hour, most of it the 2063-point employment grid.
% doPart(9): LaTeX tables for the write-up (LatexDoc/) that the table doParts do not write (the two calibration checks); solves nothing.
% Each doPart block save()s what later blocks need and the else branch load()s it back.

verbose=1; % 0: quiet, 1: one line per interest rate guess, 2: also the free entry solve

%% Parameters (Section 4, half-quarter periods)
% Preferences
Params.beta=0.99425;  % discount factor
Params.alpha=15.5;    % search disutility scale, u(1-eta)=alpha*((1-eta)^tau-1)/tau
Params.tau=0.98;      % search disutility curvature
Params.xi=0.98;       % job-finding probability is eta^xi
% Demographics
Params.rho=0.0031;    % retirement probability (40 years of active life)
Params.sigma=0.9917;  % retiree survival probability (15 years of retirement)
% Technology
Params.theta=0.19;    % capital share
Params.gamma=0.58;    % labour share
Params.delta=0.011;   % depreciation
Params.m=1950;        % entry cost, goods
% Idiosyncratic productivity s in {0,1,sbar}, eq (4.1)
Params.sbar=2.12;
Params.p=0.0037;      % probability of death (s=0) for s>0. The 2001 paper prints 0.037: a typo, Table 1 of the 1998 WP has 0.0037, and only 0.0037 gives the 0.73%-a-quarter job destruction from deaths, the 2.9% layoff rate and the 62-worker average size (the first GPU run with 0.037 gave 5.9% and 312)
Params.omega=0.973;   % persistence of the s>0 states
Params.psi=[0.5,0.5]; % entrants' distribution over {1,sbar}
% Policy: laissez-faire
Params.epsilon=0;     % firing tax per worker, units of wstar
Params.lambda=0;      % severance payment per fired worker, units of w
Params.endogsearch=1; % 0: eta=Params.eta fixed; 1: eta chosen on Grids.eta_grid (set per doPart block)
% Unemployment insurance (only in the calibrated U.S. economy, doPart(7)): replacement ratio, probability a fired worker
% is eligible, probability an eligible unemployed agent keeps eligibility. uirr=0 means no UI (and no eligibility states).
Params.uirr=0; Params.uikappa=0; Params.uiphi=0;
Params.uirr_US=0.66; Params.uikappa_US=0.5; % the U.S. system's replacement ratio and eligibility probability (Section 4)
% Exogenous search intensity: Params.eta is set in doPart(1) from the endogenous laissez-faire economy's
% average job-finding rate (the paper: "the average job-finding rate in the laissez-faire economy", 61.0%).
% The paper chose alpha=15.5 to get an average unemployment duration of one quarter, i.e. its 61.0% job-finding
% rate; with the grid interpolation layer the laissez-faire here gives 0.610 at alpha=15.5 (run of 2026-10-23;
% the grid convergence checks behind that number are recorded in NOTES.md).

% Initial guesses for the general equilibrium
Params.wstar_guess=0.7; % warm start for the free entry solve (updated after each equilibrium)
i_guess=(1.04)^(1/8)-1; % annual 4% real interest rate

%% Grids
% Establishments: e_grid is last period's employment (and the choice set for n).
% Dense at small sizes (low-shock establishments have a handful of workers) and
% sparse at large sizes. Mass at the top of the grid is reported as a check.
Grids.e_grid=[0,1:1:50,52:2:150,155:5:300,320:20:600]';
Grids.n_e=length(Grids.e_grid);
Grids.s_grid=[0;1;Params.sbar];
Grids.n_s=3;
Grids.Q=[1,0,0; Params.p,Params.omega*(1-Params.p),(1-Params.omega)*(1-Params.p); Params.p,(1-Params.omega)*(1-Params.p),Params.omega*(1-Params.p)];

% Households: assets in units of the after-tax wage (w=1 when solving). One
% grid for both households, so the exogenous and endogenous laissez-faire
% economies differ only by the search setup. Mean assets are about 40 wages in
% laissez-faire; mass at the top of the grid is reported as a check.
% The bottom of the grid is a_min>0 rather than 0: an unemployed agent with no
% severance and zero assets would have zero consumption for every choice
% (V=-Inf), which poisons the expectations; with a_min the worst case is
% c=i*a_min, finite. The borrowing constraint a'>=0 becomes a'>=a_min.
%
% GRID INTERPOLATION LAYER on assets (Grids.ngridinterp; 0 switches it off).
% Run of 2026-10-23 without it: 96% of ALL agents (employed, unemployed and
% retired) sat on one grid point at 36.4 wages, the savings policy never chose
% a'<a, and the job-finding rate moved from 0.95 to 0.67 with the asset grid
% density. Per-period saving flows here are ~0.1 wage (retirement saving over
% 320 periods, unemployed dissaving ~1 wage) against a stock of ~40 wages, so a
% pure-discretization step of 0.3-0.9 wages rounds every saving decision to
% a'=a. The interpolation layer (postGI) solves on the coarse grid, then
% searches a' on ngridinterp points between each pair of coarse points within
% +-maxaprimediff coarse points of the coarse optimum, with the value function
% interpolated; the distribution lives on the coarse grid with lottery weights.
% MEMORY: with endogenous search the fine arrays are
% (1+2*maxaprimediff+2*maxaprimediff*ngridinterp) x n_a x N_hh x n_eta doubles,
% three of them, plus the coarse return matrix n_eta x n_a^2 x N_hh. N_hh is
% 10 in laissez-faire and up to ~56 with twelve months of severance. The GPU
% (2026-10-23) fitted ~3 GB arrays and not 5.7 GB; the toolkit's default
% maxaprimediff=ceil(n_a/5) made the fine array 60+ GB. The settings below are
% 0.73 GB per fine array and 0.92 GB coarse at N_hh=56. AV2001_StationaryEqm
% prints both sizes with the N_R count (verbose>=1). The firm problem keeps
% its integer employment grid (no interpolation).
a_min=1e-3; a_max=150;
Grids.n_a=201;
Grids.a_grid=a_min+(a_max-a_min)*(linspace(0,1,Grids.n_a).^2)'; % quadratic spacing: step 0.73 wages at 36 wages, 0.004 at the bottom
Grids.ngridinterp=20;   % fine step 0.035 wages at 36 wages
Grids.maxaprimediff=3;  % coarse points each side of the coarse optimum searched on the fine grid (the toolkit notes postGI is exact at 3 without a decision variable; doPart(8c) checks 2/3/5 with one)
% Search intensity eta in [0,1): eta^xi is the job-finding probability, so the
% grid step (0.02) is roughly the resolution of an individual's job-finding
% rate; the reported rate averages over the asset distribution. The run of
% 2026-10-23 gave the same rate (+-0.01) with 51, 101 and 201 points.
n_eta=51;
Grids.eta_grid=linspace(0,0.99,n_eta)';
% Household value function: Howards with a sparse transition matrix (the toolkit default) while the household has
% at most this many states, by indexing above it. The sparse build expands the transition matrix densely, which
% grows with the square of the state count; the attrition-exact employment grid of doPart(8d) has ~860 states.
Grids.howardssparse_maxNhh=100;

% Toolkit options
vfoptions.verbose=0;
vfoptions.tolerance=1e-9;
simoptions.verbose=0;
simoptions.tolerance=1e-8;

fzerooptions_i=optimset('TolX',1e-7,'Display','iter');

%% doPart(1): the two laissez-faire economies
if doPart(1)==1
    % (a) endogenous search: sets the exogenous economies' eta
    fprintf('doPart(1a): laissez-faire economy with endogenous job-finding rate \n')
    Params.epsilon=0; Params.lambda=0; Params.endogsearch=1;
    tic;
    [i_LFendog,Out_LFendog]=AV2001_SolveGE(i_guess-0.0015,i_guess+0.0025,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,'laissez-faire, endogenous search');
    Params.wstar_guess=Out_LFendog.wstar;
    fprintf('doPart(1a): done in %.1f minutes. i=%.6f (annual %.3f%%), wstar=%.5f, w=%.5f \n',toc/60,i_LFendog,100*Out_LFendog.annualinterestrate,Out_LFendog.wstar,Out_LFendog.w)
    fprintf('  job-finding rate %.3f, unemployment rate %.2f%%, layoff rate %.2f%%, K/Y annual %.3f (measured output), average establishment size %.1f, %d reachable establishment states \n',Out_LFendog.jobfindingrate,100*Out_LFendog.unemploymentrate,100*Out_LFendog.layoffrate,Out_LFendog.KYannual,Out_LFendog.avgestablishmentsize,Out_LFendog.N_R)
    fprintf('  (paper: 0.61, 4.9%%, 2.9%%, 1.7, 62) \n')

    % (b) exogenous search at the endogenous laissez-faire average job-finding rate
    fprintf('doPart(1b): laissez-faire economy with exogenous job-finding rate \n')
    Params.jobfindingrate_exog=Out_LFendog.jobfindingrate;
    Params.eta=Params.jobfindingrate_exog^(1/Params.xi);
    Params.endogsearch=0;
    tic;
    [i_LF,Out_LF]=AV2001_SolveGE(i_LFendog-0.002,i_LFendog+0.001,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,'laissez-faire, exogenous search');
    Params.wstar_guess=Out_LF.wstar;
    fprintf('doPart(1b): done in %.1f minutes. i=%.6f (annual %.3f%%), wstar=%.5f, w=%.5f \n',toc/60,i_LF,100*Out_LF.annualinterestrate,Out_LF.wstar,Out_LF.w)
    fprintf('  unemployment rate %.2f%%, layoff rate %.2f%%, K/Y annual %.3f (measured output), average establishment size %.1f, %d reachable establishment states \n',100*Out_LF.unemploymentrate,100*Out_LF.layoffrate,Out_LF.KYannual,Out_LF.avgestablishmentsize,Out_LF.N_R)
    fprintf('  (paper: 4.9%%, 2.9%%, 1.7, 62) \n')
    save ./SavedOutput/AV2001_doPart1.mat Params i_LF Out_LF i_LFendog Out_LFendog % Grids, vfoptions and simoptions are defined above, not saved (a later run then sees their current values)
else
    % Only the Params fields doPart(1) sets are taken from the file, so parameters added to the script later are kept
    tmp=load('./SavedOutput/AV2001_doPart1.mat','Params','i_LF','Out_LF','i_LFendog','Out_LFendog');
    Params.jobfindingrate_exog=tmp.Params.jobfindingrate_exog; Params.eta=tmp.Params.eta; Params.wstar_guess=tmp.Params.wstar_guess;
    i_LF=tmp.i_LF; Out_LF=tmp.Out_LF; i_LFendog=tmp.i_LFendog; Out_LFendog=tmp.Out_LFendog;
    clear tmp
end

%% doPart(2): Table 2, firing taxes rebated as employment subsidies, exogenous job-finding rate
epsilon_vec=[0,0.67,2,4,8]; % months of wages: 0, 1, 3, 6, 12 (a half-quarter wage is 1.5 months)
ColHeaders_T2={'laissez-faire','1 month','3 months','6 months','12 months'};
RowNames_T2={'Job-finding rate','Layoff rate','Unemployment rate','Interest rate','Employment','Output','Capital','Assets','Consumption','S.D. consumption','Welfare: no search disutility','Welfare'};
if doPart(2)==1
    fprintf('doPart(2): Table 2, firing taxes rebated as employment subsidies (exogenous job-finding rate) \n')
    Params.endogsearch=0;
    Outs_T2=cell(1,length(epsilon_vec));
    Outs_T2{1}=Out_LF;
    i_T2=i_LF*ones(1,length(epsilon_vec));
    for ee=2:length(epsilon_vec)
        Params.epsilon=epsilon_vec(ee); Params.lambda=0;
        fprintf('doPart(2): epsilon=%.2f \n',Params.epsilon)
        tic;
        % Firing taxes lower the value of establishments and the interest rate; bracket around the previous solution
        [i_T2(ee),Outs_T2{ee}]=AV2001_SolveGE(i_T2(ee-1)-0.0015,i_T2(ee-1)+0.0005,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('Table 2, epsilon=%.2f',Params.epsilon));
        Params.wstar_guess=Outs_T2{ee}.wstar;
        fprintf('doPart(2): epsilon=%.2f done in %.1f minutes, i=%.6f \n',Params.epsilon,toc/60,i_T2(ee))
    end
    Params.epsilon=0;
    AV2001_PrintTable('Table 2: Firing taxes rebated as employment subsidies (exogenous job-finding rate)',ColHeaders_T2,Outs_T2,RowNames_T2,1,0);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table2.tex','w');
    AV2001_PrintTable('Table 2: Firing taxes rebated as employment subsidies (exogenous job-finding rate)',ColHeaders_T2,Outs_T2,RowNames_T2,fid,1);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart2.mat Outs_T2 i_T2 epsilon_vec
else
    load ./SavedOutput/AV2001_doPart2.mat Outs_T2 i_T2 epsilon_vec
end

%% doPart(3): Table 4, one-period unemployment insurance financed by payroll tax, exogenous job-finding rate
lambda_vec=[0,0.67,2,4,8];
ColHeaders_T4={'laissez-faire','1 month','3 months','6 months','12 months'};
RowNames_T4={'Job-finding rate','Layoff rate','Unemployment rate','Employment','Output','Capital','Consumption','S.D. consumption','Welfare: no search disutility','Welfare'};
if doPart(3)==1
    fprintf('doPart(3): Table 4, unemployment insurance with duration of benefits equal to one model period (exogenous job-finding rate) \n')
    Params.endogsearch=0;
    Outs_T4=cell(1,length(lambda_vec));
    Outs_T4{1}=Out_LF;
    i_T4=i_LF*ones(1,length(lambda_vec));
    for ll=2:length(lambda_vec)
        Params.epsilon=0; Params.lambda=lambda_vec(ll);
        fprintf('doPart(3): lambda=%.2f \n',Params.lambda)
        tic;
        [i_T4(ll),Outs_T4{ll}]=AV2001_SolveGE(i_T4(ll-1)-0.001,i_T4(ll-1)+0.001,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('Table 4, lambda=%.2f',Params.lambda));
        Params.wstar_guess=Outs_T4{ll}.wstar;
        fprintf('doPart(3): lambda=%.2f done in %.1f minutes, i=%.6f \n',Params.lambda,toc/60,i_T4(ll))
    end
    Params.lambda=0;
    AV2001_PrintTable('Table 4: Unemployment insurance with duration of benefits equal to one model period (exogenous job-finding rate)',ColHeaders_T4,Outs_T4,RowNames_T4,1,0);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table4.tex','w');
    AV2001_PrintTable('Table 4: Unemployment insurance with duration of benefits equal to one model period (exogenous job-finding rate)',ColHeaders_T4,Outs_T4,RowNames_T4,fid,1);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart3.mat Outs_T4 i_T4 lambda_vec
else
    load ./SavedOutput/AV2001_doPart3.mat Outs_T4 i_T4 lambda_vec
end

%% doPart(4): Table 1, severance payments (lambda=epsilon), endogenous job-finding rate
% Columns 2-6 of the paper's Table 1; column 1 (the calibrated U.S. economy with
% its UI system) is doPart(7). The laissez-faire column is the endogenous
% search economy of doPart(1a).
sev_vec=[0,0.67,2,4,8]; % lambda=epsilon, months of wages: 0, 1, 3, 6, 12
ColHeaders_T1={'laissez-faire','1 month','3 months','6 months','12 months'};
RowNames_T1={'Job-finding rate','Layoff rate','Unemployment rate','Employment','Output','Capital','Assets','Consumption','S.D. consumption','Non-search time','Welfare'};
if doPart(4)==1
    fprintf('doPart(4): Table 1, severance payments (endogenous job-finding rate) \n')
    Params.endogsearch=1;
    Outs_T1=cell(1,length(sev_vec));
    Outs_T1{1}=Out_LFendog;
    i_T1=i_LFendog*ones(1,length(sev_vec));
    for ss=2:length(sev_vec)
        Params.epsilon=sev_vec(ss); Params.lambda=sev_vec(ss);
        fprintf('doPart(4): lambda=epsilon=%.2f \n',Params.lambda)
        tic;
        [i_T1(ss),Outs_T1{ss}]=AV2001_SolveGE(i_T1(ss-1)-0.0015,i_T1(ss-1)+0.0005,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('Table 1, lambda=epsilon=%.2f',Params.lambda));
        Params.wstar_guess=Outs_T1{ss}.wstar;
        fprintf('doPart(4): lambda=epsilon=%.2f done in %.1f minutes, i=%.6f (annual %.3f%%), job-finding rate %.3f, u=%.4f \n',Params.lambda,toc/60,i_T1(ss),100*Outs_T1{ss}.annualinterestrate,Outs_T1{ss}.jobfindingrate,Outs_T1{ss}.unemploymentrate)
    end
    Params.epsilon=0; Params.lambda=0;
    AV2001_PrintTable('Table 1: Severance payments (endogenous job-finding rate; U.S. column is doPart(7))',ColHeaders_T1,Outs_T1,RowNames_T1,1,0);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table1.tex','w');
    AV2001_PrintTable('Table 1: Severance payments (endogenous job-finding rate; U.S. column is doPart(7))',ColHeaders_T1,Outs_T1,RowNames_T1,fid,1);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart4.mat Outs_T1 i_T1 sev_vec
else
    load ./SavedOutput/AV2001_doPart4.mat Outs_T1 i_T1 sev_vec
end

%% doPart(5): Table 3, firing taxes rebated as employment subsidies, endogenous job-finding rate
% As Table 2 (epsilon on, lambda=0, wstar<w is the employment subsidy) but with search chosen; the laissez-faire
% column is the endogenous search economy of doPart(1a).
ColHeaders_T3={'laissez-faire','1 month','3 months','6 months','12 months'};
RowNames_T3={'Job-finding rate','Layoff rate','Unemployment rate','Employment','Output','Capital','Assets','Consumption','S.D. consumption','Non-search time','Welfare'};
if doPart(5)==1
    fprintf('doPart(5): Table 3, firing taxes rebated as employment subsidies (endogenous job-finding rate) \n')
    Params.endogsearch=1;
    Outs_T3=cell(1,length(epsilon_vec));
    Outs_T3{1}=Out_LFendog;
    i_T3=i_LFendog*ones(1,length(epsilon_vec));
    for ee=2:length(epsilon_vec)
        Params.epsilon=epsilon_vec(ee); Params.lambda=0;
        fprintf('doPart(5): epsilon=%.2f \n',Params.epsilon)
        tic;
        [i_T3(ee),Outs_T3{ee}]=AV2001_SolveGE(i_T3(ee-1)-0.0015,i_T3(ee-1)+0.0005,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('Table 3, epsilon=%.2f',Params.epsilon));
        Params.wstar_guess=Outs_T3{ee}.wstar;
        fprintf('doPart(5): epsilon=%.2f done in %.1f minutes, i=%.6f (annual %.3f%%), job-finding rate %.3f, u=%.4f \n',Params.epsilon,toc/60,i_T3(ee),100*Outs_T3{ee}.annualinterestrate,Outs_T3{ee}.jobfindingrate,Outs_T3{ee}.unemploymentrate)
    end
    Params.epsilon=0;
    AV2001_PrintTable('Table 3: Firing taxes rebated as employment subsidies (endogenous job-finding rate)',ColHeaders_T3,Outs_T3,RowNames_T3,1,0);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table3.tex','w');
    AV2001_PrintTable('Table 3: Firing taxes rebated as employment subsidies (endogenous job-finding rate)',ColHeaders_T3,Outs_T3,RowNames_T3,fid,1);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart5.mat Outs_T3 i_T3
else
    load ./SavedOutput/AV2001_doPart5.mat Outs_T3 i_T3
end

%% doPart(6): Table 5, one-period unemployment insurance financed by payroll tax, endogenous job-finding rate
% As Table 4 (lambda on, epsilon=0, w<wstar is the payroll tax) but with search chosen.
ColHeaders_T5={'laissez-faire','1 month','3 months','6 months','12 months'};
RowNames_T5={'Job-finding rate','Layoff rate','Unemployment rate','Employment','Output','Capital','Consumption','S.D. consumption','Non-search time','Welfare'};
if doPart(6)==1
    fprintf('doPart(6): Table 5, unemployment insurance with duration of benefits equal to one model period (endogenous job-finding rate) \n')
    Params.endogsearch=1;
    Outs_T5=cell(1,length(lambda_vec));
    Outs_T5{1}=Out_LFendog;
    i_T5=i_LFendog*ones(1,length(lambda_vec));
    for ll=2:length(lambda_vec)
        Params.epsilon=0; Params.lambda=lambda_vec(ll);
        fprintf('doPart(6): lambda=%.2f \n',Params.lambda)
        tic;
        [i_T5(ll),Outs_T5{ll}]=AV2001_SolveGE(i_T5(ll-1)-0.001,i_T5(ll-1)+0.001,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('Table 5, lambda=%.2f',Params.lambda));
        Params.wstar_guess=Outs_T5{ll}.wstar;
        fprintf('doPart(6): lambda=%.2f done in %.1f minutes, i=%.6f (annual %.3f%%), job-finding rate %.3f, u=%.4f \n',Params.lambda,toc/60,i_T5(ll),100*Outs_T5{ll}.annualinterestrate,Outs_T5{ll}.jobfindingrate,Outs_T5{ll}.unemploymentrate)
    end
    Params.lambda=0;
    AV2001_PrintTable('Table 5: Unemployment insurance with duration of benefits equal to one model period (endogenous job-finding rate)',ColHeaders_T5,Outs_T5,RowNames_T5,1,0);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table5.tex','w');
    AV2001_PrintTable('Table 5: Unemployment insurance with duration of benefits equal to one model period (endogenous job-finding rate)',ColHeaders_T5,Outs_T5,RowNames_T5,fid,1);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart6.mat Outs_T5 i_T5
else
    load ./SavedOutput/AV2001_doPart6.mat Outs_T5 i_T5
end

%% doPart(7): the calibrated U.S. economy, first column of Table 1
% Section 4 of AV2001: a stylised unemployment insurance system. A fired worker becomes eligible with probability 0.5
% (half of laid-off workers do not qualify or do not apply); eligible unemployed receive 66% of the wage every period
% and keep eligibility with probability uiphi; benefits are financed by firing taxes (60% of them, the experience-rated
% part) and a payroll tax; no severance payments (lambda=0). Search is endogenous.
% uiphi and epsilon are calibrated inside the model: uiphi so that agents spend about 9 weeks collecting benefits
% while unemployed (periods of benefits per spell of eligibility, stock over inflow; a period is 6.5 weeks), and
% epsilon so that firing taxes pay 60% of benefits. Both are updated between general equilibrium solves:
%  - weeks: an eligible spell lasts D=1/(1-c*uiphi) periods, c the probability an eligible agent is still unemployed
%    (and not retired) next period; c is read off the equilibrium as (1-1/D)/uiphi and uiphi set to (1-1/D*)/c;
%  - share: the share is proportional to epsilon given the layoffs, so epsilon is scaled by 0.6/share.
% Starting values uiphi=0.5, epsilon=0.3.
UIweeks_target=9; UIshare_target=0.6;
if doPart(7)==1
    fprintf('doPart(7): the calibrated U.S. economy (unemployment insurance, endogenous search) \n')
    Params.endogsearch=1; Params.lambda=0;
    Params.uirr=Params.uirr_US; Params.uikappa=Params.uikappa_US;
    Params.uiphi=0.5; Params.epsilon=0.3;
    i_US=i_LFendog;
    tic;
    for it=1:10
        [i_US,Out_US]=AV2001_SolveGE(i_US-0.001,i_US+0.001,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('U.S. economy, calibration iteration %d',it));
        Params.wstar_guess=Out_US.wstar;
        fprintf('doPart(7): iteration %d: uiphi=%.4f epsilon=%.4f -> weeks collecting benefits %.2f (target %d), firing taxes pay %.3f of benefits (target %.2f); job-finding %.3f, u %.4f, layoff %.4f, i %.3f%% \n', ...
            it,Params.uiphi,Params.epsilon,Out_US.UIweeks,UIweeks_target,Out_US.firingtaxshareofUI,UIshare_target,Out_US.jobfindingrate,Out_US.unemploymentrate,Out_US.layoffrate,100*Out_US.annualinterestrate)
        if abs(Out_US.UIweeks-UIweeks_target)<0.05 && abs(Out_US.firingtaxshareofUI-UIshare_target)<0.005
            break
        end
        c=(1-1/Out_US.UIduration)/Params.uiphi;
        Dtarget=UIweeks_target/6.5;
        Params.uiphi=min(max((1-1/Dtarget)/c,0),0.99);
        Params.epsilon=Params.epsilon*UIshare_target/Out_US.firingtaxshareofUI;
    end
    if abs(Out_US.UIweeks-UIweeks_target)>=0.05 || abs(Out_US.firingtaxshareofUI-UIshare_target)>=0.005
        warning('doPart(7): the U.S. calibration did not converge in %d iterations (weeks %.2f, share %.3f)',it,Out_US.UIweeks,Out_US.firingtaxshareofUI)
    end
    uiphi_US=Params.uiphi; epsilon_US=Params.epsilon;
    fprintf('doPart(7): done in %.1f minutes. uiphi=%.4f, epsilon=%.4f; job-finding %.3f, unemployment %.4f, layoff %.4f (paper: 51.7%%, 5.7%%, 2.8%%) \n',toc/60,uiphi_US,epsilon_US,Out_US.jobfindingrate,Out_US.unemploymentrate,Out_US.layoffrate)
    Params.uirr=0; Params.uikappa=0; Params.uiphi=0; Params.epsilon=0;
    % Table 1 with the U.S. column, relative to the endogenous-search laissez-faire
    ColHeaders_T1US=[{'U.S.'},ColHeaders_T1];
    AV2001_PrintTable('Table 1: Severance payments (endogenous job-finding rate)',ColHeaders_T1US,[{Out_US},Outs_T1],RowNames_T1,1,0,Out_LFendog);
    fid=fopen('./SavedOutput/LatexInputs/AV2001_Table1.tex','w');
    AV2001_PrintTable('Table 1: Severance payments (endogenous job-finding rate)',ColHeaders_T1US,[{Out_US},Outs_T1],RowNames_T1,fid,1,Out_LFendog);
    fclose(fid);
    save ./SavedOutput/AV2001_doPart7.mat Out_US i_US uiphi_US epsilon_US
else
    load ./SavedOutput/AV2001_doPart7.mat Out_US i_US uiphi_US epsilon_US
end

%% doPart(8): diagnostics for the capital/interest-rate side
% Run of 2026-10-23: the labour-market rows of Tables 1, 2 and 4 replicate (layoff 2.8->2.0, unemployment
% 4.9->3.7 in Table 2, job-finding 61->69 in Table 1), but in Table 2 the paper's interest rate falls 4.0->3.8%
% and capital RISES (100.7) with twelve months of firing taxes, where here i moves 0.03pp and capital FALLS
% (97.7). Assets fall 6.7% in the paper, 3.5% here, and S.D. consumption rises in the paper, falls here: the
% paper's mechanism (firing taxes cut the value of establishments, savings exceed it, i falls) is weak here.
% Here assets per wage are unchanged (45.9 -> 45.9): establishment value falls about as much as the wage.
%  (a) per-wage decomposition of household assets S, capital K and establishment value A-K (eq 3.19) for the
%      Table 2 and Table 1 columns, with i, w, wstar, and the eq (3.19) check S=A.
%  (b) supply and demand of assets per wage, S/w and A/w, over the annual interest rate {3,3.5,4,4.5,5}% for
%      laissez-faire and twelve months of firing taxes (exogenous search): the slopes say whether the gap is the
%      demand shift (establishment value) or the supply elasticity (savings).
%  (c) fixed-i labour productivity, the paper's out-of-equilibrium exercise (Section 5.1.1): at the laissez-faire
%      i, firing taxes of 3, 6, 12 months lower output per worker by 0.1, 0.7, 4.3% in the paper.
%  (d) an employment grid on which attrition is exact: e_grid(k)=(1-rho)^-k, so (1-rho)e is the next point down
%      and an establishment can shrink by retirements alone. On the integer grid it cannot: it must either hire
%      rho*e workers or fire a fractional worker, so the inaction bands of the (s,S) rule are not representable.
%      (c) and the laissez-faire and twelve-month equilibria are re-solved on this grid (exogenous search; ~2000
%      employment points). The twelve-month equilibrium has ~860 reachable employment states; the household's sparse
%      Howards matrix then ran out of memory (2026-10-28), so above Grids.howardssparse_maxNhh states Howards is by indexing.
%  (c) and (d) also report S/w and A/w at the fixed i, and the change in the equilibrium interest rate they imply
%  with the supply and demand slopes from (b).
% Run of 2026-10-28 (a)-(c): S=A holds; S/w is very interest-elastic (+39% per point of annual i) against A/w
% (-12%); twelve months of firing taxes shift S/w up 1.4% and A/w down 1.2% at fixed i, implying i falls only
% ~0.05pp (paper 0.2pp); fixed-i output per worker -0.05/-0.4/-4.1% (paper -0.1/-0.7/-4.3), and -0.2/-0.7/-4.5%
% on the (d) grid. The firm side replicates; the gap is how little establishment value falls relative to w.
if doPart(8)==1
    fprintf('doPart(8): diagnostics for the capital/interest-rate side \n')
    tic;
    % (a)
    fprintf('doPart(8a): per-wage decomposition at the equilibria \n')
    fprintf('  %-22s %8s %8s %8s %8s %8s %8s %8s %9s %8s \n','column','i(ann)','w','wstar','S/w','K/w','(A-K)/w','A/w','S/A-1','Y/N')
    for T={{'Table 2',Outs_T2},{'Table 1',Outs_T1}}
        for k=1:5
            O=T{1}{2}{k};
            fprintf('  %-22s %7.3f%% %8.4f %8.4f %8.2f %8.2f %8.2f %8.2f %9.2e %8.4f \n',sprintf('%s col %d',T{1}{1},k),100*O.annualinterestrate,O.w,O.wstar,O.S/O.w,O.K/O.w,(O.A-O.K)/O.w,O.A/O.w,O.S/O.A-1,O.Ymeasured/O.N_hh)
        end
    end
    % (b)
    fprintf('doPart(8b): asset supply S/w and demand A/w over the interest rate (exogenous search) \n')
    iann_vec=[0.03,0.035,0.04,0.045,0.05];
    SD_b=zeros(2,length(iann_vec),2); % (policy, i, [S/w, A/w])
    Params.endogsearch=0; Params.lambda=0;
    eps_b=[0,8];
    for pp=1:2
        Params.epsilon=eps_b(pp);
        for ii=1:length(iann_vec)
            [~,O]=AV2001_StationaryEqm((1+iann_vec(ii))^(1/8)-1,Params,Grids,vfoptions,simoptions,0);
            SD_b(pp,ii,:)=[O.S/O.w,O.A/O.w];
            fprintf('  epsilon=%d, i=%.1f%%: S/w %.2f, A/w %.2f, (A-K)/w %.2f, K/w %.2f, w %.4f \n',eps_b(pp),100*iann_vec(ii),O.S/O.w,O.A/O.w,(O.A-O.K)/O.w,O.K/O.w,O.w)
        end
        fprintf('  epsilon=%d: d log(S/w)/d i = %.1f, d log(A/w)/d i = %.1f (per percentage point of annual i, in %%) \n',eps_b(pp), ...
            100*mean(diff(log(SD_b(pp,:,1))))/(100*0.005),100*mean(diff(log(SD_b(pp,:,2))))/(100*0.005))
    end
    Params.epsilon=0;
    % slopes of log(S/w) and log(A/w) per percentage point of annual i, from (b) at epsilon=0
    slopeS=mean(diff(log(SD_b(1,:,1))))/0.5; slopeA=mean(diff(log(SD_b(1,:,2))))/0.5;
    % (c)
    fprintf('doPart(8c): output per worker at the laissez-faire i (paper: -0.1, -0.7, -4.3%% at 3, 6, 12 months) \n')
    eps_c=[0,2,4,8];
    YN_c=zeros(1,length(eps_c)); SA_c=zeros(length(eps_c),2);
    for ee=1:length(eps_c)
        Params.epsilon=eps_c(ee);
        [~,O]=AV2001_StationaryEqm(i_LF,Params,Grids,vfoptions,simoptions,0);
        YN_c(ee)=O.Ymeasured/O.N_hh; SA_c(ee,:)=[O.S/O.w,O.A/O.w];
        fprintf('  epsilon=%d: output per worker %.2f (index), K/Y annual %.3f, S/w %.2f, A/w %.2f, implied change in the equilibrium i %+.3f pp \n',eps_c(ee),100*YN_c(ee)/YN_c(1),O.KYannual,SA_c(ee,1),SA_c(ee,2), ...
            (log(SA_c(ee,2)/SA_c(1,2))-log(SA_c(ee,1)/SA_c(1,1)))/(slopeS-slopeA))
    end
    Params.epsilon=0;
    % (d) attrition-exact employment grid
    Grids_geo=Grids;
    Kgeo=ceil(log(Grids.e_grid(end))/(-log(1-Params.rho)));
    Grids_geo.e_grid=[0;(1-Params.rho).^(-(0:1:Kgeo)')];
    Grids_geo.n_e=length(Grids_geo.e_grid);
    fprintf('doPart(8d): employment grid with exact attrition, %d points on [0,%.0f] (step %.2f workers at 116) \n',Grids_geo.n_e,Grids_geo.e_grid(end),116*(1/(1-Params.rho)-1))
    YN_d=zeros(1,length(eps_c)); SA_d=zeros(length(eps_c),2);
    for ee=1:length(eps_c)
        Params.epsilon=eps_c(ee);
        [~,O]=AV2001_StationaryEqm(i_LF,Params,Grids_geo,vfoptions,simoptions,0);
        YN_d(ee)=O.Ymeasured/O.N_hh; SA_d(ee,:)=[O.S/O.w,O.A/O.w];
        fprintf('  fixed i: epsilon=%d: output per worker %.2f (index), layoff %.4f, N_R=%d, S/w %.2f, A/w %.2f, implied change in the equilibrium i %+.3f pp \n',eps_c(ee),100*YN_d(ee)/YN_d(1),O.layoffrate,O.N_R,SA_d(ee,1),SA_d(ee,2), ...
            (log(SA_d(ee,2)/SA_d(1,2))-log(SA_d(ee,1)/SA_d(1,1)))/(slopeS-slopeA))
    end
    Params.epsilon=0;
    Outs_geo=cell(1,2); i_geo=zeros(1,2);
    for pp=1:2
        Params.epsilon=eps_b(pp);
        [i_geo(pp),Outs_geo{pp}]=AV2001_SolveGE(i_LF-0.0015,i_LF+0.0005,Params,Grids_geo,vfoptions,simoptions,verbose,fzerooptions_i,sprintf('attrition-exact e-grid, epsilon=%d',eps_b(pp)));
        Params.wstar_guess=Outs_geo{pp}.wstar;
    end
    Params.epsilon=0; Params.wstar_guess=Out_LF.wstar;
    O0=Outs_geo{1}; O8=Outs_geo{2};
    fprintf('  general equilibrium on this grid, laissez-faire vs twelve months (paper Table 2: i 4.0 -> 3.8%%, capital 100.7, assets 93.3, output 98.9, u 4.9 -> 3.7%%, layoff 2.9 -> 2.0%%) \n')
    fprintf('    i %.3f%% -> %.3f%%, capital %.1f, assets %.1f, output %.1f, u %.2f%% -> %.2f%%, layoff %.2f%% -> %.2f%%, S/w %.2f -> %.2f, (A-K)/w %.2f -> %.2f, N_R %d -> %d \n',100*O0.annualinterestrate,100*O8.annualinterestrate, ...
        100*O8.K/O0.K,100*O8.S/O0.S,100*O8.Ymeasured/O0.Ymeasured,100*O0.unemploymentrate,100*O8.unemploymentrate,100*O0.layoffrate,100*O8.layoffrate,O0.S/O0.w,O8.S/O8.w,(O0.A-O0.K)/O0.w,(O8.A-O8.K)/O8.w,O0.N_R,O8.N_R)
    fprintf('doPart(8): done in %.1f minutes \n',toc/60)
    save ./SavedOutput/AV2001_doPart8.mat SD_b iann_vec eps_b YN_c YN_d SA_c SA_d eps_c slopeS slopeA Outs_geo i_geo
end

%% doPart(9): LaTeX tables for the write-up that the table doParts do not write; solves nothing
if doPart(9)==1
    fprintf('doPart(9): LaTeX tables for the write-up \n')
    AV2001_LatexCalibration(Out_LFendog,Params,'./SavedOutput/LatexInputs/AV2001_Calibration.tex');
    AV2001_LatexCalibrationUS(Out_US,uiphi_US,epsilon_US,Params,'./SavedOutput/LatexInputs/AV2001_CalibrationUS.tex');
end

diary off
