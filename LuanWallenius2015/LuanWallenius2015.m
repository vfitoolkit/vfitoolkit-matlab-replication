% Laun & Wallenius (2015) - A Life Cycle Model of Health and Retirement: Swedish Pension Reform
% Full replication of LW2015, does all the model exercises from that paper.
%
%
% NOTE: the wage profiles (Params.w_low_i, Params.w_peak_premium_i) are our best-guess interpretations
% of LW2015; we cannot they are correct from the paper or any replication source. See inline comments for details. 
% The health-transition Markov chain (pi_z_J) is also a best-guess interpretation, but confident in it as explained in inline comments.

%% Sizes
Params.agejshifter=24;                          % first model period = age 25
N_j=80-Params.agejshifter;                      % 56 model periods, ages 25..80
Params.J=N_j;
Params.N_j=N_j;
Params.agej=1:N_j;
Params.age=Params.agej+Params.agejshifter;

n_d=5;                                          % action codes (see aprimeFn header)
n_k=100;                                        % assets
n_hist=1237;                                    % experienceasset: history index
n_a=[n_k,n_hist];                               % hist_idx last (experienceasset)
n_z=5;                                          % health state {1..5}: very bad..very good
N_i=4;                                          % PType: NCL, NCH, CL, CH
Names_i={'NCL','NCH','CL','CH'};


doPart=[0,0,0,1,1,1];
% doPart(1): pre-reform baseline
% doPart(2): regular pension reform only (post-reform PB + DI; OPB unchanged), partial eqm
% doPart(3): regular pension reform only, general eqm
% doPart(4): full reform (also post-reform OPB DC scheme), partial eqm
% doPart(5): full reform, general eqm
% doPart(6): welfare CEV (full GE vs pre-reform) -- closed form for log utility

% Make sure the subfolders to save output exist
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

%% hist_idx layout -- READ THIS BEFORE CHANGING ANYTHING ABOVE
% hist_idx is the experienceasset and encodes the agent's "employment history"
% as a single integer in {1..1237}. It tracks three pieces of information:
%   swa = first stop-work age      (model period; 0 = never stopped)
%   dia = first DI-claim age       (model period; 0 = never claimed DI)
%   pba = first PB-claim age       (model period; 0 = never claimed PB)
% All three transitions are absorbing (work is absorbing once stopped, DI and
% PB claims cannot be reversed). LW2016's replication code uses the same
% scheme; LW2015 describes it as "each possible employment history" (p.129).
%
% The unconstrained product {swa,dia,pba} has 57 x 41 x 21 = 49,077 entries,
% but only ~2.5% of those are reachable because:
%   (i)  if DI is claimed then swa = dia (you can't work while on DI), and
%   (ii) DI and PB are mutually exclusive pre-65 (auto-conversion at 65 is
%        mechanical, not a separate claim).
% We therefore compress the reachable subset into 1237 integer indices, split
% into two contiguous blocks:
%   Block A (no DI, indices 1..1197):  dia = 0; (swa, pba) cartesian.
%       swa ranges over {0,1..56} (57 values; 0 = never stopped).
%       pba ranges over {0,37..56} (21 values; 0 = never, 37 corresponds to
%       age 61, the earliest claim age).
%       Encoding: hist_idx = swa*21 + pba_idx, with pba_idx = 1 if pba=0,
%       else pba_idx = pba - 35.
%   Block B (on DI, indices 1198..1237):  dia in {1..40} (DI claimed at model
%       period 1..40 i.e. real age 25..64). swa = dia (forced), pba = 0.
%       Encoding: hist_idx = 1197 + dia.
% Inversion (decode) and encoding are inlined in both LuanWallenius2015_aprimeFn.m
% and LuanWallenius2015_ReturnFn.m (look for `n_no_di = 1197`). The agent's
% decision variable d in {1..5} dispatches the transition: see the aprimeFn
% header for the d -> (swa', dia', pba') rules.
%
% Benefit calculations are NOT precomputed into tables. They are evaluated on
% the fly by LuanWallenius2015_PB.m, _DI.m, _OPB.m at every ReturnFn call,
% which in turn call LuanWallenius2015_AP.m to evaluate average pension points
% from the (swa, dia) history under the literal top-15-best rule (with
% 3-year-pre-DI projection inside AP when dia > 0).

%% Parameters
Params.r=0.03;
Params.beta=1/(1+Params.r);

% Disutility from work, by health state (h=1 verybad..h=5 verygood)
Params.b1=3.5;
Params.b2=3.0;
Params.b3=2.5;
Params.b4=2.0;
Params.b5=1.5;

% Labor
Params.lbar=1/3;

% Taxes (paper Section 3)
Params.tau_l1=0.31;
Params.tau_l2=0.51;
Params.tau_ss=0.23;
Params.tau_c=0.25;
Params.tax_thresh=6;                            % ~219300 SEK / 36300 SEK = 6 BA

% Pension policy ages (Swedish pre-reform)
Params.agej_pb_first=61-Params.agejshifter;     % first PB-claim model age (=37)
Params.agej_di_max =64-Params.agejshifter;      % last new DI claim (=40)
Params.agej_65     =65-Params.agejshifter;      % DI->PB transfer (=41)

% Pension-points cap (7.5 BA earnings cap; BA subtracted)
Params.ppmax=6.5;

% Lump-sum transfer (paper: ~15000 SEK = ~0.4 BA)
Params.T=0.4;
% Pre-reform: T held at the paper's reported equilibrium value; we skip the GE budget-balance loop here (we do solve for it for the post-reform run below).

% Reform regime:
%   0 = pre-reform (Swedish PAYG-DB pension; AP-based PB and DI;
%       pre-reform OPB).
%   1 = regular pension reform only (NDC PB; 64%-of-pre-DI-avg DI;
%       OPB unchanged) -- maps to LW2015 Fig 6's "regular pension reform
%       only" decomposition line.
%   2 = full reform (also DC OPB: 4.5%/30% contributions, benefit =
%       capital_OPB / annuity_factor) -- maps to LW2015 Fig 5's headline
%       full-reform line.
% OPB_DI top-up is held at the pre-reform formula in all regimes (paper
% does not specify a post-reform OPB_DI rule).
Params.regime=0;

% Age weights (paper assumes equal mass per age)
Params.mewj=ones(1,N_j)/N_j;

%% Wage profiles, per-PType (w_low, w_peak_premium) scalars
% Hump-shaped quadratic peaking at age 50, in BA units per year of full-time
% work. Closed form is inlined where used (ReturnFn for current-period
% earnings; AP, OPB, OPB_DI for per-age wages):
%   w(real_age) = w_low_i + w_peak_premium_i * max(0, 1 - ((real_age-50)/30)^2)
%
% Values eyeballed from LW2015 Fig. 1 (LINDA 1968-1997 regression -- the
% regression coefficients are not published in the paper or any working-
% paper appendix we can locate).
% NOTE: we CANNOT borrow from the LW2016 sequel replication code (RED
% 14-161, Compusswe/Version 1/w.m) because LW2016 publishes only TWO
% regressed wage profiles (HS, CL); the within-education fixed-effect
% dimension that LW2015 uses for the NCL/NCH and CL/CH split is gone in
% LW2016. Fitting all four LW2015 profiles requires either the original
% regression output or LINDA access to re-estimate.
Params.w_low_i         =[2.0; 3.0; 3.0; 5.0];   % NCL, NCH, CL, CH
Params.w_peak_premium_i=[1.5; 2.0; 2.5; 5.0];

%% Occupational pension: blue vs white collar by PType
Params.collar_i=[0; 0; 1; 1];                   % 0 = blue-collar; 1 = white-collar

%% PType mass (~17% college; split 50/50 low/high within each education group)
Params.ptypedist=[0.415; 0.415; 0.085; 0.085];
PTypeDistParamNames={'ptypedist'};

%% Figure 1: life-cycle earnings by PType (LW2015 Fig 1)
% Reproduces the four hump-shaped earnings profiles in BA/year.
fig1_ages=25:70;
fig1_w=zeros(4,length(fig1_ages));
for ii=1:4
    hump_ii=max(0,1-((fig1_ages-50)/30).^2);
    fig1_w(ii,:)=Params.w_low_i(ii)+Params.w_peak_premium_i(ii)*hump_ii;
end
fig1=figure(1);
plot(fig1_ages,fig1_w(1,:),'-b','LineWidth',1.5); hold on;
plot(fig1_ages,fig1_w(2,:),'--b','LineWidth',1.5);
plot(fig1_ages,fig1_w(3,:),'-r','LineWidth',1.5);
plot(fig1_ages,fig1_w(4,:),'--r','LineWidth',1.5);
xlabel('Age');
ylabel('Yearly Earnings in Base Amounts (BA)');
legend({'Non-College Low','Non-College High','College Low','College High'},'Location','northwest');
title('Figure 1: Life cycle earnings of Swedish males by skill type');
grid on;
saveas(fig1,'SavedOutput/Graphs/LuanWallenius2015_Figure1.png');

%% Table 4: calibrated parameter values (LW2015 Table 4)
fid=fopen('SavedOutput/LatexInputs/LuanWallenius2015_Table4.tex','w');
fprintf(fid,'\\begin{tabular}{lll}\n\\hline\n');
fprintf(fid,'Parameter & Value & Explanation \\\\\n\\hline\n');
fprintf(fid,'\\multicolumn{3}{l}{\\textit{Policy parameters}} \\\\\n');
fprintf(fid,'$\\tau_{l1}$ & %.2f & Tax on labor income below threshold \\\\\n',Params.tau_l1);
fprintf(fid,'$\\tau_{l2}$ & %.2f & Tax on labor income above threshold \\\\\n',Params.tau_l2);
fprintf(fid,'$\\tau_{ss}$ & %.2f & Social security tax \\\\\n',Params.tau_ss);
fprintf(fid,'$\\tau_c$ & %.2f & Consumption tax \\\\\n',Params.tau_c);
fprintf(fid,'$\\bar{Y}_{\\text{thr}}$ & %g BA & Tax-bracket threshold \\\\\n',Params.tax_thresh);
fprintf(fid,'$T$ & %.2f BA & Lump-sum transfer (pre-reform) \\\\\n',Params.T);
fprintf(fid,'\\hline\n');
fprintf(fid,'\\multicolumn{3}{l}{\\textit{Utility parameters}} \\\\\n');
fprintf(fid,'$\\beta$ & %.4f & Discount factor \\\\\n',Params.beta);
fprintf(fid,'$r$ & %.2f & Interest rate \\\\\n',Params.r);
fprintf(fid,'$b_5$ & %.1f & Disutility from work when health very good \\\\\n',Params.b5);
fprintf(fid,'$b_4$ & %.1f & Disutility from work when health good \\\\\n',Params.b4);
fprintf(fid,'$b_3$ & %.1f & Disutility from work when health fair \\\\\n',Params.b3);
fprintf(fid,'$b_2$ & %.1f & Disutility from work when health bad \\\\\n',Params.b2);
fprintf(fid,'$b_1$ & %.1f & Disutility from work when health very bad \\\\\n',Params.b1);
fprintf(fid,'$\\bar{l}$ & 1/3 & Labor supply (full-time fraction) \\\\\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\multicolumn{3}{l}{\\textit{Health parameters}} \\\\\n');
fprintf(fid,'$\\varepsilon^l$ & 1 & Health drop from small negative shock \\\\\n');
fprintf(fid,'$\\varepsilon^h$ & 3 & Health drop from large negative shock \\\\\n');
fprintf(fid,'$\\varepsilon^i$ & 1 & Health gain from positive shock \\\\\n');
fprintf(fid,'$p^l_a$ & 0.1 $\\to$ 0.6 & Prob. small neg. shock (linear in age) \\\\\n');
fprintf(fid,'$p^l_1$ & 0.6 & Prob. small neg. shock at $h=1$ \\\\\n');
fprintf(fid,'$p^h$ & 0.01 & Prob. large neg. shock \\\\\n');
fprintf(fid,'$p^h_1$ & 0.1 & Prob. large neg. shock at $h=1$ \\\\\n');
fprintf(fid,'$p^i_{a,5}$ & 1.0 $\\to$ 0.5 & Prob. positive shock at $h=5$ \\\\\n');
fprintf(fid,'$p^i_{a,4}$ & 0.9 $\\to$ 0.5 & Prob. positive shock at $h=4$ \\\\\n');
fprintf(fid,'$p^i_{a,3}$ & 0.8 $\\to$ 0.5 & Prob. positive shock at $h=3$ \\\\\n');
fprintf(fid,'$p^i_{a,2}$ & 0.7 $\\to$ 0.5 & Prob. positive shock at $h=2$ \\\\\n');
fprintf(fid,'$p^i_{a,1}$ & 0.6 $\\to$ 0.5 & Prob. positive shock at $h=1$ \\\\\n');
fprintf(fid,'\\hline\n\\end{tabular}\n');
fclose(fid);

%% Grids
k_grid=linspace(0,50,n_k)';                     % 0 to 50 BA
hist_grid=(1:n_hist)';                          % experienceasset values are the indices themselves
a_grid=[k_grid; hist_grid];                     % hist_grid last (experienceasset)

z_grid=(1:5)';                                  % health 1..5

d_grid=(1:n_d)';                                % action codes 1..5

%% Age-dependent health transition matrix pi_z_J: 5 x 5 x N_j
% Three exogenous shocks per LW2015 Table 4:
%   small negative (Delta h = -1), prob p_l, age-varying (and elevated at h=1)
%   large negative (Delta h = -3), prob p_h, mostly age-invariant
%   positive       (Delta h = +1), prob p_i, age- and health-varying
% Reading A (independent Bernoulli draws): see comment block immediately
% following this loop.
pi_z_J=zeros(5,5,N_j);
p_pos_young=[0.6 0.7 0.8 0.9 1.0];
for jj=1:N_j
    age_frac=(jj-1)/(N_j-1);
    for hh=1:5
        if hh==1
            p_l=0.6;                                            % LW2015 p^l_1 override
            p_h=0.1;                                            % LW2015 p^h_1 override
        else
            p_l=0.1+0.5*age_frac;                               % LW2015 p^l_a default
            p_h=0.01;                                           % LW2015 p^h default
        end
        p_i=p_pos_young(hh)-(p_pos_young(hh)-0.5)*age_frac;     % LW2015 p^i_{a,h}

        % Enumerate the 2^3 = 8 joint outcomes
        for s_l=0:1
            for s_h=0:1
                for s_i=0:1
                    prob = (s_l*p_l + (1-s_l)*(1-p_l)) ...
                         * (s_h*p_h + (1-s_h)*(1-p_h)) ...
                         * (s_i*p_i + (1-s_i)*(1-p_i));
                    delta = -s_l - 3*s_h + s_i;
                    h_prime = max(1, min(5, hh + delta));
                    pi_z_J(hh, h_prime, jj) = pi_z_J(hh, h_prime, jj) + prob;
                end
            end
        end
    end
end

% READING A FOR THE HEALTH-TRANSITION TABLE -- WHY AND HOW
% LW2015 Table 4 lists three shock probabilities (p^l small negative, p^h
% large negative, p^i positive) and does NOT specify how to aggregate them
% into a one-step Markov transition. At several (h, age) cells the three
% sum to more than 1 if read as mutually-exclusive event probabilities:
%   young h=4:  p^l=0.10, p^h=0.01, p^i=0.90  -> sum 1.01
%   old h=4:    p^l=0.60, p^h=0.01, p^i=0.50  -> sum 1.11
%   old h=3:    p^l=0.60, p^h=0.01, p^i=0.50  -> sum 1.11
% So the natural mutually-exclusive reading needs ad-hoc renormalization.
%
% Reading A (above): the three shocks are independent Bernoulli draws per
% period. The joint distribution has 8 outcomes (each combination of the
% three shocks firing or not), the net change in h is the algebraic sum of
% the firing-shock values, and we cap to [1, 5] at the end. This:
%   (i)   uses every Table-4 probability at face value, no clipping,
%   (ii)  is consistent with how the LW2016 replication code combines
%         negative shocks with the (endogenous, success-probability-driven)
%         positive event -- as independent draws whose joint probabilities
%         multiply,
%   (iii) gives a richer five-state Markov chain than Reading B would: at
%         healthy states the positive shock partially "shields" the agent
%         from a coincident negative shock (e.g. at h=5 young with p^i=1,
%         a coincident small negative shock nets to 0 and stays at h=5;
%         coincident large negative nets to -2 and lands at h=3 instead of
%         h=2), and at h=1 the elevated p^l keeps the agent stuck through
%         the cap rather than moving.
%
% Important caveat: LW2016 cannot directly confirm this interpretation --
% LW2016 dropped the exogenous positive shock entirely and replaced it with
% endogenous health investment, so their code structure doesn't include the
% three-event aggregation LW2015 leaves ambiguous.
%
% Empirical validation: LW2015 p.132 reports the model's health distribution
% at real age 70 as (42, 25, 13, 9, 11) % across (very good, good, fair,
% bad, very bad). Health is purely exogenous in LW2015, so this distribution
% is just the initial dist P(h=5)=1 propagated through 45 periods of pi_z_J
% -- independent of all other model parameters. See LWtest.m. Propagating
% under Reading A gives an L1 distance of 5.6 percentage points to the
% paper's reported values; Reading B gives 34.4. Reading A is therefore not
% just a defensible interpretation but is empirically consistent with the
% paper's own model output, while Reading B is decisively ruled out.

%% Experienceasset setup (hist_idx)
vfoptions.experienceasset=1;
vfoptions.aprimeFn=@(d,hist_idx,agej,agej_pb_first,agej_di_max) LuanWallenius2015_aprimeFn(d,hist_idx,agej,agej_pb_first,agej_di_max);
simoptions.experienceasset=1;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

%% Divide-and-conquer and grid interpolation layer
vfoptions.divideandconquer=1; % This both cuts memory use and speeds it up
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Return function
DiscountFactorParamNames={'beta'};
ReturnFn=@(d,kprime,k,hist_idx,h,agej,r,b1,b2,b3,b4,b5,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_pb_first,agej_di_max,agej_65,regime) ...
    LuanWallenius2015_ReturnFn(d,kprime,k,hist_idx,h,agej,r,b1,b2,b3,b4,b5,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_pb_first,agej_di_max,agej_65,regime);

%% Initial distribution: k=0, hist_idx=1 (swa=0,dia=0,pba=0), h=verygood (h_idx=5)
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,1,5)=1;

AgeWeightParamNames={'mewj'};

%% Functions to evaluate (inline)
% hist_idx -> on_DI, on_PB inferred via the same arithmetic used inside ReturnFn
FnsToEvaluate.working=@(d,kprime,k,hist_idx,h) (d==1 || d==5);
FnsToEvaluate.onDI=@(d,kprime,k,hist_idx,h) (hist_idx>1197);
FnsToEvaluate.onPB=@(d,kprime,k,hist_idx,h) ((hist_idx<=1197) && (mod(hist_idx-1,21)+1>1));
FnsToEvaluate.assets=@(d,kprime,k,hist_idx,h) k;
FnsToEvaluate.health=@(d,kprime,k,hist_idx,h) h;
FnsToEvaluate.histidx=@(d,kprime,k,hist_idx,h) hist_idx;
% Retirement age: decoded swa (with swa=0 -> "never stopped, count as age N_j+agejshifter")
FnsToEvaluate.retire_age=@(d,kprime,k,hist_idx,h,agejshifter,N_j) ...
    ((hist_idx>1197)*(hist_idx-1197+agejshifter)) + ...
    ((hist_idx<=1197)*((floor((hist_idx-1)/21)>0)*(floor((hist_idx-1)/21)+agejshifter) + (floor((hist_idx-1)/21)==0)*(N_j+agejshifter)));
% DI claim age: dia+agejshifter if on DI, 0 otherwise. Divide by mean(onDI) at j=N_j for the conditional mean.
FnsToEvaluate.di_claim_age=@(d,kprime,k,hist_idx,h,agejshifter) ...
    (hist_idx>1197)*(hist_idx-1197+agejshifter);

%% Pre-reform
if doPart(1)==1
    % Solve value function
    vfoptions.verbose=1;
    tic;
    [V_pre,Policy_pre]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    vftime=toc
    vfoptions.verbose=0;

    % Stationary distribution
    StationaryDist_pre=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_pre,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);

    %% All-stats and life-cycle profiles
    AllStats_pre=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_pre,Policy_pre,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_pre=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_pre,Policy_pre,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    save('./SavedOutput/LuanWallenius2015_results.mat','V_pre','Policy_pre','StationaryDist_pre','AllStats_pre','AgeConditionalStats_pre','Params','-v7.3');
else
    load ./SavedOutput/LuanWallenius2015_results.mat
end
% Free up some memory
clear V_pre Policy_pre StationaryDist_pre

%% Figure 3: model employment distribution by age bin
bin_labels={'50-54','55-59','60-64','65-69','70-74'};

emp_age=100*AgeConditionalStats_pre.working.Mean;
bin_starts=[50 55 60 65 70]-Params.agejshifter;
model_emp=zeros(1,5);
for bb=1:5
    js=bin_starts(bb):(bin_starts(bb)+4);
    model_emp(bb)=mean(emp_age(js));
end

fig3=figure(3);
plot(1:5,model_emp,'-bo','LineWidth',1.5,'MarkerFaceColor','b');
set(gca,'XTick',1:5,'XTickLabel',bin_labels);
xlabel('Age');
ylabel('Fraction Employed (%)');
ylim([0 100]);
legend({'Model'},'Location','northeast');
title('Figure 3: Model employment distribution');
grid on;
saveas(fig3,'SavedOutput/Graphs/LuanWallenius2015_Figure3.png');

%% Figure 4: model DI incidence by age bin
di_bin_labels={'25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64'};

di_age=100*AgeConditionalStats_pre.onDI.Mean;
di_bin_starts=[25 30 35 40 45 50 55 60]-Params.agejshifter;
model_di=zeros(1,8);
for bb=1:8
    js=di_bin_starts(bb):(di_bin_starts(bb)+4);
    model_di(bb)=mean(di_age(js));
end

fig4=figure(4);
plot(1:8,model_di,'-bo','LineWidth',1.5,'MarkerFaceColor','b');
set(gca,'XTick',1:8,'XTickLabel',di_bin_labels);
xlabel('Age');
ylabel('Fraction on DI (%)');
ylim([0 25]);
legend({'Model'},'Location','northwest');
title('Figure 4: Model DI incidence');
grid on;
saveas(fig4,'SavedOutput/Graphs/LuanWallenius2015_Figure4.png');

%% ===================================================================
%% POST-REFORM EXPERIMENT (PARTIAL EQUILIBRIUM)
%% ===================================================================
% Switch the regime flag to 1 (post-reform NDC) and resolve. T is held at
% the pre-reform value of 0.4 BA -- this is LW2015's "partial equilibrium"
% specification (Section 6: "ignore this general equilibrium aspect and
% only consider the partial equilibrium decision problem of agents"). PE
% gives a 2.7-year retirement-age increase in the paper vs 2.5 under full
% GE. Welfare CVs would differ more substantially and require the GE
% T-clearing loop; we don't compute welfare here.
%
% Only PB and DI are switched to their post-reform formulas (NDC pension
% and 64%-of-pre-DI-avg DI). OPB and OPB_DI are held at pre-reform values.
% This corresponds to LW2015's "regular pension reform only" decomposition
% (Fig 6), which gives a slightly larger labor-supply response than the
% full reform.


% Switch regime
Params.regime=1;

%% Reform, partial eqm
if doPart(2)==1
    % Solve post-reform value function
    vfoptions.verbose=1;
    tic;
    [V_post,Policy_post]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    vftime_post=toc
    vfoptions.verbose=0;

    %% Post-reform stationary distribution
    StationaryDist_post=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_post,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);

    %% Post-reform AllStats / AgeConditionalStats
    AllStats_post=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_post,Policy_post,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_post=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_post,Policy_post,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    save('./SavedOutput/LuanWallenius2015_results_partialeqm.mat','V_post','Policy_post','StationaryDist_post','AllStats_post','AgeConditionalStats_post','Params','-v7.3');
else
    load ./SavedOutput/LuanWallenius2015_results_partialeqm.mat
end
% Free up some memory
clear V_post Policy_post StationaryDist_post

%% Pre- vs regular-only PE employment by age bin
% (Used downstream by Figs 5/5b/6 and by the PE-summary block.)
emp_age_pre =100*AgeConditionalStats_pre.working.Mean;
emp_age_post=100*AgeConditionalStats_post.working.Mean;
model_emp_pre =zeros(1,5);
model_emp_post=zeros(1,5);
for bb=1:5
    js=bin_starts(bb):(bin_starts(bb)+4);
    model_emp_pre(bb) =mean(emp_age_pre(js));
    model_emp_post(bb)=mean(emp_age_post(js));
end

%% Quick summary: average retirement age proxy (1 - mean-employed share)
%  weighted by age. Not the paper's exact statistic but a usable comparison.
mean_emp_pre =mean(AgeConditionalStats_pre.working.Mean);
mean_emp_post=mean(AgeConditionalStats_post.working.Mean);
fprintf('\nPost-reform PE summary:\n');
fprintf('  Mean employment rate (pre):  %.3f\n',mean_emp_pre);
fprintf('  Mean employment rate (post): %.3f\n',mean_emp_post);
fprintf('  Delta employment (post-pre): %+.3f\n',mean_emp_post-mean_emp_pre);
fprintf('  Employment 65-69 (pre):  %.1f%%\n',model_emp_pre(4));
fprintf('  Employment 65-69 (post): %.1f%%\n',model_emp_post(4));


%% Reform, general eqm
% Same regime switch as the PE block above, but now T is a GE price solved by
% VFI Toolkit's HeteroAgentStationaryEqm_Case1_FHorz_PType. The GeneralEqmEqn
% is inline: at equilibrium the aggregate per-agent net surplus
%   (income_tax + ss_tax + consumption_tax) - (PB + DI + OPB + OPB_DI)
% equals T (the lump-sum transfer per agent). LW2015 p.134 reports T rising
% from ~0.4 BA (pre-reform) to ~0.6 BA (post-reform GE) as pension outlays fall.
if doPart(3)==1
    % Switch to post-reform; T is the GE price (initial guess = PE value)
    Params.regime=1;
    Params.T=0.4;

    % FnsToEvaluate: per-agent budget surplus (helper in LuanWallenius2015_BudgetFn.m)
    FnsToEvaluate_GE.surplus=@(d,kprime,k,hist_idx,h,agej,r,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_65,regime) ...
        LuanWallenius2015_BudgetFn(d,kprime,k,hist_idx,h,agej,r,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_65,regime);

    % GeneralEqmEqn (inline): aggregate surplus minus T = 0
    GeneralEqmEqns.budget=@(surplus,T) surplus-T;

    GEPriceParamNames={'T'};
    heteroagentoptions.verbose=1;

    [p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,[],pi_z_J,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate_GE,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

    Params.T=p_eqm.T;
    Params.T_ge=p_eqm.T;

    % Re-solve at GE T for stats
    [V_ge,Policy_ge]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    StationaryDist_ge=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_ge,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);

    AllStats_ge=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_ge,Policy_ge,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_ge=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_ge,Policy_ge,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    save('./SavedOutput/LuanWallenius2015_results_generaleqm.mat','V_ge','Policy_ge','StationaryDist_ge','AllStats_ge','AgeConditionalStats_ge','Params','-v7.3');
else
    load ./SavedOutput/LuanWallenius2015_results_generaleqm.mat
end
% Free up some memory
clear V_ge Policy_ge StationaryDist_ge

% Restore regime and T to baseline values for any subsequent re-runs
Params.regime=0;
Params.T=0.4;

%% Figure 5b: pre- vs post-reform (PE) vs post-reform (GE) employment
emp_age_ge=100*AgeConditionalStats_ge.working.Mean;
model_emp_ge=zeros(1,5);
for bb=1:5
    js=bin_starts(bb):(bin_starts(bb)+4);
    model_emp_ge(bb)=mean(emp_age_ge(js));
end

fig5b=figure(6);
plot(1:5,model_emp_pre,'--bo','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
plot(1:5,model_emp_post,'-rs','LineWidth',1.5,'MarkerFaceColor','r');
plot(1:5,model_emp_ge,'-.gd','LineWidth',1.5,'MarkerFaceColor','g');
set(gca,'XTick',1:5,'XTickLabel',bin_labels);
xlabel('Age');
ylabel('Fraction Employed (%)');
ylim([0 100]);
legend({'Pre-reform','Post-reform (regular only, PE)','Post-reform (regular only, GE)'},'Location','northeast');
title('Figure 5b: PE vs GE for the regular-only pension reform');
grid on;
saveas(fig5b,'SavedOutput/Graphs/LuanWallenius2015_Figure5b.png');

%% Quick GE summary
mean_emp_ge=mean(AgeConditionalStats_ge.working.Mean);
fprintf('\nPost-reform GE summary:\n');
fprintf('  Equilibrium T (GE):           %.4f BA  (paper: ~0.6)\n',Params.T_ge);
fprintf('  Mean employment rate (GE):    %.3f\n',mean_emp_ge);
fprintf('  Delta employment (GE - pre):  %+.3f\n',mean_emp_ge-mean_emp_pre);
fprintf('  Delta employment (GE - PE):   %+.3f\n',mean_emp_ge-mean_emp_post);
fprintf('  Employment 65-69 (pre):  %.1f%%\n',model_emp_pre(4));
fprintf('  Employment 65-69 (PE):   %.1f%%\n',model_emp_post(4));
fprintf('  Employment 65-69 (GE):   %.1f%%\n',model_emp_ge(4));


%% ===================================================================
%% FULL REFORM (regular + occupational), partial equilibrium
%% ===================================================================
% Same setup as doPart(2) PE, but now regime=2 so PB, DI AND OPB are all on
% their post-reform formulas (NDC pension, 64%-of-pre-DI-avg DI, and DC OPB
% with 4.5%/30% contributions). T is held at the pre-reform 0.4 BA (PE per
% LW2015 Section 6). This corresponds to LW2015 Figure 5's headline
% pre-vs-full post-reform comparison; combined with doPart(2)'s
% regular-only line, it recovers the LW2015 Figure 6 decomposition. Paper
% p.133 notes the full-reform line should sit BELOW the regular-only line
% at 65-69 and 70-74 because the new DC OPB is more generous and creates
% an income effect that dampens additional labor supply.

% Switch to full-reform regime, PE
Params.regime=2;
Params.T=0.4;

if doPart(4)==1
    vfoptions.verbose=1;
    tic;
    [V_full,Policy_full]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    vftime_full=toc
    vfoptions.verbose=0;

    StationaryDist_full=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_full,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);

    AllStats_full=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_full,Policy_full,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_full=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_full,Policy_full,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    save('./SavedOutput/LuanWallenius2015_results_fullreform.mat','V_full','Policy_full','StationaryDist_full','AllStats_full','AgeConditionalStats_full','Params','-v7.3');
else
    load ./SavedOutput/LuanWallenius2015_results_fullreform.mat
end
% Free up some memory
clear V_full Policy_full StationaryDist_full

% Restore baselines
Params.regime=0;
Params.T=0.4;

%% Per-bin full-PE employment (used by Fig 6 PE companion + summary stats below)
emp_age_full=100*AgeConditionalStats_full.working.Mean;
model_emp_full=zeros(1,5);
for bb=1:5
    js=bin_starts(bb):(bin_starts(bb)+4);
    model_emp_full(bb)=mean(emp_age_full(js));
end

%% Full reform PE summary
mean_emp_full=mean(AgeConditionalStats_full.working.Mean);
fprintf('\nFull reform PE summary:\n');
fprintf('  Mean employment rate (full PE):    %.3f\n',mean_emp_full);
fprintf('  Delta employment (full - pre):     %+.3f\n',mean_emp_full-mean_emp_pre);
fprintf('  Delta employment (full - reg only):%+.3f  (paper: should be NEGATIVE)\n',mean_emp_full-mean_emp_post);
fprintf('  Employment 65-69 (pre):              %.1f%%\n',model_emp_pre(4));
fprintf('  Employment 65-69 (reg-only PE):      %.1f%%\n',model_emp_post(4));
fprintf('  Employment 65-69 (full PE):          %.1f%%\n',model_emp_full(4));


%% ===================================================================
%% FULL REFORM (regular + occupational), general equilibrium
%% ===================================================================
% Like doPart(3) but with regime=2 (full reform). T is the GE price; the
% inline GeneralEqmEqn forces aggregate net surplus = T at convergence.
% Paper's headline 2.5-year retirement-age increase is computed under this
% full-GE regime.
if doPart(5)==1
    Params.regime=2;
    Params.T=0.4;

    FnsToEvaluate_GE.surplus=@(d,kprime,k,hist_idx,h,agej,r,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_65,regime) ...
        LuanWallenius2015_BudgetFn(d,kprime,k,hist_idx,h,agej,r,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_65,regime);

    GeneralEqmEqns.budget=@(surplus,T) surplus-T;

    GEPriceParamNames={'T'};
    heteroagentoptions.verbose=1;

    [p_eqm,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,[],pi_z_J,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate_GE,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

    Params.T=p_eqm.T;
    T_full_ge=p_eqm.T;

    [V_full_ge,Policy_full_ge]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    StationaryDist_full_ge=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_full_ge,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);

    AllStats_full_ge=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_full_ge,Policy_full_ge,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_full_ge=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_full_ge,Policy_full_ge,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    save('./SavedOutput/LuanWallenius2015_results_fullreform_ge.mat','V_full_ge','Policy_full_ge','StationaryDist_full_ge','AllStats_full_ge','AgeConditionalStats_full_ge','Params','T_full_ge','-v7.3');
else
    load ./SavedOutput/LuanWallenius2015_results_fullreform_ge.mat
end

% Restore baseline regime and T
Params.regime=0;
Params.T=0.4;

%% Full reform GE summary
emp_age_full_ge=100*AgeConditionalStats_full_ge.working.Mean;
model_emp_full_ge=zeros(1,5);
for bb=1:5
    js=bin_starts(bb):(bin_starts(bb)+4);
    model_emp_full_ge(bb)=mean(emp_age_full_ge(js));
end
mean_emp_full_ge=mean(AgeConditionalStats_full_ge.working.Mean);

fprintf('\n Full reform GE summary: \n');
fprintf('  Equilibrium T (full GE):           %.4f BA  (paper: ~0.6)\n',T_full_ge);
fprintf('  Mean employment rate (full GE):    %.3f\n',mean_emp_full_ge);
fprintf('  Delta employment (full GE - pre):  %+.3f\n',mean_emp_full_ge-mean_emp_pre);
fprintf('  Delta employment (full GE - PE):   %+.3f\n',mean_emp_full_ge-mean_emp_full);
fprintf('  Employment 65-69 (full GE):        %.1f%%\n',model_emp_full_ge(4));

%% Figure 5 (LW2015 Fig 5 headline): pre-reform vs full post-reform employment
% Paper's Figure 5 is pre-reform vs full post-reform under GE. doPart(5)
% provides the full-GE solve, so this is the paper-aligned headline plot.
fig5=figure(5);
plot(1:5,model_emp_pre,'--bo','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
plot(1:5,model_emp_full_ge,'-rs','LineWidth',1.5,'MarkerFaceColor','r');
set(gca,'XTick',1:5,'XTickLabel',bin_labels);
xlabel('Age');
ylabel('Fraction Employed (%)');
ylim([0 100]);
legend({'Pre-reform','Post-reform (full, GE)'},'Location','northeast');
title('Figure 5: Employment, pre- vs full post-reform (GE)');
grid on;
saveas(fig5,'SavedOutput/Graphs/LuanWallenius2015_Figure5.png');

%% Figure 6 (LW2015 Fig 6): pre vs regular-only GE vs full-reform GE
% Paper-aligned: both post-reform lines are GE (Section 5's main calibration).
% Paper p.133 result: the full-reform line sits BELOW the regular-only line
% at 65-69 and 70-74, because the new DC OPB is more generous and creates an
% income effect that dampens additional labor supply.
fig6=figure(7);
plot(1:5,model_emp_pre,'--bo','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
plot(1:5,model_emp_ge,'-rs','LineWidth',1.5,'MarkerFaceColor','r');
plot(1:5,model_emp_full_ge,'-md','LineWidth',1.5,'MarkerFaceColor','m');
set(gca,'XTick',1:5,'XTickLabel',bin_labels);
xlabel('Age');
ylabel('Fraction Employed (%)');
ylim([0 100]);
legend({'Pre-reform','Post-reform (regular only, GE)','Post-reform (full, GE)'},'Location','northeast');
title('Figure 6: Role of occupational pensions in pension reform');
grid on;
saveas(fig6,'SavedOutput/Graphs/LuanWallenius2015_Figure6.png');


%% ===================================================================
%% Average retirement age and DI claim age across all four regimes
%% ===================================================================
% retire_age FnsToEvaluate returns the decoded swa (real age); at the
% terminal model period j=N_j the cross-sectional mean across PTypes is the
% population-average final retirement age. di_claim_age returns 0 for
% non-DI claimants, so dividing AgeConditionalStats.di_claim_age.Mean by
% AgeConditionalStats.onDI.Mean at j=N_j gives the conditional mean DI
% claim age among DI claimants.

cond_di_age=@(acs) acs.di_claim_age.Mean(N_j)/max(acs.onDI.Mean(N_j),eps);

fprintf('\n=== Average retirement age and DI claim age (real ages) ===\n');
fprintf('                            Retire age   DI claim age\n');
fprintf('  Pre-reform                   %5.2f        %5.2f\n', ...
    AgeConditionalStats_pre.retire_age.Mean(N_j),  cond_di_age(AgeConditionalStats_pre));
fprintf('  Post-reform (reg only, PE)   %5.2f        %5.2f\n', ...
    AgeConditionalStats_post.retire_age.Mean(N_j), cond_di_age(AgeConditionalStats_post));
fprintf('  Post-reform (reg only, GE)   %5.2f        %5.2f\n', ...
    AgeConditionalStats_ge.retire_age.Mean(N_j),   cond_di_age(AgeConditionalStats_ge));
fprintf('  Post-reform (full, PE)       %5.2f        %5.2f\n', ...
    AgeConditionalStats_full.retire_age.Mean(N_j), cond_di_age(AgeConditionalStats_full));
fprintf('  Post-reform (full, GE)       %5.2f        %5.2f\n', ...
    AgeConditionalStats_full_ge.retire_age.Mean(N_j), cond_di_age(AgeConditionalStats_full_ge));
fprintf('\n  LW2015 reports:  pre-reform retire age = 62.1, post-reform (full GE) = 64.6 (Delta = +2.5) \n');
fprintf('                   pre-reform DI claim age = 52.8, post-reform = 51.6 (Delta = -1.2) \n');


%% ===================================================================
%% Welfare CEV (LW2015 Section 5.3, full GE vs pre-reform)
%% ===================================================================
% Closed-form CV for log utility. The per-period utility log(c) - b*l + h has
% the labor-disutility and health terms additively separable from
% consumption, so scaling c by (1+cev) shifts each period's utility by
% log(1+cev) -- a constant that (a) does NOT change the optimal policy
% (constant shift doesn't affect arg max) and (b) shifts lifetime value by
% exactly Gamma*log(1+cev), where Gamma = (1-beta^N_j)/(1-beta) is the
% discount sum. Setting V_post + log(1+cev)*Gamma = V_pre and solving:
%   cev = exp((V_pre - V_post)/Gamma) - 1
%   CV  = -cev = 1 - exp((V_pre - V_post)/Gamma)   [paper convention: CV > 0 means reform welfare-improving]
if doPart(6)==1
    load ./SavedOutput/LuanWallenius2015_results.mat V_pre
    load ./SavedOutput/LuanWallenius2015_results_fullreform_ge.mat V_full_ge

    Gamma=(1-Params.beta^N_j)/(1-Params.beta);
    init_k=1; init_hist=1; init_h=5; init_j=1;

    CV=struct();
    CV_pop_weighted=0;
    for ii=1:N_i
        iistr=Names_i{ii};
        V_pre_init  =V_pre.(iistr)(init_k,init_hist,init_h,init_j);
        V_post_init =V_full_ge.(iistr)(init_k,init_hist,init_h,init_j);
        CV.(iistr)=100*(1-exp((V_pre_init-V_post_init)/Gamma));
        CV_pop_weighted=CV_pop_weighted+Params.ptypedist(ii)*CV.(iistr);
    end

    save('./SavedOutput/LuanWallenius2015_results_cev.mat','CV','CV_pop_weighted','Gamma');
else
    load ./SavedOutput/LuanWallenius2015_results_cev.mat
end

fprintf('\n=== Welfare CEV (LW2015 Section 5.3, full GE vs pre-reform) ===\n');
fprintf('  CV is the per-period consumption reduction in the post-reform world\n');
fprintf('  required to equate lifetime utility with the pre-reform world.\n');
fprintf('  Positive CV = post-reform is welfare-improving.\n\n');
fprintf('  PType            CV (%%)\n');
for ii=1:N_i
    fprintf('  %-14s  %+7.3f\n',Names_i{ii},CV.(Names_i{ii}));
end
fprintf('  %-14s  %+7.3f\n','Pop-weighted',CV_pop_weighted);
fprintf('\n  LW2015 (p.134, against full GE): NCL 4.3%%, CH 8.6%% (reform favours high-income types)\n');