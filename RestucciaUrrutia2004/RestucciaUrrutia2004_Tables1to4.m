% Restuccia & Urrutia (2004)
% Tables 1,2,3,4
% Figure 1,2


%% Start with Table 1
% FOR SOME REASON USING CNTRL+ENTER LEADS TO ERROR 'amibiguity about matlab
% path' for the second function to be evalutated (_kappaF). But running the
% same code in piece by piece using F9 does not lead to same error. I have
% no idea why!!!! FIXED BY RENAMING _kappaF to _FnkappaF and _F to _FnF.
% Still have no idea why.
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'agej','plowerbar', 'pupperbar','nlowerbar','nupperbar','psi0','psi1'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1) RestucciaUrrutia2004_HFn(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1); % H, Human capital
FnsToEvaluateParamNames(2).Names={'agej','w','kappa0','kappa1','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_FnkappaF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1); % kappaF (a3_val==1, is that the state s=1 for old)
FnsToEvaluateParamNames(3).Names={'agej'};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (agej==2)*(1-a3_val); % fraction of j==2 population with s==0: Table 1 (i) fraction of non-college
FnsToEvaluateParamNames(4).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % fraction college dropout
FnsToEvaluateParamNames(5).Names={'agej','gamma','g'};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,gamma,g) (agej==1)*((d1_val/z1_val)^(1/gamma)-g); % e=(bhat/b)^(1/gamma)-g for age j==1, private expenditure in early education
FnsToEvaluateParamNames(6).Names={'agej','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_FnF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1); % F (a3_val==1, is that the state s=1 for old)
% The next six are about getting the wages for non-college, college dropout, and college so as to calculate the (completed and dropout) college premium
FnsToEvaluateParamNames(7).Names={'agej','w'};
FnsToEvaluateFn_7 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w) w*a1_val*(a3_val==0)*(agej==2); % Total wages for non-college
FnsToEvaluateParamNames(8).Names={'agej'};
FnsToEvaluateFn_8 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==0)*(agej==2); % Fraction non-college
FnsToEvaluateParamNames(9).Names={'agej','w','psi0','psi1','plowerbar'};
FnsToEvaluateFn_9 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,psi0,psi1,plowerbar) w*plowerbar*a1_val*(a3_val==1)*(agej==2)*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % Total wages for college dropouts
FnsToEvaluateParamNames(10).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_10 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (a3_val==1)*(agej==2)*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % Fraction college dropout
FnsToEvaluateParamNames(11).Names={'agej','w','psi0','psi1','pupperbar'};
FnsToEvaluateFn_11 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,psi0,psi1,pupperbar) w*pupperbar*a1_val*(a3_val==1)*(agej==2)*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % Total wages for college completeds
FnsToEvaluateParamNames(12).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_12 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (a3_val==1)*(agej==2)*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % Fraction college completed
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9,FnsToEvaluateFn_10,FnsToEvaluateFn_11,FnsToEvaluateFn_12};

MomentsForTable1=gather(EvalFnOnAgentDist_AggVars_FHorz_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_gridfn, a_gridfn, z_gridfn, simoptions, AgeDependentGridParamNames));

GDP=Params.A*MomentsForTable1(1); % Y=A*H from firms production function.

% Table1.FractionOfNonCollege=MomentsForTable1(3); % fraction of agej=2 population with s=0 % SHOULD BE A CONDITIONAL MOMENT
% Table1.DropoutRate=MomentsForTable1(4); % mean of (1-q(bhat)) for population with agej=2 and s=1 % SHOULD BE A CONDITIONAL MOMENT
Table1.EarlyEducDivGDP=(MomentsForTable1(5)+Params.g)/GDP; % (e+g)/Y
Table1.CollegeDivGDP=MomentsForTable1(6)/GDP; % F/Y
Table1.PublicDivTotalCollege=MomentsForTable1(2)/MomentsForTable1(6); % kappaF/F

% Average wage of non-college
AvgWage_NoCollege=MomentsForTable1(7)/MomentsForTable1(8);
% Average wage of college dropout
AvgWage_CollegeDropout=MomentsForTable1(9)/MomentsForTable1(10);
% Average wage of college
AvgWage_CollegeCompleted=MomentsForTable1(11)/MomentsForTable1(12);

Table1.AverageDropoutPremium=AvgWage_CollegeDropout/AvgWage_NoCollege;
Table1.AverageCollegePremium=AvgWage_CollegeCompleted/AvgWage_NoCollege;

disp('Tables1to4: progress marker 1')
% The next moment is standard deviation of log earnings. Earnings is simply w*h.
% Since w is a constant and equal to one we actually just need
% stddev(log(h)), but I will use w as part of calculation so that if we wanted to change it we could without breaking codes.
% It is unclear if Restuccia & Urrutia do this at the level of households
% or workers. (And if workers, why are the 'young members of older households' not included?)
% Turns out they do an even more restricted version:
% For reasons that are unclear this is not how Restuccia & Urrutia (2004)
% actually calculate this. At bottom right of pg 1363 they instead
% calculate the 'standard deviation of log earnings' based only on old
% household earnings, and not those of all households. "We measure
% disparity of earnings in the model as the standard deviations of
% log(wh_0) across parents". So I now do this instead in the following
% lines, even though I feel this is the 'incorrect' calibration mapping of
% model to data.
FnsToEvaluateParamNames(13).Names={'w'};
FnsToEvaluateFn_13 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,w) log(w*a1_val); % H, Human capital
FnsToEvaluateParamNames(14).Names={};
FnsToEvaluateFn_14 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a3_val; % s
FnsToEvaluateParamNames(15).Names={'agej'};
FnsToEvaluateFn_15 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % agej (for the panel data)
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9,FnsToEvaluateFn_10,FnsToEvaluateFn_11,FnsToEvaluateFn_12,FnsToEvaluateFn_13,FnsToEvaluateFn_14,FnsToEvaluateFn_15};

% options.agegroupings=1:1:N_j; % for each age, this is anyway the default
% simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
% simoptions.parallel=2;
% The following are elements of Table 1 that are really age-contional moments
Table1.StdDevLogEarnings=gather(AgeConditionalStats(13).StdDeviation(2));
Table1.FractionOfNonCollege=gather(AgeConditionalStats(3).Mean(2)); % fraction of agej=2 population with s=0
Table1.DropoutRate=gather(AgeConditionalStats(4).Mean(2)/AgeConditionalStats(14).Mean(2)); % mean of (1-q(bhat)) for population with agej=2 and s=1
Table1.DropoutRate_Alt=gather(MomentsForTable1(10)/AgeConditionalStats(14).Mean(2));

disp('Tables1to4: progress marker 2')

% Restuccia & Urrutia (2004), pg 1363: "[We measure] the degree of
% intergenerational persistence as the estimated beta1 coefficient in the
% regression: log(why')=beta0+beta1*log(wh0)+epsilon
% using parent-children pairs properly weighted by the invariante
% distribution of types mew.
% Their h0 is the 'a' variable h at age j==2. Their hy' is the 'a' variable
% at age j==1 (but it can also be calculated, see formula on bottom right of their pg 1359). I will calculate it based on
% actual simulated panel data, just to illustrate how this is done (and
% because in more complex models the closed form formula approach wouldn't work)
% simoptions.numbersims=1000;
simoptions.simperiods=6;
YetMoreMomentsForTable1=SimPanelValues_FHorz_Case2(StationaryDist,Policy, FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,PhiaprimeParamNames, simoptions);
X1=YetMoreMomentsForTable1(13,1:end-1,:).*(YetMoreMomentsForTable1(15,1:end-1,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=YetMoreMomentsForTable1(13,2:end,:).*(YetMoreMomentsForTable1(15,2:end,:)==1);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
Table1.IntergenerationalCorrelationOfEarnings=regcoeff(2);

%Table 1
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table1.tex', 'w');
fprintf(FID, 'Calibration of the Benchmark Economy \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}rlcccc} \n \\hline \\hline \n');
fprintf(FID, 'Target & Data & Model & Parameter & Value \\\\ \\hline \n');
fprintf(FID, '(i) Fraction of non-college & 0.54 & %8.2f & $\\psi_0$ & %8.2f \\\\ \n ', Table1.FractionOfNonCollege, Params.psi0);
fprintf(FID, '(ii) Dropout rate & 0.50 & %8.2f & $\\psi_1$ & %8.2f \\\\ \n ', Table1.DropoutRate, Params.psi1);
fprintf(FID, '(iii) Early education/GDP & 0.044 & %8.2f & $\\gamma$ & %8.2f \\\\ \n ', Table1.EarlyEducDivGDP, Params.gamma);
fprintf(FID, '(iv) College/GDP & 0.028 & %8.2f & $f$ & %8.2f \\\\ \n ', Table1.CollegeDivGDP, Params.f);
fprintf(FID, '(v) Public/Total college & 0.64 & %8.2f & $\\kappa_0$ & %8.2f \\\\ \n ', Table1.PublicDivTotalCollege, Params.kappa0);
fprintf(FID, '(vi) Average dropout premium & 1.41 & %8.2f & $\\underbar{p}$ & %8.2f \\\\ \n ', Table1.AverageDropoutPremium, Params.plowerbar);
fprintf(FID, '(vii) Average college premium & 2.33 & %8.2f & $\\bar{p}$ & %8.2f \\\\ \n ', Table1.AverageCollegePremium, Params.pupperbar);
fprintf(FID, '(viii) std(log earnings) & 0.60 & %8.2f & $\\sigma_b$ & %8.2f \\\\ \n ', Table1.StdDevLogEarnings, Params.sigma_b);
fprintf(FID, '(ix) Intergenerational correlation of earnings & 0.40 & %8.2f & $\\rho_b$ & %8.2f \\\\ \n ', Table1.IntergenerationalCorrelationOfEarnings, Params.rho_b);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Data column is copy of original from Restuccia \\& Urrutia (2004), is not part of the replication. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Now Table 2
disp('Tables1to4: progress marker 3')

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'agej','w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w) log(w*a1_val); % Earnings
FnsToEvaluateParamNames(2).Names={'agej'};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) log(a2_val); % bhat, (log) aquired ability
FnsToEvaluateParamNames(3).Names={'agej'};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) log(z1_val); % b, (log) innate ability
FnsToEvaluateParamNames(4).Names={'agej'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % age
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4};

% First, the standard deviations in the cross section. Since the number for
% earnings is same as that in Table 1 this has presumably also been
% calcuated as conditional on agej==2. I assume this is true of all the
% numbers here.
% simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
% simoptions.parallel=2;
% The following are elements of Table 2 that are really age-contional moments
Table2.StdDevLogEarnings=gather(AgeConditionalStats(1).StdDeviation(2));
Table2.StdDevLogAcquiredAbility=gather(AgeConditionalStats(2).StdDeviation(2));
Table2.StdDevLogInnateAbility=gather(AgeConditionalStats(3).StdDeviation(2));

% Second, the intergenerational correlations.
PanelForCalculatingIntergenCorr=SimPanelValues_FHorz_Case2(StationaryDist,Policy, FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,PhiaprimeParamNames, simoptions);
% Earnings
X1=PanelForCalculatingIntergenCorr(1,1:end-1,:).*(PanelForCalculatingIntergenCorr(4,1:end-1,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(1,2:end,:).*(PanelForCalculatingIntergenCorr(4,2:end,:)==1);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
Table2.IntergenerationalCorrelationOfEarnings=regcoeff(2);
% Acquired Ability (only measured when j==2, so have to do slightly differently)
X1=PanelForCalculatingIntergenCorr(2,1:end-2,:).*(PanelForCalculatingIntergenCorr(4,1:end-2,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(2,3:end,:).*(PanelForCalculatingIntergenCorr(4,3:end,:)==2);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
Table2.IntergenerationalCorrelationOfAcquiredAbility=regcoeff(2); % corr(X1,Y,'rows','pairwise')
% Innate Ability
X1=PanelForCalculatingIntergenCorr(3,1:end-1,:).*(PanelForCalculatingIntergenCorr(4,1:end-1,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(3,2:end,:).*(PanelForCalculatingIntergenCorr(4,2:end,:)==1);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
Table2.IntergenerationalCorrelationOfInnateAbility=regcoeff(2); % corr(X1,Y,'rows','pairwise')

%Table 2
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table2.tex', 'w');
fprintf(FID, 'Disparity and Persistence in the Benchmark Economy \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & Innate & Acquired &  \\\\ \n');
fprintf(FID, ' & Ability & Ability & Earnings \\\\ \\hline \n');
fprintf(FID, 'Cross-sectional disparity: & %8.2f & %8.2f & %8.2f \\\\ \n ', Table2.StdDevLogInnateAbility,Table2.StdDevLogAcquiredAbility, Table2.StdDevLogEarnings);
fprintf(FID, '$\\quad$ std(log x) &  &  & \\\\ \n ');
fprintf(FID, 'Intergenerational correlation: & %8.2f & %8.2f & %8.2f \\\\ \n ', Table2.IntergenerationalCorrelationOfInnateAbility,Table2.IntergenerationalCorrelationOfAcquiredAbility, Table2.IntergenerationalCorrelationOfEarnings);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: In model notation these columns are: $b$, $\\hat{b}$, and $wh$. I follow Restuccia \\& Urrutia (2004) in reporting as cross-sectional the number conditional on being an elderly household, not cross-sectional over the whole model economy. \\\\ \n');
fprintf(FID, 'Restuccia \\& Urrutia (2004) explain calculation of intergenerational correlation of earnings at bottom of pg 1363. I assume the intergeneration correlations of (log) innate and acquired ability are calculated by the analagous regressions (with modification for acquired ability as is only observed for old). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Now for Table 3, Figure 1 & Figure 2
disp('Tables1to4: progress marker 4')

% These could be calculated from the StationaryDist directly, but I have
% just done it from a simulated panel instead.

% Many of these are same as for Table 2, but now without log (actually
% given we are just interested in terciles I could still use log if I
% wanted)
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,w) w*a1_val; % Earnings
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a2_val; % bhat, aquired ability.
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) z1_val; % b, innate ability
FnsToEvaluateParamNames(4).Names={'agej'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % age
FnsToEvaluateParamNames(5).Names={'agej','gamma','g'};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,gamma,g) (agej==1)*((d1_val/z1_val)^(1/gamma)-g); % e=(bhat/b)^(1/gamma)-g for age j==1, private expenditure in early education
FnsToEvaluateParamNames(6).Names={'agej'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==1)*(agej==2); % college enrollment decision
FnsToEvaluateParamNames(7).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_7 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % indicator for college completion
FnsToEvaluateParamNames(8).Names={};
FnsToEvaluateFn_8 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) d1_val; % bhat, aquired ability (again, but in a way that is easily paired with the young parent HH)
FnsToEvaluateParamNames(9).Names={};
FnsToEvaluateFn_9 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a1_val; % Human capital
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9};

PanelForCalculatingTable3=SimPanelValues_FHorz_Case2(StationaryDist,Policy, FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,PhiaprimeParamNames, simoptions);

% Table 3:
% For Panel A and B we work with young HHs
Indicator_YoungHH=shiftdim((PanelForCalculatingTable3(4,:,:)==1),1);
% Calculate the tercile cutoffs for young parent earnings
Earnings=shiftdim(PanelForCalculatingTable3(1,:,:),1);
Earnings=Earnings(Indicator_YoungHH);
YoungEarningsTercileCutoffs=prctile(Earnings,[33,66]);
% Calculate the tercile cutoffs for child innate ability
ChildInnateAbility=shiftdim(PanelForCalculatingTable3(3,:,:),1);
ChildInnateAbility=ChildInnateAbility(Indicator_YoungHH);
ChildInnateAbilityTercileCutoffs=prctile(ChildInnateAbility,[33,66]);
% Get the 'Expenditures in early education' (note, is just private expenditures)
EarlyEducExpend=shiftdim(PanelForCalculatingTable3(5,:,:),1);
EarlyEducExpend=EarlyEducExpend(Indicator_YoungHH);
% Get the 'Childs aquired ability'
ChildAquiredAbility=shiftdim(PanelForCalculatingTable3(8,:,:),1);
ChildAquiredAbility=ChildAquiredAbility(Indicator_YoungHH);
% Set up indicators based on the terciles
Indicator_LowEarningsTercile=(Earnings<YoungEarningsTercileCutoffs(1));
Indicator_MidEarningsTercile=(Earnings>=YoungEarningsTercileCutoffs(1)).*(Earnings<YoungEarningsTercileCutoffs(2));
Indicator_HighEarningsTercile=(Earnings>=YoungEarningsTercileCutoffs(2));
Indicator_LowChildInnateAbilityTercile=(ChildInnateAbility<ChildInnateAbilityTercileCutoffs(1));
Indicator_MidChildInnateAbilityTercile=(ChildInnateAbility>=ChildInnateAbilityTercileCutoffs(1)).*(ChildInnateAbility<ChildInnateAbilityTercileCutoffs(2));
Indicator_HighChildInnateAbilityTercile=(ChildInnateAbility>=ChildInnateAbilityTercileCutoffs(2));

Table3.PanelA.i11=mean(EarlyEducExpend(logical(Indicator_LowEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelA.i12=mean(EarlyEducExpend(logical(Indicator_MidEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelA.i13=mean(EarlyEducExpend(logical(Indicator_HighEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelA.i21=mean(EarlyEducExpend(logical(Indicator_LowEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelA.i22=mean(EarlyEducExpend(logical(Indicator_MidEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelA.i23=mean(EarlyEducExpend(logical(Indicator_HighEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelA.i31=mean(EarlyEducExpend(logical(Indicator_LowEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));
Table3.PanelA.i32=mean(EarlyEducExpend(logical(Indicator_MidEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));
Table3.PanelA.i33=mean(EarlyEducExpend(logical(Indicator_HighEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));

Table3.PanelB.i11=mean(ChildAquiredAbility(logical(Indicator_LowEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelB.i12=mean(ChildAquiredAbility(logical(Indicator_MidEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelB.i13=mean(ChildAquiredAbility(logical(Indicator_HighEarningsTercile.*Indicator_LowChildInnateAbilityTercile)));
Table3.PanelB.i21=mean(ChildAquiredAbility(logical(Indicator_LowEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelB.i22=mean(ChildAquiredAbility(logical(Indicator_MidEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelB.i23=mean(ChildAquiredAbility(logical(Indicator_HighEarningsTercile.*Indicator_MidChildInnateAbilityTercile)));
Table3.PanelB.i31=mean(ChildAquiredAbility(logical(Indicator_LowEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));
Table3.PanelB.i32=mean(ChildAquiredAbility(logical(Indicator_MidEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));
Table3.PanelB.i33=mean(ChildAquiredAbility(logical(Indicator_HighEarningsTercile.*Indicator_HighChildInnateAbilityTercile)));

% Note: the following overwrite the ChildInnateAbility (the period 2 state instead of the period 1 decision)

% For Panel C we work with old HHs
Indicator_OldHH=shiftdim((PanelForCalculatingTable3(4,:,:)==2),1);
% Calculate the tercile cutoffs for old parent earnings
Earnings=shiftdim(PanelForCalculatingTable3(1,:,:),1); % redundant as already calculated above
Earnings=Earnings(Indicator_OldHH);
OldEarningsTercileCutoffs=prctile(Earnings,[33,66]);
% Calculate the tercile cutoffs for child innate ability (now based on Old HHs rather than young HHs)
ChildAcquiredAbility=shiftdim(PanelForCalculatingTable3(2,:,:),1);
ChildAcquiredAbility=ChildAcquiredAbility(Indicator_OldHH);
ChildAcquiredAbilityTercileCutoffs=prctile(ChildAcquiredAbility,[33,66]);
% Get the 'College enrollment rate' (note, prior to taking mean is actually a binary for enrollment)
CollegeEnrollment=shiftdim(PanelForCalculatingTable3(6,:,:),1);
CollegeEnrollment=CollegeEnrollment(Indicator_OldHH);
% Set up indicators based on the terciles
Indicator_LowEarningsTercile=(Earnings<OldEarningsTercileCutoffs(1));
Indicator_MidEarningsTercile=(Earnings>=OldEarningsTercileCutoffs(1)).*(Earnings<OldEarningsTercileCutoffs(2));
Indicator_HighEarningsTercile=(Earnings>=OldEarningsTercileCutoffs(2));
Indicator_LowChildAcquiredAbilityTercile=(ChildAcquiredAbility<ChildAcquiredAbilityTercileCutoffs(1));
Indicator_MidChildAcquiredAbilityTercile=(ChildAcquiredAbility>=ChildAcquiredAbilityTercileCutoffs(1)).*(ChildAcquiredAbility<ChildAcquiredAbilityTercileCutoffs(2));
Indicator_HighChildAcquiredAbilityTercile=(ChildAcquiredAbility>=ChildAcquiredAbilityTercileCutoffs(2));

Table3.PanelC.i11=mean(CollegeEnrollment(logical(Indicator_LowEarningsTercile.*Indicator_LowChildAcquiredAbilityTercile)));
Table3.PanelC.i12=mean(CollegeEnrollment(logical(Indicator_MidEarningsTercile.*Indicator_LowChildAcquiredAbilityTercile)));
Table3.PanelC.i13=mean(CollegeEnrollment(logical(Indicator_HighEarningsTercile.*Indicator_LowChildAcquiredAbilityTercile)));
Table3.PanelC.i21=mean(CollegeEnrollment(logical(Indicator_LowEarningsTercile.*Indicator_MidChildAcquiredAbilityTercile)));
Table3.PanelC.i22=mean(CollegeEnrollment(logical(Indicator_MidEarningsTercile.*Indicator_MidChildAcquiredAbilityTercile)));
Table3.PanelC.i23=mean(CollegeEnrollment(logical(Indicator_HighEarningsTercile.*Indicator_MidChildAcquiredAbilityTercile)));
Table3.PanelC.i31=mean(CollegeEnrollment(logical(Indicator_LowEarningsTercile.*Indicator_HighChildAcquiredAbilityTercile)));
Table3.PanelC.i32=mean(CollegeEnrollment(logical(Indicator_MidEarningsTercile.*Indicator_HighChildAcquiredAbilityTercile)));
Table3.PanelC.i33=mean(CollegeEnrollment(logical(Indicator_HighEarningsTercile.*Indicator_HighChildAcquiredAbilityTercile)));

Table3matrix=zeros(3,3,3); % Third dimension is panel
Table3matrix(:,:,1)=[Table3.PanelA.i11,Table3.PanelA.i12,Table3.PanelA.i13; Table3.PanelA.i21,Table3.PanelA.i22,Table3.PanelA.i23; Table3.PanelA.i31,Table3.PanelA.i32,Table3.PanelA.i33];
Table3matrix(:,:,2)=[Table3.PanelB.i11,Table3.PanelB.i12,Table3.PanelB.i13; Table3.PanelB.i21,Table3.PanelB.i22,Table3.PanelB.i23; Table3.PanelB.i31,Table3.PanelB.i32,Table3.PanelB.i33];
Table3matrix(:,:,3)=[Table3.PanelC.i11,Table3.PanelC.i12,Table3.PanelC.i13; Table3.PanelC.i21,Table3.PanelC.i22,Table3.PanelC.i23; Table3.PanelC.i31,Table3.PanelC.i32,Table3.PanelC.i33];

%Table 3
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table3.tex', 'w');
fprintf(FID, 'Decision Rules by Parents Earnings and Childs Ability \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, '\\multicolumn{4}{c}{Panel A: Expenditures in Early Education} \\\\ \n');
fprintf(FID, ' & \\multicolumn{3}{c}{Young Parent Earnings Tercile} \\\\ \n');
fprintf(FID, 'Child Innate Ability: & I & II & III \\\\ \n \\hline ');
fprintf(FID, 'Low: & %8.4f & %8.4f & %8.4f \\\\ \n ', Table3matrix(1,:,1));
fprintf(FID, 'Medium: & %8.4f & %8.4f & %8.4f \\\\ \n ', Table3matrix(2,:,1));
fprintf(FID, 'High: & %8.4f & %8.4f & %8.4f \\\\ \n ', Table3matrix(3,:,1));
fprintf(FID, '\\hline \n ');
fprintf(FID, '\\multicolumn{4}{c}{Panel B: Childs Acquired Ability} \\\\ \n');
fprintf(FID, ' & \\multicolumn{3}{c}{Young Parent Earnings Tercile} \\\\ \n');
fprintf(FID, 'Child Innate Ability: & I & II & III \\\\ \n \\hline ');
fprintf(FID, 'Low: & %8.3f & %8.3f & %8.3f \\\\ \n ', Table3matrix(1,:,2));
fprintf(FID, 'Medium: & %8.3f & %8.3f & %8.3f \\\\ \n ', Table3matrix(2,:,2));
fprintf(FID, 'High: & %8.3f & %8.3f & %8.3f \\\\ \n ', Table3matrix(3,:,2));
fprintf(FID, '\\hline \n ');
fprintf(FID, '\\multicolumn{4}{c}{Panel C: College Enrollment Rate (percentage)} \\\\ \n');
fprintf(FID, ' & \\multicolumn{3}{c}{Old Parent Earnings Tercile} \\\\ \n');
fprintf(FID, 'Child Acquired Ability: & I & II & III \\\\ \n \\hline ');
fprintf(FID, 'Low: & %8.2f & %8.2f & %8.2f \\\\ \n ', 100*Table3matrix(1,:,3));
fprintf(FID, 'Medium: & %8.2f & %8.2f & %8.2f \\\\ \n ', 100*Table3matrix(2,:,3));
fprintf(FID, 'High: & %8.2f & %8.2f & %8.2f \\\\ \n ', 100*Table3matrix(3,:,3));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Earnings, Innate ability, and aquired ability are $wh$, $b$ and $bhat$. The three panels report $e$, $bhat$ and $s$ respectively. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


% For Figure 1
Indicator_OldHH=shiftdim((PanelForCalculatingTable3(4,:,:)==2),1);
% Calculate the quintile cutoffs for old parent earnings
Earnings=shiftdim(PanelForCalculatingTable3(1,:,:),1); % redundant as already calculated above
Earnings=Earnings(Indicator_OldHH);
OldEarningsTercileCutoffs=prctile(Earnings,[20,40,60,80]);
% The three variables that we wish to graph
ChildInnateAbility=shiftdim(PanelForCalculatingTable3(3,:,:),1);
ChildInnateAbility=ChildInnateAbility(Indicator_OldHH);
ChildAcquiredAbility=shiftdim(PanelForCalculatingTable3(2,:,:),1);
ChildAcquiredAbility=ChildAcquiredAbility(Indicator_OldHH);
HumanCapital=shiftdim(PanelForCalculatingTable3(9,:,:),1);
HumanCapital=HumanCapital(Indicator_OldHH);
% Set up indicators based on the quintiles
Indicator_EarningsQuintile1=(Earnings<OldEarningsTercileCutoffs(1));
Indicator_EarningsQuintile2=(Earnings>=OldEarningsTercileCutoffs(1)).*(Earnings<OldEarningsTercileCutoffs(2));
Indicator_EarningsQuintile3=(Earnings>=OldEarningsTercileCutoffs(2)).*(Earnings<OldEarningsTercileCutoffs(3));
Indicator_EarningsQuintile4=(Earnings>=OldEarningsTercileCutoffs(3)).*(Earnings<OldEarningsTercileCutoffs(4));
Indicator_EarningsQuintile5=(Earnings>=OldEarningsTercileCutoffs(4));

Figure1.InnateAbility_raw=[mean(ChildInnateAbility(logical(Indicator_EarningsQuintile1)));...
    mean(ChildInnateAbility(logical(Indicator_EarningsQuintile2)));...
    mean(ChildInnateAbility(logical(Indicator_EarningsQuintile3)));...
    mean(ChildInnateAbility(logical(Indicator_EarningsQuintile4)));...
    mean(ChildInnateAbility(logical(Indicator_EarningsQuintile5)))];
% Normalize relative to bottom quintile
Figure1.InnateAbility=(Figure1.InnateAbility_raw)./(Figure1.InnateAbility_raw(1));

Figure1.AcquiredAbility_raw=[mean(ChildAcquiredAbility(logical(Indicator_EarningsQuintile1)));...
    mean(ChildAcquiredAbility(logical(Indicator_EarningsQuintile2)));...
    mean(ChildAcquiredAbility(logical(Indicator_EarningsQuintile3)));...
    mean(ChildAcquiredAbility(logical(Indicator_EarningsQuintile4)));...
    mean(ChildAcquiredAbility(logical(Indicator_EarningsQuintile5)))];
% Normalize relative to bottom quintile
Figure1.AcquiredAbility=(Figure1.AcquiredAbility_raw)./(Figure1.AcquiredAbility_raw(1));

Figure1.HumanCapital_raw=[mean(HumanCapital(logical(Indicator_EarningsQuintile1)));...
    mean(HumanCapital(logical(Indicator_EarningsQuintile2)));...
    mean(HumanCapital(logical(Indicator_EarningsQuintile3)));...
    mean(HumanCapital(logical(Indicator_EarningsQuintile4)));...
    mean(HumanCapital(logical(Indicator_EarningsQuintile5)))];
% Normalize relative to bottom quintile
Figure1.HumanCapital=(Figure1.HumanCapital_raw)./(Figure1.HumanCapital_raw(1));


% Figure 2
% Is also off of old HH earnings quintiles. So just need to get the
% variables themselves and can reuse all the indicators from Figure 1
CollegeEnrollment=shiftdim(PanelForCalculatingTable3(6,:,:),1);
CollegeEnrollment=CollegeEnrollment(Indicator_OldHH);
CollegeCompletion=shiftdim(PanelForCalculatingTable3(7,:,:),1);
CollegeCompletion=CollegeCompletion(Indicator_OldHH);

Figure2.CollegeEnrollment=100*[mean(CollegeEnrollment(logical(Indicator_EarningsQuintile1)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile2)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile3)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile4)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile5)))];

Figure2.CollegeCompletion=100*[mean(CollegeCompletion(logical(Indicator_EarningsQuintile1)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile2)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile3)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile4)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile5)))];


% Fig 1
figure(1);
bar([Figure1.InnateAbility, Figure1.AcquiredAbility, Figure1.HumanCapital])
xlabel('Parent Earnings Quintile')
ylabel('Fraction of Bottom Quintile')
legend('Innate Ability','Acquired Ability','Human Capital', 'Location','northwest')
saveas(gcf,'./SavedOutput/Graphs/RestucciaUrrutia2004_Figure1.png')

% Fig 2 (only bottom panel, top panel is from data, not model)
figure(2)
bar([Figure2.CollegeEnrollment, Figure2.CollegeCompletion])
xlabel('Parent Earnings Quintile')
ylim([0,100])
ylabel('Percent of Quintile')
legend('Enrolled','Completed', 'Location','northwest')
saveas(gcf,'./SavedOutput/Graphs/RestucciaUrrutia2004_Figure2.png')


%% Stats for Table 4 (table 4 requires comparison to baseline, so calculate the baseline here).
% The early education expenditure decision is 'e': 5
% Change in early education is 'bhat': 8
% The college enrollment decision is 's': 6
% Private expenditures on college (F-kappaF?):
FnsToEvaluateParamNames(10).Names={'agej','w','kappa0','kappa1','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_10 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_kappaF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1); % kappaF (a3_val==1, is that the state s=1 for old)
FnsToEvaluateParamNames(11).Names={'agej','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_11 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_F(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1); % F (a3_val==1, is that the state s=1 for old)
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9,FnsToEvaluateFn_10,FnsToEvaluateFn_11};

% simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
% simoptions.parallel=2;

%% Table 4
% Uniform lump-sum transfer of 30 percent of the total college tuition cost
% to either yound or old parents, keeping decision rules constant.

% Since this is partial equilibrium and there is no consideration of
% changes in the decision rules the only change is being made
% to the StationaryDist. The following lines implement the necessary change
% to the StationaryDist. The remaining lines then use this modified 
% StationaryDist to calculate the model statistics reported in Table 4.

% Implementing this is actually very awkward. 'Assets' are not part of the
% model, so the only way to approximate a lump-sum transfer appears to be
% to give the lump-sum in the form of human capital (so that wh is changed
% to capture 'transfer of 30 percent of the total college tuition cost').

% Note that the transfer to old households of one-third of the costs of
% college actually has nothing to do with aleviating borrowing constraints,
% as is made clear in the formulation of the model necessitated by the use
% of the VFI toolkit this decision (whether to send child to college) is
% taken by younger household. I do not understand how exactly Restuccia &
% Urrutia (2004) implement the transfers to old HHs in such a way that the
% effect is non-zero. (MAYBE EMAIL THEM TO ASK???)
% I instead just report a version of Table 4 in which the transfer is made
% to all households and the first column reports changes in
% early-education decisions (h) and the second column reports changes in
% college attendence decisions (s). However it should be remembered that
% all the effect comes from transfers to young HHs, and there is zero
% effect on these two decisions from the transfers to old HHs. (Since old
% HHs make neither of the decisions, and since reoptimization of policy
% rules is ruled out by RU2004 as part of the exercise. RU2004 description
% of this as coming from transfers to old comes from their assigning
% (misleadingly) the decision of s to old HHs. Worth observing that their
% setup of having 'intervening' decisions made between periods is very
% common in the Economics literature. Suggests a conflict between the 
% mathematics of when these decisions are made, and how much of Economics
% interprets their timing. Note how in RU2004 formulation, the state space
% of the 'intermediate old HH' is the state space of a young HH, and not 
% the state space of an old HH.)

sizeOfTransfer=0.3*Params.f*Params.nupperbar; % Not certain, but seems the most plausible definition of 'total college tuition cost'
StatDistj001old=StationaryDist.j001;
StationaryDist.j001=zeros(size(StationaryDist.j001),'gpuArray');
h_grid=daz_gridstructure.a_grid.j001(1:n_h); % (in j=1) h is the first endog state
for h_c=1:n_h
    h_val=daz_gridstructure.a_grid.j001(1:n_h);
    h_val=h_val+sizeOfTransfer/Params.w; %So wh will be increased by 'sizeOfTransfer'
    [~, hplustransfer_c]=min(abs(h_grid-h_val));
    StationaryDist.j001(hplustransfer_c,:,:,:,:)=StationaryDist.j001(hplustransfer_c,:,:,:,:)+StatDistj001old(h_c,:,:,:,:); % Near top of grid multiple h_c values will hit the max, but hopefully this is irrelevant as in stationary dist noone is up there anyway (ie. need to set h_grid to ensure this else would generate a bias/error in model results for Table 4).
end

% Use this 'StationaryDist including transfers' to calculate moments for Table 4
% simoptions.parallel=3;
AgeConditionalStats_withTransfer=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
% simoptions.parallel=2;
% The early education expenditure decision is 'e': 5
% Change in early education is 'bhat': 8
% The college enrollment decision is 's': 6
% Private expenditures on college (F-kappaF?): 11-10

Pre_e=gather(AgeConditionalStats(5).Mean(1));
Post_e=gather(AgeConditionalStats_withTransfer(5).Mean(1));
Pre_bhat=gather(AgeConditionalStats(8).Mean(1));
Post_bhat=gather(AgeConditionalStats_withTransfer(8).Mean(1));
Pre_s=gather(AgeConditionalStats(6).Mean(1));
Post_s=gather(AgeConditionalStats_withTransfer(6).Mean(1));
Pre_FminuskappaF=gather(AgeConditionalStats(11).Mean(1))-gather(AgeConditionalStats(10).Mean(1));
Post_FminuskappaF=gather(AgeConditionalStats_withTransfer(11).Mean(1))-gather(AgeConditionalStats_withTransfer(10).Mean(1));

Table4.PercentChange_e=(Post_e-Pre_e)/Pre_e;
Table4.Change_e_PercentOfTransfer=(Post_e-Pre_e)/sizeOfTransfer; % As percent of transfer
Table4.PercentChange_bhat=(Post_bhat-Pre_bhat)/Pre_bhat;
Table4.PercentChange_s=(Post_s-Pre_s)/Pre_s;
Table4.PercentChange_FminuskappaF=(Post_FminuskappaF-Pre_FminuskappaF)/Pre_FminuskappaF;
Table4.Change_FminuskappaF_PercentOfTransfer=(Post_FminuskappaF-Pre_FminuskappaF)/sizeOfTransfer; % As percent of transfer


%Table 4
FID = fopen('./SavedOutput/LatexInputs/RestucciaUrrutia2004_Table4.tex', 'w');
fprintf(FID, 'Uniform Transfer Experiment \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, ' & Young & Old \\\\ \n');
fprintf(FID, 'Transfer to & Parents & Parents \\\\ \n');
fprintf(FID, 'Increase in education expenditures: & %8.2f & %8.2f \\\\ \n ', 100*Table4.Change_e_PercentOfTransfer, 100*Table4.Change_FminuskappaF_PercentOfTransfer);
fprintf(FID, '$\\quad$ (percent of transfer) &  &  \\\\ \n ');
fprintf(FID, 'Parents changing education decisions: & %8.2f & %8.2f \\\\ \n ', 100*Table4.PercentChange_bhat, 100*Table4.PercentChange_s);
fprintf(FID, '$\\quad$ (percent of parents in age group) &  & \\\\ \n ');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Clockwise from top-left, these are percent change in e, F-kappaF, s, and bhat. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


% %% Take a look at college completion probability fn: not in paper, just out of interest
% 
% figure(3)
% bhatgrid=daz_gridstructure.a_grid.j002(n_a(1,2)+1:n_a(1,2)+n_a(2,2));
% qbhat=min(Params.psi0*(1+bhatgrid).^Params.psi1,1);
% plot(qbhat)
% 
% %% Do people go to the 'top' of the bhat and h grids?: not in paper, just out of interest
% 
% figure(4)
% pop_h=sum(sum(sum(sum(StationaryDist.j002,4),5),2),3);
% pop_bhat=sum(sum(sum(sum(StationaryDist.j002,4),5),1),3);
% subplot(2,1,1); plot(cumsum(pop_bhat))
% subplot(2,1,2); plot(cumsum(pop_h))