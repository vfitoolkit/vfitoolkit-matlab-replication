% HV2000 Figure 7: saving rates by age x income quartile, Model 4 (a_underbar=0).
% Compares qualitatively with Attanasio (1994) Table 2.10. The paper omits the
% lowest quartile (always negative); we plot Q2, Q3, Q4 and the Total Sample.
%
% Uses per-agent SavingRate=Savings/Income (HV2000_SavingRateFn_M4) and takes the
% age-conditional MEDIAN per income quartile (the paper's preferred summary).
%
% Workspace expected (after Model 4 with borrowFlag=0, BEFORE the simoptions restore).

fprintf('\n========== HV2000 Figure 7 (Model 4, Attanasio comparison) ==========\n');

simoptions_Fig=simoptions4;
simoptions_Fig.nquantiles=100;
simoptions_Fig.conditionalrestrictions=struct();
AllStatsFig=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_M4,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions_Fig);
incQcuts=AllStatsFig.Income.QuantileCutoffs;
inc25=incQcuts(26);
inc50=incQcuts(51);
inc75=incQcuts(76);

simoptions_Fig.conditionalrestrictions=struct();
% GPU arrayfun cannot compile closures that capture outer-workspace vars
% ("Use of functional workspace is not supported"); bake the quartile cutoffs
% as numeric literals via str2func+sprintf.
simoptions_Fig.conditionalrestrictions.Q2= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     '(HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)>%.16g) ' ...
     '&& (HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)<=%.16g)'], inc25, inc50));
simoptions_Fig.conditionalrestrictions.Q3= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     '(HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)>%.16g) ' ...
     '&& (HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)<=%.16g)'], inc50, inc75));
simoptions_Fig.conditionalrestrictions.Q4= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     'HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)>%.16g'], inc75));

AgeCondFig=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_M4,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions_Fig);

age=Params.agejshifter+(1:N_j);
sav_total=AgeCondFig.SavingRate.Median;
sav_q2   =AgeCondFig.Q2.SavingRate.Median;
sav_q3   =AgeCondFig.Q3.SavingRate.Median;
sav_q4   =AgeCondFig.Q4.SavingRate.Median;

% HV2000 Fig 7 covers ages ~21-75; restrict plot to that window for comparability with Attanasio Table 2.10.
ageRange=21:75;
ageIdx=ageRange-Params.agejshifter;

if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

fig7=figure(7);
plot(ageRange,100*sav_total(ageIdx),'-k','LineWidth',1); hold on;
plot(ageRange,100*sav_q2(ageIdx),'-b','LineWidth',1);
plot(ageRange,100*sav_q3(ageIdx),'-g','LineWidth',1);
plot(ageRange,100*sav_q4(ageIdx),'-r','LineWidth',1);
xlabel('Age group'); ylabel('Saving rate (%)');
title('Saving rates by income quartile (Model 4, $\underline{a}=0$) - cf Attanasio (1994) Table 2.10','Interpreter','latex');
legend({'Total Sample','Second Quartile','Third Quartile','Fourth Quartile'},'Location','best');
hold off;
saveas(fig7,'./SavedOutput/Graphs/HV2000_Figure7.png');

fprintf('Wrote Figure 7 to ./SavedOutput/Graphs/\n');
