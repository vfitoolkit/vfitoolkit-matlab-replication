% HV2000 Figure 5: Model 2 with NO social security (haveSS=0, a_underbar=0).
% Same structure as Figure 3 but for the no-SS scenario.
% Panel a: saving rates by age × income group (Lowest 10%, 50-60%, Upper 10%).
% Panel b: age-income distribution by quantile.
%
% Called from inside HV2000_Model2.m when haveSS=0 (during HV2000_Table7 / Tables89 baseline-noSS pass).

fprintf('\n========== HV2000 Figure 5 (Model 2, NO social security, a_underbar=0) ==========\n');

simoptions_Fig=simoptions2;
simoptions_Fig.nquantiles=100;
simoptions_Fig.conditionalrestrictions=struct();
AllStatsFig=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_Fig);
incQcuts=AllStatsFig.Income.QuantileCutoffs;
inc10=incQcuts(11);
inc50=incQcuts(51);
inc60=incQcuts(61);
inc90=incQcuts(91);

simoptions_Fig.conditionalrestrictions=struct();
% GPU arrayfun cannot compile closures that capture outer-workspace vars
% ("Use of functional workspace is not supported"); bake the percentile cutoffs
% as numeric literals via str2func+sprintf.
simoptions_Fig.conditionalrestrictions.Low10= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     'HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)<=%.16g'], inc10));
simoptions_Fig.conditionalrestrictions.Mid50_60= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     '(HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)>=%.16g) ' ...
     '&& (HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)<=%.16g)'], inc50, inc60));
simoptions_Fig.conditionalrestrictions.High10= str2func(sprintf( ...
    ['@(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ' ...
     'HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)>=%.16g'], inc90));

AgeCondFig=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_Fig);

age=Params.agejshifter+(1:N_j);
sav_low10 =AgeCondFig.Low10.Savings.Mean   ./ AgeCondFig.Low10.Income.Mean;
sav_mid   =AgeCondFig.Mid50_60.Savings.Mean./ AgeCondFig.Mid50_60.Income.Mean;
sav_high10=AgeCondFig.High10.Savings.Mean  ./ AgeCondFig.High10.Income.Mean;

if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

%% Single figure with two stacked panels (matches HV2000 Fig 5)
fig5=figure(5);

%% Panel a: saving rates by age x income group
subplot(2,1,1);
plot(age,sav_low10,'-b','LineWidth',1); hold on;
plot(age,sav_mid,'-g','LineWidth',1);
plot(age,sav_high10,'-r','LineWidth',1);
xlabel('Age'); ylabel('Saving rate');
title('(a) Saving rates (Model 2, no SS, $\underline{a}=0$)','Interpreter','latex');
legend({'Lowest 10%','50%-60%','Upper 10%'},'Location','best');
yline(0,'k-','HandleVisibility','off');  % HV2000 draw a zero line on this panel; HandleVisibility off keeps it out of the legend (this legend() is called before the yline, so without it MATLAB appends a spurious 'data1' entry)
ylim([-3 0.5]); yticks(-3:0.5:0.5);  % match HV2000's axis
xlim([20 100]); xticks(20:10:100);
hold off;

%% Panel b: age-conditional income at percentiles + mean
simoptions_FigB=simoptions2;
simoptions_FigB.nquantiles=100;
AgeCondB=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_FigB);
incQj=AgeCondB.Income.QuantileCutoffs;
inc10j  =incQj(11,:);
inc25j  =incQj(26,:);
inc50j  =incQj(51,:);
incMeanj=AgeCondB.Income.Mean;

%% Panel b: age-income distribution
subplot(2,1,2);
plot(age,inc10j,'-c','LineWidth',1); hold on;
plot(age,inc25j,'-b','LineWidth',1);
plot(age,inc50j,'-g','LineWidth',1);
plot(age,incMeanj,'-r','LineWidth',1);
xlabel('Age'); ylabel('Income');
title('(b) Age-income distribution (Model 2, no SS, $\underline{a}=0$)','Interpreter','latex');
legend({'10% Quantile','25% Quantile','50% Quantile','Mean'},'Location','best');
ylim([0 2.5]); yticks(0:0.5:2.5);  % match HV2000's axis
xlim([20 100]); xticks(20:10:100);
hold off;

saveas(fig5,'./SavedOutput/Graphs/HV2000_Figure5.png');

fprintf('Wrote Figure 5 to ./SavedOutput/Graphs/HV2000_Figure5.png\n');
