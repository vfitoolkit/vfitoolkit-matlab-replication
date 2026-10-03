% HV2000 Figure 6: alternative bequest treatments in Model 2.
% Two panels: (a) a_underbar=0; (b) a_underbar=-w.
% Two lines per panel:
%   "Equal Transfers"     - baseline Model 2 (single population-wide T), from the doPart(1) baseline on disk
%   "Different Transfers" - PType variant with bequests by ability type, from Results.Model2_ptype_abar{0,w}
%
% Both saving-rate vectors are computed via simoptions.conditionalrestrictions on income bins
% so the numbers are directly comparable.

fprintf('\n========== HV2000 Figure 6 (alternative bequest treatments, Model 2) ==========\n');

HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];
HV2000_multiples_labels={'0.25','0.5','0.75','1','1.5','2','3','4','7','10'};

% The "Equal Transfers" baseline is Model 2 as solved in doPart(1). Do NOT read it
% from the live Results: HV2000_Table7 (haveSS=0) and HV2000_Tables89 (within-age
% averaging) both overwrite Results.Model2_abar{0,w}, and when doPart(1)==0 those
% fields do not exist at all. Read the doPart(1) baseline from disk instead.
if ~exist('./SavedOutput/HV2000_Results_Part1.mat','file')
    error('HV2000_Figure6 needs ./SavedOutput/HV2000_Results_Part1.mat (the doPart(1) baseline); run doPart(1) to completion first')
end
BaseResults=load('./SavedOutput/HV2000_Results_Part1.mat','Results');
BaseResults=BaseResults.Results;

fig6=figure(6);
% Panel a: a_underbar=0
subplot(2,1,1);
% HV2000 Fig 6 plots the multiples as evenly spaced CATEGORIES with a linear axis,
% not on a log scale, so plot against 1:10 and label the ticks with the multiples.
xcat=1:length(HV2000_multiples);
plot(xcat,100*BaseResults.Model2_abar0.savingRates,'-+b','LineWidth',1); hold on;
plot(xcat,100*Results.Model2_ptype_abar0.savingRates,'-xr','LineWidth',1);
yline(0,'k-');
xlim([1 length(HV2000_multiples)]); xticks(xcat); xticklabels(HV2000_multiples_labels);
ylim([-5 35]); yticks(-5:5:35);  % match HV2000's axis
xlabel('Income Multiple'); ylabel('Saving Rate (%)');
title('(a) $\underline{a}=0$','Interpreter','latex');
% northwest, not 'best': the series rise to the top right, so 'best' put the
% legend box straight over the 7x and 10x points and made them look missing.
legend({'Equal Transfers','Different Transfers'},'Location','northwest');
hold off;
% Panel b: a_underbar=-w
subplot(2,1,2);
plot(xcat,100*BaseResults.Model2_abarw.savingRates,'-+b','LineWidth',1); hold on;
plot(xcat,100*Results.Model2_ptype_abarw.savingRates,'-xr','LineWidth',1);
yline(0,'k-');
xlim([1 length(HV2000_multiples)]); xticks(xcat); xticklabels(HV2000_multiples_labels);
ylim([-20 40]); yticks(-20:10:40);  % match HV2000's axis
xlabel('Income Multiple'); ylabel('Saving Rate (%)');
title('(b) $\underline{a}=-w$','Interpreter','latex');
legend({'Equal Transfers','Different Transfers'},'Location','northwest');
hold off;

saveas(fig6,'./SavedOutput/Graphs/HV2000_Figure6.png');
