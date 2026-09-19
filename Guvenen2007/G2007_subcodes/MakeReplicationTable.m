% MakeReplicationTable.m
% Stage (d): assemble the baseline replication table + HIP-vs-RIP comparison from the
% saved life-cycle profiles. Pure post-processing (CPU): loads RunBaseline_results.mat
% (HIP, lambda in {0,0.4,0.62,1}, each calibrated to W/Y=4 -> delta*~0.95) and
% RunRIP_results.mat (RIP, W/Y=4 -> delta*~0.95). No solve, no GPU.
%
% Reports the ACCURATE solution (per the 2026-09-14 decision): W/Y=4 with the
% grid-converged delta*, NOT the paper's delta=0.966. The var(log c) rises therefore
% sit above the paper's Fig-8 numbers by construction (see NOTES.md / tex App A.4);
% the paper's numbers are shown alongside for reference, not as a target we retune to.
%
% Outputs: console+diary table, a LaTeX table fragment (SavedOutput/replication_table.tex),
% and the consumption-inequality figures (SavedOutput/Graphs/): Fig 7 (RIP) and Fig 8
% (HIP by lambda), matching Guvenen's separate Figures 7 and 8.

clearvars -except doPart
if exist('MakeReplicationTable_diary.txt','file'); delete('MakeReplicationTable_diary.txt'); end
if ~exist('SavedOutput/Graphs','dir'); mkdir('SavedOutput/Graphs'); end
diary('MakeReplicationTable_diary.txt');
fprintf('=== MakeReplicationTable.m run %s ===\n',char(datetime('now')));

Jwork=40; agevec=25:95; ageswork=25:65; % 25..65 inclusive is the working-life span for var(log c) rise

%% Load
B=load('RunBaseline_results.mat');
R=load('RunRIP_results.mat');
lambdavec=gather(B.lambdavec);
nlam=numel(lambdavec);

% paper Fig-8 reference rises (dashes where the paper reports no number)
paperFig8=[0.32, NaN, 0.21, 0.14]; % lambda = 0, 0.4, 0.62, 1
paperRIP=0.26;                     % paper Fig 7 RIP rise (~26 log points)

%% Extract HIP profiles per lambda
HIP=struct();
for i=1:nlam
    HIP(i).lambda=gather(B.results(i).lambda);
    HIP(i).delta=gather(B.results(i).delta);
    HIP(i).WYtot=gather(B.results(i).WYtot);
    HIP(i).varlogc=gather(B.results(i).varlogc);
    HIP(i).varlogy=gather(B.results(i).varlogy);
    HIP(i).meanc=gather(B.results(i).meanc);
    HIP(i).meany=gather(B.results(i).meany);
    HIP(i).meana=gather(B.results(i).meana);
    HIP(i).rise=HIP(i).varlogc(Jwork+1)-HIP(i).varlogc(1);
    HIP(i).riseY=HIP(i).varlogy(Jwork+1)-HIP(i).varlogy(1);
end

%% Extract RIP profiles
RIP.delta=gather(R.Params.delta);
RIP.varlogc=gather(R.AgeConditionalStats.logc.Variance);
RIP.varlogy=gather(R.AgeConditionalStats.logy.Variance);
RIP.meanc=gather(R.AgeConditionalStats.cons.Mean);
RIP.meany=gather(R.AgeConditionalStats.income.Mean);
RIP.meana=gather(R.AgeConditionalStats.assets.Mean);
RIP.rise=RIP.varlogc(Jwork+1)-RIP.varlogc(1);
% RIP W/Y (total income) from its own profiles
rmewj=ones(1,R.N_j)/R.N_j; rr=gather(R.Params.r);
RIP.WYtot=sum(rmewj.*RIP.meana)/sum(rmewj.*(RIP.meany+rr*RIP.meana));

%% Console table
fprintf('\n=== Baseline replication: within-cohort consumption-inequality rise (ages 25->65) ===\n');
fprintf('Model      lambda   delta*    W/Y    var(logc)_25  var(logc)_65   RISE    paper Fig8\n');
for i=1:nlam
    pf=paperFig8(i); pfs='   --'; if ~isnan(pf); pfs=sprintf('%6.2f',pf); end
    fprintf('HIP        %5.2f   %6.4f  %5.2f      %6.3f        %6.3f     %6.3f   %s\n',...
        HIP(i).lambda,HIP(i).delta,HIP(i).WYtot,HIP(i).varlogc(1),HIP(i).varlogc(Jwork+1),HIP(i).rise,pfs);
end
fprintf('RIP          --    %6.4f  %5.2f      %6.3f        %6.3f     %6.3f   %6.2f\n',...
    RIP.delta,RIP.WYtot,RIP.varlogc(1),RIP.varlogc(Jwork+1),RIP.rise,paperRIP);
fprintf('US data (Deaton-Paxson / Guvenen Fig 1): rise ~= 0.21 log points\n');

%% var(log y) cross-check (should track the analytic HIP decomposition regardless of delta)
fprintf('\n=== var(log y) rise (income inequality, ages 25->65), cross-check ===\n');
for i=1:nlam
    fprintf('HIP lambda=%.2f: var(logy) rise = %.3f\n',HIP(i).lambda,HIP(i).riseY);
end
fprintf('RIP: var(logy) rise = %.3f\n',RIP.varlogy(Jwork+1)-RIP.varlogy(1));

%% Carroll-Summers co-movement (aggregate, all-sample baseline): does mean consumption
%% growth parallel mean income growth over the life cycle? (Section III.C general fact.)
%% NB: the DISCRIMINATING Figure-11 test (growth DIFFERS by education under HIP, not under
%% RIP) needs the education-group calibrations -> RunEducationGroups.m (separate).
% windows: 25->55 (paper's) and 25->64 (end of working life); both within working life so
% income growth is clean (25->65 would cross into the pension drop and go negative)
fprintf('\n=== Carroll-Summers co-movement (all-sample baseline), total growth ===\n');
fprintf('Model         cons 25-55   inc 25-55    cons 25-64   inc 25-64\n');
for i=1:nlam
    cg55=HIP(i).meanc(31)/HIP(i).meanc(1)-1; yg55=HIP(i).meany(31)/HIP(i).meany(1)-1;
    cg64=HIP(i).meanc(40)/HIP(i).meanc(1)-1; yg64=HIP(i).meany(40)/HIP(i).meany(1)-1;
    fprintf('HIP lam=%.2f     %6.1f%%      %6.1f%%       %6.1f%%     %6.1f%%\n',...
        HIP(i).lambda,100*cg55,100*yg55,100*cg64,100*yg64);
end
cg55=RIP.meanc(31)/RIP.meanc(1)-1; yg55=RIP.meany(31)/RIP.meany(1)-1;
cg64=RIP.meanc(40)/RIP.meanc(1)-1; yg64=RIP.meany(40)/RIP.meany(1)-1;
fprintf('RIP             %6.1f%%      %6.1f%%       %6.1f%%     %6.1f%%\n',100*cg55,100*yg55,100*cg64,100*yg64);
fprintf('(Carroll-Summers: consumption growth should PARALLEL income growth over the life cycle)\n');
fprintf('(the discriminating Fig-11 test -- growth DIFFERS by education under HIP, not RIP -- needs RunEducationGroups.m)\n');

%% LaTeX table fragment
fid=fopen('SavedOutput/replication_table.tex','w');
fprintf(fid,'%% Auto-generated by MakeReplicationTable.m -- do not edit by hand\n');
fprintf(fid,'\\begin{tabular}{lccccc}\n\\toprule\n');
fprintf(fid,'Model & $\\lambda$ & $\\delta^*$ & $W/Y$ & rise in var$(\\log c)$ & paper \\\\\n');
fprintf(fid,'\\midrule\n');
for i=1:nlam
    pf=paperFig8(i); pfs='--'; if ~isnan(pf); pfs=sprintf('%.2f',pf); end
    fprintf(fid,'HIP & %.2f & %.4f & %.2f & %.3f & %s \\\\\n',...
        HIP(i).lambda,HIP(i).delta,HIP(i).WYtot,HIP(i).rise,pfs);
end
fprintf(fid,'RIP & -- & %.4f & %.2f & %.3f & %.2f \\\\\n',RIP.delta,RIP.WYtot,RIP.rise,paperRIP);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('\nWrote SavedOutput/replication_table.tex\n');

%% Table 2 (income-inequality decomposition) -- analytic, all-sample HIP params.
% This is a pure income-process object (no solved model needed); the same numbers are
% checked in tests/test_kalman.m. Written as a LaTeX fragment (table2_decomposition.tex)
% shown against the paper's original in the tex. Uses UNCONDITIONAL Table-1 row-2 values
% (sab per NOTES.md); experience t=age-25.
t2.rho=0.821; t2.s2a=0.022; t2.s2eps=0.047; t2.s2eta=0.029; t2.s2b=0.00038; t2.sab=-0.00045;
t2ages=[30 35 45 55 65];
fid=fopen('SavedOutput/table2_decomposition.tex','w');
fprintf(fid,'%% Auto-generated by MakeReplicationTable.m -- do not edit by hand (overwritten each run)\n');
fprintf(fid,'\\begin{tabular}{lcccc}\n\\toprule\n');
fprintf(fid,'Age & (1) $\\sigma^2_\\alpha+\\sigma^2_\\varepsilon$ & (2) $\\frac{1-\\rho^{2t+1}}{1-\\rho^2}\\sigma^2_\\eta$ & (3) $2\\sigma_{\\alpha\\beta}t+\\sigma^2_\\beta t^2$ & (4) $\\frac{(3)}{(1)+(2)+(3)}$ \\\\\n');
fprintf(fid,'\\midrule\n');
for aa=t2ages
    t=aa-25;
    c1=t2.s2a+t2.s2eps;
    c2=t2.s2eta*(1-t2.rho^(2*t+1))/(1-t2.rho^2);
    c3=2*t2.sab*t+t2.s2b*t^2;
    c4=c3/(c1+c2+c3);
    fprintf(fid,'%d & %.3f & %.3f & %.3f & %.3f \\\\\n',aa,c1,c2,c3,c4);
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('Wrote SavedOutput/table2_decomposition.tex\n');

% US-data reference used in BOTH Fig 7 and Fig 8 (Guvenen's plots overlay the
% empirical consumption-inequality profile). We prefer the recomputed CEX line
% (Guvenen2007_DataWork.m -> Guvenen2007_DataWork.mat: Krueger-Perri nondurables,
% Deaton-Paxson age-cohort decomposition, 95% bootstrap band). If that .mat is
% not present yet, fall back to a linear stand-in (rise ~0.21 from ~0.15 at 25).
% Either way we do NOT reproduce the Storesletten et al. (2004) line from Fig 7.
haveData=exist('SavedOutput/Guvenen2007_DataWork.mat','file');
if haveData
    DW=load('SavedOutput/Guvenen2007_DataWork.mat'); DL=DW.DataLine;
    usAge=DL.ages(:)'; usData=DL.varlogc(:)'; usLo=DL.ciLow(:)'; usHi=DL.ciHigh(:)';
    usLbl='US data (CEX, DP)';
else
    usAge=ageswork; usData=0.15+(0.21/40)*(ageswork-25); usLo=[]; usHi=[];
    usLbl='US data (approx)';
    fprintf('NOTE: Guvenen2007_DataWork.mat not found -- using linear US-data stand-in (run doPart(1))\n');
end

%% Figure 7: RIP model -- var(log c) over the life cycle + US data reference
fig=figure(1); clf; hold on;
if haveData; fill([usAge,fliplr(usAge)],[usLo,fliplr(usHi)],[0.85 0.85 0.85],'EdgeColor','none','HandleVisibility','off'); end
plot(ageswork,RIP.varlogc(1:Jwork+1),'k--','LineWidth',1.5);
plot(usAge,usData,'ko','MarkerSize',4,'MarkerFaceColor','k');
xlabel('age'); ylabel('var(log c)');
legend({'RIP',usLbl},'Location','northwest');
title('Fig 7: consumption inequality, RIP model (accurate solution, W/Y=4)');
saveas(fig,'SavedOutput/Graphs/Replication_Fig7_RIP.png');
fprintf('Wrote SavedOutput/Graphs/Replication_Fig7_RIP.png\n');

%% Figure 8: HIP model by lambda -- var(log c) over the life cycle + US data reference
fig=figure(2); clf; hold on;
if haveData; fill([usAge,fliplr(usAge)],[usLo,fliplr(usHi)],[0.85 0.85 0.85],'EdgeColor','none','HandleVisibility','off'); end
cols=lines(nlam);
for i=1:nlam
    plot(ageswork,HIP(i).varlogc(1:Jwork+1),'-','Color',cols(i,:),'LineWidth',1.5);
end
plot(usAge,usData,'ko','MarkerSize',4,'MarkerFaceColor','k');
xlabel('age'); ylabel('var(log c)');
legend([arrayfun(@(i) sprintf('HIP \\lambda=%.2f',lambdavec(i)),1:nlam,'UniformOutput',false),...
    {usLbl}],'Location','northwest');
title('Fig 8: consumption inequality, HIP model (accurate solution, W/Y=4)');
saveas(fig,'SavedOutput/Graphs/Replication_Fig8_HIP.png');
fprintf('Wrote SavedOutput/Graphs/Replication_Fig8_HIP.png\n');

%% Figure: HIP baseline (lambda=0.62) vs RIP -- consumption & income profiles
ib=find(abs(lambdavec-0.62)<1e-9,1);
fig=figure(3); clf; hold on;
plot(agevec,HIP(ib).meanc,'b-','LineWidth',1.5);
plot(agevec,RIP.meanc,'r-','LineWidth',1.5);
plot(agevec,HIP(ib).meany,'b--','LineWidth',1);
plot(agevec,RIP.meany,'r--','LineWidth',1);
xlabel('age'); legend('HIP mean c','RIP mean c','HIP mean y','RIP mean y','Location','best');
title('HIP (\lambda=0.62) vs RIP: consumption and income');
saveas(fig,'SavedOutput/Graphs/Replication_meancy_HIPvsRIP.png');
fprintf('Wrote SavedOutput/Graphs/Replication_meancy_HIPvsRIP.png\n');

%% Concavity of var(log c) over ages 25-55 (Guvenen's "second fact"; tex Appendix C).
% t^2 coefficient of a quadratic fit; >0 = non-concave (convex), <0 = concave.
fprintf('\n=== Concavity of var(log c), ages 25-55 (t^2 coeff: >0 non-concave/convex) ===\n');
prof={}; nm={};
if haveData
    prof{end+1}=DL.varlogc(:);    nm{end+1}='Data household';
    prof{end+1}=DL.varlogc_pe(:); nm{end+1}='Data per-equiv';
end
prof{end+1}=RIP.varlogc(:); nm{end+1}='RIP';
for i=1:nlam; prof{end+1}=HIP(i).varlogc(:); nm{end+1}=sprintf('HIP lam=%.2f',lambdavec(i)); end
for k=1:numel(prof)
    v=prof{k}; p=polyfit((0:30)',v(1:31),2);
    fprintf('%-16s t2=%+.5f  (%s)\n',nm{k},p(1),ternary(p(1)>0,'non-concave','concave'));
end

fprintf('\n=== MakeReplicationTable done ===\n');
diary off;

function s=ternary(c,a,b); if c; s=a; else; s=b; end; end
