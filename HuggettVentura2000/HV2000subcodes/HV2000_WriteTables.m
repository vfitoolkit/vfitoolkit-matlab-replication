% Writes HV2000 Tables 2, 3, 4 to ./SavedOutput/LatexInputs/ as standalone
% .tex files (one tabular environment per file).
%
% Table 2: Model parameters (HV2000 Table 2)
% Table 3: Earnings process parameters (HV2000 Table 3)
% Table 4: Descriptive statistics across the 4 models (HV2000 Table 4)
%
% Requires Results.Model{1..4}.AllStats and Params in the workspace.



%% Table 2: model parameters (shared across models)
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table2.tex','w');
fprintf(FID,'\\begin{tabular}{ccccccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,'$\\beta$ & $\\sigma$ & $A$ & $\\alpha$ & $\\delta$ & $g$ & $N$ & $R$ & $s_j$ & $n$ & $\\underline{\\hat{a}}$ \\\\ \n');
fprintf(FID,'\\hline\n');
fprintf(FID,'%.3f & %.2f & %.5f & %.2f & %.2f & %.3f & %d & %d & US 1990 & %.3f & 0 \\\\ \n',...
    Params.beta,Params.sigma,Params.A,Params.alpha,Params.delta,Params.g,Params.J,Params.Jr,Params.n);
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
fclose(FID);

%% Table 3: earnings process parameters (calibration choices)
% Note: this table is independent of abar, so arbirarily doing the following using the *_abar0 version
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table3.tex','w');
fprintf(FID,'\\begin{tabular}{lcccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,'Model & $\\sigma^{2}_{y_{1}}$ & $\\gamma$ & $\\sigma^{2}_{\\varepsilon_{1}}$ & $\\sigma^{2}_{\\varepsilon_{2}}$ \\\\ \n');
fprintf(FID,'\\hline\n');
fprintf(FID,'Model 1 & --- & --- & --- & --- \\\\ \n');
fprintf(FID,'Model 2 & %1.2f & --- & --- & --- \\\\ \n',Results.Model2_abar0.Params.sigma2_y1);
fprintf(FID,'Model 3 & %1.2f & %1.3f & %1.2f & --- \\\\ \n',Results.Model3_abar0.Params.sigma2_y1,Results.Model3_abar0.Params.gamma,Results.Model3_abar0.Params.sigma2_eps1);
fprintf(FID,'Model 4 & %1.2f & %1.3f & %1.2f & %1.2f \\\\ \n',Results.Model4_abar0.Params.sigma2_y1,Results.Model4_abar0.Params.gamma,Results.Model4_abar0.Params.sigma2_eps1,Results.Model4_abar0.Params.sigma2_eps2);
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
fclose(FID);

%% Table 4: descriptive statistics
% Top X% income share = 1 - LorenzCurve(P), where P = length-X (the
% LorenzCurve is the share of income going to the bottom P percentiles).
%
% S/Y is NET national saving, SnetY, not the aggregate of agent-level saving.
% AllStats.Savings.Mean aggregates (1+g)*aprime-a over the whole stationary
% distribution, which counts the assets of the agents who die this period --
% those become bequests, not next-period capital -- and so overstates national
% saving by about 50% (0.134 against HV2000's 0.090). Subtracting the
% growth-adjusted assets of the deceased fixes it:
%   SnetY = [ Savings - (1+g)*AccBeq/(1+r(1-tau)) ] / Y
% where AccBeq/(1+r*(1-tau)) undoes the return factor in the AccBeq definition
% to recover the undiscounted dead assets. This reproduces the steady-state
% identity (n+g+n*g)*K/Y to all printed digits in all eight rows, which is what
% HV2000's column satisfies. NB the agent-level measure is the correct one for
% the saving rates of Tables 5-9 (HV2000 p.380 is a group's total saving over
% its total income), so those are left alone.
% Transfer wealth (Kotlikoff-Summers decomposition) not yet computed.
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table4.tex','w');
fprintf(FID,'\\begin{tabular}{lccccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,' & & &                    & Trans  & Income & \\multicolumn{4}{l}{Percentage of income received} \\\\ \n');
fprintf(FID,' & & &                    & wealth & Gini   & \\multicolumn{4}{l}{by the top} \\\\  \\cline{7-10} \n');
fprintf(FID,' & K/Y & S/Y & $r$ (\\%%) & (\\%%) &        & Top 1\\%% & Top 5\\%% & Top 10\\%% & Top 20\\%% \\\\ \n');
fprintf(FID,'\\hline\n');
R=Results.Model1_abar0; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 1, $\\underline{a}=0$  & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model1_abarw; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 1, $\\underline{a}=-w$ & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model2_abar0; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 2, $\\underline{a}=0$  & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model2_abarw; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 2, $\\underline{a}=-w$ & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model3_abar0; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 3, $\\underline{a}=0$  & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model3_abarw; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 3, $\\underline{a}=-w$ & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model4_abar0; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 4, $\\underline{a}=0$  & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));

R=Results.Model4_abarw; A=R.AllStats; P=R.Params; Y=P.A*A.K.Mean^P.alpha*A.L.Mean^(1-P.alpha); L=length(A.Income.LorenzCurve); SnetY=(A.Savings.Mean-(1+P.g)*A.AccBeq.Mean/(1+P.r*(1-P.tau)))/Y;
fprintf(FID,'Model 4, $\\underline{a}=-w$ & %.2f & %.3f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ \n',A.K.Mean/Y,SnetY,100*P.r,R.TransferWealthPct,A.Income.Gini,100*(1-A.Income.LorenzCurve(round(L*0.99))),100*(1-A.Income.LorenzCurve(round(L*0.95))),100*(1-A.Income.LorenzCurve(round(L*0.90))),100*(1-A.Income.LorenzCurve(round(L*0.80))));
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
% Note sits OUTSIDE the tabular: a p{} multicolumn wider than the table's own
% natural width pushes all the slack into the final column, visibly detaching it.
fprintf(FID,'\\\\[4pt]\n');
fprintf(FID,'\\parbox{15cm}{\\footnotesize Note: with the lower limit of assets $\\underline{a}=-w$, Models 2--4 all have some households at the borrowing corner with negative income (earnings plus $ra$ plus $T$ can fall below zero). The income Gini and the Lorenz shares in that row are therefore computed over a variable that is not everywhere positive, and are not the standard Lorenz-based measure; they are reported because the affected mass is small and each sits under $0.01$ above the corresponding $\\underline{a}=0$ row, but they should be read with that caveat. Model 1 is unaffected because its earnings are deterministic.}\n');
fclose(FID);

%% Table 5: saving rates at multiples of mean income (HV2000 Table 5)
% Layout matches HV2000 Table 5: income multiples are ROWS, models are COLUMNS.
% US data column is the "US^a" column of HV2000 Table 5 (averages from Kuznets 1953 and Projector 1968).
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table5.tex','w');
% Model 1 is reported at a_underbar=0 only, as in HV2000 Table 5 (its Model 1 block
% has a single column; only Models 2-4 get both borrowing limits).
fprintf(FID,'\\begin{tabular}{lcccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,' & & Model 1 & \\multicolumn{2}{c}{Model 2} & \\multicolumn{2}{c}{Model 3} & \\multicolumn{2}{c}{Model 4} \\\\ \n');
fprintf(FID,'\\cline{3-3} \\cline{4-5} \\cline{6-7} \\cline{8-9} \n');
fprintf(FID,'Income multiple & US & $\\underline{a}=0$ & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ \\\\ \n');
fprintf(FID,'\\hline\n');
s1_0=100*Results.Model1_abar0.savingRates;
s2_0=100*Results.Model2_abar0.savingRates; s2_w=100*Results.Model2_abarw.savingRates;
s3_0=100*Results.Model3_abar0.savingRates; s3_w=100*Results.Model3_abarw.savingRates;
s4_0=100*Results.Model4_abar0.savingRates; s4_w=100*Results.Model4_abarw.savingRates;
HV2000_multiples_labels={'0.25','0.50','0.75','1.0','1.5','2.0','3.0','4.0','7.0','10.0'};
HV2000_us=[-19.3,-1.3,4.8,7.9,13.0,16.5,22.4,27.1,37.3,39.2];
for ii=1:length(HV2000_multiples_labels)
    fprintf(FID,'%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\ \n', ...
        HV2000_multiples_labels{ii},HV2000_us(ii),s1_0(ii),s2_0(ii),s2_w(ii),s3_0(ii),s3_w(ii),s4_0(ii),s4_w(ii));
end
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
% Note sits OUTSIDE the tabular, as for Table 4.
fprintf(FID,'\\\\[4pt]\n');
fprintf(FID,'\\parbox{15cm}{\\footnotesize Note: Model 1 has deterministic earnings, so no household reaches an income above roughly $1.5$ times the mean and the income bands at the $2$ multiple and above are empty. An empty band gives $0/0$, reported here as NaN rather than as a saving rate. HV2000 report no Model 1 entries at these multiples either.}\n');
fclose(FID);

%% Figure 2: HV2000 age-earnings profile (ybar_j)
% Plots Params.ybar_j (set in HuggettVentura2000_USdata.m from SSA 1990 medians x BLS LFP).
% Compare against HV2000 Fig 2 (page 374).
fig2=figure(2);
plot(Params.agejshifter+(1:Params.J),Params.ybar_j,'-k+','LineWidth',1);
xlabel('Age'); ylabel('Earnings (ratio to overall mean)');
title('HV2000 Figure 2: age-earnings profile');
xlim([20 100]); ylim([0 1.05*max(Params.ybar_j)]);
saveas(fig2,'./SavedOutput/Graphs/HV2000_Figure2.png');
fprintf('Wrote Figure 2 to ./SavedOutput/Graphs/HV2000_Figure2.png\n');
