% Replicates the results (and more) of Aiyagari (1994)
% Note: The algorithm used for calculating the general equilibrium is not what you would want to use normally for solving this model. It finds the general equilibrium using a discrete grid on interest rates, rather
% than just solving the fixed-point problem on interest rates directly by using optimization. This is done for robustness reasons and because these codes were written as part of a paper on
% a convergent algorithm for solving Bewley-Huggett-Aiyagari-models.

%% Set some basic variables

n_k=2^10;%2^9;
n_z=27; %21;
n_p=451; %151

%Parameters
Params.beta=0.96; %Model period is one-sixth of a year
Params.alpha=0.36;
Params.delta=0.08;

mu_vec=[1,3,5]; % {1,3,5}
sigma_vec=[0.2,0.4]; % {0.2,0.4}
rho_vec=[0,0.3,0.6,0.9]; % {0,0.3,0.6,0.9}

Params.q=3; %Hyperparameter for Tauchen method to approximate AR(1). Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3

%% Some Toolkit options
vfoptions.lowmemory=1;
simoptions.ncores=feature('numcores'); % Number of CPU cores
simoptions.parallel=4; % Use sparse matrices

%% Solve the model a bunch of times
Table1=zeros(2,3,2,4);
Table2=zeros(2,3,2,4);
Table3=zeros(6,3,2,4); % My own creation. Inequality stats.
TimeTable=zeros(1,3,2,4);
for mu_c=1:3
    Params.mu=mu_vec(mu_c);
    for sigma_c=1:2
        Params.sigma=sigma_vec(sigma_c);
        for rho_c=1:4
            Params.rho=rho_vec(rho_c);
            fprintf('Current iteration mu_c=%d, sigma_c=%d, rho_c=%d \n', mu_c,sigma_c,rho_c')
            tic;
            OutputVector=Aiyagari1994_Fn(n_k,n_z,n_p,Params, vfoptions, simoptions);
            time1=toc
            % OutputVector=[sqrt(s_variance), s_corr, p_eqm*100, aggsavingsrate, EarningsGini, IncomeGini, WealthGini, EarningsParetoCoeff,IncomeParetoCoeff, WealthParetoCoeff];
            Table1(:,mu_c,sigma_c,rho_c)=[OutputVector(1); OutputVector(2)];
            Table2(:,mu_c,sigma_c,rho_c)=[OutputVector(3); OutputVector(4)];
            Table3(:,mu_c,sigma_c,rho_c)=[OutputVector(5); OutputVector(6); OutputVector(7); OutputVector(8); OutputVector(9); OutputVector(10)];            

            TimeTable(1,mu_c,sigma_c,rho_c)=time1;
        end
    end
end

save ./SavedOutput/Aiyagari1994Tables.mat Table1 Table2 Table3 TimeTable

%% Save output as Latex tables summarizing the results

FilenameString=['./SavedOutput/LatexInputs/Aiyagari1994_Table1.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Markov Chain Approximation to the Labour Endowment Shock \\\\ \n');
fprintf(FID, 'Markov Chain $\\sigma$/Markov Chain $\\rho$ \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|cccc} \\hline \\hline \n');
fprintf(FID, ' $\\sigma$/$\\rho$ & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \\hline \n', rho_vec(1), rho_vec(2), rho_vec(3), rho_vec(4));
fprintf(FID, ' %1.1f & %1.2f/%1.2f & %1.2f/%1.2f & %1.2f/%1.2f & %1.2f/%1.2f \\\\ \n', sigma_vec(1), Table1(1,1,1,1), Table1(2,1,1,1), Table1(1,1,1,2), Table1(2,1,1,2), Table1(1,1,1,3), Table1(2,1,1,3), Table1(1,1,1,4), Table1(2,1,1,4));
fprintf(FID, ' %1.1f & %1.2f/%1.2f & %1.2f/%1.2f & %1.2f/%1.2f & %1.2f/%1.2f \\\\ \n', sigma_vec(2), Table1(1,1,2,1), Table1(2,1,2,1), Table1(1,1,2,2), Table1(2,1,2,2), Table1(1,1,2,3), Table1(2,1,2,3), Table1(1,1,2,4), Table1(2,1,2,4));
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 1 of Aiyagari (1994) using grid sizes $ n_k=%d $, $ n_z=%d $, $ n_p=%d $ \n', n_k, n_z, n_p);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

FilenameString=['./SavedOutput/LatexInputs/Aiyagari1994_Table2.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccc} \\hline \\hline \n');
fprintf(FID, '\\multicolumn{4}{l}{A. Net Return to Capital in \\%%/Aggregate savings rate in \\%% ($\\sigma=0.2$)} \\\\ \n');
fprintf(FID, ' $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(1), Table2(1,1,1,1),Table2(2,1,1,1), Table2(1,2,1,1),Table2(2,2,1,1), Table2(1,3,1,1),Table2(2,3,1,1));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(2), Table2(1,1,1,2),Table2(2,1,1,2), Table2(1,2,1,2),Table2(2,2,1,2), Table2(1,3,1,2),Table2(2,3,1,2));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(3), Table2(1,1,1,3),Table2(2,1,1,3), Table2(1,2,1,3),Table2(2,2,1,3), Table2(1,3,1,3),Table2(2,3,1,3));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(4), Table2(1,1,1,4),Table2(2,1,1,4), Table2(1,2,1,4),Table2(2,2,1,4), Table2(1,3,1,4),Table2(2,3,1,4));
fprintf(FID, '\\multicolumn{4}{l}{B. Net Return to Capital in \\%%/Aggregate savings rate in \\%% ($\\sigma=0.4$)} \\\\ \n');
fprintf(FID, ' $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(1), Table2(1,1,2,1),Table2(2,1,2,1), Table2(1,2,2,1),Table2(2,2,2,1), Table2(1,3,2,1),Table2(2,3,2,1));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(2), Table2(1,1,2,2),Table2(2,1,2,2), Table2(1,2,2,2),Table2(2,2,2,2), Table2(1,3,2,2),Table2(2,3,2,2));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(3), Table2(1,1,2,3),Table2(2,1,2,3), Table2(1,2,2,3),Table2(2,2,2,3), Table2(1,3,2,3),Table2(2,3,2,3));
fprintf(FID, ' %1.1f & %1.4f/%1.2f & %1.4f/%1.2f & %1.4f/%1.2f \\\\ \n', rho_vec(4), Table2(1,1,2,4),Table2(2,1,2,4), Table2(1,2,2,4),Table2(2,2,2,4), Table2(1,3,2,4),Table2(2,3,2,4));
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 2 of Aiyagari (1994) using grid sizes $ n_k=%d $, $ n_z=%d $, $ n_p=%d $ \n', n_k, n_z, n_p');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

FilenameString=['./SavedOutput/LatexInputs/Aiyagari1994_TableGini.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Gini Coefficients for Earnings, Income, and Wealth \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccc} \\hline \\hline \n');
fprintf(FID, '\\multicolumn{4}{l}{A. Earnings Gini/Income Gini/Wealth Gini ($\\sigma=0.2$)} \\\\ \n');
fprintf(FID, '  $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(1), Table3(1,1,1,1),Table3(2,1,1,1),Table3(3,1,1,1), Table3(1,2,1,1),Table3(2,2,1,1),Table3(3,2,1,1), Table3(1,3,1,1),Table3(2,3,1,1),Table3(3,3,1,1));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(2), Table3(1,1,1,2),Table3(2,1,1,2),Table3(3,1,1,2), Table3(1,2,1,2),Table3(2,2,1,2),Table3(3,2,1,2), Table3(1,3,1,2),Table3(2,3,1,2),Table3(3,3,1,2));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(3), Table3(1,1,1,3),Table3(2,1,1,3),Table3(3,1,1,3), Table3(1,2,1,3),Table3(2,2,1,3),Table3(3,2,1,3), Table3(1,3,1,3),Table3(2,3,1,3),Table3(3,3,1,3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(4), Table3(1,1,1,4),Table3(2,1,1,4),Table3(3,1,1,4), Table3(1,2,1,4),Table3(2,2,1,4),Table3(3,2,1,4), Table3(1,3,1,4),Table3(2,3,1,4),Table3(3,3,1,4));
fprintf(FID, '\\multicolumn{4}{l}{B. Earnings Gini/Income Gini/Wealth Gini ($\\sigma=0.4$)} \\\\ \n');
fprintf(FID, '  $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(1), Table3(1,1,2,1),Table3(2,1,2,1),Table3(3,1,2,1), Table3(1,2,2,1),Table3(2,2,2,1),Table3(3,2,2,1), Table3(1,3,2,1),Table3(2,3,2,1),Table3(3,3,2,1));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(2), Table3(1,1,2,2),Table3(2,1,2,2),Table3(3,1,2,2), Table3(1,2,2,2),Table3(2,2,2,2),Table3(3,2,2,2), Table3(1,3,2,2),Table3(2,3,2,2),Table3(3,3,2,2));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(3), Table3(1,1,2,3),Table3(2,1,2,3),Table3(3,1,2,3), Table3(1,2,2,3),Table3(2,2,2,3),Table3(3,2,2,3), Table3(1,3,2,3),Table3(2,3,2,3),Table3(3,3,2,3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(4), Table3(1,1,2,4),Table3(2,1,2,4),Table3(3,1,2,4), Table3(1,2,2,4),Table3(2,2,2,4),Table3(3,2,2,4), Table3(1,3,2,4),Table3(2,3,2,4),Table3(3,3,2,4));
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Aiyagari (1994) reports a few Gini coefficients, but no Table. \n');
fprintf(FID, 'Uses grid sizes $ n_k=%d $, $ n_z=%d $, $ n_p=%d $ \n', n_k, n_z, n_p);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

FilenameString=['./SavedOutput/LatexInputs/Aiyagari1994_TableInvParetoCoeff.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, '(Inverse) Pareto Coefficients for Earnings, Income, and Wealth \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccc} \\hline \\hline \n');
fprintf(FID, '\\multicolumn{4}{l}{A. Earnings Pareto Coeff/Income Pareto Coeff/Wealth Pareto Coeff ($\\sigma=0.2$)} \\\\ \n');
fprintf(FID, ' $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(1), Table3(4,1,1,1),Table3(5,1,1,1),Table3(6,1,1,1), Table3(4,2,1,1),Table3(5,2,1,1),Table3(6,2,1,1), Table3(4,3,1,1),Table3(5,3,1,1),Table3(6,3,1,1));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(2), Table3(4,1,1,2),Table3(5,1,1,2),Table3(6,1,1,2), Table3(4,2,1,2),Table3(5,2,1,2),Table3(6,2,1,2), Table3(4,3,1,2),Table3(5,3,1,2),Table3(6,3,1,2));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(3), Table3(4,1,1,3),Table3(5,1,1,3),Table3(6,1,1,3), Table3(4,2,1,3),Table3(5,2,1,3),Table3(6,2,1,3), Table3(4,3,1,3),Table3(5,3,1,3),Table3(6,3,1,3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(4), Table3(4,1,1,4),Table3(5,1,1,4),Table3(6,1,1,4), Table3(4,2,1,4),Table3(5,2,1,4),Table3(6,2,1,4), Table3(4,3,1,4),Table3(5,3,1,4),Table3(6,3,1,4));
fprintf(FID, '\\multicolumn{4}{l}{B. Earnings Pareto Coeff/Income Pareto Coeff/Wealth Pareto Coeff ($\\sigma=0.4$)} \\\\ \n');
fprintf(FID, ' $\\rho$/$\\mu$ & %d & %d & %d \\\\ \\hline \n', mu_vec(1), mu_vec(2), mu_vec(3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(1), Table3(4,1,2,1),Table3(5,1,2,1),Table3(6,1,2,1), Table3(4,2,2,1),Table3(5,2,2,1),Table3(6,2,2,1), Table3(4,3,2,1),Table3(5,3,2,1),Table3(6,3,2,1));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(2), Table3(4,1,2,2),Table3(5,1,2,2),Table3(6,1,2,2), Table3(4,2,2,2),Table3(5,2,2,2),Table3(6,2,2,2), Table3(4,3,2,2),Table3(5,3,2,2),Table3(6,3,2,2));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(3), Table3(4,1,2,3),Table3(5,1,2,3),Table3(6,1,2,3), Table3(4,2,2,3),Table3(5,2,2,3),Table3(6,2,2,3), Table3(4,3,2,3),Table3(5,3,2,3),Table3(6,3,2,3));
fprintf(FID, ' %1.1f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f & %1.2f/%1.2f/%1.2f \\\\ \n', rho_vec(4), Table3(4,1,2,4),Table3(5,1,2,4),Table3(6,1,2,4), Table3(4,2,2,4),Table3(5,2,2,4),Table3(6,2,2,4), Table3(4,3,2,4),Table3(5,3,2,4),Table3(6,3,2,4));
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Aiyagari (1994) does not report Inverted Pareto coefficients. \n');
%fprintf(FID, 'Inverted Pareto coefficients, b, are calculated from top earnings/income/wealth shares as b=1/[log(S1\\%%/S0.1\\%%)/log(10)]. \\\\ \n');
fprintf(FID, 'Uses grid sizes $ n_k=%d $, $ n_z=%d $, $ n_p=%d $ \n', n_k, n_z, n_p);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

mu_ind=2;
FilenameString=['./SavedOutput/LatexInputs/Aiyagari1994_TableRminusG.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccccc} \\hline \\hline \n');
fprintf(FID, '\\multicolumn{4}{c}{$\\sigma=0.2$} \\\\ \n');
fprintf(FID, '  &  & \\multicolumn{2}{c}{Gini Coeffs}& \\multicolumn{2}{c}{Inverted Pareto Coeff.} \\\\ \n');
fprintf(FID, ' $\\rho$ & Net Return to Capital in \\%% & Income & Wealth & Income & Wealth \\\\ \n');
rr=1; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,1,rr), Table3(1,mu_ind,1,rr), Table3(3,mu_ind,1,rr), Table3(4,mu_ind,1,rr), Table3(6,mu_ind,1,rr) );
rr=2; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,1,rr), Table3(1,mu_ind,1,rr), Table3(3,mu_ind,1,rr), Table3(4,mu_ind,1,rr), Table3(6,mu_ind,1,rr) );
rr=3; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,1,rr), Table3(1,mu_ind,1,rr), Table3(3,mu_ind,1,rr), Table3(4,mu_ind,1,rr), Table3(6,mu_ind,1,rr) );
rr=4; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,1,rr), Table3(1,mu_ind,1,rr), Table3(3,mu_ind,1,rr), Table3(4,mu_ind,1,rr), Table3(6,mu_ind,1,rr) );
fprintf(FID, '\\hline \\multicolumn{4}{c}{$\\sigma=0.4$} \\\\ \n');
fprintf(FID, '  &  & \\multicolumn{2}{c}{Gini Coeffs}& \\multicolumn{2}{c}{Inverted Pareto Coeff.} \\\\ \n');
fprintf(FID, ' $\\rho$ & Net Return to Capital in \\%% & Income & Wealth & Income & Wealth \\\\ \n');
rr=1; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,2,rr), Table3(1,mu_ind,2,rr), Table3(3,mu_ind,2,rr), Table3(4,mu_ind,2,rr), Table3(6,mu_ind,2,rr) );
rr=2; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,2,rr), Table3(1,mu_ind,2,rr), Table3(3,mu_ind,2,rr), Table3(4,mu_ind,2,rr), Table3(6,mu_ind,2,rr) );
rr=3; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,2,rr), Table3(1,mu_ind,2,rr), Table3(3,mu_ind,2,rr), Table3(4,mu_ind,2,rr), Table3(6,mu_ind,2,rr) );
rr=4; fprintf(FID, ' %1.1f & %1.4f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', rho_vec(rr), Table2(1,mu_ind,2,rr), Table3(1,mu_ind,2,rr), Table3(3,mu_ind,2,rr), Table3(4,mu_ind,2,rr), Table3(6,mu_ind,2,rr) );
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Aiyagari (1994) with $\\mu$=%d using grid sizes $ n_k=%d $, $ n_z=%d $, $ n_p=%d $ \n', mu_vec(mu_ind), n_k, n_z, n_p');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);








