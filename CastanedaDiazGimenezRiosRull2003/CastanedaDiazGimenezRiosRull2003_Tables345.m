% Recreates Tables 3, 4, and 5 of CastanedaDiazGimenezRiosRull2003.
% These Tables simply report the calibrated model parameters.


FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table3.tex', 'w');
fprintf(FID, 'Parameter Values for the Benchmark Model Economy \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llr} \n \\hline \\hline \n');
fprintf(FID, '\\multicolumn{3}{l}{\\sl Preferences} \\\\ \n');
fprintf(FID, '\\quad Time discount factor  \\hspace*{6cm}       & $\\beta$ & %8.3f\\\\ \n', Params.beta);
fprintf(FID, '\\quad Curvature of consumption                  & $\\sigma_1$ & %8.3f\\\\ \n', Params.sigma1);
fprintf(FID, '\\quad Curvature of leisure                      & $\\sigma_2$ & %8.3f\\\\ \n', Params.sigma2);
fprintf(FID, '\\quad Relative share of consumption and leisure & $\\chi$ & %8.3f\\\\ \n', Params.chi);
fprintf(FID, '\\quad Endowment of productive time           & $\\ell$ & %8.3f\\\\ \n', Params.elle);
fprintf(FID, '\\multicolumn{3}{l}{\\sl Age and endowment process} \\\\ \n');
fprintf(FID, '\\quad Common probability of retiring                   & $p_{e\\varrho}$ & %8.3f\\\\ \n', Params.p_eg);
fprintf(FID, '\\quad Common probability of dying                      & $1-p_{\\varrho\\varrho}$ & %8.3f\\\\ \n', 1-Params.p_gg);
fprintf(FID, '\\quad Life cycle earnings profile               & $\\phi_1$ & %8.3f\\\\ \n', Params.phi1);
fprintf(FID, '\\quad Intergenerational persistence of earnings & $\\phi_2$ & %8.3f\\\\ \n', Params.phi2);
fprintf(FID, '\\multicolumn{3}{l}{\\sl Technology} \\\\ \n');
fprintf(FID, '\\quad Capital share of income                      & $\\theta$ & %8.3f\\\\ \n', Params.theta);
fprintf(FID, '\\quad Capital depreciation rate                 & $\\delta$ & %8.3f\\\\ \n', Params.delta);
fprintf(FID, '\\multicolumn{3}{l}{\\sl Fiscal policy} \\\\ \n');
fprintf(FID, '\\quad Government consumption                    & $G$      & %8.3f\\\\ \n', Params.G);
fprintf(FID, '\\quad Normalized Retirement pensions            & $\\omega$ & %8.3f\\\\ \n', Params.omega);
fprintf(FID, '\\quad Income tax function parameters            & $a_0$    & %8.3f\\\\ \n', Params.a0);
fprintf(FID, '\\quad                                           & $a_1$    & %8.3f\\\\ \n', Params.a1);
fprintf(FID, '\\quad                                           & $a_2$    & %8.3f\\\\ \n', Params.a2);
fprintf(FID, '\\quad                                           & $a_3$    & %8.3f\\\\ \n', Params.a3);
fprintf(FID, '\\quad Estate tax function parameters:             \\\\ \n');
fprintf(FID, '\\quad \\quad Tax-exempt level                   & $\\underbar{z}$    & %8.3f\\\\ \n', Params.zlowerbar);
fprintf(FID, '\\quad \\quad Marginal tax rate                  & $\\tau_{E}$        & %8.3f\\\\ \n', Params.tauE);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Some minor changes to the precise description of the parameters are made from the original.');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table for the transition matix
Gamma_ee=Gamma(1:Params.J,1:Params.J)./(1-Params.p_eg);
FID = fopen('./SavedOutput/LatexInputs/CastanedaDiazGimenezRiosRull2003_Table4and5_TransMatrix.tex', 'w');
fprintf(FID, 'Transition Probabilities of the Process on the Endowment of Efficiency Labor  \\\\ \n');
fprintf(FID, 'Units for Working-Age Households That Remain at Working Age One Period \\\\ \n');
fprintf(FID, 'Later $\\Gamma_{\\mathcal{E}\\mathcal{E}}$ (\\%%) \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}crrrrrr} \n \\hline \\hline \n');
fprintf(FID, ' & & & \\multicolumn{4}{c}{{$\\Gamma_{\\cal{E} \\cal{E}}$ {\\footnotesize(\\%%)} From $s$ To $s''$}} \\\\ \n');
fprintf(FID, ' \\cline{4-7} & $e(s)$ & $\\gamma^*_s$ {\\footnotesize(\\%%)} &  $s''=1$ & $s''=2$ & $s''=3$ & $s''=4$ \\\\ \n \\hline \n');
fprintf(FID, ' $s = 1$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', Params.e1, 100*gammastarfull(1)/sum(gammastarfull(1:Params.J)),100*Gamma_ee(1,1),100*Gamma_ee(1,2),100*Gamma_ee(1,3),100*Gamma_ee(1,4));
fprintf(FID, ' $s = 2$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', Params.e2, 100*gammastarfull(2)/sum(gammastarfull(1:Params.J)),100*Gamma_ee(2,1),100*Gamma_ee(2,2),100*Gamma_ee(2,3),100*Gamma_ee(2,4));
fprintf(FID, ' $s = 3$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', Params.e3, 100*gammastarfull(3)/sum(gammastarfull(1:Params.J)),100*Gamma_ee(3,1),100*Gamma_ee(3,2),100*Gamma_ee(3,3),100*Gamma_ee(3,4));
fprintf(FID, ' $s = 4$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', Params.e4, 100*gammastarfull(4)/sum(gammastarfull(1:Params.J)),100*Gamma_ee(4,1),100*Gamma_ee(4,2),100*Gamma_ee(4,3),100*Gamma_ee(4,4));
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fclose(FID);
