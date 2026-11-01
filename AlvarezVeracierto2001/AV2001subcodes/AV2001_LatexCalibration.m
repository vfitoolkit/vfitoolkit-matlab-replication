function AV2001_LatexCalibration(Out_LFendog,Params,filename)
% Writes the laissez-faire calibration check (endogenous search economy against the targets the paper states in
% its Section 4 text) as a LaTeX tabular. The paper has no table of these, so its numbers are shown alongside.

g=@(x) double(gather(x));
O=Out_LFendog;
fid=fopen(filename,'w');
fprintf(fid,'\\begin{tabular}{lcc} \\hline \n');
fprintf(fid,' & Model & Paper \\\\ \\hline \n');
fprintf(fid,'Job-finding rate (per half-quarter) & %.1f\\%% & 61.0\\%% \\\\ \n',100*g(O.jobfindingrate));
fprintf(fid,'Unemployment rate & %.1f\\%% & 4.9\\%% \\\\ \n',100*g(O.unemploymentrate));
fprintf(fid,'Layoff rate (per half-quarter) & %.1f\\%% & 2.9\\%% \\\\ \n',100*g(O.layoffrate));
fprintf(fid,'Annual real interest rate & %.1f\\%% & 4.0\\%% \\\\ \n',100*g(O.annualinterestrate));
fprintf(fid,'Capital to (measured) output, annual & %.2f & 1.7 \\\\ \n',g(O.KYannual));
fprintf(fid,'Average establishment size & %.1f & 62 \\\\ \n',g(O.avgestablishmentsize));
fprintf(fid,'Search intensity $\\alpha$ (set, not targeted) & %.1f & 15.5 \\\\ \n',Params.alpha);
fprintf(fid,'\\hline \n\\end{tabular} \n');
fclose(fid);

end
