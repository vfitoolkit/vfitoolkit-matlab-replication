function AV2001_LatexCalibrationUS(Out_US,uiphi_US,epsilon_US,Params,filename)
% Writes the calibration of the U.S. unemployment insurance system (doPart(7)) as a LaTeX tabular: the two
% parameters calibrated inside the model and the targets from AV2001's Section 4 text, which the paper reports
% only in words, so its numbers are shown alongside.

g=@(x) double(gather(x));
O=Out_US;
fid=fopen(filename,'w');
fprintf(fid,'\\begin{tabular}{lcc} \\hline \n');
fprintf(fid,' & Model & Paper \\\\ \\hline \n');
fprintf(fid,'\\multicolumn{3}{l}{\\emph{Targets}} \\\\ \n');
fprintf(fid,'Weeks collecting benefits while unemployed & %.1f & about 9 \\\\ \n',g(O.UIweeks));
fprintf(fid,'Share of benefits paid by firing taxes & %.1f\\%% & 60\\%% \\\\ \n',100*g(O.firingtaxshareofUI));
fprintf(fid,'\\multicolumn{3}{l}{\\emph{Set directly}} \\\\ \n');
fprintf(fid,'Replacement ratio & %.2f & 0.66 \\\\ \n',Params.uirr_US);
fprintf(fid,'Probability a laid-off worker is eligible & %.2f & 0.50 \\\\ \n',Params.uikappa_US);
fprintf(fid,'\\multicolumn{3}{l}{\\emph{Calibrated}} \\\\ \n');
fprintf(fid,'Probability of keeping eligibility, $\\phi$ & %.3f & \\\\ \n',uiphi_US);
fprintf(fid,'Firing tax, $\\epsilon$ (wages) & %.3f & \\\\ \n',epsilon_US);
fprintf(fid,'\\hline \n\\end{tabular} \n');
fclose(fid);

end
