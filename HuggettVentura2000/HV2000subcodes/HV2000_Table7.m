% HV2000 Table 7: saving rates with no social security.
% Re-runs Models 2/3/4 in GE with haveSS=0, theta=0, b_common=0.
% Both borrowing limits (a_underbar in {0, -w}). Table 7 omits Model 1 (HV2000 p. 388).
%
% Runtime: 6 extra GE solves on top of Tables 2-5.
%
% NOTE: this overwrites Results.Model{2,3,4}_abar{0,w} (and their savingRates).
% Tables 2-5 are already on disk by the time this script runs, so the loss is harmless.
% Baseline Params (theta, b_common, haveSS, GE prices/eqns) are restored at the end.

fprintf('\n========== HV2000 Table 7 (no social security) ==========\n');

% Save baseline GE setup and Params we will modify
GEPriceParamNames_orig=GEPriceParamNames;
GeneralEqmEqns_orig=GeneralEqmEqns;
theta_orig=Params.theta;
b_common_orig=Params.b_common;
haveSS_orig=Params.haveSS;

% Switch to no-SS configuration
Params.haveSS=0;
Params.theta=0;
Params.b_common=0;
GEPriceParamNames=setdiff(GEPriceParamNames,{'theta','b_common'},'stable');
GeneralEqmEqns=rmfield(GeneralEqmEqns,{'SSBalance','bCommon'});

HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];
rates_noSS=struct();

% Each of the 6 GE solves (Models 2-4 x 2 a_underbar) is checkpointed to a separate
% .mat under ./SavedOutput/Intermediates/ so an interrupted run can resume. At the
% current grid these do not reliably fit one 24h job. To force a full re-solve, set
% doPart3_force=true OR delete the mats.
doPart3_force=false;
if ~exist('./SavedOutput/Intermediates','dir'); mkdir('./SavedOutput/Intermediates'); end

for borrowFlag=0:1
    Params.borrowFlag=borrowFlag;
    if borrowFlag==0
        abarstr='abar0';
    else
        abarstr='abarw';
    end

    % --- Model 2: load if cached, else solve+save ---
    unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Table7_%s_Model2.mat',abarstr);
    if ~doPart3_force && exist(unit_file,'file')
        S=load(unit_file); sr2=S.sr;
        fprintf('Table7 %s Model2: LOAD %s\n',abarstr,unit_file);
    else
        HV2000_Model2;
        sr=Results.(['Model2_',abarstr]).savingRates; sr2=sr;
        borrowFlag_stamp=borrowFlag; amax_stamp=Params.amax;
        save(unit_file,'sr','borrowFlag_stamp','amax_stamp');
    end

    % --- Model 3: load if cached, else solve+save ---
    unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Table7_%s_Model3.mat',abarstr);
    if ~doPart3_force && exist(unit_file,'file')
        S=load(unit_file); sr3=S.sr;
        fprintf('Table7 %s Model3: LOAD %s\n',abarstr,unit_file);
    else
        HV2000_Model3;
        sr=Results.(['Model3_',abarstr]).savingRates; sr3=sr;
        borrowFlag_stamp=borrowFlag; amax_stamp=Params.amax;
        save(unit_file,'sr','borrowFlag_stamp','amax_stamp');
    end
    clear V StationaryDist Policy

    % --- Model 4: load if cached, else solve+save ---
    unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Table7_%s_Model4.mat',abarstr);
    if ~doPart3_force && exist(unit_file,'file')
        S=load(unit_file); sr4=S.sr;
        fprintf('Table7 %s Model4: LOAD %s\n',abarstr,unit_file);
    else
        HV2000_Model4;
        sr=Results.(['Model4_',abarstr]).savingRates; sr4=sr;
        borrowFlag_stamp=borrowFlag; amax_stamp=Params.amax;
        save(unit_file,'sr','borrowFlag_stamp','amax_stamp');
    end
    clear V StationaryDist Policy

    rates_noSS.(['Model2_',abarstr])=sr2;
    rates_noSS.(['Model3_',abarstr])=sr3;
    rates_noSS.(['Model4_',abarstr])=sr4;
end

% Restore baseline GE setup and Params
GEPriceParamNames=GEPriceParamNames_orig;
GeneralEqmEqns=GeneralEqmEqns_orig;
Params.theta=theta_orig;
Params.b_common=b_common_orig;
Params.haveSS=haveSS_orig;

%% Write Table 7
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table7.tex','w');
fprintf(FID,'\\begin{tabular}{lccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,'& & \\multicolumn{2}{c}{Model 2} & \\multicolumn{2}{c}{Model 3} & \\multicolumn{2}{c}{Model 4} \\\\ \\cline{3-4} \\cline{5-6} \\cline{7-8} \n');
fprintf(FID,'Income multiple & US & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ \\\\ \n');
fprintf(FID,'\\hline\n');
us=[-19.3, -1.3, 4.8, 7.9, 13.0, 16.5, 22.4, 27.1, 37.3, 39.2];
multlabels={'0.25','0.50','0.75','1.0','1.5','2.0','3.0','4.0','7.0','10.0'};
m2_0=100*rates_noSS.Model2_abar0;
m2_w=100*rates_noSS.Model2_abarw;
m3_0=100*rates_noSS.Model3_abar0;
m3_w=100*rates_noSS.Model3_abarw;
m4_0=100*rates_noSS.Model4_abar0;
m4_w=100*rates_noSS.Model4_abarw;
for ii=1:length(HV2000_multiples)
    fprintf(FID,'%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\ \n', ...
        multlabels{ii},us(ii),m2_0(ii),m2_w(ii),m3_0(ii),m3_w(ii),m4_0(ii),m4_w(ii));
end
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
% Note sits OUTSIDE the tabular: a p{} multicolumn wider than the table's own
% natural width pushes all the slack into the final column, visibly detaching it.
fprintf(FID,'\\\\[4pt]\n');
fprintf(FID,'\\parbox{15cm}{\\footnotesize Note: the asset grid ceiling is $a_{\\max}=%g$. Whether it binds, and for how much mass, is reported by the CheckTopOfGrid diagnostic at the end of the main script and discussed in the computational appendix; where it binds, the largest income multiples are the entries affected.}\n',Params.amax);
fclose(FID);

