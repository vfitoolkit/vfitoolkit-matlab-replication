% HV2000 Tables 8 and 9 (Appendix B) - "equal saving rates within age groups" counterfactual.
% Table 8: baseline (SS on). Models 2-4 x {a_underbar in 0, -w}.
% Table 9: no social security. Same models.
%
% Mechanism: each agent's individual saving RATE is replaced with their age group's own saving
% rate sbar(j) (age group j's total saving over its total income), so agent i's counterfactual
% saving is sbar(j)*y_i. The bin saving rate = sum(sbar(j)*income | in bin) / sum(income | in bin),
% i.e. an income-weighted average of the age-group saving rates of the agents in the bin.
% This isolates how much of the saving-rate-by-income pattern is driven by cross-age vs within-age variation.
% NB it is the RATE that is equalised within an age group, not the saving level: equalising the
% level makes the bin ratio decay like 1/(income multiple) and inverts the profile. See the comment
% at the Params.savingRateAgeMean assignment in Models 2-4.
%
% Runtime: 12 extra GE solves (Models 2-4 x 2 a_underbar x 2 SS scenarios). The within-age
% calculation reuses the workspace that each Model script leaves behind (StationaryDist + Policy)
% via the Params.computeWithinAgeMean flag added to Models 2-4.

fprintf('\n========== HV2000 Tables 8 and 9 (within-age-group averaging) ==========\n');

% Each of the 12 GE solves (Models 2-4 x 2 a_underbar x 2 SS scenarios) is
% checkpointed to a separate .mat under ./SavedOutput/Intermediates/ so an
% interrupted run can resume. To force a full re-solve, set doPart4_force=true
% OR delete the mats.
doPart4_force=false;
if ~exist('./SavedOutput/Intermediates','dir'); mkdir('./SavedOutput/Intermediates'); end

% Save baseline state we will modify
Params_haveSS_orig=Params.haveSS;
Params_theta_orig=Params.theta;
Params_bcommon_orig=Params.b_common;
GEPriceParamNames_orig=GEPriceParamNames;
GeneralEqmEqns_orig=GeneralEqmEqns;

% Activate within-age-group block in Models 2-4 (Fig 3/4/7 are already guarded off by doPart1figures==0 outside doPart(1); Fig 5 will be re-rendered by Model 2 haveSS=0 pass but is identical to the Table 7 render).
Params.computeWithinAgeMean=1;

HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];
rates_T8=struct(); % baseline (with SS)
rates_T9=struct(); % no SS

for pass=1:2
    if pass==1
        Params.haveSS=1;
        Params.theta=Params_theta_orig;
        Params.b_common=Params_bcommon_orig;
        GEPriceParamNames=GEPriceParamNames_orig;
        GeneralEqmEqns=GeneralEqmEqns_orig;
        fprintf('\n--- Pass 1/2: baseline (SS on) ---\n');
    elseif pass==2
        Params.haveSS=0;
        Params.theta=0;
        Params.b_common=0;
        GEPriceParamNames=setdiff(GEPriceParamNames_orig,{'theta','b_common'},'stable');
        GeneralEqmEqns=rmfield(GeneralEqmEqns_orig,{'SSBalance','bCommon'});
        fprintf('\n--- Pass 2/2: no social security ---\n');
    end

    for borrowFlag=0:1
        Params.borrowFlag=borrowFlag;
        if borrowFlag==0
            abarstr='abar0';
        else
            abarstr='abarw';
        end

        % --- Model 2: load if cached, else solve+save ---
        unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Tables89_pass%d_%s_Model2.mat',pass,abarstr);
        if ~doPart4_force && exist(unit_file,'file')
            S=load(unit_file); srWA2=S.srWA;
            fprintf('Tables89 pass%d %s Model2: LOAD %s\n',pass,abarstr,unit_file);
        else
            HV2000_Model2;
            srWA=Results.(['Model2_',abarstr]).savingRatesWithinAge; srWA2=srWA;
            pass_stamp=pass; borrowFlag_stamp=borrowFlag; haveSS_stamp=Params.haveSS;
            save(unit_file,'srWA','pass_stamp','borrowFlag_stamp','haveSS_stamp');
        end

        % --- Model 3: load if cached, else solve+save ---
        unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Tables89_pass%d_%s_Model3.mat',pass,abarstr);
        if ~doPart4_force && exist(unit_file,'file')
            S=load(unit_file); srWA3=S.srWA;
            fprintf('Tables89 pass%d %s Model3: LOAD %s\n',pass,abarstr,unit_file);
        else
            HV2000_Model3;
            srWA=Results.(['Model3_',abarstr]).savingRatesWithinAge; srWA3=srWA;
            pass_stamp=pass; borrowFlag_stamp=borrowFlag; haveSS_stamp=Params.haveSS;
            save(unit_file,'srWA','pass_stamp','borrowFlag_stamp','haveSS_stamp');
        end
        clear V StationaryDist Policy

        % --- Model 4: load if cached, else solve+save ---
        unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Tables89_pass%d_%s_Model4.mat',pass,abarstr);
        if ~doPart4_force && exist(unit_file,'file')
            S=load(unit_file); srWA4=S.srWA;
            fprintf('Tables89 pass%d %s Model4: LOAD %s\n',pass,abarstr,unit_file);
        else
            HV2000_Model4;
            srWA=Results.(['Model4_',abarstr]).savingRatesWithinAge; srWA4=srWA;
            pass_stamp=pass; borrowFlag_stamp=borrowFlag; haveSS_stamp=Params.haveSS;
            save(unit_file,'srWA','pass_stamp','borrowFlag_stamp','haveSS_stamp');
        end
        clear V StationaryDist Policy

        if pass==1
            rates_T8.(['Model2_',abarstr])=srWA2;
            rates_T8.(['Model3_',abarstr])=srWA3;
            rates_T8.(['Model4_',abarstr])=srWA4;
        elseif pass==2
            rates_T9.(['Model2_',abarstr])=srWA2;
            rates_T9.(['Model3_',abarstr])=srWA3;
            rates_T9.(['Model4_',abarstr])=srWA4;
        end
    end
end

% Restore baseline state
Params.haveSS=Params_haveSS_orig;
Params.theta=Params_theta_orig;
Params.b_common=Params_bcommon_orig;
Params.computeWithinAgeMean=0;
GEPriceParamNames=GEPriceParamNames_orig;
GeneralEqmEqns=GeneralEqmEqns_orig;

%% Write Table 8
us=[-19.3, -1.3, 4.8, 7.9, 13.0, 16.5, 22.4, 27.1, 37.3, 39.2];
multlabels={'0.25','0.50','0.75','1.0','1.5','2.0','3.0','4.0','7.0','10.0'};

FID=fopen('./SavedOutput/LatexInputs/HV2000_Table8.tex','w');
fprintf(FID,'\\begin{tabular}{lccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,'& & \\multicolumn{2}{c}{Model 2} & \\multicolumn{2}{c}{Model 3} & \\multicolumn{2}{c}{Model 4} \\\\ \\cline{3-4} \\cline{5-6} \\cline{7-8} \n');
fprintf(FID,'Income multiple & US & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ \\\\ \n');
fprintf(FID,'\\hline\n');
m2_0=100*rates_T8.Model2_abar0; m2_w=100*rates_T8.Model2_abarw;
m3_0=100*rates_T8.Model3_abar0; m3_w=100*rates_T8.Model3_abarw;
m4_0=100*rates_T8.Model4_abar0; m4_w=100*rates_T8.Model4_abarw;
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

%% Write Table 9
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table9.tex','w');
fprintf(FID,'\\begin{tabular}{lccccccc}\n');
fprintf(FID,'\\hline\n');
fprintf(FID,'& & \\multicolumn{2}{c}{Model 2} & \\multicolumn{2}{c}{Model 3} & \\multicolumn{2}{c}{Model 4} \\\\ \\cline{3-4} \\cline{5-6} \\cline{7-8} \n');
fprintf(FID,'Income multiple & US & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ & $\\underline{a}=0$ & $\\underline{a}=-w$ \\\\ \n');
fprintf(FID,'\\hline\n');
m2_0=100*rates_T9.Model2_abar0; m2_w=100*rates_T9.Model2_abarw;
m3_0=100*rates_T9.Model3_abar0; m3_w=100*rates_T9.Model3_abarw;
m4_0=100*rates_T9.Model4_abar0; m4_w=100*rates_T9.Model4_abarw;
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
