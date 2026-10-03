% HV2000 Table 6: saving rates in Model 4 at three temp-shock variances.
% Columns: US (HV2000 Table 5 baseline) | "Including all agents" sigma2_eps2 in {0.01, 0.04, 0.09}
% | "Excluding agents w/ temp shocks" (restricted to middle e) for the same three variances.
% Model 4 only, a_underbar=0 only.
%
% Runtime: 3 extra GE solves on top of Tables 2-5.

fprintf('\n========== HV2000 Table 6 (Model 4 temp-shock-variance sensitivity) ==========\n');

Params.borrowFlag=0;
HV2000_Table6_variances=[0.01, 0.04, 0.09];
HV2000_multiples=[0.25,0.5,0.75,1,1.5,2,3,4,7,10];
ratesAll=zeros(length(HV2000_multiples),3);
ratesExcl=zeros(length(HV2000_multiples),3);

% Each of the 3 GE solves is checkpointed to a separate .mat under
% ./SavedOutput/Intermediates/ so an interrupted run can resume. At the current
% grid a Model 4 solve is ~8h, so all three do not fit one 24h job. To force a
% full re-solve, set doPart2_force=true OR delete the mats.
doPart2_force=false;
if ~exist('./SavedOutput/Intermediates','dir'); mkdir('./SavedOutput/Intermediates'); end

for kk=1:length(HV2000_Table6_variances)
    unit_file=sprintf('./SavedOutput/Intermediates/HV2000_Table6_var%d.mat',kk);
    if ~doPart2_force && exist(unit_file,'file')
        S=load(unit_file); ratesAll(:,kk)=S.rAll; ratesExcl(:,kk)=S.rExcl;
        fprintf('Table6 variance %d (sigma2_eps2=%g): LOAD %s\n',kk,HV2000_Table6_variances(kk),unit_file);
        continue
    end
    Params.sigma2_eps2=HV2000_Table6_variances(kk);
    HV2000_Model4;  % solves Model 4, populates Results.Model4_abar0.savingRates (standard income bins)
    ratesAll(:,kk)=Results.Model4_abar0.savingRates;

    % "Excluding agents w/ temp shocks": restrict each income bin to the middle e value.
    % Income-bin spec is the same as HV2000 Table 5 (p. 380):
    %   "multiples are calculated by taking a 10% band around each income multiple
    %    and then dividing total saving of agents in the band by total income of
    %    agents in the band. Income is defined as earnings after social security
    %    taxes plus interest income and transfers."
    % Here we add an extra (e == emid) clause to drop agents with nonzero transitory shock.
    % StationaryDist, Policy, FnsToEvaluate_M4, e_grid, z_grid_J etc. are still in workspace from the just-finished HV2000_Model4 run.
    Params.IncomeMean=Results.Model4_abar0.AllStats.Income.Mean;
    emid=e_grid(2); % middle e value (n_e=3 so index 2 is middle)
    % simoptions4 (Model 4 variant with experienceassetze/n_e/e_grid) is in
    % scope from the just-finished HV2000_Model4 run.
    % The bin saving rate is total saving in the band over total income in the band,
    % which is a ratio of two AGGREGATES: the bin mass cancels, so no conditional
    % restriction is needed. Folding the bin indicator (here including the extra
    % e==emid clause that drops agents with a nonzero transitory shock) into the
    % functions and taking one AggVars pass avoids AllStats entirely -- AllStats
    % sorts the whole grid for medians/Lorenz/quantiles we never read, and that
    % sort is what exhausts memory once the asset grid grows. An empty bin gives
    % 0/0=NaN, which is what the old RestrictedSampleMass>0 guard produced.
    % GPU arrayfun cannot compile closures that capture outer-workspace vars, so
    % bake the bin edges and emid in as numeric literals via str2func+sprintf.
    binsig='@(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS) ';
    binincfn='HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)';
    FnsToEvaluate_bins=struct();
    for ii=1:length(HV2000_multiples)
        loVal=0.9*HV2000_multiples(ii)*Params.IncomeMean;
        hiVal=1.1*HV2000_multiples(ii)*Params.IncomeMean;
        inbin=sprintf('((%s>=%.16g)*(%s<=%.16g)*(e==%.16g))',binincfn,loVal,binincfn,hiVal,emid);
        FnsToEvaluate_bins.(sprintf('Sav%d',ii))=str2func([binsig '((1+g)*aprime-a)*' inbin]);
        FnsToEvaluate_bins.(sprintf('Inc%d',ii))=str2func([binsig binincfn '*' inbin]);
    end
    AggBinsExcl=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_bins,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions4);
    for ii=1:length(HV2000_multiples)
        ratesExcl(ii,kk)=gather(AggBinsExcl.(sprintf('Sav%d',ii)).Mean/AggBinsExcl.(sprintf('Inc%d',ii)).Mean);
    end

    % V, Policy and StationaryDist are no longer needed, and memory available is limited, so clear them
    clear V Policy StationaryDist

    % Checkpoint this variance. amax_stamp is recorded so a leftover file from a
    % different asset grid is identifiable rather than silently loaded.
    rAll=ratesAll(:,kk); rExcl=ratesExcl(:,kk);
    sigma2_eps2_stamp=HV2000_Table6_variances(kk); amax_stamp=Params.amax;
    save(unit_file,'rAll','rExcl','sigma2_eps2_stamp','amax_stamp');
end

% Restore baseline sigma2_eps2 (so any later inspection of Params reads the baseline value)
Params.sigma2_eps2=0.01;

%% Write Table 6
FID=fopen('./SavedOutput/LatexInputs/HV2000_Table6.tex','w');
fprintf(FID,'\\begin{tabular}{lccccccc} \n');
fprintf(FID,'\\hline\n');
fprintf(FID,'& & \\multicolumn{3}{c}{Including all agents} & \\multicolumn{3}{c}{Excluding agents w/ temp shocks} \\\\ \\cline{3-5} \\cline{6-8} \n');
fprintf(FID,'& & \\multicolumn{3}{c}{Variance ($\\sigma_{\\epsilon,2}^2$)} & \\multicolumn{3}{c}{Variance ($\\sigma_{\\epsilon,2}^2$)} \\\\ \\cline{3-5} \\cline{6-8} \n');
fprintf(FID,'Income multiple & US & 0.01 & 0.04 & 0.09 & 0.01 & 0.04 & 0.09 \\\\ \n');
fprintf(FID,'\\hline\n');
us=[-19.3, -1.3, 4.8, 7.9, 13.0, 16.5, 22.4, 27.1, 37.3, 39.2];
multlabels={'0.25','0.50','0.75','1.0','1.5','2.0','3.0','4.0','7.0','10.0'};
for ii=1:length(HV2000_multiples)
    fprintf(FID,'%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\ \n', ...
        multlabels{ii},us(ii), ...
        100*ratesAll(ii,1),100*ratesAll(ii,2),100*ratesAll(ii,3), ...
        100*ratesExcl(ii,1),100*ratesExcl(ii,2),100*ratesExcl(ii,3));
end
fprintf(FID,'\\hline\n');
fprintf(FID,'\\end{tabular}\n');
% Note sits OUTSIDE the tabular: a p{} multicolumn wider than the table's own
% natural width pushes all the slack into the final column, visibly detaching it.
fprintf(FID,'\\\\[4pt]\n');
fprintf(FID,'\\parbox{15cm}{\\footnotesize Note: the asset grid ceiling is $a_{\\max}=%g$. Whether it binds, and for how much mass, is reported by the CheckTopOfGrid diagnostic at the end of the main script and discussed in the computational appendix; where it binds, the largest income multiples are the entries affected.}\n',Params.amax);
fclose(FID);

