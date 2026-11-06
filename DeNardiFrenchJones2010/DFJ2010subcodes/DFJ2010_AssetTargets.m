% DFJ2010_AssetTargets: the median-asset GMM targets by cohort x PI quintile x year from DataWork.Targets.MomentList, as
% TargetMoments.AgeConditionalStats.cohortC.aliveqX.assets.Median, with the weighting and covariance matrices, tkorder,
% FnsToEvaluate_gmm and simoptions_gmm. Script shared by DeNardiFrenchJones2010.m and DFJ2010_EndogMedex.m; needs
% DataWork, cohortentryage, ncohorts, N_j, Params.agejshifter, simoptions in the workspace.
% [Moved out of DeNardiFrenchJones2010.m on 2026-09-29.]

%% GMM targets from DFJ2010_DataWork.m, by cohort (no pooling across cohorts)
% DFJ2010 match median assets by cohort p, permanent-income quintile q and calendar year v (1998,...,2006), pooling men
% and women: 102 moment conditions in the HRS rebuild (101 in DFJ2010). Cohort c enters the model at cohortentryage(c),
% so its year-v cell is the model's median assets at age cohortentryage(c)+(v-1996), under the restriction aliveq_q
% (alive and in quintile q, both genders), for cohort c. The toolkit's cohort machinery (estimoptions.cohortagejshifter
% with the N_j-dimensional jequaloneDist above) computes each cohort's distribution from its own entry age with the same
% decision rules, and the ptype weights within a cohort are the ptypes' 1996 masses in that cohort.

% The GMM only needs the median of assets under the aliveq restrictions
simoptions_gmm=simoptions;
simoptions_gmm.conditionalrestrictions=struct();
for q=1:5
    simoptions_gmm.conditionalrestrictions.(['aliveq',num2str(q)])=simoptions.conditionalrestrictions.(['aliveq',num2str(q)]);
end
simoptions_gmm.whichstats=[0,1,0,0,0,0,0]; % median only
FnsToEvaluate_gmm.assets=@(aprime,a,h,zeta,xi) a;

MomentList=DataWork.Targets.MomentList; % one row per DFJ moment: cohort, PIq, year, avgage, median, se, n (row = index into Targets.VCV_medians)
Jdata=height(MomentList);
% Targets in the toolkit's format: TargetMoments.AgeConditionalStats.cohortC.aliveqX.assets.Median, an N_j vector with
% NaN at untargeted ages. The toolkit stacks the moments cohort by cohort, then in the order the aliveqX fields were created,
% then by age; tkorder records that order so the weighting and covariance matrices can be put into it.
TargetMoments=struct();
tkorder=zeros(0,1);
for c=1:ncohorts
    for q=1:5
        rows=find(MomentList.cohort==c & MomentList.PIq==q);
        if ~isempty(rows)
            ages=cohortentryage(c)+(MomentList.year(rows)-1996); % model age of each cell
            v=nan(N_j,1);
            v(ages-Params.agejshifter)=MomentList.median(rows);
            TargetMoments.AgeConditionalStats.(['cohort',num2str(c)]).(['aliveq',num2str(q)]).assets.Median=v;
            [~,ord]=sort(ages);
            tkorder=[tkorder; rows(ord)]; %#ok<AGROW>
        end
    end
end
if length(tkorder)~=Jdata || length(unique(tkorder))~=Jdata
    error('the DFJ moments do not map one-to-one onto (cohort, quintile, age) targets')
end
WeightingMatrix=diag(1./MomentList.se(tkorder).^2); % diagonal: inverse bootstrap variances of the cell medians
CoVarMatrixDataMoments=DataWork.Targets.VCV_medians(tkorder,tkorder);
fprintf('GMM targets: %d cohort-year-quintile medians over %d cohorts \n',Jdata,ncohorts)

