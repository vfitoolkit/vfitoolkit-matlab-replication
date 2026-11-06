% DFJ2010_EndogMedex3: GMM estimation of the endogenous medical-expense model (stage 4 of EndogenousMedex_proposal.md), with
% the autocorrelation targets in the definition ACdefinition set by the main script: 'toolkit' (doPart(7): AutoCorr targets, data
% rows ac1/ac2 of DFJ2010_MedexTargets, the pair population's own means and sds) or 'DFJ' (doPart(8): CustomModelStats targets
% through DFJ2010_EM_CustomStats, data rows ac1dfj/ac2dfj, the full cell's means and sds). Saves ..._AC<definition>.mat, and
% ..._both.mat with the two side by side when both parts ran.
% the four preference parameters nu, omega, beta and the utility-floor level eMin (DFJ2010 Table A1 column 5), starting
% from DFJ2010's published estimates, on the targets of DFJ2010_EMTargets (asset medians by cohort x quintile x year,
% mean out-of-pocket medex and the first and second autocorrelations of log medex by cohort x PI half x year). The
% medical-needs shifter coefficients (mxvec) are held at DFJ2010's estimates: making them estimable needs a
% ParametrizeParamsFn from mxvec to the per-ptype age profiles, not done yet. Called from DeNardiFrenchJones2010.m as
% doPart(7); needs the endogenous model (DFJ2010_EndogMedex1) and its targets (DFJ2010_EndogMedex2) in the workspace.
%
% Cost: one objective evaluation = the value function (109 s for the 10 ptypes on the RTX 4000) + five cohorts' distributions,
% restricted life-cycle profiles and autocorrelations (the autocorrelation command pushes the measures forward with the
% Tan two-step since 2026-10-04, so it needs no transition matrix: memory is not a constraint). With 4 parameters
% lsqnonlin needs about 5 evaluations per iteration, and the standard errors a further 33.

if ~exist('ACdefinition','var')
    ACdefinition='toolkit'; % set by DeNardiFrenchJones2010.m: 'toolkit' (part 7) or 'DFJ' (part 8)
end
DFJ2010_EMTargets % the targets, weighting and covariance matrices in this definition (TargetMoments_EM, WeightingMatrix_EM, CoVarMatrixDataMoments_EM, FnsToEvaluate_EM, simoptions_EM, okEM, ...)
estimoptions_EM=struct();
estimoptions_EM.verbose=1;
estimoptions_EM.cohortagejshifter=1; % jequaloneDist has an N_j dimension: the entry ages are the slices with mass (the five cohorts)
estimoptions_EM.eedefault=4; % finite-difference step for the standard errors: the largest, as for the exogenous model (the asset medians sit on a grid)
EstimParamNames_EM={'nu','omega','beta','eMin'};
estimoptions_EM.constrainpositive={'nu','omega','eMin'};
estimoptions_EM.constrain0to1={'beta'};
if strcmp(ACdefinition,'DFJ')
    estimoptions_EM.CustomModelStats=@DFJ2010_EM_CustomStats; % the DFJ-definition autocorrelations, by cohort (TargetMoments_EM.CustomModelStats.cohortC.AC1_h1, ...)
end
vfoptions_EM=vfoptions;
vfoptions_EM.verbose=0; % the value fn is solved at every objective evaluation; keep the diary readable
vfoptions_EM.verboseparams=0;
fprintf('\n===== GMM estimation of the endogenous medical-expense model, autocorrelation targets in the %s definition: estimating %s on %d moments =====\n',ACdefinition,strjoin(EstimParamNames_EM,', '),sum(okEM))
% Starting point: DFJ2010's published estimates (Table A1 column 5), which DFJ2010_EndogMedex1 puts in Params. [The part-8 estimates
% of 2026-10-31 (SavedOutput/DFJ2010_EndogMedex_Estimation_ACDFJ.mat) came from a warm start at the point where lsqnonlin, started
% from these values, had stopped in the run of 2026-10-29, which the 24 h time limit then cut while it computed the standard errors.]
Params_EMstart=Params;
fprintf('   starting from DFJ2010 Table A1 column 5: nu %.4f, omega %.4f, beta %.4f, eMin %.1f\n',Params_EMstart.nu,Params_EMstart.omega,Params_EMstart.beta,Params_EMstart.eMin)
tic;
[EstimParams_EM,EstimParamsConfInts_EM,estsummary_EM]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames_EM,TargetMoments_EM,WeightingMatrix_EM,CoVarMatrixDataMoments_EM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_EMstart,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_EM,estimoptions_EM,vfoptions_EM,simoptions_EM);
estimtime_EM=toc;
fprintf('GMM estimation of the endogenous model (%s definition) took %.1f minutes\n',ACdefinition,estimtime_EM/60)
fprintf('Estimated parameters (DFJ2010 Table A1 column 5 in brackets: nu 2.15 (.055), omega 3.19 (.054), beta .99 (.01), eMin 202 (23.4)), standard errors in brackets:\n')
for pp=1:length(EstimParamNames_EM)
    fprintf('   %-6s %12.4f  (%.4f)\n',EstimParamNames_EM{pp},EstimParams_EM.(EstimParamNames_EM{pp}),estsummary_EM.EstimParamsStdDev.(EstimParamNames_EM{pp}))
end
fprintf('Objective at the estimates: %.1f  [%d moments, %d parameters; DFJ2010: 667.0 on 242 df]\n',estsummary_EM.objectivefnval,sum(okEM),length(EstimParamNames_EM))
contrib_est=((estsummary_EM.currentmomentvec(:)-targetvec_EM(okEM))./sevec_EM(okEM)).^2; % by moment type at the estimates
typelist=unique(EMmomenttype(okEM)); EMtype_ok=EMmomenttype(okEM);
for tt=1:length(typelist)
    ind=strcmp(EMtype_ok,typelist{tt});
    fprintf('   %-40s %4d moments, objective %7.1f, mean |model-data|/se %.2f\n',typelist{tt},sum(ind),sum(contrib_est(ind)),mean(sqrt(contrib_est(ind))))
end
save(['./SavedOutput/DFJ2010_EndogMedex_Estimation_AC',ACdefinition,'.mat'],'EstimParams_EM','EstimParamsConfInts_EM','estsummary_EM','EstimParamNames_EM','estimtime_EM','targetvec_EM','sevec_EM','EMmomenttype','okEM','ACdefinition')
% Keep the result in the workspace: when parts 7 and 8 both run, the second prints the two side by side
if ~exist('Results_EM','var')
    Results_EM=struct();
end
Results_EM.(ACdefinition).EstimParams=EstimParams_EM; Results_EM.(ACdefinition).estsummary=estsummary_EM; Results_EM.(ACdefinition).estimtime=estimtime_EM;
Results_EM.(ACdefinition).targetvec=targetvec_EM; Results_EM.(ACdefinition).sevec=sevec_EM; Results_EM.(ACdefinition).EMmomenttype=EMmomenttype; Results_EM.(ACdefinition).okEM=okEM;
if isfield(Results_EM,'toolkit') && isfield(Results_EM,'DFJ')
    fprintf('\nThe two estimations side by side (standard errors in brackets): toolkit definition | DFJ2010 definition | DFJ2010 Table A1 column 5\n')
    dfjA1=[2.15,3.19,0.99,202]; dfjse=[0.055,0.054,0.01,23.4];
    for pp=1:length(EstimParamNames_EM)
        fprintf('   %-6s %10.4f (%.4f) | %10.4f (%.4f) | %10.4f (%.4f)\n',EstimParamNames_EM{pp},Results_EM.toolkit.EstimParams.(EstimParamNames_EM{pp}),Results_EM.toolkit.estsummary.EstimParamsStdDev.(EstimParamNames_EM{pp}),Results_EM.DFJ.EstimParams.(EstimParamNames_EM{pp}),Results_EM.DFJ.estsummary.EstimParamsStdDev.(EstimParamNames_EM{pp}),dfjA1(pp),dfjse(pp))
    end
    fprintf('   objective %.1f on %d moments | %.1f on %d moments | 667.0 on 249\n',Results_EM.toolkit.estsummary.objectivefnval,sum(Results_EM.toolkit.okEM),Results_EM.DFJ.estsummary.objectivefnval,sum(Results_EM.DFJ.okEM))
    save ./SavedOutput/DFJ2010_EndogMedex_Estimation_both.mat Results_EM EstimParamNames_EM
end
