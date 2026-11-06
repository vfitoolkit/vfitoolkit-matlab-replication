% DFJ2010_EndogMedex2: the endogenous medical-expense model, part 2 (stage 1 of the proposal): the initial distribution
% and cohorts (rebuilt on the endogenous shock grids), the asset targets, and the comparison with the data at the
% published parameters: median assets by cohort x PI quintile x year (DFJ2010 Figure 12), mean and 90th percentile of
% out-of-pocket medical expenses by cohort x PI half x wave (Figure A4; 2-year averages as the data, toolkit windows), the autocorrelations of log medical expenses
% (of the log of the 2-year average, toolkit PType autocorrelation command, horizons 2 and 4 among those alive at both window starts), and DFJ2010's full GMM
% objective at the published parameters over the endogenous model's targets (DFJ2010_EMTargets: asset medians, mean
% medex, 90th percentiles, AC1, AC2), hand-assembled and then by the estimation command itself at the published parameters (one
% evaluation per autocorrelation definition: estimoptions.skipestimation=1 with skipjacobian=1). Called from DeNardiFrenchJones2010.m as doPart(6), after
% DFJ2010_EndogMedex1. No estimation (that is DFJ2010_EndogMedex3, doPart(7)).

%% Initial distribution, cohorts and quintile restrictions (shared subcode)
DFJ2010_InitialDistCohorts
% PI half of each ptype for the medical-expense moments (DFJ2010 appendix A.2: Q=2 halves of the PI distribution; the
% quintile ptypes sit at PI percentiles 10,30,50,70,90, and (q-1)/Q<=I<q/Q puts the 50th percentile in the upper half)
Params.half.femaleq1=1; Params.half.femaleq2=1; Params.half.femaleq3=2; Params.half.femaleq4=2; Params.half.femaleq5=2;
Params.half.maleq1=1;   Params.half.maleq2=1;   Params.half.maleq3=2;   Params.half.maleq4=2;   Params.half.maleq5=2;
simoptions.conditionalrestrictions.aliveh1=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==1);
simoptions.conditionalrestrictions.aliveh2=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==2);

%% GMM targets: median assets by cohort, quintile and year (shared subcode)
DFJ2010_AssetTargets
MedexCells=DataWork.Targets.MedexCells;
MedexList=DataWork.Targets.MedexList;

%% Model vs data, cohort by cohort (as the estimator does it)
FnsToEvaluate_val.assets=@(aprime,a,h,zeta,xi) a;
FnsToEvaluate_val.oopmedex=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    max(DFJ2010EM_OOPMedex(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6),125); % floored at 125 per period: DFJ2010 bottom-code the TWO-YEAR average at minmedex=250 (gmmcrit51.g); the toolkit averages the raw function, so a per-period floor of half that is used (equal when both years are above 250; see the write-up)
FnsToEvaluate_val.totalmedex=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_TotalMedex(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6);
FnsToEvaluate_val.consumption=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_Consumption(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6);
FnsToEvaluate_val.transfers=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_Transfers(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6);
simoptions_val=simoptions;
simoptions_val.whichstats=[1,1,0,0,0,1,0]; % mean, median, quantile cutoffs (the same for every window length)
simoptions_val.nquantiles=10; % QuantileCutoffs has nquantiles+1 rows incl. min and max, so row 10 is the 90th percentile
simoptions_val.agemass_withCohort=1;
% Medical expenses in the data are 2-year averages (a wave reports the two years before it, and DFJ2010's gmmcrit51.g averages
% adjacent simulated years before every medical-expense moment): simoptions.ntimeperiods=[1,2] gives the single-period
% statistics (the asset medians; the annual medex of the cohort-1 print) and the 2-period averages (oopmedex_2periodAverage:
% the means and 90th percentiles the data measure, by window start age). The windows need the z transitions, and whichcombos
% limits the (expensive, exact) window statistics to out-of-pocket medex on the two PI-half pages.
simoptions_val.ntimeperiods=[1,2];
simoptions_val.pi_z=pi_z_J;
simoptions_val.windowmethod='tdigest1'; % the window statistics in one enumeration in blocks folded into summaries (memory bounded by simoptions.windowblocksize): the Mean exact, the 90th percentile from the t-digest (approximate, the mean of the centroid holding its rank); 'tdigest2' makes the quantiles exact with a second enumeration at about twice the time, and the default 'list' keeps every path and ran out of memory in the estimator (2026-10-15). Runtime is the binding constraint (user, 2026-10-20).
fnval=fieldnames(FnsToEvaluate_val); rnval=fieldnames(simoptions_val.conditionalrestrictions);
simoptions_val.whichcombos=zeros(2*length(fnval),N_j,1+length(rnval)); % rows: the M=1 block then the M=2 block; pages: unrestricted then the restrictions
simoptions_val.whichcombos(1:length(fnval),:,:)=1; % every single-period statistic
% (the 2-period average of oopmedex on the PI-half pages is switched on per cohort in the loop, at the window start ages the moments use: the cohort's age in 1997, 1999, ..., 2005)
iow=length(fnval)+find(strcmp(fnval,'oopmedex')); pagesw=1+[find(strcmp(rnval,'aliveh1')),find(strcmp(rnval,'aliveh2'))];
% Autocorrelations of the log of the 2-period average of out-of-pocket medex (DFJ2010 appendix A.2 eqns A4-A5; gmmcrit51.g logs
% the averaged series): between the windows of consecutive waves (horizon 2) and of waves two apart (horizon 4), among those
% alive at both window start ages, from the toolkit's distribution-based autocorrelation command (no panel simulation).
% whichcombos (set per cohort in the loop) limits it to the window start ages the moments use.
FnsToEvaluate_ac.oopmedex=FnsToEvaluate_val.oopmedex; % floored at 125 per period (half of DFJ2010's 250 on the two-year average), so every window is positive and its log is defined
simoptions_ac=simoptions;
simoptions_ac.conditionalrestrictions=struct(); % the PI halves only: pages unrestricted, aliveh1, aliveh2
simoptions_ac.conditionalrestrictions.aliveh1=simoptions.conditionalrestrictions.aliveh1;
simoptions_ac.conditionalrestrictions.aliveh2=simoptions.conditionalrestrictions.aliveh2;
simoptions_ac.ntimeperiods=2;
simoptions_ac.logwhichcombos=1; % the statistics of log(oopmedex_2periodAverage), reported as oopmedex_2periodAverage_log
simoptions_ac.timehorizons=[2,4];
simoptions_ac.agemass_withCohort=1;
haveautocorr=1;

%% The endogenous model's GMM targets in the toolkit's format (shared with the estimation, DFJ2010_EndogMedex3)
DFJ2010_EMTargets
ModelMoments_EM=struct(); % the model's values in the same slots, filled cohort by cohort below

years=1996:2:2006;
modelmedian=nan(Jdata,1);
ModelMedex.Mean=nan(5,2,6); ModelMedex.P90=nan(5,2,6); ModelMedex.AC1=nan(5,2,6); ModelMedex.AC2=nan(5,2,6); ModelMedex.AC1dfj=nan(5,2,6); ModelMedex.AC2dfj=nan(5,2,6);
ModelCMS=struct(); % the DFJ-definition autocorrelations by cohort, in the slots of the CustomModelStats targets (ACdefinition=DFJ)
Params_c=Params;
tic; cohorttimes=zeros(ncohorts,3); % per cohort: agent distribution, LifeCycleProfiles (incl. the 2-period windows: means, 90th percentiles), AutoCorr (the log windows)
for c=1:ncohorts
    jequaloneDist_c=struct();
    cohortmass_c=zeros(1,length(Names_i));
    for ii=1:length(Names_i)
        temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c)-Params.agejshifter);
        cohortmass_c(ii)=gather(sum(temp(:)));
        if cohortmass_c(ii)>0
            jequaloneDist_c.(Names_i{ii})=temp/cohortmass_c(ii);
        else
            jequaloneDist_c.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
        end
    end
    Params_c.ptypemasses=Params.ptypemasses(:)'.*cohortmass_c;
    Params_c.ptypemasses=Params_c.ptypemasses/sum(Params_c.ptypemasses);
    Params_c.mewj=Params.mewj(cohortentryage(c)-Params.agejshifter)*ones(1,N_j);
    simoptions_val.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
    simoptions_ac.jequaloneDistAge=simoptions_val.jequaloneDistAge;
    t_c=tic;
    starts=cohortj0(c)+[1,3,5,7,9]; starts=starts(starts<=N_j-1); % the window start ages of this cohort's waves 1998..2006
    simoptions_val.whichcombos(iow,:,pagesw)=0; simoptions_val.whichcombos(iow,starts,pagesw)=1;
    StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_val);
    cohorttimes(c,1)=toc(t_c); t_c=tic;
    ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy,FnsToEvaluate_val,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_val);
    cohorttimes(c,2)=toc(t_c);
    rows=find(MomentList.cohort==c);
    for rr=rows'
        modelmedian(rr)=gather(ACS_c.(['aliveq',num2str(MomentList.PIq(rr))]).assets.Median(cohortentryage(c)+MomentList.year(rr)-1996-Params.agejshifter));
    end
    % Autocorrelations: the window start ages the moments use are the cohort's age in 1997, 1999, 2001 and 2003 (horizon 2 to the
    % next wave's window, horizon 4 to the one after; the 1998 and 2000 waves' partners back to 1996 are not model objects)
    simoptions_ac.whichcombos=zeros(1,N_j,3);
    origins=cohortj0(c)+[1,3,5,7]; origins=origins(origins<=N_j-3);
    simoptions_ac.whichcombos(1,origins,2:3)=1;
    t_c=tic;
    AC_c=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_c,Policy,FnsToEvaluate_ac,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,simoptions_ac);
    cohorttimes(c,3)=toc(t_c);
    fprintf('   cohort %d (enters at %d): agent distribution %.1f s, LifeCycleProfiles with the 2-period windows %.1f s, AutoCorr with the log windows %.1f s\n',c,cohortentryage(c),cohorttimes(c,1),cohorttimes(c,2),cohorttimes(c,3))
    % DFJ2010's own definition of the autocorrelation (full-cell means and sds; the formula of DFJ2010_EM_CustomStats), from the same outputs
    cstr=['cohort',num2str(c)];
    for hh=1:2
        R=AC_c.(['aliveh',num2str(hh)]).oopmedex_2periodAverage_log;
        mR=gather(R.Mean); sR=gather(R.StdDeviation);
        for K=[2,4]
            num=gather(R.(['AutoCovariance_horizon',num2str(K)]))+(gather(R.(['PairMean_j_horizon',num2str(K)]))-mR(1:N_j-K)).*(gather(R.(['PairMean_jplusk_horizon',num2str(K)]))-mR(1+K:N_j));
            ModelCMS.(cstr).(['AC',num2str(K/2),'_h',num2str(hh)])=num./(sR(1:N_j-K).*sR(1+K:N_j));
        end
    end
    % Medical expenses: a wave in year y covers years y-1 and y, so its model analogue is the 2-period average starting at the
    % cohort's age in year y-1
    for hh=1:2
        hstr=['aliveh',num2str(hh)];
        for y=2:6
            s=cohortj0(c)+(years(y)-1996)-1; % start age of the window (y-1, y)
            if s<=N_j-1
                ModelMedex.Mean(c,hh,y)=gather(ACS_c.(hstr).oopmedex_2periodAverage.Mean(s));
                ModelMedex.P90(c,hh,y)=gather(ACS_c.(hstr).oopmedex_2periodAverage.QuantileCutoffs(10,s));
                if s-2>=cohortj0(c) % AC1: the windows of waves y-1 and y
                    ModelMedex.AC1(c,hh,y)=gather(AC_c.(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon2(s-2));
                    ModelMedex.AC1dfj(c,hh,y)=ModelCMS.(cstr).(['AC1_h',num2str(hh)])(s-2);
                end
                if s-4>=cohortj0(c) % AC2: the windows of waves y-2 and y
                    ModelMedex.AC2(c,hh,y)=gather(AC_c.(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon4(s-4));
                    ModelMedex.AC2dfj(c,hh,y)=ModelCMS.(cstr).(['AC2_h',num2str(hh)])(s-4);
                end
            end
        end
    end
    % the model's values of this cohort's targets, in the targets' slots (the rows and statistics DFJ2010_EMTargets built; the
    % parser vectorises a QuantileCutoffs matrix the same way on both sides)
    cstr=['cohort',num2str(c)];
    if isfield(TargetMoments_EM.AgeConditionalStats,cstr)
        rnames=fieldnames(TargetMoments_EM.AgeConditionalStats.(cstr)); % the restrictions
        for rr=1:length(rnames)
            fnames=fieldnames(TargetMoments_EM.AgeConditionalStats.(cstr).(rnames{rr})); % the rows: assets, oopmedex_2periodAverage
            for ff=1:length(fnames)
                snames=fieldnames(TargetMoments_EM.AgeConditionalStats.(cstr).(rnames{rr}).(fnames{ff}));
                for ss=1:length(snames)
                    ModelMoments_EM.AgeConditionalStats.(cstr).(rnames{rr}).(fnames{ff}).(snames{ss})=gather(ACS_c.(rnames{rr}).(fnames{ff}).(snames{ss}));
                end
            end
        end
    end
    if isfield(TargetMoments_EM,'AutoCorr') && isfield(TargetMoments_EM.AutoCorr,cstr)
        rnames=fieldnames(TargetMoments_EM.AutoCorr.(cstr));
        for rr=1:length(rnames)
            fnames=fieldnames(TargetMoments_EM.AutoCorr.(cstr).(rnames{rr})); % the row: oopmedex_2periodAverage_log
            for ff=1:length(fnames)
                snames=fieldnames(TargetMoments_EM.AutoCorr.(cstr).(rnames{rr}).(fnames{ff}));
                for ss=1:length(snames)
                    ModelMoments_EM.AutoCorr.(cstr).(rnames{rr}).(fnames{ff}).(snames{ss})=gather(AC_c.(rnames{rr}).(fnames{ff}).(snames{ss})(:));
                end
            end
        end
    end
    if c==1
        fprintf('Cohort 1 (enters at 74): mean out-of-pocket medex by age 74..84, lower half %s, upper half %s (thousands); mean total medical goods upper half %s; mean transfers upper half %s\n',mat2str(round(gather(ACS_c.aliveh1.oopmedex.Mean(5:15))/1000,1)),mat2str(round(gather(ACS_c.aliveh2.oopmedex.Mean(5:15))/1000,1)),mat2str(round(gather(ACS_c.aliveh2.totalmedex.Mean(5:15))/1000,1)),mat2str(round(gather(ACS_c.aliveh2.transfers.Mean(5:15))/1000,2)))
    end
end
cohortlooptime=toc;
fprintf('Five-cohort loop (distributions, medians, medical-expense statistics) took %.1f seconds: agent distributions %.1f, LifeCycleProfiles %.1f, AutoCorr %.1f\n',cohortlooptime,sum(cohorttimes(:,1)),sum(cohorttimes(:,2)),sum(cohorttimes(:,3)))

%% Assets: the objective at the published parameters (DFJ2010: 667 on 242 df for their full criterion) and Figure 12
gap=modelmedian(tkorder)-MomentList.median(tkorder);
fprintf('\nAsset-median part of the GMM objective at DFJ2010''s published endogenous-model parameters (diagonal W): %.1f over %d moments (exogenous benchmark at its published parameters: 115.5)\n',gap'*WeightingMatrix*gap,Jdata)
fprintf('   by quintile: %s;  by cohort: %s\n',mat2str(round(arrayfun(@(q) sum((gap(MomentList.PIq(tkorder)==q)).^2./MomentList.se(tkorder(MomentList.PIq(tkorder)==q)).^2),1:5),1)),mat2str(round(arrayfun(@(c) sum((gap(MomentList.cohort(tkorder)==c)).^2./MomentList.se(tkorder(MomentList.cohort(tkorder)==c)).^2),1:5),1)))
% DFJ2010 Figure 12: as Figure 5, the cohorts aged 74 and 84 in 1996, data (solid) and the endogenous model (dashed), 1996-2006;
% both lines start from the 1996 data median, the initial distribution
figure_c=figure_c+1;
figure(figure_c);
hold on
for c=[1,3]
    for q=1:5
        rows=find(MomentList.cohort==c & MomentList.PIq==q);
        [~,yy]=sort(MomentList.year(rows)); rows=rows(yy);
        ages=[cohortentryage(c); cohortentryage(c)+MomentList.year(rows(:))-1996];
        plot(ages,[DataWork.Targets.MedianAssets(c,q,1); MomentList.median(rows(:))]/1000,'k-',ages,[DataWork.Targets.MedianAssets(c,q,1); modelmedian(rows(:))]/1000,'k--')
    end
end
hold off
xlim([70,102]); xlabel('Age'); ylabel('Assets (000s of 1998 dollars)')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_EM_Fig12_assets.png')

%% The full GMM objective of the endogenous model at the published parameters (DFJ2010: 667.0 on 242 df over 249 moments)
% The model's values are stacked by the same parser as the targets (DFJ2010_EMTargets), so the two vectors line up entry by
% entry. Conventions that differ from DFJ2010's: the model AC1/AC2 are the toolkit's pair-based correlations of the log of the 2-year
% average (centred and scaled by the pair population's means and sds; DFJ: by the full cell's); the 90th percentile is a level
% moment (DFJ: an indicator moment); the 20 autocorrelations with the 1996 wave are not targeted; diagonal weighting from the bootstrap standard errors.
[modelvec_EM,~]=SetupTargetMoments_FHorz_withCohorts(ModelMoments_EM,FnsToEvaluate_EM,N_j,cohortj0,1,simoptions_EM,Names_i);
contrib_EM=((modelvec_EM-targetvec_EM)./sevec_EM).^2;
fprintf('\nFull GMM objective of the endogenous model at DFJ2010''s published parameters (diagonal W): %.1f over %d moments (DFJ2010: 667.0 over 249)\n',sum(contrib_EM(okEM)),sum(okEM))
typelist=unique(EMmomenttype(okEM));
for tt=1:length(typelist)
    ind=okEM & strcmp(EMmomenttype,typelist{tt});
    fprintf('   %-32s %4d moments, objective %7.1f, mean |model-data|/se %.2f\n',typelist{tt},sum(ind),sum(contrib_EM(ind)),mean(sqrt(contrib_EM(ind))))
end
if any(isnan(modelvec_EM(okEM)))
    fprintf('   WARNING: %d targeted moments have a NaN model value (the objective above omits them)\n',sum(isnan(modelvec_EM(okEM))))
end

% The same objective with the autocorrelations in DFJ2010's definition (data rows ac1dfj/ac2dfj, model values from the same push)
ACdefinition='DFJ'; DFJ2010_EMTargets
ModelMoments_EM_DFJ=struct(); ModelMoments_EM_DFJ.AgeConditionalStats=ModelMoments_EM.AgeConditionalStats;
cnames=fieldnames(TargetMoments_EM.CustomModelStats); % the model values in exactly the targets' slots (same cohorts, same fields, same order: the parser stacks the fields of the structure it is given)
for cc2=1:length(cnames)
    fnames=fieldnames(TargetMoments_EM.CustomModelStats.(cnames{cc2}));
    for ff=1:length(fnames)
        ModelMoments_EM_DFJ.CustomModelStats.(cnames{cc2}).(fnames{ff})=ModelCMS.(cnames{cc2}).(fnames{ff});
    end
end
[modelvec_EM_DFJ,~]=SetupTargetMoments_FHorz_withCohorts(ModelMoments_EM_DFJ,FnsToEvaluate_EM,N_j,cohortj0,1,simoptions_EM,Names_i);
contrib_EM_DFJ=((modelvec_EM_DFJ-targetvec_EM)./sevec_EM).^2;
fprintf('\nThe same with the autocorrelations in DFJ2010''s definition (full-cell means and sds, CustomModelStats targets): %.1f over %d moments\n',sum(contrib_EM_DFJ(okEM)),sum(okEM))
typelist=unique(EMmomenttype(okEM));
for tt=1:length(typelist)
    ind=okEM & strcmp(EMmomenttype,typelist{tt});
    fprintf('   %-40s %4d moments, objective %7.1f, mean |model-data|/se %.2f\n',typelist{tt},sum(ind),sum(contrib_EM_DFJ(ind)),mean(sqrt(contrib_EM_DFJ(ind))))
end

%% The estimator's own objective at the published parameters (estimoptions.skipestimation=1 with skipjacobian=1), both definitions
% The hand-assembled objectives above stack this part's model moments through the targets' parser. Here the estimation command itself
% (as DFJ2010_EndogMedex3 calls it, same targets, weights, options and CustomModelStats) evaluates its objective once at the published
% parameters: skipestimation=1 takes the input parameters as the estimates, and skipjacobian=1 (toolkit, 2026-10-23) stops it from
% then computing the Jacobian for the standard errors (1+2*4 evaluations at each of four step sizes, about nine hours here), so it is
% one evaluation (value fn plus the five cohorts, about 15 minutes) per definition. Its objective should equal the hand-assembled one
% up to roundoff (the same model, the same moments, the same stacking), and its per-cohort timing lines are those of one estimation step.
EstimatorEval_EM=struct();
for dd=1:2
    if dd==1
        ACdefinition='toolkit'; handcontrib=contrib_EM;
    else
        ACdefinition='DFJ'; handcontrib=contrib_EM_DFJ;
    end
    DFJ2010_EMTargets % TargetMoments_EM, WeightingMatrix_EM, CoVarMatrixDataMoments_EM, FnsToEvaluate_EM, simoptions_EM, okEM, EMmomenttype in this definition
    estimoptions_ev=struct();
    estimoptions_ev.verbose=1;
    estimoptions_ev.cohortagejshifter=1;
    estimoptions_ev.skipestimation=1; % the input parameters are the estimates
    estimoptions_ev.skipjacobian=1; % one evaluation, no Jacobian (so no standard errors; those are DFJ2010_EndogMedex3's)
    estimoptions_ev.constrainpositive={'nu','omega','eMin'}; % as the estimation (irrelevant to a single evaluation, which turns the constraints off)
    estimoptions_ev.constrain0to1={'beta'};
    if strcmp(ACdefinition,'DFJ')
        estimoptions_ev.CustomModelStats=@DFJ2010_EM_CustomStats;
    end
    vfoptions_ev=vfoptions; vfoptions_ev.verbose=0; vfoptions_ev.verboseparams=0;
    fprintf('\nThe estimation command''s objective at the published parameters, %s definition (skipestimation=1, skipjacobian=1: one evaluation)\n',ACdefinition)
    evaltimer=tic;
    [~,~,estsummary_ev]=EstimateLifeCycleModel_PType_MethodOfMoments({'nu','omega','beta','eMin'},TargetMoments_EM,WeightingMatrix_EM,CoVarMatrixDataMoments_EM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,[],FnsToEvaluate_EM,estimoptions_ev,vfoptions_ev,simoptions_EM);
    evaltime=toc(evaltimer);
    contrib_ev=((estsummary_ev.currentmomentvec(:)-estsummary_ev.targetmomentvec(:))./sevec_EM(okEM)).^2; % the estimator's moments are the targeted (non-NaN) entries in the parser's order, so EMmomenttype(okEM) labels them
    fprintf('   objective %.1f over %d moments in %.1f minutes; hand-assembled above %.1f, difference %.3e (should be roundoff); targets vs the stacked data targets, max abs difference %.3e (should be zero)\n',estsummary_ev.objectivefnval,length(estsummary_ev.currentmomentvec),evaltime/60,sum(handcontrib(okEM)),estsummary_ev.objectivefnval-sum(handcontrib(okEM)),max(abs(estsummary_ev.targetmomentvec(:)-targetvec_EM(okEM))))
    typelist=unique(EMmomenttype(okEM)); labels=EMmomenttype(okEM); handok=handcontrib(okEM);
    for tt=1:length(typelist)
        ind=strcmp(labels,typelist{tt});
        fprintf('   %-40s %4d moments, objective %7.1f (hand-assembled %7.1f)\n',typelist{tt},sum(ind),sum(contrib_ev(ind)),sum(handok(ind)))
    end
    EstimatorEval_EM.(ACdefinition).objective=estsummary_ev.objectivefnval; EstimatorEval_EM.(ACdefinition).currentmomentvec=estsummary_ev.currentmomentvec; EstimatorEval_EM.(ACdefinition).targetmomentvec=estsummary_ev.targetmomentvec;
    EstimatorEval_EM.(ACdefinition).contrib=contrib_ev; EstimatorEval_EM.(ACdefinition).labels=labels; EstimatorEval_EM.(ACdefinition).evaltime=evaltime; EstimatorEval_EM.(ACdefinition).handobjective=sum(handcontrib(okEM));
end
ACdefinition='toolkit'; DFJ2010_EMTargets % back to the toolkit definition (DFJ2010_EndogMedex3 sets its own)

%% Medical expenses: mean and 90th percentile by cohort x PI half x year, data vs model (DFJ2010 Figure A4)
fprintf('\nOut-of-pocket medical expenses by cohort, PI half and year (thousands): data mean / model mean | data p90 / model p90\n')
for c=1:5
    for hh=1:2
        fprintf('   cohort %d half %d:',c,hh)
        for y=2:6
            fprintf('  %4.1f/%4.1f | %4.1f/%4.1f',MedexCells.Mean(c,hh,y)/1000,ModelMedex.Mean(c,hh,y)/1000,MedexCells.P90(c,hh,y)/1000,ModelMedex.P90(c,hh,y)/1000)
        end
        fprintf('\n')
    end
end
okm=~isnan(MedexCells.Mean(:,:,2:6)) & ~isnan(ModelMedex.Mean(:,:,2:6));
dm=MedexCells.Mean(:,:,2:6); mm=ModelMedex.Mean(:,:,2:6); dp=MedexCells.P90(:,:,2:6); mp=ModelMedex.P90(:,:,2:6);
fprintf('   over the %d cells with data: mean medex data %.0f vs model %.0f; 90th pct data %.0f vs model %.0f; corr of cell means %.3f\n',sum(okm(:)),mean(dm(okm)),mean(mm(okm)),mean(dp(okm)),mean(mp(okm)),corr(dm(okm),mm(okm)))
if haveautocorr
    fprintf('   autocorrelations of log medex, cohort 1 (rows: PI half; columns: waves), toolkit definition: data AC1 %s / model %s, data AC2 %s / model %s\n',mat2str(round(squeeze(MedexCells.AC1(1,:,2:6)),2)),mat2str(round(squeeze(ModelMedex.AC1(1,:,2:6)),2)),mat2str(round(squeeze(MedexCells.AC2(1,:,3:6)),2)),mat2str(round(squeeze(ModelMedex.AC2(1,:,3:6)),2)))
    fprintf('   DFJ2010 definition: data AC1 %s / model %s, data AC2 %s / model %s\n',mat2str(round(squeeze(MedexCells.AC1dfj(1,:,2:6)),2)),mat2str(round(squeeze(ModelMedex.AC1dfj(1,:,2:6)),2)),mat2str(round(squeeze(MedexCells.AC2dfj(1,:,3:6)),2)),mat2str(round(squeeze(ModelMedex.AC2dfj(1,:,3:6)),2)))
end
figure_c=figure_c+1;
figure(figure_c);
for hh=1:2
    subplot(2,2,hh)
    hold on
    for c=1:5
        ages=cohortentryage(c)+years(2:6)-1996;
        plot(ages,squeeze(MedexCells.Mean(c,hh,2:6))/1000,'k-',ages,squeeze(ModelMedex.Mean(c,hh,2:6))/1000,'b--')
    end
    hold off
    title(['Mean OOP medex, PI half ',num2str(hh)]); xlabel('age'); ylabel('thousands of 1998 dollars') % (data solid, model dashed: said in the write-up caption; a long title overlaps the neighbouring subplot)
    subplot(2,2,2+hh)
    hold on
    for c=1:5
        ages=cohortentryage(c)+years(2:6)-1996;
        plot(ages,squeeze(MedexCells.P90(c,hh,2:6))/1000,'k-',ages,squeeze(ModelMedex.P90(c,hh,2:6))/1000,'b--')
    end
    hold off
    title(['90th pct of OOP medex, PI half ',num2str(hh)]); xlabel('age'); ylabel('thousands of 1998 dollars')
end
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_EM_FigA4_medex.png')

save ./SavedOutput/DFJ2010_EndogMedex_Stage1.mat Params modelmedian MomentList tkorder ModelMedex MedexCells cohortentryage vftime cohortlooptime haveautocorr targetvec_EM sevec_EM modelvec_EM EMmomenttype okEM ModelCMS modelvec_EM_DFJ contrib_EM_DFJ EstimatorEval_EM
