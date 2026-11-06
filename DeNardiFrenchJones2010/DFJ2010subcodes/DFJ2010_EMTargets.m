% DFJ2010_EMTargets: the GMM targets of the ENDOGENOUS medical-expense model (DFJ2010 Sec. IX; online appendix A.2) in
% the toolkit's cohort format, with the weighting and covariance matrices: the median-asset targets of the exogenous model
% (DFJ2010_AssetTargets) plus, by cohort (5) x PI half (2) x wave (1998-2006), the mean and the 90th percentile of
% out-of-pocket medical expenses (A2, A3) and the first and second autocorrelations of log out-of-pocket medical expenses
% (A4, A5). A wave in year y reports the two years y-1 and y, and DFJ2010's criterion (gmmcrit51.g line 215) averages
% adjacent simulated years before every medical-expense moment, so the model objects are the toolkit's 2-period averages
% (simoptions.ntimeperiods=2, rows oopmedex_2periodAverage) starting at the cohort's age in year y-1; the autocorrelation
% with the previous wave is horizon 2 (two model periods) between those windows, with the wave before that horizon 4, of
% the log of the window (oopmedex_2periodAverage_log, simoptions.logwhichcombos), among those alive at both window starts.
% The 1998 first autocorrelation and the 2000 second autocorrelation pair a wave with the 1996 observation, which in
% DFJ2010's simulation is each person's DATA value (mxsim96); the model has no such object, so those 20 moments are not
% targeted. The 90th percentile is a level moment of the window (QuantileCutoffs row 10 with nquantiles=10; DFJ2010's is
% an indicator moment, the same difference as for the asset medians, see the write-up); set usep90=0 to leave it out.
% Params.agejshifter, simoptions and FnsToEvaluate_val (oopmedex) in the workspace.
%
% Produces: TargetMoments_EM, WeightingMatrix_EM, CoVarMatrixDataMoments_EM (diagonal, from the bootstrap standard
% errors of every moment; the asset medians' bootstrap covariance across cells is not used here), FnsToEvaluate_EM,
% simoptions_EM, and targetvec_EM/sevec_EM/EMmomenttype (the stacked targets, their standard errors and a type label per
% entry, in the toolkit's stacking order, NaN where a cell has no data).
%
% Ordering: the toolkit stacks the moments cohort by cohort (AllStats, AgeConditionalStats, AutoCorr), within those in
% the order the fields were created, within a statistic by age. Rather than reproduce that by hand (as tkorder does for
% the asset medians), a twin structure holding the standard errors in the same slots is run through the same parser
% (SetupTargetMoments_FHorz_withCohorts), which guarantees the weighting matrix lines up with the targets.

if ~exist('ACdefinition','var')
    ACdefinition='toolkit'; % 'toolkit': the autocorrelations are AutoCorr targets (data rows ac1/ac2, the pair population's own means and sds);
                            % 'DFJ': CustomModelStats targets (DFJ2010_EM_CustomStats; data rows ac1dfj/ac2dfj, the full cell's means and sds)
end
if strcmp(ACdefinition,'DFJ')
    acsuffix='dfj';
else
    acsuffix='';
end
years_EM=DataWork.Targets.years; % 1996:2:2006; year index 1 is the initial condition and never a moment
cohortj0=cohortentryage-Params.agejshifter; % model age index at which each cohort enters

%% Functions and restrictions the targets need
FnsToEvaluate_EM.assets=@(aprime,a,h,zeta,xi) a;
FnsToEvaluate_EM.oopmedex=FnsToEvaluate_val.oopmedex;     % floored at 125 per period (DFJ2010 bottom-code the two-year average at 250, gmmcrit51.g; half of it per period is the toolkit-window analogue), so every 2-period window is positive and its log defined
simoptions_EM=simoptions;
simoptions_EM.windowmethod='tdigest1'; % the 2-period window statistics (means exact, 90th percentiles from the t-digest, approximate) in one enumeration in blocks folded into summaries; 'tdigest2' would make the percentiles exact at about twice the time, the default 'list' ran out of memory at the first objective evaluation (2026-10-15). Runtime is the binding constraint (user, 2026-10-20).
simoptions_EM.nquantiles=10; % QuantileCutoffs has nquantiles+1 rows incl. min and max, so row 10 is the 90th percentile
usep90=1; % =1 targets the 90th percentiles of the 2-period averages as level moments; =0 leaves them out (227 moments as before 2026-10-14)
simoptions_EM.conditionalrestrictions=struct();
for q=1:5
    simoptions_EM.conditionalrestrictions.(['aliveq',num2str(q)])=simoptions.conditionalrestrictions.(['aliveq',num2str(q)]);
end
simoptions_EM.conditionalrestrictions.aliveh1=simoptions.conditionalrestrictions.aliveh1;
simoptions_EM.conditionalrestrictions.aliveh2=simoptions.conditionalrestrictions.aliveh2;
% (the autocorrelation command uses the Tan two-step push, one per-age policy step at a time, so no transition matrix is
% formed and no memory option is needed)

%% Targets and their standard errors, same structures
TargetMoments_EM=struct();
SEMoments_EM=struct();
% (i) the asset medians, as in the exogenous model (TargetMoments and MomentList from DFJ2010_AssetTargets)
TargetMoments_EM.AgeConditionalStats=TargetMoments.AgeConditionalStats;
for c=1:ncohorts
    for q=1:5
        rows=find(MomentList.cohort==c & MomentList.PIq==q);
        if ~isempty(rows)
            v=nan(N_j,1);
            v(cohortentryage(c)+(MomentList.year(rows)-1996)-Params.agejshifter)=MomentList.se(rows);
            SEMoments_EM.AgeConditionalStats.(['cohort',num2str(c)]).(['aliveq',num2str(q)]).assets.Median=v;
        end
    end
end
% (ii) mean and 90th percentile of out-of-pocket medical expenses by cohort x PI half x wave, and (iii) the autocorrelations
MedexList=DataWork.Targets.MedexList;
nMean=0; nP90=0; nAC1=0; nAC2=0;
for c=1:ncohorts
    cstr=['cohort',num2str(c)];
    for hh=1:2
        hstr=['aliveh',num2str(hh)];
        vmean=nan(N_j,1); smean=nan(N_j,1);
        vp90=nan(simoptions_EM.nquantiles+1,N_j); sp90=nan(simoptions_EM.nquantiles+1,N_j); % a matrix target: the age is its column
        vac1=nan(N_j-2,1); sac1=nan(N_j-2,1);
        vac2=nan(N_j-4,1); sac2=nan(N_j-4,1);
        for y=2:length(years_EM)
            s=cohortj0(c)+(years_EM(y)-1996)-1; % start age of the 2-period window (years y-1, y) of this cohort's wave y
            if s>N_j-1 % the window runs past the last age
                continue
            end
            r=find(MedexList.cohort==c & MedexList.half==hh & MedexList.year==years_EM(y) & strcmp(MedexList.stat,'mean'));
            if ~isempty(r)
                vmean(s)=MedexList.value(r); smean(s)=MedexList.se(r); nMean=nMean+1;
            end
            r=find(MedexList.cohort==c & MedexList.half==hh & MedexList.year==years_EM(y) & strcmp(MedexList.stat,'p90'));
            if ~isempty(r) && usep90==1
                vp90(10,s)=MedexList.value(r); sp90(10,s)=MedexList.se(r); nP90=nP90+1;
            end
            r=find(MedexList.cohort==c & MedexList.half==hh & MedexList.year==years_EM(y) & strcmp(MedexList.stat,['ac1',acsuffix]));
            if ~isempty(r) && s-2>=cohortj0(c) % AC1: the windows of waves y-1 and y, horizon 2 from start age s-2 (not when wave y-1 is 1996)
                vac1(s-2)=MedexList.value(r); sac1(s-2)=MedexList.se(r); nAC1=nAC1+1;
            end
            r=find(MedexList.cohort==c & MedexList.half==hh & MedexList.year==years_EM(y) & strcmp(MedexList.stat,['ac2',acsuffix]));
            if ~isempty(r) && s-4>=cohortj0(c) % AC2: the windows of waves y-2 and y, horizon 4 from start age s-4 (not when wave y-2 is 1996)
                vac2(s-4)=MedexList.value(r); sac2(s-4)=MedexList.se(r); nAC2=nAC2+1;
            end
        end
        if any(~isnan(vmean))
            TargetMoments_EM.AgeConditionalStats.(cstr).(hstr).oopmedex_2periodAverage.Mean=vmean;
            SEMoments_EM.AgeConditionalStats.(cstr).(hstr).oopmedex_2periodAverage.Mean=smean;
        end
        if any(~isnan(vp90),'all')
            TargetMoments_EM.AgeConditionalStats.(cstr).(hstr).oopmedex_2periodAverage.QuantileCutoffs=vp90;
            SEMoments_EM.AgeConditionalStats.(cstr).(hstr).oopmedex_2periodAverage.QuantileCutoffs=sp90;
        end
        if strcmp(ACdefinition,'DFJ') % CustomModelStats targets, the fields DFJ2010_EM_CustomStats returns (index = window start age)
            if any(~isnan(vac1))
                TargetMoments_EM.CustomModelStats.(cstr).(['AC1_h',num2str(hh)])=vac1;
                SEMoments_EM.CustomModelStats.(cstr).(['AC1_h',num2str(hh)])=sac1;
            end
            if any(~isnan(vac2))
                TargetMoments_EM.CustomModelStats.(cstr).(['AC2_h',num2str(hh)])=vac2;
                SEMoments_EM.CustomModelStats.(cstr).(['AC2_h',num2str(hh)])=sac2;
            end
        else % AutoCorr targets (the toolkit's definition)
            if any(~isnan(vac1))
                TargetMoments_EM.AutoCorr.(cstr).(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon2=vac1;
                SEMoments_EM.AutoCorr.(cstr).(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon2=sac1;
            end
            if any(~isnan(vac2))
                TargetMoments_EM.AutoCorr.(cstr).(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon4=vac2;
                SEMoments_EM.AutoCorr.(cstr).(hstr).oopmedex_2periodAverage_log.AutoCorrelation_horizon4=sac2;
            end
        end
    end
end

%% Stack them in the toolkit's order (the same parser the estimator uses), and the weighting matrix over the non-NaN entries
[targetvec_EM,cohortmoments_EM]=SetupTargetMoments_FHorz_withCohorts(TargetMoments_EM,FnsToEvaluate_EM,N_j,cohortj0,1,simoptions_EM,Names_i);
[sevec_EM,~]=SetupTargetMoments_FHorz_withCohorts(SEMoments_EM,FnsToEvaluate_EM,N_j,cohortj0,1,simoptions_EM,Names_i);
if length(sevec_EM)~=length(targetvec_EM) || any(isnan(sevec_EM)~=isnan(targetvec_EM))
    error('DFJ2010_EMTargets: the standard errors do not line up with the targets')
end
% A type label per stacked entry, for reporting the objective by moment type (the parser's names: AllStats first, then
% AgeConditionalStats, then AutoCorr, cohort by cohort; a vector target's block is N_j, N_j-2 or N_j-4 long, a QuantileCutoffs one (nquantiles+1)*N_j)
EMmomenttype=cell(length(targetvec_EM),1);
for c=1:ncohorts
    pos=cohortmoments_EM.offset(c);
    names=cohortmoments_EM.acsmomentnames{c}; sizes=cohortmoments_EM.acscummomentsizes{c};
    for mm=1:size(names,1)*(sizes(end)>0)
        if mm==1, i1=pos+1; else, i1=pos+sizes(mm-1)+1; end
        i2=pos+sizes(mm);
        EMmomenttype(i1:i2)={[names{mm,2},'.',names{mm,3}]}; % assets.Median, oopmedex_2periodAverage.Mean or .QuantileCutoffs
    end
    pos=pos+sizes(end);
    names=cohortmoments_EM.acrmomentnames{c}; sizes=cohortmoments_EM.acrcummomentsizes{c}; % (placeholders of size 0 when the cohort has no AutoCorr targets, e.g. ACdefinition=DFJ)
    for mm=1:size(names,1)*(sizes(end)>0)
        if mm==1, i1=pos+1; else, i1=pos+sizes(mm-1)+1; end
        i2=pos+sizes(mm);
        EMmomenttype(i1:i2)={[names{mm,2},'.',names{mm,4}]}; % oopmedex_2periodAverage_log.AutoCorrelation_horizon2 or _horizon4
    end
    pos=pos+sizes(end);
    names=cohortmoments_EM.cmsmomentnames{c}; sizes=cohortmoments_EM.cmscummomentsizes{c}; % the DFJ-definition autocorrelations (ACdefinition=DFJ)
    for mm=1:length(names)*(sizes(end)>0)
        if mm==1, i1=pos+1; else, i1=pos+sizes(mm-1)+1; end
        i2=pos+sizes(mm);
        EMmomenttype(i1:i2)={['oopmedex_2periodAverage_log.',names{mm}(1:3),'dfj']}; % AC1dfj or AC2dfj
    end
end
okEM=~isnan(targetvec_EM);
WeightingMatrix_EM=diag(1./sevec_EM(okEM).^2);
CoVarMatrixDataMoments_EM=diag(sevec_EM(okEM).^2);
fprintf('Endogenous-model GMM targets: %d moments = %d asset medians + %d mean medex + %d 90th pct medex + %d AC1 + %d AC2, all medex moments of 2-year averages, autocorrelations in the %s definition (DFJ2010: 249 = 101 + 148), over %d cohorts\n',sum(okEM),Jdata,nMean,nP90,nAC1,nAC2,ACdefinition,ncohorts)
