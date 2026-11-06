% DFJ2010_EndogMedex5: the endogenous medical-expense model's figures at the estimates (DFJ2010 Figures 12, 13, 14 and A4).
% Called from DeNardiFrenchJones2010.m as doPart(10), on its own: panel A of Figure 14 is the EXOGENOUS model, so this part
% starts with the exogenous model in the workspace (solved at the current parameters, DFJ2010 Table 3 column 1, as parts 2-5) and
% only then sets up the endogenous model (DFJ2010_EndogMedex1), which overwrites it. The endogenous model is solved at the
% estimates of doPart(7) or doPart(8), read from SavedOutput/DFJ2010_EndogMedex_Estimation_AC<EMfigdefinition>.mat (EMfigdefinition
% is set by the main script: 'DFJ', the estimates with the autocorrelations in DFJ2010's definition, or 'toolkit').
%   Figure 12: median assets by cohort and PI quintile, data and model (as doPart(6) draws it at the published parameters)
%   Figure A4: mean and 90th percentile of out-of-pocket medical expenses by cohort and PI half, data and model (2-year averages)
%   Figure 13: median assets by PI quintile, cohort aged 74 in 1996, baseline vs no medical expenses (DFJ2010 Sec. IX.C: 'setting
%              the medical needs parameter mt to zero')
%   Figure 14: the same cohort, baseline vs 50 percent of the floor: A, the exogenous model (consumption floor cfloor halved); B,
%              the endogenous model (footnote 8: 'we index the estimated utility floor by the consumption level that provides the
%              floor when m=1. To run the counterfactual, we cut that consumption level in half', so eMin halved)
% The counterfactuals follow doPart(5) (DFJ2010 Sec. VIII): the policy is solved with mortality risk, and the cohort's distribution
% is iterated with the health transitions conditional on survival (no attrition), from age 74 to 100.
% Cost: two exogenous value fns (about 20 s each), the endogenous value fn at the published parameters inside DFJ2010_EndogMedex1
% and at the estimates (about 110 s each), the five cohorts' distributions and life-cycle profiles (about 10 minutes), and three
% endogenous value fns for the counterfactuals.

if ~exist('EMfigdefinition','var')
    EMfigdefinition='DFJ'; % set by DeNardiFrenchJones2010.m
end
figtimer=tic;

%% Figure 14, panel A: the EXOGENOUS model, baseline vs 50 percent of the consumption floor (before the endogenous model overwrites it)
% the health transitions conditional on survival (as doPart(5))
pi_z_J_noattrition=pi_z_J;
for ii=1:length(Names_i)
    pih=pi_h_J.(Names_i{ii});
    for jj=1:N_j
        for hh=1:2
            pih(hh,:,jj)=[pih(hh,1:2,jj)/sum(pih(hh,1:2,jj)),0,0];
        end
        pi_z_J_noattrition.(Names_i{ii})(:,:,jj)=kron(pi_zeta,pih(:,:,jj));
    end
end
% the cohort entering at 74
c_cf=1;
jequaloneDist_cf=struct();
cohortmass_cf=zeros(1,length(Names_i));
for ii=1:length(Names_i)
    temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c_cf)-Params.agejshifter);
    cohortmass_cf(ii)=gather(sum(temp(:)));
    if cohortmass_cf(ii)>0
        jequaloneDist_cf.(Names_i{ii})=temp/cohortmass_cf(ii);
    else
        jequaloneDist_cf.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
    end
end
simoptions_cf=simoptions_gmm; % aliveq1..5 restrictions, median only
simoptions_cf.jequaloneDistAge=cohortentryage(c_cf)-Params.agejshifter;
simoptions_cf.agemass_withCohort=1;
vfoptions_cf=vfoptions;
vfoptions_cf.verbose=0;
vfoptions_cf.verboseparams=0;
fprintf('\nDFJ2010 Figure 14 panel A: the exogenous model at Table 3 column %d (cfloor %.0f), baseline vs 50 percent of the consumption floor\n',Table3column(1),Params.cfloor)
cfnames_exog={'baseline','cfloor50'};
CF_exog=struct();
for cc=1:length(cfnames_exog)
    Params_cf=Params;
    Params_cf.ptypemasses=Params.ptypemasses(:)'.*cohortmass_cf;
    Params_cf.ptypemasses=Params_cf.ptypemasses/sum(Params_cf.ptypemasses);
    Params_cf.mewj=Params.mewj(cohortentryage(c_cf)-Params.agejshifter)*ones(1,N_j);
    if strcmp(cfnames_exog{cc},'cfloor50')
        Params_cf.cfloor=0.5*Params.cfloor;
    end
    [~,Policy_cf]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_cf,DiscountFactorParamNames,vfoptions_cf);
    StationaryDist_cf=StationaryDist_Case1_FHorz_PType(jequaloneDist_cf,AgeWeightsParamNames,PTypeDistParamNames,Policy_cf,n_d,n_a,n_z,N_j,Names_i,pi_z_J_noattrition,Params_cf,simoptions_cf);
    ACS_cf=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_cf,Policy_cf,FnsToEvaluate_gmm,Params_cf,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cf);
    CF_exog.(cfnames_exog{cc})=zeros(5,N_j);
    for q=1:5
        CF_exog.(cfnames_exog{cc})(q,:)=gather(ACS_cf.(['aliveq',num2str(q)]).assets.Median);
    end
    fprintf('   exogenous %-9s median assets at ages 74/90/100 by quintile (thousands): %s\n',cfnames_exog{cc},mat2str(round(CF_exog.(cfnames_exog{cc})(:,[74,90,100]-Params.agejshifter)/1000)))
end
clear Policy_cf StationaryDist_cf ACS_cf

%% The endogenous model (overwrites the exogenous model's shock process, grids, ReturnFn, V and Policy; solves at the published parameters)
DFJ2010_EndogMedex1
% at the estimates
Est_EM=load(['./SavedOutput/DFJ2010_EndogMedex_Estimation_AC',EMfigdefinition,'.mat'],'EstimParams_EM','EstimParamNames_EM');
for pp=1:length(Est_EM.EstimParamNames_EM)
    Params.(Est_EM.EstimParamNames_EM{pp})=Est_EM.EstimParams_EM.(Est_EM.EstimParamNames_EM{pp});
end
fprintf('\nEndogenous model at the estimates with the autocorrelations in the %s definition: nu %.4f, omega %.4f, beta %.4f, eMin %.1f (published: 2.1461, 3.1910, 0.9905, 202.1)\n',EMfigdefinition,Params.nu,Params.omega,Params.beta,Params.eMin)
tic;
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
fprintf('Value fn at the estimates took %.1f seconds\n',toc)

%% Initial distribution and cohorts on the endogenous shock grids, the PI-half restrictions, the asset and medical-expense data (as doPart(6))
DFJ2010_InitialDistCohorts
Params.half.femaleq1=1; Params.half.femaleq2=1; Params.half.femaleq3=2; Params.half.femaleq4=2; Params.half.femaleq5=2;
Params.half.maleq1=1;   Params.half.maleq2=1;   Params.half.maleq3=2;   Params.half.maleq4=2;   Params.half.maleq5=2;
simoptions.conditionalrestrictions.aliveh1=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==1);
simoptions.conditionalrestrictions.aliveh2=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==2);
DFJ2010_AssetTargets
MedexCells=DataWork.Targets.MedexCells;
cohortj0=cohortentryage-Params.agejshifter;

%% Figures 12 and A4 at the estimates: each cohort from its 1996 distribution (as doPart(6), without the autocorrelations)
% (assets, and out-of-pocket medical expenses floored at 125 per period: the same functions as doPart(6) and the estimation)
FnsToEvaluate_fig.assets=@(aprime,a,h,zeta,xi) a;
FnsToEvaluate_fig.oopmedex=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    max(DFJ2010EM_OOPMedex(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6),125); % floored at 125 per period: DFJ2010 bottom-code the TWO-YEAR average at minmedex=250 (gmmcrit51.g); the toolkit averages the raw function, so a per-period floor of half that is used (equal when both years are above 250; see the write-up)
simoptions_fig=simoptions;
simoptions_fig.whichstats=[1,1,0,0,0,1,0]; % mean, median, quantile cutoffs
simoptions_fig.nquantiles=10; % QuantileCutoffs has nquantiles+1 rows incl. min and max, so row 10 is the 90th percentile
simoptions_fig.agemass_withCohort=1;
simoptions_fig.ntimeperiods=[1,2]; % the single-period statistics (asset medians) and the 2-period averages (the medical-expense data are 2-year averages)
simoptions_fig.pi_z=pi_z_J;
simoptions_fig.windowmethod='tdigest1'; % as the estimation
fnfig=fieldnames(FnsToEvaluate_fig); rnfig=fieldnames(simoptions_fig.conditionalrestrictions);
simoptions_fig.whichcombos=zeros(2*length(fnfig),N_j,1+length(rnfig)); % rows: the M=1 block then the M=2 block; pages: unrestricted then the restrictions
simoptions_fig.whichcombos(1:length(fnfig),:,:)=1; % every single-period statistic
iowf=length(fnfig)+find(strcmp(fnfig,'oopmedex')); pageswf=1+[find(strcmp(rnfig,'aliveh1')),find(strcmp(rnfig,'aliveh2'))]; % the 2-period average of oopmedex on the PI-half pages, set per cohort below

years=1996:2:2006;
modelmedian=nan(Jdata,1);
ModelMedex=struct(); ModelMedex.Mean=nan(5,2,6); ModelMedex.P90=nan(5,2,6);
Params_c=Params;
tic;
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
    simoptions_fig.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
    starts=cohortj0(c)+[1,3,5,7,9]; starts=starts(starts<=N_j-1); % the window start ages of this cohort's waves 1998..2006
    simoptions_fig.whichcombos(iowf,:,pageswf)=0; simoptions_fig.whichcombos(iowf,starts,pageswf)=1;
    StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_fig);
    ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy,FnsToEvaluate_fig,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_fig);
    rows=find(MomentList.cohort==c);
    for rr=rows'
        modelmedian(rr)=gather(ACS_c.(['aliveq',num2str(MomentList.PIq(rr))]).assets.Median(cohortentryage(c)+MomentList.year(rr)-1996-Params.agejshifter));
    end
    % a wave in year y covers years y-1 and y: the 2-period average starting at the cohort's age in year y-1
    for hh=1:2
        hstr=['aliveh',num2str(hh)];
        for y=2:6
            s=cohortj0(c)+(years(y)-1996)-1;
            if s<=N_j-1
                ModelMedex.Mean(c,hh,y)=gather(ACS_c.(hstr).oopmedex_2periodAverage.Mean(s));
                ModelMedex.P90(c,hh,y)=gather(ACS_c.(hstr).oopmedex_2periodAverage.QuantileCutoffs(10,s));
            end
        end
    end
end
clear StationaryDist_c ACS_c
fprintf('Five cohorts at the estimates (distributions and life-cycle profiles) took %.1f seconds\n',toc)

% Assets: the asset-median block of the objective at the estimates should equal the one the estimation reports (doPart(7): 151.4,
% doPart(8): 152.4), a check that this is the same model at the same parameters
gap=modelmedian(tkorder)-MomentList.median(tkorder);
if strcmp(EMfigdefinition,'DFJ')
    estassetblock=152.4; % DFJ2010_EndogMedex_Estimation_ACDFJ, 2026-10-31
else
    estassetblock=151.4; % DFJ2010_EndogMedex_Estimation_ACtoolkit, 2026-10-28
end
fprintf('Asset-median block of the objective at the estimates (diagonal W): %.1f over %d moments (the %s-definition estimation reports %.1f)\n',gap'*WeightingMatrix*gap,Jdata,EMfigdefinition,estassetblock)
okm=~isnan(MedexCells.Mean(:,:,2:6)) & ~isnan(ModelMedex.Mean(:,:,2:6));
dm=MedexCells.Mean(:,:,2:6); mm=ModelMedex.Mean(:,:,2:6); dp=MedexCells.P90(:,:,2:6); mp=ModelMedex.P90(:,:,2:6);
fprintf('Medical expenses over the %d cells with data: mean data %.0f vs model %.0f; 90th pct data %.0f vs model %.0f; corr of cell means %.3f (published parameters: model mean 4362, 90th pct 9793, corr 0.853)\n',sum(okm(:)),mean(dm(okm)),mean(mm(okm)),mean(dp(okm)),mean(mp(okm)),corr(dm(okm),mm(okm)))

% DFJ2010 Figure 12 at the estimates: the cohorts aged 74 and 84 in 1996, data (solid) and the endogenous model (dashed),
% 1996-2006; both lines start from the 1996 data median, the initial distribution
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
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_EM_Fig12_assets_estimates.png')

% DFJ2010 Figure A4 at the estimates (data solid, model dashed)
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
    title(['Mean OOP medex, PI half ',num2str(hh)]); xlabel('age'); ylabel('thousands of 1998 dollars')
    subplot(2,2,2+hh)
    hold on
    for c=1:5
        ages=cohortentryage(c)+years(2:6)-1996;
        plot(ages,squeeze(MedexCells.P90(c,hh,2:6))/1000,'k-',ages,squeeze(ModelMedex.P90(c,hh,2:6))/1000,'b--')
    end
    hold off
    title(['90th pct of OOP medex, PI half ',num2str(hh)]); xlabel('age'); ylabel('thousands of 1998 dollars')
end
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_EM_FigA4_medex_estimates.png')

%% Figures 13 and 14 panel B: the endogenous model's counterfactuals at the estimates, cohort aged 74 in 1996, no attrition (as doPart(5))
pi_z_J_noattrition=pi_z_J; % rebuilt on the endogenous model's zeta transitions
for ii=1:length(Names_i)
    pih=pi_h_J.(Names_i{ii});
    for jj=1:N_j
        for hh=1:2
            pih(hh,:,jj)=[pih(hh,1:2,jj)/sum(pih(hh,1:2,jj)),0,0];
        end
        pi_z_J_noattrition.(Names_i{ii})(:,:,jj)=kron(pi_zeta,pih(:,:,jj));
    end
end
jequaloneDist_cf=struct();
cohortmass_cf=zeros(1,length(Names_i));
for ii=1:length(Names_i)
    temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c_cf)-Params.agejshifter);
    cohortmass_cf(ii)=gather(sum(temp(:)));
    if cohortmass_cf(ii)>0
        jequaloneDist_cf.(Names_i{ii})=temp/cohortmass_cf(ii);
    else
        jequaloneDist_cf.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp); % placeholder, gets weight zero
    end
end
simoptions_cf=simoptions_gmm; % aliveq1..5 restrictions, median only (rebuilt by DFJ2010_AssetTargets above)
simoptions_cf.jequaloneDistAge=cohortentryage(c_cf)-Params.agejshifter;
simoptions_cf.agemass_withCohort=1;
vfoptions_cf=vfoptions;
vfoptions_cf.verbose=0;
vfoptions_cf.verboseparams=0;
cfnames_EM={'baseline','nomedex','eMin50'};
CF_EM=struct();
tic;
for cc=1:length(cfnames_EM)
    Params_cf=Params;
    Params_cf.ptypemasses=Params.ptypemasses(:)'.*cohortmass_cf;
    Params_cf.ptypemasses=Params_cf.ptypemasses/sum(Params_cf.ptypemasses);
    Params_cf.mewj=Params.mewj(cohortentryage(c_cf)-Params.agejshifter)*ones(1,N_j);
    switch cfnames_EM{cc}
        case 'nomedex' % the medical-needs shifter mu set to zero: ln mu=-100 and no shocks, so medical goods (and out-of-pocket expenses) are zero
            for ii=1:length(Names_i)
                Params_cf.mu_coeff_healthgood.(Names_i{ii})=-100*ones(size(Params.mu_coeff_healthgood.(Names_i{ii})));
                Params_cf.mu_coeff_healthbad.(Names_i{ii})=-100*ones(size(Params.mu_coeff_healthbad.(Names_i{ii})));
                Params_cf.s_coeff_healthgood.(Names_i{ii})=zeros(size(Params.s_coeff_healthgood.(Names_i{ii})));
                Params_cf.s_coeff_healthbad.(Names_i{ii})=zeros(size(Params.s_coeff_healthbad.(Names_i{ii})));
            end
        case 'eMin50' % the utility floor indexed by half the consumption level that provides it when mu=1 (DFJ2010 footnote 8)
            Params_cf.eMin=0.5*Params.eMin;
    end
    [~,Policy_cf]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_cf,DiscountFactorParamNames,vfoptions_cf);
    StationaryDist_cf=StationaryDist_Case1_FHorz_PType(jequaloneDist_cf,AgeWeightsParamNames,PTypeDistParamNames,Policy_cf,n_d,n_a,n_z,N_j,Names_i,pi_z_J_noattrition,Params_cf,simoptions_cf);
    ACS_cf=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_cf,Policy_cf,FnsToEvaluate_gmm,Params_cf,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cf);
    CF_EM.(cfnames_EM{cc})=zeros(5,N_j);
    for q=1:5
        CF_EM.(cfnames_EM{cc})(q,:)=gather(ACS_cf.(['aliveq',num2str(q)]).assets.Median);
    end
    fprintf('   endogenous %-9s median assets at ages 74/90/100 by quintile (thousands): %s (%.0f s so far)\n',cfnames_EM{cc},mat2str(round(CF_EM.(cfnames_EM{cc})(:,[74,90,100]-Params.agejshifter)/1000)),toc)
end
clear Policy_cf StationaryDist_cf ACS_cf

ages_cf=74:100;
jj_cf=ages_cf-Params.agejshifter;
% DFJ2010 Figure 13: baseline (dashed) vs no medical expenses (solid), endogenous model
figure_c=figure_c+1;
figure(figure_c);
hold on
for q=1:5
    plot(ages_cf,CF_EM.baseline(q,jj_cf)/1000,'--','Color',[0 0 0.4+0.12*q],'LineWidth',1)
    plot(ages_cf,CF_EM.nomedex(q,jj_cf)/1000,'-','Color',[0 0 0.4+0.12*q],'LineWidth',1)
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Endogenous medical expenses: baseline (dashed) vs no medical expenses (solid)')
legend({'q1 baseline','q1','q2 baseline','q2','q3 baseline','q3','q4 baseline','q4','q5 baseline','q5'},'Location','northeast','NumColumns',1)
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_EM_Fig13_nomedex.png')
% DFJ2010 Figure 14: baseline (dashed) vs 50 percent of the floor (solid); A exogenous, B endogenous. One colour per quintile (shared by
% the two panels, legend in panel B), a common y-axis, the y-label on panel A only, and the line styles explained in an overall title.
figure_c=figure_c+1;
figure(figure_c);
set(figure(figure_c),'Position',[100,100,1200,500]);
CFpanels={CF_exog.baseline,CF_exog.cfloor50; CF_EM.baseline,CF_EM.eMin50};
paneltitles={'A: exogenous medical expenses','B: endogenous medical expenses'};
qcolors=lines(5);
ymax14=20*ceil(max([CFpanels{1,1}(:,jj_cf);CFpanels{1,2}(:,jj_cf);CFpanels{2,1}(:,jj_cf);CFpanels{2,2}(:,jj_cf)],[],'all')/20000);
hq=gobjects(1,5);
for pp=1:2
    subplot(1,2,pp)
    hold on
    for q=1:5
        plot(ages_cf,CFpanels{pp,1}(q,jj_cf)/1000,'--','Color',qcolors(q,:),'LineWidth',1.2)
        hq(q)=plot(ages_cf,CFpanels{pp,2}(q,jj_cf)/1000,'-','Color',qcolors(q,:),'LineWidth',1.2);
    end
    hold off
    xlim([ages_cf(1),ages_cf(end)]); ylim([0,ymax14]); box on
    xlabel('Age'); title(paneltitles{pp})
    if pp==1
        ylabel('Median assets (thousands of 1998 dollars)')
    else
        legend(hq,{'PI quintile 1','PI quintile 2','PI quintile 3','PI quintile 4','PI quintile 5'},'Location','northeast')
    end
end
sgtitle('Median assets, cohort aged 74 in 1996: baseline (dashed) and 50 percent of the floor (solid)')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig14_cfloor50.png')
fprintf('\nDFJ2010 Sec. IX.C: eliminating medical expenses has effects similar to the exogenous model (Figure 9); halving the floor raises the savings of all quintiles, by less in the endogenous model\n')
fprintf('   here, top quintile at 90: endogenous baseline %.0fk, no medex %.0fk, half floor %.0fk; exogenous baseline %.0fk, half floor %.0fk\n',CF_EM.baseline(5,90-Params.agejshifter)/1000,CF_EM.nomedex(5,90-Params.agejshifter)/1000,CF_EM.eMin50(5,90-Params.agejshifter)/1000,CF_exog.baseline(5,90-Params.agejshifter)/1000,CF_exog.cfloor50(5,90-Params.agejshifter)/1000)
fprintf('   bottom quintile at 90: endogenous baseline %.0fk, half floor %.0fk; exogenous baseline %.0fk, half floor %.0fk\n',CF_EM.baseline(1,90-Params.agejshifter)/1000,CF_EM.eMin50(1,90-Params.agejshifter)/1000,CF_exog.baseline(1,90-Params.agejshifter)/1000,CF_exog.cfloor50(1,90-Params.agejshifter)/1000)
figtime=toc(figtimer);
fprintf('doPart(10), the endogenous-model figures at the estimates, took %.1f minutes\n',figtime/60)

save ./SavedOutput/DFJ2010_EndogMedex_Figures.mat EMfigdefinition Est_EM modelmedian MomentList tkorder ModelMedex MedexCells cohortentryage CF_exog CF_EM cfnames_exog cfnames_EM figtime
