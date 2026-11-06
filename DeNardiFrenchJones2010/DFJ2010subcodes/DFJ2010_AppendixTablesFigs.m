% DFJ2010_AppendixTablesFigs: DFJ2010 Tables 1 and 2, Figure 8 (both panels), and appendix Figures A1 and A2. Called from
% DeNardiFrenchJones2010.m as doPart(11), on its own (like doPart(10)): it starts with the EXOGENOUS model in the workspace and
% then sets up the endogenous model (DFJ2010_EndogMedex1), which overwrites it.
%   Table 1:   the medical-expense shock process (persistence and innovation variances of the persistent and transitory
%              parts), directly from Params: the exogenous model's values (DFJ2010_FirstStageParams) and the endogenous model's
%              (DFJ2010_EndogMedex1, DFJ2010's mxvec)
%   Table 2:   life expectancy at 70 by PI quintile, sex and health, from the model's health-and-survival chain pi_h_J (no solve)
%   Figure A1: median assets by cohort and PI quintile in the data (the target cells, DFJ2010_AssetTargets)
%   Figure 8:  the CDF of assets held one period before death, data and model: A, the bequest model (Table 3 column 3, the one
%              DFJ2010 focus on); B, the benchmark (column 1). Both at DFJ2010's published values, as doPart(3) draws panel B.
%   Figure A2: the CDF of out-of-pocket medical expenses, data and model: A, the exogenous model (column 1), out-of-pocket
%              expenses net of government transfers, max(m-b,0); B, the endogenous model at the estimates of doPart(8) (or 7,
%              EMfigdefinition), whose out-of-pocket expenses are already net of the transfer. The model's 2-year averages (as
%              the data), pooled over the cohort x PI half x wave cells of the medical-expense moments with the data's cell
%              counts as weights, so model and data have the same composition; both bottom-coded at $250 on the 2-year average, as DFJ2010
%              (the HRS extract floors the data; the model's pooled 2-year averages are floored here).
% Writes SavedOutput/LatexInputs/DFJ2010_Table1.tex and DFJ2010_Table2.tex (tabular environments the write-up inputs).
% Needs SavedOutput/DFJ2010_DataWork.mat with Targets.MedexObs (the person-wave medical expenses, added 2026-11-01).

if ~exist('EMfigdefinition','var')
    EMfigdefinition='DFJ'; % set by DeNardiFrenchJones2010.m
end
apptimer=tic;
if ~exist('./SavedOutput/LatexInputs','dir')
    mkdir('./SavedOutput/LatexInputs')
end
years=1996:2:2006;
cellmin=10; % the medical-expense cells with fewer observations are not moments (DFJ2010_MedexTargets)
qcolors=lines(5);
quintilenames={'Bottom','Second','Third','Fourth','Top'};

%% DFJ2010 Table 1: the medical-expense shock process, directly from Params (exogenous model now; the endogenous model's column below)
Table1=struct();
Table1.rownames={'rho_m','sigma_mepsilon^2','sigma_mxi^2'};
Table1.exog=[Params.rho_m; Params.sigma_mepsilon^2; Params.sigma_mxi^2];
Table1.dfj=[0.922; 0.050; 0.665]; Table1.dfjse=[0.010; 0.008; 0.014]; % DFJ2010 Table 1
fprintf('\nDFJ2010 Table 1 (exogenous model, from Params): rho_m %.3f, var of persistent innovations %.3f, var of transitory part %.3f (DFJ2010: .922, .050, .665); persistent share of the variance %.3f (DFJ2010 text: .335)\n',Table1.exog(1),Table1.exog(2),Table1.exog(3),Table1.exog(2)/(1-Table1.exog(1)^2)/Params.normalizepsi^2)

%% DFJ2010 Table 2: life expectancy at 70 from the model's health-and-survival chain (no solve)
% pi_h_J.(ptype)(:,:,jj) is the model's transition from period jj to jj+1 over health states (good, bad, death, dead); the
% period of age 70 is jj=1. Life expectancy at 70 in health h is the expected number of further periods alive (good or bad
% health), sum over t>=1 of the probability of being alive at 70+t, plus one half for the year of death (deaths are spread
% over the year). The aggregates weight by the 1996 sample's composition (DFJ2010: 'the gender and health distributions
% observed in each permanent income quintile', and so on): the ptype masses and each ptype's health shares, all 1996 ages.
LE=zeros(5,2,2); % PI quintile x health at 70 (1 good, 2 bad) x sex (1 female, 2 male)
Wcomp=zeros(5,2,2); % the same cells' 1996 sample weights
for ii=1:length(Names_i)
    q=Params.quintile.(Names_i{ii});
    g=1+strncmp(Names_i{ii},'male',4);
    P=gather(pi_h_J.(Names_i{ii}));
    for h0=1:2
        dist=zeros(1,n_z(1)); dist(h0)=1;
        le=0.5;
        for jj=1:N_j-1
            dist=dist*P(:,:,jj);
            le=le+sum(dist(1:2));
        end
        LE(q,h0,g)=le;
    end
    D96=DataWork.InitialDist.InitDist.(Names_i{ii}); % n_a x 2 (health) x 33 (1996 age)
    hshare=squeeze(sum(sum(D96,1),3)); hshare=hshare(:)'/sum(hshare);
    Wcomp(q,:,g)=Params.ptypemasses(ii)*hshare;
end
LEall=sum(sum(LE.*Wcomp,3),2)./sum(sum(Wcomp,3),2); % by quintile
LEsex=squeeze(sum(sum(LE.*Wcomp,1),2)./sum(sum(Wcomp,1),2)); % female, male
LEhealth=squeeze(sum(sum(LE.*Wcomp,1),3)./sum(sum(Wcomp,1),3)); % good, bad
Table2dfj=[7.6,5.9,12.8,10.9,11.1; 8.4,6.6,13.8,12.0,12.4; 9.3,7.4,14.7,13.2,13.1; 10.5,8.4,15.7,14.2,14.4; 11.3,9.3,16.7,15.1,14.7];
Table2dfj_sex=[14.3,9.7]; Table2dfj_health=[14.4,11.6]; % women, men; healthy, unhealthy
Table2=[squeeze(LE(:,1,2)),squeeze(LE(:,2,2)),squeeze(LE(:,1,1)),squeeze(LE(:,2,1)),LEall]; % columns as DFJ2010: healthy male, unhealthy male, healthy female, unhealthy female, all
fprintf('\nDFJ2010 Table 2: life expectancy at 70 [DFJ2010 in brackets], rows PI quintile, columns healthy male, unhealthy male, healthy female, unhealthy female, all\n')
for q=1:5
    fprintf('   %-7s',quintilenames{q}); fprintf('  %5.1f [%4.1f]',[Table2(q,:);Table2dfj(q,:)]); fprintf('\n')
end
fprintf('   men %.1f [9.7], women %.1f [14.3]; unhealthy %.1f [11.6], healthy %.1f [14.4]\n',LEsex(2),LEsex(1),LEhealth(2),LEhealth(1))
fid=fopen('./SavedOutput/LatexInputs/DFJ2010_Table2.tex','w');
fprintf(fid,'\\begin{tabular}{lccccc}\n\\toprule\n');
fprintf(fid,'PI quintile & Healthy male & Unhealthy male & Healthy female & Unhealthy female & All \\\\\n\\midrule\n');
for q=1:5
    fprintf(fid,'%s',quintilenames{q}); fprintf(fid,' & %.1f',Table2(q,:)); fprintf(fid,' \\\\\n'); % replication only: the write-up shows DFJ2010's Table 2 above it
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'Men & %.1f & & & & \\\\\nWomen & %.1f & & & & \\\\\n',LEsex(2),LEsex(1));
fprintf(fid,'Unhealthy & %.1f & & & & \\\\\nHealthy & %.1f & & & & \\\\\n',LEhealth(2),LEhealth(1));
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% DFJ2010 appendix Figure A1: median assets by cohort and PI quintile, data (the target cells, each line one cohort)
figure_c=figure_c+1;
figure(figure_c);
set(figure(figure_c),'Position',[100,100,900,550]);
hold on
hq=gobjects(1,5);
for c=1:ncohorts
    for q=1:5
        rows=find(MomentList.cohort==c & MomentList.PIq==q);
        hq(q)=plot(MomentList.avgage(rows),MomentList.median(rows)/1000,'-o','Color',qcolors(q,:),'MarkerSize',3,'LineWidth',1.2);
    end
end
hold off
box on
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Median assets by cohort and PI quintile, data (each line one cohort, 1998-2006)')
legend(hq,{'PI quintile 1','PI quintile 2','PI quintile 3','PI quintile 4','PI quintile 5'},'Location','northwest')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_FigA1_data.png')

%% The data for Figures 8 and A2
AssetsBeforeDeath=DataWork.Targets.AssetsBeforeDeath;
MedexObs=DataWork.Targets.MedexObs;
MedexCnt=DataWork.Targets.MedexCells.Cnt; % data observations per cohort x PI half x year cell
dataA2=[];
for c=1:ncohorts
    for hh=1:2
        for y=2:6
            if MedexCnt(c,hh,y)>=cellmin
                ind=(MedexObs.cohort==c) & (MedexObs.half==hh) & MedexObs.obs(:,y);
                dataA2=[dataA2; MedexObs.M(ind,y)]; %#ok<AGROW>
            end
        end
    end
end
xgA2=linspace(0,60000,601); % the CDFs of Figure A2 are evaluated on this grid (dollars)
FdataA2=mean(dataA2(:)'<=xgA2(:),2)';
fprintf('\nFigure A2 data: %d person-waves in the medical-expense cells; quartiles %s, 90th percentile %.0f\n',numel(dataA2),mat2str(round(quantile(dataA2,[0.25,0.5,0.75]))),quantile(dataA2,0.9))

%% Figure 8 (panels A and B) and Figure A2 panel A: the exogenous model at DFJ2010 Table 3 columns 3 and 1
% PI halves for the medical-expense cells (as DFJ2010_EndogMedex2)
Params.half.femaleq1=1; Params.half.femaleq2=1; Params.half.femaleq3=2; Params.half.femaleq4=2; Params.half.femaleq5=2;
Params.half.maleq1=1;   Params.half.maleq2=1;   Params.half.maleq3=2;   Params.half.maleq4=2;   Params.half.maleq5=2;
FnsToEvaluate_A2x.netmedex=@(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi,cfloor,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6) ...
    max(DeNardiFrenchJones2010_MedicalExpenses(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi) ...
    -DeNardiFrenchJones2010_GovTransfers(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi,cfloor,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6),0);
simoptions_A2=simoptions;
simoptions_A2.conditionalrestrictions=struct();
simoptions_A2.conditionalrestrictions.aliveh1=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==1);
simoptions_A2.conditionalrestrictions.aliveh2=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==2);
simoptions_A2.whichstats=[1,0,0,0,0,1,0]; % mean and quantile cutoffs
simoptions_A2.nquantiles=100; % percentiles: QuantileCutoffs has 101 rows (min, 1st..99th percentile, max)
simoptions_A2.agemass_withCohort=1;
simoptions_A2.ntimeperiods=[1,2]; % the 2-year averages (as the data); the single-period block is cheap and kept as in doPart(6)
simoptions_A2.pi_z=pi_z_J;
simoptions_A2.windowmethod='tdigest1';
simoptions_A2.whichcombos=zeros(2,N_j,3); % rows: the M=1 then the M=2 block of the one function; pages: unrestricted, aliveh1, aliveh2
simoptions_A2.whichcombos(1,:,:)=1;
vfoptions_8=vfoptions;
vfoptions_8.verbose=0;
vfoptions_8.verboseparams=0;
simoptions_8=simoptions_gmm;
simoptions_8.agemass_withCohort=1;
cohortj0=cohortentryage-Params.agejshifter;
fig8cols=[3,1]; % panel A: the bequest model DFJ2010 focus on (Table 3 column 3); panel B: the benchmark (column 1)
deathdist8=zeros(n_a,2);
ptsA2=cell(1,2); wA2=cell(1,2); % the model's pooled 2-year medical expenses (percentile midpoints) and their weights, by panel of Figure A2
for pp=1:2
    col=fig8cols(pp);
    Params_8=Params;
    Params_8.upsilon=Table3.upsilon(col);
    Params_8.beta=Table3.beta(col);
    Params_8.delta=Table3.delta(col);
    Params_8.k=Table3.k(col);
    Params_8.cfloor=exp(Table3.logcfloor(col));
    if Table3.phi0(col)==1
        Params_8.theta=0;
    else
        Params_8.theta=(((1/0.98)*(1-Table3.phi0(col))/Table3.phi0(col))^Params_8.upsilon)/((1/0.98)*Params_8.beta); % as the main script
    end
    [~,Policy_8]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params_8,DiscountFactorParamNames,vfoptions_8);
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
        Params_c=Params_8;
        Params_c.ptypemasses=Params.ptypemasses(:)'.*cohortmass_c;
        Params_c.ptypemasses=Params_c.ptypemasses/sum(Params_c.ptypemasses);
        Params_c.mewj=Params.mewj(cohortentryage(c)-Params.agejshifter)*ones(1,N_j);
        simoptions_8.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
        StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy_8,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_8);
        % Figure 8 (as doPart(3)): the assets of those in the death state (h index 3), which are the assets held one period
        % before death, from the year after entry to 2006, cohorts pooled with their 1996 shares
        jentry=cohortentryage(c)-Params.agejshifter;
        jend=min(jentry+10,N_j);
        cohortshare_c=sum(Params.ptypemasses(:)'.*cohortmass_c);
        Dh_all=zeros(n_a,n_z(1),N_j);
        for ii=1:length(Names_i)
            Dh=gather(reshape(sum(reshape(StationaryDist_c.(Names_i{ii}),[n_a,n_z(1),n_z(2)*n_e,N_j]),3),[n_a,n_z(1),N_j]));
            Dh_all=Dh_all+Params_c.ptypemasses(ii)*Dh;
        end
        deathdist8(:,pp)=deathdist8(:,pp)+cohortshare_c*sum(Dh_all(:,3,jentry+1:jend),3);
        if col==1 % Figure A2 panel A: the benchmark's out-of-pocket medical expenses net of transfers, 2-year averages
            simoptions_A2.jequaloneDistAge=simoptions_8.jequaloneDistAge;
            starts=cohortj0(c)+[1,3,5,7,9]; starts=starts(starts<=N_j-1); % the window start ages of the cohort's waves 1998..2006
            simoptions_A2.whichcombos(2,:,2:3)=0; simoptions_A2.whichcombos(2,starts,2:3)=1;
            ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy_8,FnsToEvaluate_A2x,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_A2);
            for hh=1:2
                for y=2:6
                    s=cohortj0(c)+(years(y)-1996)-1; % a wave in year y covers years y-1 and y
                    if s<=N_j-1 && MedexCnt(c,hh,y)>=cellmin
                        Qc=gather(ACS_c.(['aliveh',num2str(hh)]).netmedex_2periodAverage.QuantileCutoffs(:,s));
                        ptsA2{1}=[ptsA2{1}; (Qc(1:end-1)+Qc(2:end))/2]; % percentile midpoints, each 1/100 of the cell
                        wA2{1}=[wA2{1}; MedexCnt(c,hh,y)/100*ones(numel(Qc)-1,1)];
                    end
                end
            end
        end
    end
    cdf8=cumsum(deathdist8(:,pp))/sum(deathdist8(:,pp));
    fprintf('Figure 8 panel %s (Table 3 column %d): model quartiles of assets one period before death %s thousand, share at the zero grid point %.3f (data quartiles %s thousand, share a<=0 %.3f)\n',char('A'+pp-1),col,mat2str(round(a_grid(arrayfun(@(p) find(cdf8>=p,1),[0.25,0.5,0.75]))/1000)),cdf8(1),mat2str(round(quantile(AssetsBeforeDeath,[0.25,0.5,0.75])/1000)),mean(AssetsBeforeDeath<=0))
end
clear Policy_8 StationaryDist_c ACS_c

% DFJ2010 Figure 8: data (solid) and model (dashed); A with bequest motives, B without
abd_sorted=sort(AssetsBeforeDeath);
figure_c=figure_c+1;
figure(figure_c);
set(figure(figure_c),'Position',[100,100,1100,450]);
paneltitles8={'A: model with bequest motives (Table 3 column 3)','B: model without bequest motives (Table 3 column 1)'};
for pp=1:2
    subplot(1,2,pp)
    plot(abd_sorted/1000,(1:numel(abd_sorted))/numel(abd_sorted),'k-','LineWidth',1.2)
    hold on
    plot(a_grid/1000,cumsum(deathdist8(:,pp))/sum(deathdist8(:,pp)),'b--','LineWidth',1.2)
    hold off
    box on
    xlim([0,500]); ylim([0,1]); xlabel('Assets (thousands of 1998 dollars)'); title(paneltitles8{pp})
    if pp==1
        ylabel('CDF')
    else
        legend({'data (HRS)','model'},'Location','southeast')
    end
end
sgtitle('Assets held one period before death')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_Fig8_AB_assetsbeforedeath.png')

%% Figure A2 panel B and Table 1's endogenous column: the endogenous model (overwrites the exogenous model) at the estimates
DFJ2010_EndogMedex1
Table1.endog=[Params.rho_m; Params.sigma_mepsilon^2; Params.sigma_mxi^2];
Est_EM=load(['./SavedOutput/DFJ2010_EndogMedex_Estimation_AC',EMfigdefinition,'.mat'],'EstimParams_EM','EstimParamNames_EM');
for pp=1:length(Est_EM.EstimParamNames_EM)
    Params.(Est_EM.EstimParamNames_EM{pp})=Est_EM.EstimParams_EM.(Est_EM.EstimParamNames_EM{pp});
end
fprintf('\nEndogenous model at the estimates with the autocorrelations in the %s definition: nu %.4f, omega %.4f, beta %.4f, eMin %.1f\n',EMfigdefinition,Params.nu,Params.omega,Params.beta,Params.eMin)
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
DFJ2010_InitialDistCohorts % on the endogenous model's shock grids
Params.half.femaleq1=1; Params.half.femaleq2=1; Params.half.femaleq3=2; Params.half.femaleq4=2; Params.half.femaleq5=2;
Params.half.maleq1=1;   Params.half.maleq2=1;   Params.half.maleq3=2;   Params.half.maleq4=2;   Params.half.maleq5=2;
cohortj0=cohortentryage-Params.agejshifter;
% out-of-pocket medical expenses, already net of the transfer (DFJ2010EM_OOPMedex), without doPart(6)'s floor at 125
FnsToEvaluate_A2e.oopmedex=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_OOPMedex(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6); % NOT floored (doPart(6) floors at 125 per period for the logs of the autocorrelation moments)
simoptions_A2e=simoptions;
simoptions_A2e.conditionalrestrictions=struct();
simoptions_A2e.conditionalrestrictions.aliveh1=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==1);
simoptions_A2e.conditionalrestrictions.aliveh2=@(aprime,a,h,zeta,xi,half) ((h==0) || (h==1)) && (half==2);
simoptions_A2e.whichstats=[1,0,0,0,0,1,0];
simoptions_A2e.nquantiles=100;
simoptions_A2e.agemass_withCohort=1;
simoptions_A2e.ntimeperiods=[1,2];
simoptions_A2e.pi_z=pi_z_J;
simoptions_A2e.windowmethod='tdigest1';
simoptions_A2e.whichcombos=zeros(2,N_j,3);
simoptions_A2e.whichcombos(1,:,:)=1;
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
    Params_c=Params;
    Params_c.ptypemasses=Params.ptypemasses(:)'.*cohortmass_c;
    Params_c.ptypemasses=Params_c.ptypemasses/sum(Params_c.ptypemasses);
    Params_c.mewj=Params.mewj(cohortentryage(c)-Params.agejshifter)*ones(1,N_j);
    simoptions_A2e.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
    starts=cohortj0(c)+[1,3,5,7,9]; starts=starts(starts<=N_j-1);
    simoptions_A2e.whichcombos(2,:,2:3)=0; simoptions_A2e.whichcombos(2,starts,2:3)=1;
    StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_A2e);
    ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy,FnsToEvaluate_A2e,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_A2e);
    for hh=1:2
        for y=2:6
            s=cohortj0(c)+(years(y)-1996)-1;
            if s<=N_j-1 && MedexCnt(c,hh,y)>=cellmin
                Qc=gather(ACS_c.(['aliveh',num2str(hh)]).oopmedex_2periodAverage.QuantileCutoffs(:,s));
                ptsA2{2}=[ptsA2{2}; (Qc(1:end-1)+Qc(2:end))/2];
                wA2{2}=[wA2{2}; MedexCnt(c,hh,y)/100*ones(numel(Qc)-1,1)];
            end
        end
    end
end
clear StationaryDist_c ACS_c

%% DFJ2010 Table 1 (written now that both models' values are known)
fprintf('\nDFJ2010 Table 1 from Params: exogenous model %s, endogenous model %s (rho_m, var persistent innovations, var transitory; DFJ2010: .922 .050 .665)\n',mat2str(round(Table1.exog',3)),mat2str(round(Table1.endog',3)))
Table1labels={'$\\rho_m$ & Autocorrelation, persistent part','$\\sigma^2_{\\epsilon}$ & Innovation variance, persistent part','$\\sigma^2_{\\xi}$ & Innovation variance, transitory part'};
fid=fopen('./SavedOutput/LatexInputs/DFJ2010_Table1.tex','w');
fprintf(fid,'\\begin{tabular}{llcc}\n\\toprule\n');
fprintf(fid,'Parameter & Variable & Exogenous model & Endogenous model \\\\\n\\midrule\n');
for rr=1:3 % replication only: the write-up shows DFJ2010's Table 1 above it
    fprintf(fid,[Table1labels{rr},' & %.3f & %.3f \\\\\n'],Table1.exog(rr),Table1.endog(rr));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% DFJ2010 appendix Figure A2: CDF of out-of-pocket medical expenses (2-year averages), data (solid) and model (dashed)
minmedex=250; % DFJ2010 bottom-code the 2-year average at 250 (gmmcrit51.g), and so does the HRS extract for the data
FmodelA2=zeros(2,numel(xgA2));
for pp=1:2
    FmodelA2(pp,:)=sum(wA2{pp}.*(max(ptsA2{pp},minmedex)<=xgA2),1)/sum(wA2{pp});
end
figure_c=figure_c+1;
figure(figure_c);
set(figure(figure_c),'Position',[100,100,1100,450]);
paneltitlesA2={'A: exogenous medical expenses (net of transfers)','B: endogenous medical expenses (at the estimates)'};
for pp=1:2
    subplot(1,2,pp)
    plot(xgA2/1000,FdataA2,'k-','LineWidth',1.2)
    hold on
    plot(xgA2/1000,FmodelA2(pp,:),'b--','LineWidth',1.2)
    hold off
    box on
    xlim([0,xgA2(end)/1000]); ylim([0,1]); xlabel('Out-of-pocket medical expenses (thousands of 1998 dollars)'); title(paneltitlesA2{pp})
    if pp==1
        ylabel('CDF')
    else
        legend({'data (HRS)','model'},'Location','southeast')
    end
end
sgtitle('Out-of-pocket medical expenses, 2-year averages, 1998-2006 waves')
saveas(figure(figure_c),'./SavedOutput/Graphs/DFJ2010_FigA2_medexCDF.png')
for pp=1:2
    fprintf('Figure A2 panel %s: model CDF at 1k/5k/10k/20k/40k %s, data %s\n',char('A'+pp-1),mat2str(round(interp1(xgA2,FmodelA2(pp,:),[1000,5000,10000,20000,40000]),3)),mat2str(round(interp1(xgA2,FdataA2,[1000,5000,10000,20000,40000]),3)))
end

apptime=toc(apptimer);
fprintf('doPart(11), Tables 1 and 2 and Figures 8, A1 and A2, took %.1f minutes\n',apptime/60)
save ./SavedOutput/DFJ2010_AppendixTablesFigs.mat Table1 LE Wcomp Table2 LEall LEsex LEhealth deathdist8 fig8cols ptsA2 wA2 xgA2 FdataA2 FmodelA2 EMfigdefinition apptime
