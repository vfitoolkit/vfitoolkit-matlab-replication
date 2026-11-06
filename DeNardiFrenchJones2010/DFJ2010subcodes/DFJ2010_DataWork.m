% DFJ2010_DataWork: build every data object the De Nardi, French & Jones (2010) replication needs,
% from the public RAND HRS Longitudinal File.
%
% Steps
%   1. DFJ2010_HRSextract: AHEAD singles panel 1994-2006 (rebuilds wealth_singles/dataprep2.dta)
%   2. DFJ-style arrays (GAUSS getdata mirror): permanent-income quintiles, cohorts, top-coding
%   3. DFJ2010_Targets: median assets by cohort x PI quintile x year, the 101 moment conditions,
%      their covariance matrix and weighting matrices; Figures 1 and A1 of the paper
%   4. DFJ2010_InitialDist: the 1996 joint distribution of (age, sex, PI, health, assets) on the
%      model's grids, as jequaloneDist-style arrays by permanent type
%   5. If DFJ's own wealth_singles/data/wlthmat13b.out is present, run steps 2-3 on it too and
%      print the two sets of targets side by side
% Output: SavedOutput/DFJ2010_HRSpanel.mat, SavedOutput/DFJ2010_DataWork.mat, figures in
% SavedOutput/Graphs/, and the diary DFJ2010_DataWork_Diary.txt.
%
% Lives in DFJ2010subcodes/ (moved 2026-11-01), but all its paths are relative to the DeNardiFrenchJones2010 folder: run it
% as doPart(1) of DeNardiFrenchJones2010.m, or on its own with the DeNardiFrenchJones2010 folder as the current folder
% (addpath('./DFJ2010subcodes'); DFJ2010_DataWork). Not with run('DFJ2010subcodes/DFJ2010_DataWork.m'), which changes the
% current folder to DFJ2010subcodes/ and so breaks the relative paths.
% Needs the RAND HRS Longitudinal File 2022 (V1) Stata file in ../../DataSets/HRS/HRS Longitudinal 2022/
% (see HRS_rebuild_plan.md). CPU only, no toolkit calls, runs in about a minute.

% Record the full run (incl. any error) to DFJ2010_DataWork_Diary.txt, fresh each run.
% (Not when called from DeNardiFrenchJones2010.m as its doPart(1): a nested diary would cut the main diary; the main
%  script sets calledfromDFJmain=1 before calling this script.)
if ~exist('calledfromDFJmain','var')
    diary off
    if exist('./DFJ2010_DataWork_Diary.txt','file'), delete('./DFJ2010_DataWork_Diary.txt'); end
    diary ./DFJ2010_DataWork_Diary.txt
end

close all
addpath('./DFJ2010subcodes');
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end

headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end

%% Model conventions that the data objects must match (keep in line with DeNardiFrenchJones2010.m)
DFJ2010_AssetGrid; % sets n_a and a_grid (shared subcode)
Names_i={'femaleq1','femaleq2','femaleq3','femaleq4','femaleq5','maleq1','maleq2','maleq3','maleq4','maleq5'};

%% 1. Extract the AHEAD singles panel from the RAND HRS file
DFJ2010_HRSextract
% leaves: panel (person-wave table), person (person table), waves, DFJyearlist

%% 2. DFJ-style arrays (GAUSS getdata(), with ageshft=2 so the first data column is 1996)
% GAUSS settings (singlsav15_wave8.gau): bornage=70, dieage=102, asstmax=2e6, momyr1=1996, momyr2=2006
bornage=70;
dieage=102;
asstmax=2000000;
cohortedges=[69,74,79,84,89,dieage]; % cohorts_a, defined on 1994 ages (=1996 age - ageshft)
years=[1996,1998,2000,2002,2004,2006];

N=height(person);
nW=numel(waves);
[~,loc]=ismember(panel.HHIDPN,person.HHIDPN);
idx=sub2ind([N,nW],loc,panel.RANDwave-waves(1)+1);
Aw=nan(N,nW);    Aw(idx)=panel.assets;
ALIVEw=false(N,nW); ALIVEw(idx)=panel.alive;
AGEw=nan(N,nW);  AGEw(idx)=panel.age;
Hw=nan(N,nW);    Hw(idx)=panel.badhealth;
Mw=nan(N,nW);    Mw(idx)=panel.medex;
mcols=find(ismember(DFJyearlist,years)); % the 1996..2006 columns

% DFJ's wlthmat13b.out holds only the people interviewed in 1996 (the start of the moment period;
% those who died between 1993 and 1995 are not in it), and getdata() then drops anyone whose 1996
% age is outside the model's ages. Mirror both, so N is the same denominator as GAUSS totobs.
age1996=person.age1996;
inmodel=ALIVEw(:,mcols(1)) & (age1996>=bornage & age1996<=dieage);
fprintf('\nStep 2: %d of %d persons are interviewed in 1996 with age in [%d,%d] and enter the moment calculations (DFJ: 2,688)\n',sum(inmodel),N,bornage,dieage)
person=person(inmodel,:);
Aw=Aw(inmodel,:); ALIVEw=ALIVEw(inmodel,:); AGEw=AGEw(inmodel,:); Hw=Hw(inmodel,:); Mw=Mw(inmodel,:);
age1996=age1996(inmodel);
N=height(person);

% top-code assets, observation indicator (alive and assets observed)
A_in=min(Aw(:,mcols),asstmax);
obs_in=ALIVEw(:,mcols) & ~isnan(Aw(:,mcols));
AGE_in=AGEw(:,mcols);
years_in=years;
% permanent-income quintile (PIquant2: quintile i is PI in ((i-1)/5, i/5]); PI is a percentile rank in (0,1]
PIq_in=ceil(5*person.PI-1e-9);
PIq_in(isnan(person.PI))=0;
% cohort (getchrt on 1996 ages minus ageshft=2, edges cohortedges; outside -> 0)
age1994=age1996-2;
cohort_in=zeros(N,1);
for c=1:5
    cohort_in(age1994>cohortedges(c) & age1994<=cohortedges(c+1))=c;
end
fprintf('   PI quintile counts: %s;  cohort counts (0=none): %s\n',mat2str(histcounts(PIq_in,0.5:1:5.5)),mat2str(histcounts(cohort_in,-0.5:1:5.5)))
fprintf('   cohorts (1996 ages): 1: 72-76, 2: 77-81, 3: 82-86, 4: 87-91, 5: 92-102 (paper Sec. VI.A)\n')
person.PIq=PIq_in;
person.cohort=cohort_in;

%% 3. Targets: medians, moment conditions, covariance and weighting matrices
DFJ2010_Targets
nY_pool=numel(years_in);
Targets.MedianAssets=MedA;          % 5 cohorts x 5 PI quintiles x 6 years (NaN where fewer than 10 obs)
Targets.CellCounts=Cnt;
Targets.AvgAge=AvgAge;
Targets.MedianAssets_survivors=MedA_surv;
Targets.MomentList=MomentList;      % one row per moment condition (cohort, PIq, year, avgage, median, n)
Targets.mmtmtx=mmtmtx;              % N x J individual moment functions
Targets.VCV=VCV;
Targets.W=W;
Targets.Wdiag=Wdiag;
Targets.VCV_medians=VCV_medians;    % bootstrap covariance of the cell medians (levels), for the toolkit's estimator
Targets.SE_medians=SE_medians;
Targets.years=years;
Targets.cohortedges=cohortedges;
Targets.N=N;
Targets.J=J;
% the same targets indexed by (PI quintile, age) for use with the model's age-conditional stats:
% each cohort-year cell is placed at the rounded average age of its members
Targets.ByAge=sortrows(table(MomentList.PIq,round(MomentList.avgage),MomentList.median,MomentList.se,MomentList.n,MomentList.cohort,MomentList.year,(1:J)',...
    'VariableNames',{'PIq','age','medianassets','se','n','cohort','year','momentindex'}),{'PIq','age'});
[~,ia]=unique(Targets.ByAge(:,{'PIq','age'}),'rows');
if numel(ia)<J
    warning('DFJ2010_DataWork: %d (PI quintile, age) pairs receive more than one cohort-year moment; the model script must resolve these',J-numel(ia))
end

fprintf('\nMedian assets by cohort (rows, 1996 age 74/79/84/89/94+) and year (cols 1996..2006), PI quintile 5 (top):\n')
disp(round(squeeze(MedA(:,5,:))))
fprintf('PI quintile 1 (bottom):\n')
disp(round(squeeze(MedA(:,1,:))))
fprintf('Moments per cohort: %s; per PI quintile: %s\n',mat2str(histcounts(MomentList.cohort,0.5:1:5.5)),mat2str(histcounts(MomentList.PIq,0.5:1:5.5)))

% Figure 1 of the paper: median assets by PI quintile for the cohorts aged 74 and 84 in 1996
fig1=figure(1);
hold on
for q=1:5 % the cohort aged 74 first, then 84, so the legend (filled column by column) has the blue 74s on the left and the red 84s on the right
    plot(squeeze(AvgAge(1,q,:)),squeeze(MedA(1,q,:))/1000,'-o','Color',[0 0 0.5+0.1*q],'LineWidth',1)
end
for q=1:5
    plot(squeeze(AvgAge(3,q,:)),squeeze(MedA(3,q,:))/1000,'--s','Color',[0.5+0.1*q 0 0],'LineWidth',1)
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Median assets by cohort and PI quintile: data (solid: aged 74 in 1996, dashed: aged 84)')
legend({'q1 (74)','q2 (74)','q3 (74)','q4 (74)','q5 (74)','q1 (84)','q2 (84)','q3 (84)','q4 (84)','q5 (84)'},'Location','northwest','NumColumns',2)
saveas(fig1,'./SavedOutput/Graphs/DFJ2010_DataWork_Figure1.png')

% Figure A1 of the appendix: all five cohorts
fig2=figure(2);
hold on
for c=1:5
    for q=1:5
        plot(squeeze(AvgAge(c,q,:)),squeeze(MedA(c,q,:))/1000,'-','Color',[0 0 0.4+0.12*q],'LineWidth',1)
    end
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Median assets by cohort and PI quintile, all five cohorts (unbalanced panel)')
saveas(fig2,'./SavedOutput/Graphs/DFJ2010_DataWork_FigureA1.png')

% Figure 2 of the paper: median assets by birth cohort, everyone alive (unbalanced panel) vs the survivors to 2006
% (balanced panel), all PI quintiles pooled. Paper: 74-year-olds in 1996 have median assets $60,000 (everyone) vs
% $84,000 (survivors to 2006), falling to $44,000 by 2006 for the survivors.
MedA_pool=nan(5,nY_pool); MedA_pool_surv=nan(5,nY_pool); AvgAge_pool=nan(5,nY_pool);
for c=1:5
    for y=1:nY_pool
        ind=(cohort_in==c) & obs_in(:,y);
        if sum(ind)>=quantmin
            MedA_pool(c,y)=median(A_in(ind,y));
            AvgAge_pool(c,y)=mean(AGE_in(ind,y));
        end
        ind=ind & survivor;
        if sum(ind)>=quantmin
            MedA_pool_surv(c,y)=median(A_in(ind,y));
        end
    end
end
fprintf('\nMedian assets by cohort, all PI quintiles pooled, everyone alive (rows 1-5) then survivors to 2006 (DFJ2010 Figure 2: 74-year-olds 60k vs 84k in 1996):\n')
disp(round(MedA_pool/1000)); disp(round(MedA_pool_surv/1000))
fig3=figure(3);
hold on
for c=1:5
    plot(AvgAge_pool(c,:),MedA_pool(c,:)/1000,'-','Color',[0 0 0.8],'LineWidth',1)
    plot(AvgAge_pool(c,:),MedA_pool_surv(c,:)/1000,'--','Color',[0.8 0 0],'LineWidth',1)
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Median assets by birth cohort: everyone alive (solid) vs survivors to 2006 (dashed)')
legend({'everyone alive','survivors to 2006'},'Location','northeast')
saveas(fig3,'./SavedOutput/Graphs/DFJ2010_DataWork_Figure2.png')
% the same for the top PI quintile only (not in the paper)
fig4=figure(4);
hold on
for c=1:5
    plot(squeeze(AvgAge(c,5,:)),squeeze(MedA(c,5,:))/1000,'-','Color',[0 0 0.8],'LineWidth',1)
    plot(squeeze(AvgAge(c,5,:)),squeeze(MedA_surv(c,5,:))/1000,'--','Color',[0.8 0 0],'LineWidth',1)
end
hold off
xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
title('Top PI quintile: everyone alive (solid) vs survivors to 2006 (dashed)')
saveas(fig4,'./SavedOutput/Graphs/DFJ2010_DataWork_FigureA2.png')

% Figure 8 of the paper (data part): assets held one period before death, i.e. at the last interview before the
% person died (paper footnote 5: reported estates are too badly imputed to use). Decedents are those in the 1996
% sample who die by 2006; their last interview is one of the 1996-2004 waves.
died=~isnan(person.deathyear) & person.deathyear<=2006;
AssetsBeforeDeath=nan(N,1);
for ii=find(died)'
    lastalive=find(ALIVEw(ii,:),1,'last');
    if lastalive>=mcols(1) && ~isnan(Aw(ii,lastalive))
        AssetsBeforeDeath(ii)=min(Aw(ii,lastalive),asstmax);
    end
end
AssetsBeforeDeath=AssetsBeforeDeath(~isnan(AssetsBeforeDeath));
fprintf('\nAssets one period before death: %d decedents by 2006 with assets observed at their last interview (of %d who died); quartiles %s thousand, fraction with a<=0: %.3f\n',numel(AssetsBeforeDeath),sum(died),mat2str(round(quantile(AssetsBeforeDeath,[0.25,0.5,0.75])/1000)),mean(AssetsBeforeDeath<=0))
Targets.MedianAssets_pooled=MedA_pool;             % 5 cohorts x 6 years, all PI quintiles, everyone alive
Targets.MedianAssets_pooled_survivors=MedA_pool_surv;
Targets.AvgAge_pooled=AvgAge_pool;
Targets.AssetsBeforeDeath=AssetsBeforeDeath;       % decedents' assets at the last interview before death (1998 dollars, top-coded)

%% 3b. Medical-expense targets for the endogenous medical-expense model (DFJ2010 Sec. IX; online appendix A.2)
M_in=Mw(:,mcols);
mobs_in=ALIVEw(:,mcols) & ~isnan(Mw(:,mcols));
PI_in=person.PI;
DFJ2010_MedexTargets
Targets.MedexList=MedexList;     % one row per medical-expense moment (cohort, PI half, year, stat, value, n, se)
Targets.MedexCells=MedexCells;   % 5 cohorts x 2 halves x 6 years arrays: Mean, P90, AC1, AC2, counts
% the person-wave values behind those moments, for DFJ2010 appendix Figure A2 (the CDF of medical expenses, data vs model):
% annualised two-year out-of-pocket medical expenses by person (rows) and year 1996..2006 (columns), bottom-coded at 250 by the HRS extract (as DFJ2010),
% with the observed flag, the PI half and the cohort of each person (the cells of the moments are cohort x half x year)
Targets.MedexObs.M=M_in;
Targets.MedexObs.obs=mobs_in;
Targets.MedexObs.half=half_in;
Targets.MedexObs.cohort=cohort_in;
Targets.MedexObs.years=years;

%% 4. Initial distribution (1996) on the model's grids
c96=mcols(1);
got96=ALIVEw(:,c96) & ~isnan(Aw(:,c96)) & ~isnan(Mw(:,c96)) & ~isnan(Hw(:,c96)) & PIq_in>0;
A96=min(Aw(:,c96),asstmax);
H96=Hw(:,c96);
AGE96=AGEw(:,c96);
AGE96(isnan(AGE96))=age1996(isnan(AGE96));
fprintf('\nStep 4: ')
DFJ2010_InitialDist
InitialDist.InitDist=InitDist;      % .(ptype) n_a x 2 x 33 over (a, health, 1996 age)
InitialDist.AgeDist=AgeDist;
InitialDist.ptypemasses=ptypemasses;
InitialDist.InitCounts=InitCounts;
InitialDist.InitObs=InitObs;
InitialDist.a_grid=a_grid;
InitialDist.n_a=n_a;
InitialDist.Names_i=Names_i;
% Note for the model script: the toolkit's jequaloneDist lives at a single age j=1 (age 70), while
% DFJ start each person at their 1996 age (72-102). InitDist keeps the age dimension so the model
% script can either (a) start the age-70 distribution from the youngest cohort, or (b) treat each
% cohort as its own ptype with its own starting age.

%% 5. Same calculations on DFJ's own extract, if it is here
wlthfile='./wealth_singles/data/wlthmat13b.out';
if exist(wlthfile,'file')
    fprintf('\nStep 5: comparing with DFJ''s wlthmat13b.out\n')
    Wd=readmatrix(wlthfile,'FileType','text','Delimiter','\t');
    % column layout (getdat25_wave8.g getdata(), ashft2=1): 1 HHID, 2 weight, 4/5 PI (weighted/unweighted),
    % 7:12 marital status 1996..2006 (0 dead/absent, 1 man, 2 woman), 21:26 assets, 28:33 asset-missing
    % (0 = observed), 35:40 medex, 42:47 medex-missing, 63:68 health (1 bad), 70:75 health-missing, 77:82 age
    keepD=~(Wd(:,77)>dieage | Wd(:,77)<bornage);
    Wd=Wd(keepD,:);
    ND=size(Wd,1);
    HHID_D=Wd(:,1);
    PI_D=Wd(:,5);
    MS_D=Wd(:,7:12).*(Wd(:,77:82)<=dieage);
    A_D=Wd(:,21:26);
    A_D=min(A_D,asstmax);
    obs_D=(Wd(:,28:33)==0) & (MS_D>0);
    AGE_D=Wd(:,77:82);
    AGE_D(~obs_D)=NaN;
    PIq_D=ceil(5*PI_D-1e-9); PIq_D(isnan(PI_D))=0;
    coh_D=zeros(ND,1);
    a94=Wd(:,77)-2;
    for c=1:5
        coh_D(a94>cohortedges(c) & a94<=cohortedges(c+1))=c;
    end
    fprintf('   DFJ extract: %d persons in [70,102] (rows in file %d); PI quintiles %s; cohorts %s\n',ND,numel(keepD),mat2str(histcounts(PIq_D,0.5:1:5.5)),mat2str(histcounts(coh_D,-0.5:1:5.5)))
    % keep the HRS results, run DFJ2010_Targets on DFJ's arrays
    HRS_MedA=MedA; HRS_Cnt=Cnt; HRS_AvgAge=AvgAge; HRS_MomentList=MomentList; HRS_J=J; HRS_VCV=VCV; HRS_VCVm=VCV_medians; HRS_SEm=SE_medians;
    A_in=A_D; obs_in=obs_D; AGE_in=AGE_D; cohort_in=coh_D; PIq_in=PIq_D;
    DFJ2010_Targets
    DFJcheck.MedianAssets=MedA; DFJcheck.CellCounts=Cnt; DFJcheck.AvgAge=AvgAge; DFJcheck.MomentList=MomentList; DFJcheck.VCV=VCV; DFJcheck.J=J; DFJcheck.N=ND;
    fprintf('   moment conditions: HRS rebuild %d, DFJ extract %d (paper: 101)\n',HRS_J,J)
    % the medical-expense moments on DFJ's extract (columns 35:40 medex, 42:47 medex-missing flags; the medex column is
    % already the annualised two-year amount)
    M_D=Wd(:,35:40); mobs_D=(Wd(:,42:47)==0) & (MS_D>0) & ~isnan(M_D);
    HRS_MedexList=MedexList; HRS_MedexCells=MedexCells;
    M_in=M_D; mobs_in=mobs_D; PI_in=PI_D; cohort_in=coh_D;
    DFJ2010_MedexTargets
    DFJcheck.MedexList=MedexList; DFJcheck.MedexCells=MedexCells;
    fprintf('   medical-expense moments: HRS rebuild %d, DFJ extract %d (paper: 249-101=148); cohort 1 upper-half means, HRS vs DFJ: %s vs %s\n',height(HRS_MedexList),height(MedexList),mat2str(round(squeeze(HRS_MedexCells.Mean(1,2,2:end))'/1000,1)),mat2str(round(squeeze(MedexCells.Mean(1,2,2:end))'/1000,1)))
    MedexList=HRS_MedexList; MedexCells=HRS_MedexCells;
    both=~isnan(HRS_MedA) & ~isnan(MedA);
    d=HRS_MedA(both)-MedA(both);
    big=both & MedA>10000;
    r=HRS_MedA(big)./MedA(big);
    fprintf('   cell medians defined in both: %d; correlation %.3f; median |diff| $%.0f; among cells with DFJ median > $10k (%d cells) ratio HRS/DFJ median %.3f, IQR [%.3f, %.3f]\n',sum(both,'all'),corr(HRS_MedA(both),MedA(both)),median(abs(d)),sum(big,'all'),median(r),quantile(r,0.25),quantile(r,0.75))
    for q=1:5
        bq=both; bq(:,[1:q-1,q+1:5],:)=false;
        fprintf('      PI quintile %d: %2d cells, mean HRS %8.0f vs DFJ %8.0f, median |diff| $%6.0f\n',q,sum(bq,'all'),mean(HRS_MedA(bq)),mean(MedA(bq)),median(abs(HRS_MedA(bq)-MedA(bq))))
    end
    fprintf('   cell counts: HRS total %d vs DFJ %d over the 1998-2006 cells\n',sum(HRS_Cnt(:,:,2:end),'all'),sum(Cnt(:,:,2:end),'all'))
    fprintf('   side by side, cohort 1 (aged 74 in 1996), rows PI q1..q5, cols 1996..2006; HRS then DFJ:\n')
    disp(round([squeeze(HRS_MedA(1,:,:)),nan(5,1),squeeze(MedA(1,:,:))]))
    fprintf('   cohort 3 (aged 84 in 1996):\n')
    disp(round([squeeze(HRS_MedA(3,:,:)),nan(5,1),squeeze(MedA(3,:,:))]))
    % person-level overlap on HHID
    HHIDnum=floor(person.HHIDPN/1000);
    fprintf('   persons: HRS rebuild %d, DFJ %d, in both %d\n',Targets.N,ND,sum(ismember(HHIDnum,HHID_D)))
    % restore the HRS results in the workspace
    MedA=HRS_MedA; Cnt=HRS_Cnt; AvgAge=HRS_AvgAge; MomentList=HRS_MomentList; J=HRS_J; VCV=HRS_VCV; VCV_medians=HRS_VCVm; SE_medians=HRS_SEm;
    % Figure: HRS vs DFJ medians, cohort 1 and 3, top and bottom quintile
    fig4=figure(4);
    hold on
    for c=[1,3]
        for q=[1,5]
            plot(squeeze(HRS_AvgAge(c,q,:)),squeeze(HRS_MedA(c,q,:))/1000,'-o','Color',[0 0 0.8])
            plot(squeeze(DFJcheck.AvgAge(c,q,:)),squeeze(DFJcheck.MedianAssets(c,q,:))/1000,'--s','Color',[0.8 0 0])
        end
    end
    hold off
    xlabel('Age'); ylabel('Median assets (thousands of 1998 dollars)')
    title('Cohorts 74 and 84, PI quintiles 1 and 5: HRS rebuild (solid) vs DFJ extract (dashed)')
    saveas(fig4,'./SavedOutput/Graphs/DFJ2010_DataWork_Figure_HRSvsDFJ.png')
else
    DFJcheck=[];
    fprintf('\nStep 5 skipped: %s not found\n',wlthfile)
end

%% Save
save('./SavedOutput/DFJ2010_DataWork.mat','Targets','InitialDist','DFJcheck','person','years','Names_i','a_grid','n_a','-v7.3')
fprintf('\nSaved ./SavedOutput/DFJ2010_DataWork.mat\n')

%% End of run
if headlessFigures==1
    set(0,'DefaultFigureVisible','on');
end
if ~exist('calledfromDFJmain','var')
    diary off
end
