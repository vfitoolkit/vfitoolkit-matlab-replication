%% Empirical part of Huggett, Ventura & Yaron (2011): the earnings age profiles of Figure 1
% These age profiles (mean earnings, variance of log earnings, and mean/median earnings, for ages
% 23-60, time-effects view) are the calibration targets of HVY2011.
%
% HVY2011 Section IIA: PSID 1969-2004 family files, male household heads, ages 23 to 60.
% Statistics are computed for five-year age bins, by age bin and panel year. Males over age 30
% must work 520-5820 hours and earn at least $1,500 (1968 prices); males age 30 and below must
% work 260-5820 hours and earn at least $1,000 (1968 prices). The SEO subsample is excluded.
% Age effects are then extracted with either cohort dummies (cohort-effects view) or time
% dummies (time-effects view), the method of Deaton & Paxson (1994).
%
% This follows the authors' replication code, data/hvy_FIGS/Master_Figures.m, line by line,
% using the PSID extracts that come with it (data/hvy_FIGS/out69.txt,...,out102.txt; one file per
% panel year, each following people for four years). There are no PSID data for 1997, 1999,
% 2001 and 2003, so the panel years are 1969-1996, 1998, 2000, 2002 (31 panel years).
%
% Three things the authors' code does that differ from the description in the paper, which we
% follow so as to reproduce their numbers:
% (i) the age bin labelled 'age j' is ages j to j+4 (the paper says 'centered at age j');
% (ii) the mean earnings regression is on the level of mean earnings (the paper says log);
% (iii) earnings are deflated by the CPI of the j-th panel year counted from 1969, so the 1998,
%      2000 and 2002 panel years use the CPI of 1997, 1998 and 1999 (as the missing years are skipped).
% The minimum-wage restriction (earnings at least 520 hours at the federal minimum wage) is in the
% authors' code but not mentioned in the paper.
%
% The authors' code also saved the age-by-year cell statistics as data/hvy_FIGS/dp_all.mat, and
% their calibration targets are in data/hvy-programfiles/mean-time-aer08, var-time-aer08 and
% skew-time-aer08. We check our numbers against both.

if ~exist('figure_c','var')
    figure_c=0; % Normally this will be created by HuggettVenturaYaron2011.m, which calls this HVY2011_DataWork.m script. This is just so that you can run HVY2011_DataWork independently if you want.
end
if ~exist('CreateFigures','var')
    CreateFigures=1; % Normally this will be created by HuggettVenturaYaron2011.m, which calls this HVY2011_DataWork.m script. This is just so that you can run HVY2011_DataWork independently if you want.
end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end % so HVY2011_DataWork can also be run independently

datadir='./data/hvy_FIGS/';

%% Settings (as in Master_Figures.m)
cwid=1;      % width of age cells
coverlap=5;  % size of the age bins in years
cmin=23;
cmax=60;
numcells=floor((cmax-cmin)/cwid+1); % 38 age cells
cellages=(cmin:cwid:cmax)'; % lower edge of each age bin

panelyears=[69:1:96,98,100,102]; % two-digit PSID panel years (no data for 1997, 1999, 2001, 2003)
npanels=length(panelyears); % 31

% Earnings and hours criteria
earnup=1800000;  % maximum earnings (1968 prices)
earnlb=1500;     % minimum earnings, over age 30 (1968 prices)
earnlbY=1000;    % minimum earnings, age 30 and below (1968 prices)
maxhours=5820;
minhours=520;    % minimum hours, over age 30
minhoursY=260;   % minimum hours, age 30 and below

% CPI (equal to 1 in 1969) and the federal minimum wage (1969 onwards)
cpi_def=readmatrix([datadir,'CPI_BLS.xlw'],'FileType','spreadsheet','Range','D5:D41'); % 1969-2005
minwage=readmatrix([datadir,'MinimumWage.xlsx'],'Range','B1:B35'); % 1969-2003

%% Cell statistics: by age bin (rows) and panel year (columns)
% Column order of the outYY.txt files (from the authors' readme):
% famid68 persnr sex(3) educ(4) educftr seq1(6) seq2(7) seq3 seq4 age1(10) age2 age3 age4 mar1(14)...mar4
% hrs1(18) hrs2 hrs3 hrs4 inc1(22) inc2 inc3 inc4 food1..food4 lstjobten1..4 wge1(34)..wge4(37) wealth...
% (we only use the first year of each four-year panel: age1, hrs1, inc1)
CellMean=zeros(numcells,npanels);      % mean of earnings (level)
CellVarLog=zeros(numcells,npanels);    % variance of log earnings
CellSkew=zeros(numcells,npanels);      % mean/median of earnings (level)
CellGini=zeros(numcells,npanels);      % Gini of earnings (level)
Cell9010=zeros(numcells,npanels);      % 90-10 quantile ratio of earnings (for the HVY2011 Section V.C sensitivity on the moments to match)
Cell7525=zeros(numcells,npanels);      % 75-25 quantile ratio of earnings
CellN=zeros(numcells,npanels);         % number of observations
logearnings_byage=cell(numcells,1);    % log earnings pooled across panel years, by age bin

for jj=1:npanels
    iy=panelyears(jj);
    TEMP=readmatrix([datadir,'out',num2str(iy),'.txt'],'FileType','text');
    jy=iy-68; % index of the calendar year for the minimum wage

    % Real earnings bounds, deflated with the CPI of the jj-th panel year (as in the authors' code, see note (iii) above)
    realearnup=earnup*cpi_def(jj);
    realearnlb=earnlb*cpi_def(jj);
    realearnlbY=earnlbY*cpi_def(jj);

    age=TEMP(:,10);
    earn=TEMP(:,22);
    hrs=TEMP(:,18);
    if iy==98 || iy==100 || iy==102
        head=(TEMP(:,7)==1); % the sequence-number variable moves for these years
    else
        head=(TEMP(:,6)==1);
    end
    % Selection that does not depend on the age bin
    select_all=(earn~=999999) & (earn~=99999) & (earn~=9999999) & (earn>0) ...
        & (TEMP(:,1)>=1 & TEMP(:,1)<=2930) ... % exclude the SEO sample
        & (earn>=minhours*minwage(jy)) ... % at least 520 hours at the federal minimum wage
        & (earn<realearnup) ...
        & (((age<=30).*earn>=realearnlbY)+((age>30).*earn>=realearnlb)) ... % minimum earnings, depending on age
        & (hrs<=maxhours & ((age>30).*(hrs>=minhours)+(age<=30).*(hrs>=minhoursY))) ...
        & head & (TEMP(:,4)<=4) & (TEMP(:,3)==1); % head of household, education code, male

    for kk=1:numcells
        cl=cellages(kk);
        select=select_all & (age>=cl & age<cl+coverlap);
        REARNINGS=earn(select)/cpi_def(jj); % real earnings
        rearnings=log(REARNINGS);
        if isempty(rearnings)
            error('HVY2011_DataWork: empty age-year cell (age %i, year %i)',cl,iy)
        end
        nobs=length(REARNINGS);
        CellN(kk,jj)=nobs;
        CellMean(kk,jj)=mean(REARNINGS);
        CellVarLog(kk,jj)=var(rearnings); % [std()^2 in the authors' code, i.e. normalized by N-1]
        CellSkew(kk,jj)=mean(REARNINGS)/median(REARNINGS);
        % Gini (gini1.m in the authors' code): 0.5*sum_i sum_j |x_i-x_j| / (N^2 * mean)
        sortedE=sort(REARNINGS);
        CellGini(kk,jj)=(2*sum((1:nobs)'.*sortedE)/(nobs*sum(sortedE))-(nobs+1)/nobs); % same as the double sum formula
        % Quantile ratios: the p-th quantile is the smallest earnings with at least a fraction p of the cell at or below it
        Cell9010(kk,jj)=sortedE(ceil(0.9*nobs))/sortedE(ceil(0.1*nobs));
        Cell7525(kk,jj)=sortedE(ceil(0.75*nobs))/sortedE(ceil(0.25*nobs));
        logearnings_byage{kk}=[logearnings_byage{kk}; rearnings];
    end
end

% Statistics by age bin, pooled across panel years (used to set the level of the age effects)
PooledVarLog=zeros(numcells,1);
PooledMean=zeros(numcells,1);
PooledGini=zeros(numcells,1);
PooledSkew=zeros(numcells,1);
Pooled9010=zeros(numcells,1);
Pooled7525=zeros(numcells,1);
for kk=1:numcells
    E=exp(logearnings_byage{kk});
    nobs=length(E);
    PooledVarLog(kk)=var(logearnings_byage{kk});
    PooledMean(kk)=mean(E);
    sortedE=sort(E);
    PooledGini(kk)=(2*sum((1:nobs)'.*sortedE)/(nobs*sum(sortedE))-(nobs+1)/nobs);
    PooledSkew(kk)=mean(E)/median(E);
    Pooled9010(kk)=sortedE(ceil(0.9*nobs))/sortedE(ceil(0.1*nobs));
    Pooled7525(kk)=sortedE(ceil(0.75*nobs))/sortedE(ceil(0.25*nobs));
end

fprintf('HVY2011 data: %i age-year cells, observations per cell: min %i, median %.1f (HVY2011: at least 100 in each age-year bin) \n',numel(CellN),min(CellN(:)),median(CellN(:)))

%% Check against the cell statistics saved by the authors (data/hvy_FIGS/dp_all.mat)
dpall=load([datadir,'dp_all.mat']);
fprintf('Max abs difference from the authors'' cell statistics (dp_all.mat): mean %8.2e, var(log) %8.2e, mean/median %8.2e, Gini %8.2e \n', ...
    max(max(abs(CellMean-dpall.REARNINGS_stat(:,:,1)))), max(max(abs(CellVarLog-dpall.rearnings_stat(:,:,2)))), ...
    max(max(abs(CellSkew-dpall.REARNINGS_stat(:,:,3)))), max(max(abs(CellGini-dpall.REARNINGS_stat(:,:,5)))))
fprintf('   and pooled by age: var(log) %8.2e, mean %8.2e, Gini %8.2e, mean/median %8.2e \n', ...
    max(abs(PooledVarLog-dpall.mn_age_earns(:,3))), max(abs(PooledMean-dpall.mn_age_earns(:,6))), max(abs(PooledGini-dpall.mn_age_earns(:,7))), max(abs(PooledSkew-dpall.mn_age_earns(:,8))))

%% Deaton-Paxson regressions: age effects with cohort dummies, and with time dummies
H=numcells;
T=102-69+1; % 34 calendar years, 1969-2002
% Cohort dummies: for the year t column, age h is cohort t+H-h
ncohorts=H+T-1;
cdum=zeros(H*T,ncohorts);
block=fliplr(eye(H));
for tt=1:T
    cdum((tt-1)*H+1:tt*H,tt:H+tt-1)=block;
end
% Remove the missing years (1997, 1999, 2001 are years 29, 31, 33 of 34) and the two cohorts only seen in them (as in the authors' code)
cdum(:,[end-3,end-1])=[];
cdum([H*(T-6)+1:H*(T-5), H*(T-4)+1:H*(T-3), H*(T-2)+1:H*(T-1)],:)=[];
T=T-3; % 31 panel years
adum=kron(ones(T,1),eye(H)); % age dummies
X=[adum, cdum(:,1:end-1)]; % cohort-effects view
tdum=kron(eye(T),ones(H,1));
X_td=[adum, tdum(:,1:end-1)]; % time-effects view

i1=find(cellages==40-floor(coverlap/2)); % the age bin at which the profiles are anchored (labelled 38, i.e. ages 38-42, centred on 40)

% Variance of log earnings: anchored to the pooled value at age 40
beta_c=X\CellVarLog(:);
beta_t=X_td\CellVarLog(:);
Fig1_VarLog_cohort=beta_c(1:H)-beta_c(i1)+PooledVarLog(i1);
Fig1_VarLog_time=beta_t(1:H)-beta_t(i1)+PooledVarLog(i1);

% Mean earnings (level): anchored to the pooled mean at age 40, then both scaled so the time-effects view is 100 at age 60
beta_c=X\CellMean(:);
beta_t=X_td\CellMean(:);
Fig1_Mean_cohort=beta_c(1:H)-beta_c(i1)+PooledMean(i1);
Fig1_Mean_time=beta_t(1:H)-beta_t(i1)+PooledMean(i1);
Fig1_Mean_cohort=100*Fig1_Mean_cohort/Fig1_Mean_time(end);
Fig1_Mean_time=100*Fig1_Mean_time/Fig1_Mean_time(end);

% Gini
beta_c=X\CellGini(:);
beta_t=X_td\CellGini(:);
Fig1_Gini_cohort=beta_c(1:H)-beta_c(i1)+PooledGini(i1);
Fig1_Gini_time=beta_t(1:H)-beta_t(i1)+PooledGini(i1);

% Skewness (mean/median)
beta_c=X\CellSkew(:);
beta_t=X_td\CellSkew(:);
Fig1_Skew_cohort=beta_c(1:H)-beta_c(i1)+PooledSkew(i1);
Fig1_Skew_time=beta_t(1:H)-beta_t(i1)+PooledSkew(i1);

% Quantile ratios 90-10 and 75-25 (not in HVY2011 Figure 1; the targets of the HVY2011 Section V.C sensitivity on the moments
% to match), computed in the same way as the other statistics: anchored to the pooled value at age 40
beta_c=X\Cell9010(:);
beta_t=X_td\Cell9010(:);
Data_P9010_cohort=beta_c(1:H)-beta_c(i1)+Pooled9010(i1);
Data_P9010_time=beta_t(1:H)-beta_t(i1)+Pooled9010(i1);
beta_c=X\Cell7525(:);
beta_t=X_td\Cell7525(:);
Data_P7525_cohort=beta_c(1:H)-beta_c(i1)+Pooled7525(i1);
Data_P7525_time=beta_t(1:H)-beta_t(i1)+Pooled7525(i1);
fprintf('Quantile ratios (time effects), ages 23, 40, 60: 90-10 %6.3f %6.3f %6.3f, 75-25 %6.3f %6.3f %6.3f \n',Data_P9010_time([1,18,38]),Data_P7525_time([1,18,38]))

agedata=cellages; % 23 to 60

%% Check against the authors' calibration targets (time-effects view)
AuthorTarget_Mean=readmatrix('./data/hvy-programfiles/mean-time-aer08','FileType','text');
AuthorTarget_VarLog=readmatrix('./data/hvy-programfiles/var-time-aer08','FileType','text');
AuthorTarget_Skew=readmatrix('./data/hvy-programfiles/skew-time-aer08','FileType','text');
fprintf('Max abs difference from the authors'' calibration targets: mean %8.2e, var(log) %8.2e, mean/median %8.2e (the target files are rounded to 3-4 decimals) \n', ...
    max(abs(Fig1_Mean_time-AuthorTarget_Mean)), max(abs(Fig1_VarLog_time-AuthorTarget_VarLog)), max(abs(Fig1_Skew_time-AuthorTarget_Skew)))

fprintf('HVY2011 Figure 1: rise from age 23 to 60 in the variance of log earnings: cohort effects %6.3f, time effects %6.3f (HVY2011: about 0.4 and 0.2) \n',Fig1_VarLog_cohort(end)-Fig1_VarLog_cohort(1),Fig1_VarLog_time(end)-Fig1_VarLog_time(1))

%% HVY2011 Figure 1
if CreateFigures==1
    figure_c=figure_c+1;
    fig=figure(figure_c);
    fig.Position=[100 100 900 700];
    subplot(2,2,1);
    plot(agedata,Fig1_Mean_cohort,'--b',agedata,Fig1_Mean_time,'-r','LineWidth',1.5);
    xlim([23,60]); ylim([0,140]);
    xlabel('Age');
    title('Panel A: Mean earnings');
    legend({'Cohort effects','Time effects'},'Location','southeast');
    subplot(2,2,2);
    plot(agedata,Fig1_VarLog_cohort,'--b',agedata,Fig1_VarLog_time,'-r','LineWidth',1.5);
    xlim([23,60]); ylim([0,0.9]);
    xlabel('Age');
    title('Panel B: Variance of log earnings');
    subplot(2,2,3);
    plot(agedata,Fig1_Gini_cohort,'--b',agedata,Fig1_Gini_time,'-r','LineWidth',1.5);
    xlim([23,60]); ylim([0,0.6]);
    xlabel('Age');
    title('Panel C: Earnings Gini');
    subplot(2,2,4);
    plot(agedata,Fig1_Skew_cohort,'--b',agedata,Fig1_Skew_time,'-r','LineWidth',1.5);
    xlim([23,60]); ylim([0.8,1.6]);
    xlabel('Age');
    title('Panel D: Earnings skewness (mean/median)');
    saveas(fig,'./SavedOutput/Graphs/HVY2011_Fig1.png');
end
