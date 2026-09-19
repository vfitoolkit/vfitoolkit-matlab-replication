% Guvenen2007_DataWork.m
% doPart(1) of Guvenen2007.m -- reproduce Guvenen's (2007) empirical data line:
% the age profile of the cross-sectional variance of log consumption with 95%
% confidence bands, shown as "U.S. Data" in his Figures 7 and 8.
%
% Guvenen takes this profile from Krueger & Perri (2006), who build it from the
% CEX Interview survey. We reproduce it from scratch:
%   (1) import the CEX Interview (FMLI) files, 1980-2003 (BLS releases no
%       1982/83 microdata), from the shared archive DataSets/CEX;
%   (2) build the Krueger-Perri nondurable-consumption CU-year panel via
%       Kaplan2012_CEXConsumptionPanel -- exactly the KP sample and consumption
%       definition; our BASELINE consumption measure is household-total nondurables,
%       following Guvenen (2007) ("household nondurables ... following Deaton-Paxson")
%       (the OECD per-equivalent version is also computed, as a sensitivity);
%   (3) form the weighted variance of log nondurable consumption in each
%       (age, year) cell;
%   (4) Deaton-Paxson (1994) age-cohort decomposition: WLS of the cell variance
%       on a full set of age dummies plus cohort dummies, then read off the
%       cohort-averaged age profile (following HVY2006_DataWork.m in this repo);
%   (5) 95% confidence bands by a within-year block bootstrap over CUs.
%
% Deflation is NOT needed: the WITHIN-(age,year)-cell variance of log
% consumption is invariant to a year-specific price deflator (it only shifts
% every log consumption in a cell by the same constant).
%
% Sample (Krueger-Perri 2006): ages 25-65 (Guvenen's plotted range), complete
% income reporters, no hours/earner screens (whole population, not just workers).
% Consumption = KP nondurable categories, HOUSEHOLD total (baseline); the OECD
% per-equivalent alternative is reported as a sensitivity (it roughly halves the rise).
%
% CPU only. Per-year CEX imports are cached under Data/CEXintermediates/ (lean,
% FMLI only); reruns skip the import. Output: SavedOutput/Guvenen2007_DataWork.mat
% (age grid, variance age profile, 95% band) + a diagnostic figure. This .mat is
% loaded by MakeReplicationTable.m to draw the "US data" line in Figs 7 and 8.

clearvars -except doPart
if exist('Guvenen2007_DataWork_diary.txt','file'); delete('Guvenen2007_DataWork_diary.txt'); end
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end
if ~exist('SavedOutput/Graphs','dir'); mkdir('SavedOutput/Graphs'); end
if ~exist('Data/CEXintermediates','dir'); mkdir('Data/CEXintermediates'); end
diary('Guvenen2007_DataWork_diary.txt');
fprintf('=== Guvenen2007_DataWork.m run %s ===\n',char(datetime('now')));

%% ------------------------------------------------------------------ config
years      = [1980 1981 1984:2003];   % Krueger-Perri (2006) CEX period; 1982/83 not released
envy=getenv('GUV_YEARS'); if ~isempty(envy); years=str2num(envy); end %#ok<ST2NM> % testing override
ageRange   = [25 65];                  % Guvenen Fig 7/8 x-axis
nBoot      = 200;                      % bootstrap replications for the 95% band
rng(20070101,'twister');               % reproducible bootstrap

cexDir='../../DataSets/CEX/';
if ~exist(cexDir,'dir')
    cexDir='/home/kmarshallbanana/Dropbox/Matlab_Codes/DataSets/CEX/';
end
if ~exist(cexDir,'dir'); error('Cannot find the shared CEX archive (DataSets/CEX)'); end

% path to the KP/Kaplan CEX panel builder (reused, not reinvented) and to
% ImportCEXdata (toolkit DataEtc, assumed already on the MATLAB path)
kaplanSub='../../Kaplan2012/Kaplan2012subcodes';
if exist(kaplanSub,'dir'); addpath(kaplanSub); end
if isempty(which('Kaplan2012_CEXConsumptionPanel'))
    error('Kaplan2012_CEXConsumptionPanel not on path (expected in Kaplan2012/Kaplan2012subcodes)');
end
if isempty(which('ImportCEXdata'))
    error('ImportCEXdata not on path (it lives in VFIToolkit-matlab/DataEtc)');
end
klmInter='../../KovacsLowMoran2021/Data/CEXintermediates'; % reuse KLM's cached year structs where present

%% ----------------------------------------------- import/cache FMLI by year
% Cache a LEAN struct per year (FMLI only) so reruns are fast and memory is small.
for yy=years
    yystr=num2str(yy);
    leanfile=['Data/CEXintermediates/CEXfmli_',yystr,'.mat'];
    if exist(leanfile,'file'); continue; end
    klmfile=[klmInter,'/CEXdata_',yystr,'.mat'];
    if exist(klmfile,'file')
        fprintf('Year %s: reusing KLM cached struct, stripping to FMLI\n',yystr);
        S=load(klmfile,'yearData'); yearData=S.yearData; clear S;
    else
        fprintf('Year %s: importing CEX Interview from archive\n',yystr);
        yearData=ImportCEXdata(cexDir,yystr,yystr);
    end
    if ~isfield(yearData,'FMLI')
        error('Year %s import produced no FMLI field',yystr);
    end
    fmli=yearData.FMLI;
    save(leanfile,'fmli','-v7.3');
    clear yearData fmli
end

%% ------------------------------------------------- assemble full FMLI stack
% Merge every year's FMLI quarters into one struct so the panel builder can
% de-collide NEWIDs and apply the Kaplan year-assignment rule over the whole
% stack (its intended usage).
CEXdata=struct(); CEXdata.FMLI=struct();
for yy=years
    S=load(['Data/CEXintermediates/CEXfmli_',num2str(yy),'.mat'],'fmli');
    qn=fieldnames(S.fmli);
    for k=1:numel(qn); CEXdata.FMLI.(qn{k})=S.fmli.(qn{k}); end
    clear S
end
fprintf('Assembled FMLI stack: %d quarter-files\n',numel(fieldnames(CEXdata.FMLI)));

%% ----------------------------------------------------- build the KP panel
opts=struct();
opts.ageRange=ageRange;
opts.applyHoursScreen=false;               % whole population, not just workers
opts.applyEarnerScreens=false;
opts.applyIncompleteReporterScreen=true;   % Krueger-Perri keep complete income reporters
opts.verbose=true;
panel=Kaplan2012_CEXConsumptionPanel(CEXdata,opts);
clear CEXdata

% keep CU-years with a valid (positive) consumption and an age in the plotted range.
% BASELINE consumption = household-total nondurables: the measure Guvenen (2007) uses
% ("household nondurables ... following Deaton-Paxson (1994)", data from Krueger-Perri
% 2006). We ALSO carry the OECD per-equivalent version as a sensitivity -- it gives a
% much flatter profile (see the Data appendix in the write-up).
keep = isfinite(panel.nondurable) & panel.nondurable>0 & panel.fam_size>=1 ...
       & isfinite(panel.nondurable_pe) & panel.nondurable_pe>0 ...
       & panel.age_head>=ageRange(1) & panel.age_head<=ageRange(2) ...
       & isfinite(panel.weight) & panel.weight>0;
P=panel(keep,:);
fprintf('Panel: %d CU-years after selection (ages %d-%d, complete reporters)\n',height(P),ageRange(1),ageRange(2));

age    = double(P.age_head);
yr     = double(P.year);
w      = double(P.weight);
fam    = double(P.fam_size);
xhh    = log(double(P.nondurable));      % BASELINE: household-total nondurables
xpe    = log(double(P.nondurable_pe));   % SENSITIVITY 1: OECD per-equivalent
% SENSITIVITY 2: US Census equivalence scale (Dalaker-Naifeh 1998) -- the scale
% Krueger-Perri (2006) actually use. We approximate it by the size ratios of the
% official poverty thresholds (weighted-average thresholds, ~stable over time),
% normalised to 1 for a one-person household:
censusRatio=[1.00 1.27 1.56 2.01 2.38 2.70 3.07 3.41 4.02]; % family size 1..9+
cscale=censusRatio(min(max(round(fam),1),9))';
xcen=log(double(P.nondurable)./cscale);  % household nondurable / census adult-equivalents
cohort = yr-age;                         % birth year

% index maps
Jages  = ageRange(1):ageRange(2); J=numel(Jages);
Uyears = unique(yr(:))';            Ty=numel(Uyears);
Ucoh   = unique(cohort(:))';        Nc=numel(Ucoh);
aj = age-ageRange(1)+1;                                    % 1..J
[~,ty] = ismember(yr,Uyears);                             % 1..Ty
[~,cc] = ismember(cohort,Ucoh);                           % 1..Nc
fprintf('Cells: %d ages x %d years, %d birth cohorts\n',J,Ty,Nc);

%% -------------------------------- Deaton-Paxson age profiles (point estimates)
ageprofile     = local_dp_ageprofile(aj,ty,cc,xhh,w,J,Ty,Nc); % household total (baseline)
ageprofile_pe  = local_dp_ageprofile(aj,ty,cc,xpe,w,J,Ty,Nc); % OECD per-equivalent (sens 1)
ageprofile_cen = local_dp_ageprofile(aj,ty,cc,xcen,w,J,Ty,Nc);% Census adult-equiv  (sens 2, KP's scale)

% Diagnostic: raw pooled (cohort-UNADJUSTED) household-total profile.
rawprofile=nan(J,1);
for j=1:J
    m=(aj==j); ww=w(m); xx=xhh(m);
    mu=sum(ww.*xx)/sum(ww); rawprofile(j)=sum(ww.*(xx-mu).^2)/sum(ww);
end

% robust rise summary: the endpoints are noisy (age 25 high, age 65 dips), so we report
% the linear-fit 40-year rise rather than the 25->65 endpoint difference.
riseFit=@(p) polyval(polyfit(Jages(:),p(:),1),Jages(end))-polyval(polyfit(Jages(:),p(:),1),Jages(1));
fprintf('Household-total  : linfit rise=%.3f, trough-peak=%.3f\n',riseFit(ageprofile),max(ageprofile)-min(ageprofile));
fprintf('Census adult-eq  : linfit rise=%.3f, trough-peak=%.3f\n',riseFit(ageprofile_cen),max(ageprofile_cen)-min(ageprofile_cen));
fprintf('OECD per-equiv   : linfit rise=%.3f, trough-peak=%.3f\n',riseFit(ageprofile_pe),max(ageprofile_pe)-min(ageprofile_pe));

%% ------------------------------------------- within-year bootstrap 95% bands
% Resample CUs with replacement WITHIN each year; recompute BOTH profiles each rep.
yearRows=cell(Ty,1); for t=1:Ty; yearRows{t}=find(ty==t); end
BOOT=nan(J,nBoot); BOOT_pe=nan(J,nBoot); BOOT_cen=nan(J,nBoot);
for b=1:nBoot
    idx=cell(Ty,1);
    for t=1:Ty
        rows=yearRows{t}; n=numel(rows);
        idx{t}=rows(randi(n,n,1));
    end
    ii=vertcat(idx{:});
    BOOT(:,b)    =local_dp_ageprofile(aj(ii),ty(ii),cc(ii),xhh(ii),w(ii),J,Ty,Nc);
    BOOT_pe(:,b) =local_dp_ageprofile(aj(ii),ty(ii),cc(ii),xpe(ii),w(ii),J,Ty,Nc);
    BOOT_cen(:,b)=local_dp_ageprofile(aj(ii),ty(ii),cc(ii),xcen(ii),w(ii),J,Ty,Nc);
    if mod(b,50)==0; fprintf('  bootstrap %d/%d\n',b,nBoot); end
end
ciLow    =prctile(BOOT,2.5 ,2);     ciHigh    =prctile(BOOT,97.5,2);
ciLow_pe =prctile(BOOT_pe,2.5 ,2);  ciHigh_pe =prctile(BOOT_pe,97.5,2);
ciLow_cen=prctile(BOOT_cen,2.5 ,2); ciHigh_cen=prctile(BOOT_cen,97.5,2);

%% ---------------- third fact: consumption growth 25->55 by education (Fig 11)
% Cohort-adjusted mean log REAL consumption profile by education group, then the
% total growth from age 25 to 55. College = educbucket 4 (college+), high school =
% educbucket 2. Deflate to real with annual CPI-U (1982-84=100) so the level growth is
% real; the college-HS GAP is robust to deflation anyway (both face the same CPI).
cpiYr=1980:2003;
cpi  =[82.4 90.9 96.5 99.6 103.9 107.6 109.6 113.6 118.3 124.0 130.7 136.2 140.3 ...
       144.5 148.2 152.4 156.9 160.5 163.0 166.6 172.2 177.1 179.9 184.0];
[~,ic]=ismember(yr,cpiYr); cpiRow=cpi(ic)'/100;
xhh_real=xhh-log(cpiRow); xpe_real=xpe-log(cpiRow);
edu=double(P.educbucket);
grow2555=@(xr,g) local_dp_growth(aj(edu==g),ty(edu==g),cc(edu==g),xr(edu==g),w(edu==g),J,Ty,Nc);
eduColGrowth_hh=grow2555(xhh_real,4); eduHSGrowth_hh=grow2555(xhh_real,2);
eduColGrowth_pe=grow2555(xpe_real,4); eduHSGrowth_pe=grow2555(xpe_real,2);
fprintf('\n--- data consumption growth 25->55 by education (cohort-adj, real) ---\n');
fprintf('household: college %+.0f%%, HS %+.0f%%, gap %.0fpp\n',100*eduColGrowth_hh,100*eduHSGrowth_hh,100*(eduColGrowth_hh-eduHSGrowth_hh));
fprintf('per-equiv: college %+.0f%%, HS %+.0f%%, gap %.0fpp\n',100*eduColGrowth_pe,100*eduHSGrowth_pe,100*(eduColGrowth_pe-eduHSGrowth_pe));

%% -------------------------------------------------------------------- save
DataLine.ages   = Jages(:);          % 25..65
DataLine.varlogc= ageprofile(:);     % BASELINE: household-total, cohort-averaged var(log c)
DataLine.ciLow  = ciLow(:);
DataLine.ciHigh = ciHigh(:);
DataLine.rawvarlogc = rawprofile(:); % raw pooled household-total, diagnostic
DataLine.rise   = riseFit(ageprofile);              % linear-fit 40yr rise (baseline)
DataLine.risePeakTrough = max(ageprofile)-min(ageprofile);
DataLine.varlogc_pe = ageprofile_pe(:);             % SENSITIVITY 1: OECD per-equivalent
DataLine.ciLow_pe = ciLow_pe(:);
DataLine.ciHigh_pe= ciHigh_pe(:);
DataLine.rise_pe  = riseFit(ageprofile_pe);
DataLine.varlogc_cen = ageprofile_cen(:);           % SENSITIVITY 2: Census adult-equiv (KP's scale)
DataLine.ciLow_cen = ciLow_cen(:);
DataLine.ciHigh_cen= ciHigh_cen(:);
DataLine.rise_cen  = riseFit(ageprofile_cen);
DataLine.eduColGrowth_hh=eduColGrowth_hh;  % Fig-11 data: real cons growth 25->55, college, household
DataLine.eduHSGrowth_hh =eduHSGrowth_hh;   % high school, household
DataLine.eduColGrowth_pe=eduColGrowth_pe;  % college, per-equivalent
DataLine.eduHSGrowth_pe =eduHSGrowth_pe;   % high school, per-equivalent
DataLine.scale  = 'household-total nondurables (baseline); Census adult-equiv (*_cen, KP scale) and OECD per-equivalent (*_pe) as sensitivities';
DataLine.source = 'CEX Interview 1980-2003 (excl 1982/83), Krueger-Perri/Deaton-Paxson household nondurables, Deaton-Paxson age-cohort decomposition, within-year bootstrap 95% band';
DataLine.nCUyears = height(P);
save('SavedOutput/Guvenen2007_DataWork.mat','DataLine');
fprintf('\nWrote SavedOutput/Guvenen2007_DataWork.mat\n');
fprintf('BASELINE household-total var(log c): linfit rise = %.3f (Guvenen/KP report ~0.21)\n',DataLine.rise);
fprintf('SENSITIVITY Census adult-equiv (KP scale) rise = %.3f\n',DataLine.rise_cen);
fprintf('SENSITIVITY OECD per-equivalent rise = %.3f\n',DataLine.rise_pe);

%% ------------------------------------------------------------------ figures
% (1) baseline data line + band -- this is what MakeReplicationTable overlays on Figs 7/8
fig=figure(1); clf; hold on;
fill([Jages,fliplr(Jages)],[ciLow',fliplr(ciHigh')],[0.85 0.85 0.85],'EdgeColor','none');
plot(Jages,ageprofile,'k-','LineWidth',1.8);
xlabel('age'); ylabel('var(log c)');
legend('95% bootstrap band','CEX data (household nondurables)','Location','northwest');
title('Empirical consumption inequality (CEX household nondurables, Deaton-Paxson)');
saveas(fig,'SavedOutput/Graphs/Guvenen2007_DataLine.png');
fprintf('Wrote SavedOutput/Graphs/Guvenen2007_DataLine.png\n');

% (2) scale comparison for the Data appendix: household total vs Census adult-equiv vs OECD
fig=figure(2); clf; hold on;
fill([Jages,fliplr(Jages)],[ciLow',fliplr(ciHigh')],[0.80 0.82 0.92],'EdgeColor','none','HandleVisibility','off');
fill([Jages,fliplr(Jages)],[ciLow_cen',fliplr(ciHigh_cen')],[0.82 0.92 0.82],'EdgeColor','none','HandleVisibility','off');
fill([Jages,fliplr(Jages)],[ciLow_pe',fliplr(ciHigh_pe')],[0.94 0.86 0.80],'EdgeColor','none','HandleVisibility','off');
plot(Jages,ageprofile,'b-','LineWidth',1.8);
plot(Jages,ageprofile_cen,'-','Color',[0 0.5 0],'LineWidth',1.8);
plot(Jages,ageprofile_pe,'r-','LineWidth',1.8);
xlabel('age'); ylabel('var(log c)');
legend('Household total (Guvenen/Deaton-Paxson)','Census adult-equiv (Krueger-Perri)','OECD per-equivalent','Location','northwest');
title('CEX consumption inequality under three equivalence scales (95% bands)');
saveas(fig,'SavedOutput/Graphs/Guvenen2007_DataLine_scales.png');
fprintf('Wrote SavedOutput/Graphs/Guvenen2007_DataLine_scales.png\n');

fprintf('\n=== Guvenen2007_DataWork done ===\n');
diary off;

%% ===================================================================== subfn
function [ageprofile,cellvar] = local_dp_ageprofile(aj,ty,cc,x,w,J,Ty,Nc)
% Weighted cell variances of x=log c in each (age,year) cell, then WLS on
% age + cohort dummies; return the cohort-averaged age profile (length J).
lin = aj + (ty-1)*J;                       % (age,year) cell linear index, 1..J*Ty
sw  = accumarray(lin,w,          [J*Ty,1]);
swx = accumarray(lin,w.*x,       [J*Ty,1]);
swxx= accumarray(lin,w.*x.*x,    [J*Ty,1]);
cnt = accumarray(lin,ones(size(w)),[J*Ty,1]);
mu  = swx./sw;
v   = swxx./sw - mu.^2;                     % weighted variance per cell
% cohort index per cell: cohort = year-age; recover from the (age,year) grid
[ajg,~]=ndgrid(1:J,1:Ty);
% map each cell to its cohort bucket using the same Ucoh ordering as the caller:
% cohort value = (age) - (year); we need indices matching cc. Rebuild via lookup.
% Caller passed cc per observation; reconstruct cell->cohort from observations.
cohortOfCell=accumarray(lin,cc,[J*Ty,1],@max); % all obs in a cell share one cohort
valid = cnt>=2 & sw>0 & isfinite(v) & cohortOfCell>0;   % need >=2 obs for a variance
y = v(valid); ncell=cnt(valid);
ajc = ajg(:); ajc=ajc(valid);
coc = cohortOfCell(valid);
% design: intercept + age deviations (2..J) + cohort deviations (2..Nc)
n=numel(y);
Xage=zeros(n,J-1);  for j=2:J;  Xage(:,j-1)=(ajc==j); end
Xcoh=zeros(n,Nc-1); for c=2:Nc; Xcoh(:,c-1)=(coc==c); end
X=[ones(n,1),Xage,Xcoh];
bhat=lscov(X,y,ncell);                      % WLS, weight by cell sample size
b0=bhat(1); ba=[0;bhat(2:J)]; bc=[0;bhat(J+1:J+Nc-1)];
% cohort-average: weight cohort deviations by total obs per cohort
Nperc=accumarray(coc,ncell,[Nc,1]);
gbar = sum(Nperc.*bc)/sum(Nperc);
ageprofile = b0 + ba + gbar;                % length J
cellvar = v;
end

function g = local_dp_growth(aj,ty,cc,x,w,J,Ty,Nc)
% Cohort-averaged MEAN log-consumption age profile (Deaton-Paxson: WLS of cell mean on
% age + cohort dummies), then total growth from age 25 (index 1) to 55 (index 31).
lin = aj + (ty-1)*J;
sw  = accumarray(lin,w,          [J*Ty,1]);
swx = accumarray(lin,w.*x,       [J*Ty,1]);
cnt = accumarray(lin,ones(size(w)),[J*Ty,1]);
mu  = swx./sw;                              % weighted cell mean of log c
cohortOfCell=accumarray(lin,cc,[J*Ty,1],@max);
[ajg,~]=ndgrid(1:J,1:Ty);
valid = cnt>=1 & sw>0 & isfinite(mu) & cohortOfCell>0;
y = mu(valid); ncell=cnt(valid);
ajc = ajg(:); ajc=ajc(valid); coc = cohortOfCell(valid);
n=numel(y);
Xage=zeros(n,J-1);  for j=2:J;  Xage(:,j-1)=(ajc==j); end
Xcoh=zeros(n,Nc-1); for c=2:Nc; Xcoh(:,c-1)=(coc==c); end
bhat=lscov([ones(n,1),Xage,Xcoh],y,ncell);
b0=bhat(1); ba=[0;bhat(2:J)]; bc=[0;bhat(J+1:J+Nc-1)];
Nperc=accumarray(coc,ncell,[Nc,1]); gbar=sum(Nperc.*bc)/sum(Nperc);
prof = b0 + ba + gbar;                      % cohort-averaged mean log c, length J
g = exp(prof(31)-prof(1))-1;                % real growth age 25 -> 55
end
