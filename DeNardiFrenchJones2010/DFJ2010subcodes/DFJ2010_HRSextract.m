% DFJ2010_HRSextract: build the DFJ2010 AHEAD singles panel from the RAND HRS Longitudinal File
%
% Rebuilds (the content of) wealth_singles/dataprep2.dta from public RAND HRS data, following
% DFJ2010 JPE Appendix B.1-B.3 (DFJ2010_SupplementalAppendix.pdf) and the notes in
% HRS_rebuild_plan.md. Output is a long person-wave table plus person-level variables,
% saved to ./SavedOutput/DFJ2010_HRSpanel.mat. Later steps (wlthmat layout, first-stage
% profiles, moments, initial distribution) build on this file.
%
% Requires: DFJ2010_ReadStata114.m (same folder) and the RAND HRS Longitudinal File 2022 (V1)
% Stata version in HRSdatadir (set below).
%
% Conventions (DFJ2010): RAND wave 2 (AHEAD 1993) is DFJ's "1994", RAND wave 3 (AHEAD 1995)
% is DFJ's "1996", then 1998,...,2006 = RAND waves 4..8. All dollars converted to 1998 dollars.

%% Where the data are
HRSdatadir='../../DataSets/HRS/HRS Longitudinal 2022/'; % relative to the replication folder (cwd)
HRSfile=[HRSdatadir,'randhrs1992_2022v1.dta'];
if ~exist(HRSfile,'file')
    error('DFJ2010_HRSextract: cannot find %s (cwd must be the DeNardiFrenchJones2010 replication folder)',HRSfile)
end
if ~exist('SavedOutput','dir'); mkdir('SavedOutput'); end

waves=2:8; % RAND waves
DFJyear=[1994,1996,1998,2000,2002,2004,2006]; % DFJ's calendar-year labels for those waves
nW=numel(waves);

%% CPI-U (annual average, 1982-84=100) to put everything in DFJ's '1998 dollars'
% source: BLS CPI-U all items, US city average
% Matching against dataprep2.dta shows DFJ deflate every dollar amount reported in a wave by the
% CPI of the calendar year BEFORE the interview (the year the income questions refer to), with
% the 1998 wave as base, i.e. base = CPI(1997). Their 1998-wave amounts are exactly nominal, and
% every other wave differs from a CPI(interview year) deflation by the 1997/1998 CPI ratio.
CPIyear=1992:2008;
CPIval=[140.3,144.5,148.2,152.4,156.9,160.5,163.0,166.6,172.2,177.1,179.9,184.0,188.9,195.3,201.6,207.3,215.3];
CPI1998=CPIval(CPIyear==1997); % the base DFJ actually use for '1998 dollars' (see note above)

%% Variables to pull
personvars={'HHIDPN','HHID','PN','RAGENDER','HACOHORT','RABYEAR','RADYEAR','RADMONTH',...
    'REOOPMD','REOOPMDF','REXITWV','REMSTAT','REMPART'};
wavestems_R={'IWSTAT','MSTAT','MPART','SHLT','AGEY_E','IWENDY','IWENDM','OOPMD','IEARN','ISRET','ISSDI','IPENA','IGXFR','IUNWC','WTRESP','WTHH','GOVMD','HIGOV','NRSHOM','HOSP','PRPRM1','PRPRM2','PRPRM3','MRPREM'};
wavestems_H={'ATOTB','ITOT','ICAP','HHRES'};
wantvars=personvars;
for w=waves
    for s=1:numel(wavestems_R)
        wantvars{end+1}=sprintf('R%d%s',w,wavestems_R{s}); %#ok<SAGROW>
    end
    for s=1:numel(wavestems_H)
        wantvars{end+1}=sprintf('H%d%s',w,wavestems_H{s}); %#ok<SAGROW>
    end
end

%% Read
fprintf('DFJ2010_HRSextract: reading %d variables from the RAND HRS Longitudinal File...\n',numel(wantvars))
tic
[T,~,readinfo]=DFJ2010_ReadStata114(HRSfile,wantvars,1);
fprintf('   done in %.0f seconds\n',toc)
if ~isempty(readinfo.notfound)
    fprintf('   WARNING: variables not found: %s\n',strjoin(readinfo.notfound,' '))
end
N0=height(T);

%% Wide arrays (person x wave), NaN where not applicable
% (a variable that does not exist for some wave, e.g. premiums in AHEAD 1993, gives a NaN column)
getw=@(stem,pref) cell2mat(arrayfun(@(w) DFJ2010_HRSextract_col(T,sprintf('%s%d%s',pref,w,stem)),waves,'UniformOutput',false));
iwstat=getw('IWSTAT','R');   % 1 alive&responded, 4 alive nonresp, 5 died this wave, 6 died earlier, 7 dropped, 0 inap
mstat=getw('MSTAT','R');     % 1 married, 2 married sp absent, 3 partnered, 4 sep, 5 div, 6 sep/div, 7 widowed, 8 never married
mpart=getw('MPART','R');     % 1 partnered
shlt=getw('SHLT','R');       % 1 excellent .. 5 poor
age=getw('AGEY_E','R');
iwy=getw('IWENDY','R');
iwm=getw('IWENDM','R');
oopmd=getw('OOPMD','R');     % wave 2: previous 12 months; waves 3+: previous 2 years
earn=getw('IEARN','R');
isret=getw('ISRET','R');
issdi=getw('ISSDI','R');
ipena=getw('IPENA','R');
igxfr=getw('IGXFR','R');
iunwc=getw('IUNWC','R');
wtresp=getw('WTRESP','R');
wthh=getw('WTHH','R');
govmd=getw('GOVMD','R');
higov=getw('HIGOV','R');
nrshom=getw('NRSHOM','R');
hosp=getw('HOSP','R');
% insurance premiums (monthly, nominal): private plans 1-3 and Medicare/Medicaid HMO; not asked in AHEAD 1993 (wave 2)
prem=sum(cat(3,getw('PRPRM1','R'),getw('PRPRM2','R'),getw('PRPRM3','R'),getw('MRPREM','R')),3,'omitnan');
atotb=getw('ATOTB','H');     % household net worth incl. secondary residence
itot=getw('ITOT','H');
icap=getw('ICAP','H');
hhres=getw('HHRES','H');

alive=(iwstat==1);           % responded, alive
living=(iwstat==1 | iwstat==4); % alive, whether or not they responded this wave
diedthiswave=(iwstat==5);
everresp=any(alive,2);

%% Sample selection (Appendix B.1)
fprintf('\nSample selection:\n')
fprintf('   RAND file: %d persons\n',N0)
keep=(T.HACOHORT==1 | T.HACOHORT==0); % AHEAD entry cohort (born 1923 or earlier) and their spouses; 0 is the small HRS/AHEAD overlap sample
fprintf('   AHEAD cohort (HACOHORT=1 or 0): %d (codebook: 8337+108)\n',sum(keep))
keep=keep & alive(:,1); % interviewed in AHEAD 1993 (DFJ "wave 1")
fprintf('   ...and interviewed in 1993: %d (paper: 8,222 interviewed in AHEAD 1993)\n',sum(keep))
% Single throughout: never married/partnered in any wave in which they are observed alive.
% DFJ: "We drop all individuals who were either married or cohabiting at any point in the
% analysis (so we include individuals who were never married with those who were divorced or
% widowed by wave 1)".
evermarried=any(alive & (mstat==1 | mstat==2 | mstat==3 | mpart==1),2);
evermarried=evermarried | (T.REMSTAT==1 | T.REMSTAT==2 | T.REMSTAT==3 | T.REMPART==1); % married at death per exit interview
keep=keep & ~evermarried;
fprintf('   ...and single in every observed wave (DFJ: 3,498): %d\n',sum(keep))
% Retired: drop missing 1993 labor income, and anyone with labor income > $3,000 in any wave
earn1998=earn; % deflate wave-by-wave below; earnings refer to the previous calendar year
for ww=1:nW
    cpi=nan(N0,1);
    yy=iwy(:,ww)-1; % income questions refer to the calendar year before the interview
    ok=~isnan(yy);
    [~,loc]=ismember(yy(ok),CPIyear);
    cpi(ok)=CPIval(max(loc,1));
    earn1998(:,ww)=earn(:,ww).*(CPI1998./cpi);
end
missearn1=alive(:,1) & isnan(earn(:,1));
keep=keep & ~missearn1;
fprintf('   ...and 1993 labor income observed: %d\n',sum(keep))
highearn=any(alive & earn1998>3000,2);
keep=keep & ~highearn;
fprintf('   ...and labor income <= $3,000 (1998$) in every wave (DFJ: 3,259; 592 men, 2,667 women): %d\n',sum(keep))
fprintf('      men %d, women %d; alive in 2006: %d responding, %d living (DFJ: 884)\n',sum(keep&T.RAGENDER==1),sum(keep&T.RAGENDER==2),sum(keep&alive(:,end)),sum(keep&living(:,end)))

T=T(keep,:);
iwstat=iwstat(keep,:); mstat=mstat(keep,:); shlt=shlt(keep,:); age=age(keep,:); iwy=iwy(keep,:); iwm=iwm(keep,:);
oopmd=oopmd(keep,:); earn=earn(keep,:); isret=isret(keep,:); issdi=issdi(keep,:); ipena=ipena(keep,:); igxfr=igxfr(keep,:); iunwc=iunwc(keep,:);
wtresp=wtresp(keep,:); wthh=wthh(keep,:); govmd=govmd(keep,:); higov=higov(keep,:); nrshom=nrshom(keep,:); hosp=hosp(keep,:); prem=prem(keep,:);
atotb=atotb(keep,:); itot=itot(keep,:); icap=icap(keep,:); hhres=hhres(keep,:);
alive=alive(keep,:); living=living(keep,:); diedthiswave=diedthiswave(keep,:);
N=height(T);

%% Deflators per person-wave: DFJ use the CPI of the year before the interview for everything
cpi_iw=nan(N,nW); cpi_inc=nan(N,nW);
for ww=1:nW
    yy=iwy(:,ww)-1; ok=~isnan(yy);
    [~,loc]=ismember(yy(ok),CPIyear); cpi_inc(ok,ww)=CPIval(max(loc,1));
end
cpi_iw=cpi_inc; % assets and medical expenses: same deflator as income (matches dataprep2.dta)
% people not interviewed in a wave have no interview year; fill with the wave's modal year so
% that exit-interview amounts (which sit in the wave of death) can still be deflated
for ww=1:nW
    modeyear=mode(iwy(~isnan(iwy(:,ww)),ww))-1;
    cpi_iw(isnan(cpi_iw(:,ww)),ww)=CPIval(CPIyear==modeyear);
end

%% Variables in DFJ form (1998 dollars)
female=(T.RAGENDER==2);
% health: bad = fair or poor (DFJ 'heal'=1)
badhealth=double(shlt>=4); badhealth(isnan(shlt))=NaN; badhealth(~alive)=NaN;
% net worth (DFJ 'assets'); RAND HwATOTB = all assets incl. secondary residence less all debts
assets=atotb.*(CPI1998./cpi_iw); assets(~alive)=NaN;
% non-asset income (DFJ 'anninc'): SS retirement + SSI/SSDI + pensions & annuities + other gov transfers (+ unemp/workers comp)
% (components that are not asked in a wave come through as NaN, so sum omitting NaN; treat the
%  wave as unobserved only if Social Security retirement income itself is missing)
nonassetinc=sum(cat(3,isret,issdi,ipena,igxfr,iunwc),3,'omitnan').*(CPI1998./cpi_inc);
nonassetinc(isnan(isret))=NaN; nonassetinc(~alive)=NaN;
earnings=earn.*(CPI1998./cpi_inc); earnings(~alive)=NaN;
% out-of-pocket medical expenses, annualised (DFJ 'medcost' is average over the previous two years)
medex=oopmd.*(CPI1998./cpi_iw);
medex(:,1)=medex(:,1); % wave 2 (1993) already covers 12 months
medex(:,2:end)=medex(:,2:end)/2; % waves 3+ cover the previous two years
% DFJ's medical expenses include insurance premia (paper Section V), which RAND's OOPMD excludes:
% add 12 x the monthly premiums reported at the interview (annual, 1998$)
premiums=12*prem.*(CPI1998./cpi_iw);
medex=medex+premiums;
medex(~alive)=NaN;
% last year of life: exit interview total, placed in the wave of death, annualised by time since previous interview
exitoop=T.REOOPMD.*NaN;
exitmonths=T.REOOPMD.*NaN;
for ii=1:N
    wd=find(diedthiswave(ii,:),1);
    if ~isempty(wd) && ~isnan(T.REOOPMD(ii))
        % months from previous interview to death
        if wd>1 && ~isnan(iwy(ii,wd-1))
            monthsprev=12*(T.RADYEAR(ii)-iwy(ii,wd-1))+(T.RADMONTH(ii)-iwm(ii,wd-1));
        else
            monthsprev=24;
        end
        if isnan(monthsprev) || monthsprev<6
            monthsprev=12;
        end
        monthsprev=min(monthsprev,36);
        exitmonths(ii)=monthsprev;
        % DFJ treat the exit-interview total as 'medical expenses during the last year of life' and do
        % not annualise it (checked against dataprep2.dta: annualising by months since the previous
        % interview gives death-wave medex about 0.58x theirs, the raw total matches)
        exitoop(ii)=T.REOOPMD(ii)*(CPI1998/cpi_iw(ii,wd));
        medex(ii,wd)=exitoop(ii);
    end
end
medexobs=~isnan(medex);
% DFJ bottom-code medical expenses at $250 (GAUSS minmedex=250; the 250s are visible in wlthmat13b.out)
medex(medexobs)=max(medex(medexobs),250);
% Medicaid / gov insurance
medicaid=double(govmd==1); medicaid(isnan(govmd))=NaN;

%% Permanent income: average non-asset income over the waves the person is observed, then percentile rank
% AHEAD 1993 income is incomplete (RAND doc: only SS, SSI and a few items asked), and DFJ's
% incprof.do sets 1994/1996 income to missing, so use 1995-2006 (RAND waves 3-8).
PIwaves=1:nW; % all waves (columns = RAND waves 2..8); see the income check printed below for how 1993 compares
avginc=mean(nonassetinc(:,PIwaves),2,'omitnan');
avginc(all(isnan(nonassetinc(:,PIwaves)),2))=NaN;
% percentile rank in (0,1], ties get the average rank
[~,ord]=sort(avginc);
PI=nan(N,1);
rk=tiedrank(avginc(~isnan(avginc)));
PI(~isnan(avginc))=rk/sum(~isnan(avginc));
fprintf('\nPermanent income: %d of %d have an income observation (others get PI=NaN)\n',sum(~isnan(PI)),N)

%% Ages by DFJ calendar year (DFJ index ages at 1996 so 70..102 in the model)
% age at each interview; for people not interviewed, age is still needed to place them (e.g.
% year of death), so fill forward from birth year and DFJ calendar year
ageDFJ=age;
for ww=1:nW
    miss=isnan(ageDFJ(:,ww));
    ageDFJ(miss,ww)=DFJyear(ww)-T.RABYEAR(miss);
end

%% Cohort in DFJ's terms: from the 1996 age (GAUSS cohorts_a=69|74|79|84|89|102 applied to agedat(:,1)=1996 age)
age1996=ageDFJ(:,2);
cohortedges=[69,74,79,84,89,102]; % cohort c: cohortedges(c) < age1996 <= cohortedges(c+1)
cohort=nan(N,1);
for c=1:5
    cohort(age1996>cohortedges(c) & age1996<=cohortedges(c+1))=c;
end

%% Assemble outputs
person=table(T.HHIDPN,T.HHID,T.PN,female,T.RABYEAR,T.RADYEAR,T.RADMONTH,PI,avginc,cohort,age1996,exitoop,T.REXITWV,...
    'VariableNames',{'HHIDPN','HHID','PN','female','birthyear','deathyear','deathmonth','PI','avgnonassetinc','cohort','age1996','exitoop','exitwave'});
% long table: one row per person-wave
[pp,ww]=ndgrid(1:N,1:nW);
pp=pp(:); ww=ww(:);
panel=table(T.HHIDPN(pp),waves(ww)',DFJyear(ww)',iwy(sub2ind([N,nW],pp,ww)),ageDFJ(sub2ind([N,nW],pp,ww)),female(pp),...
    alive(sub2ind([N,nW],pp,ww)),diedthiswave(sub2ind([N,nW],pp,ww)),iwstat(sub2ind([N,nW],pp,ww)),mstat(sub2ind([N,nW],pp,ww)),...
    badhealth(sub2ind([N,nW],pp,ww)),assets(sub2ind([N,nW],pp,ww)),nonassetinc(sub2ind([N,nW],pp,ww)),earnings(sub2ind([N,nW],pp,ww)),...
    medex(sub2ind([N,nW],pp,ww)),medicaid(sub2ind([N,nW],pp,ww)),wtresp(sub2ind([N,nW],pp,ww)),nrshom(sub2ind([N,nW],pp,ww)),hosp(sub2ind([N,nW],pp,ww)),...
    premiums(sub2ind([N,nW],pp,ww)),living(sub2ind([N,nW],pp,ww)),...
    isret(sub2ind([N,nW],pp,ww)),issdi(sub2ind([N,nW],pp,ww)),ipena(sub2ind([N,nW],pp,ww)),igxfr(sub2ind([N,nW],pp,ww)),iunwc(sub2ind([N,nW],pp,ww)),itot(sub2ind([N,nW],pp,ww)),icap(sub2ind([N,nW],pp,ww)),...
    'VariableNames',{'HHIDPN','RANDwave','DFJyear','interviewyear','age','female','alive','diedthiswave','iwstat','mstat',...
    'badhealth','assets','nonassetinc','earnings','medex','medicaid','weight','nursinghome','hospital','premiums','living',...
    'ISRET','ISSDI','IPENA','IGXFR','IUNWC','HITOT','HICAP'}); % last seven are nominal RAND components, for reference
panel=sortrows(panel,{'HHIDPN','RANDwave'});

%% Summary against DFJ's published numbers
fprintf('\nChecks against DFJ2010 (Appendix B / Section V):\n')
fprintf('   persons %d (DFJ 3,259), men %d (592), women %d (2,667), living 2006 %d (884)\n',N,sum(~female),sum(female),sum(living(:,end)))
nexit=sum(~isnan(exitoop));
fprintf('   exit-interview medex used for %d of %d deaths in 1996-2006; mean %.0f, median %.0f (1998$/yr)\n',nexit,sum(sum(diedthiswave(:,2:end))),mean(exitoop,'omitnan'),median(exitoop,'omitnan'))
for ww=1:nW
    a=alive(:,ww);
    fprintf('   %d: alive %4d | median assets %9.0f (n=%4d) | mean medex %7.0f (n=%4d) | mean nonasset inc %6.0f | bad health %.2f\n',...
        DFJyear(ww),sum(a),median(assets(a,ww),'omitnan'),sum(~isnan(assets(a,ww))),mean(medex(a,ww),'omitnan'),sum(~isnan(medex(a,ww))),...
        mean(nonassetinc(a,ww),'omitnan'),mean(badhealth(a,ww),'omitnan'))
end
fprintf('   paper: mean medical expenses $3,712 (sd $13,429) in 1998 dollars, sample pooled over 1996-2006\n')
m=medex(:,2:end); m=m(alive(:,2:end));
fprintf('   here : mean %.0f (sd %.0f), n=%d person-waves; of which premiums mean %.0f\n',mean(m,'omitnan'),std(m,'omitnan'),sum(~isnan(m)),mean(premiums(alive&[false(N,1),true(N,nW-1)]),'omitnan'))
fprintf('   cohorts by 1996 age: %s\n',mat2str(histcounts(cohort,0.5:1:5.5)))

%% Save
DFJyearlist=DFJyear; %#ok<NASGU>
save('./SavedOutput/DFJ2010_HRSpanel.mat','panel','person','waves','DFJyearlist','CPIyear','CPIval','-v7.3')
fprintf('\nSaved ./SavedOutput/DFJ2010_HRSpanel.mat (%d persons x %d waves)\n',N,nW)
