% DFJ2010_MedexTargets: the medical-expense moments of the endogenous medical-expense model (DFJ2010 online appendix
% A.2, eqns A2-A5, and gmmcrit51.g): by cohort (5) x permanent-income HALF (2) x year (1998-2006), among those alive
% with medical expenses observed, using the annualised two-year out-of-pocket medical expenses bottom-coded at
% minmedex=250 (singlsav51_wave8.gau line 161): the mean (A3), the 90th percentile (A2), and the first and second
% autocorrelations of log medical expenses (A4, A5), in TWO definitions. DFJ2010 (gmmcrit51.g simcrrl, appendix eqns A4-A5):
% E[R_it R_i,t-k]/(sigma_t sigma_t-k) over the people observed at both dates, with R_it = ln m_it - (FULL cell mean of ln m
% at t) and sigma_t the full cell's standard deviation at t: rows 'ac1dfj', 'ac2dfj' (MedexCells.AC1dfj, AC2dfj). Toolkit
% definition (EvalFnOnAgentDist_AutoCorrTransProbs_FHorz): the ordinary correlation over the people observed at both dates,
% with that pair population's own means and standard deviations: rows 'ac1', 'ac2' (MedexCells.AC1, AC2). Each has its
% own bootstrap standard error (the resample's full cells for the DFJ definition). The estimation is run on each
% (DFJ2010_EndogMedex3: the toolkit definition through the AutoCorr targets, the DFJ definition through CustomModelStats).
% Cells (or cell pairs) with fewer than cellmin=10 observations are dropped.
%
% Script on workspace arrays (so DFJ2010_DataWork can run it on the HRS rebuild and on DFJ's own extract):
%   M_in      N x nY  annualised two-year out-of-pocket medical expenses by moment year (NaN if not observed)
%   mobs_in   N x nY  logical: alive and medex observed
%   PI_in     N x 1   permanent-income percentile rank in (0,1] (NaN if missing)
%   cohort_in N x 1   cohort 1..5 (0 none)
%   years_in  1 x nY
% Produces:
%   MedexList table: cohort, half, year, stat ('mean','p90','ac1','ac2'), value, n, se (bootstrap, B=500, same seed)
%   MedexCells: struct of 5 x 2 x nY arrays Mean, P90, Cnt and 5 x 2 x nY AC1, AC2 (NaN where dropped)
minmedex=250;
cellmin=10;
N=size(M_in,1);
nY=size(M_in,2);
Mb=max(M_in,minmedex); % bottom-coding (gmmcrit51 loadsims)
half_in=1+(PI_in>=0.5); % DFJ: quantile q of Q=2 is (q-1)/Q <= I < q/Q, so I=0.5 exactly is in the upper half
half_in(isnan(PI_in))=0;
Bboot=500;
rng(20100118,'twister');
draws=randi(N,N,Bboot); % the same bootstrap draws for every statistic (columns are resamples)
MedexCells.Mean=nan(5,2,nY); MedexCells.P90=nan(5,2,nY); MedexCells.Cnt=zeros(5,2,nY);
MedexCells.AC1=nan(5,2,nY); MedexCells.AC2=nan(5,2,nY); MedexCells.AC1dfj=nan(5,2,nY); MedexCells.AC2dfj=nan(5,2,nY); MedexCells.CntAC1=zeros(5,2,nY); MedexCells.CntAC2=zeros(5,2,nY);
MedexList=cell(0,7);
for c=1:5
    for hh=1:2
        ind0=(cohort_in==c & half_in==hh);
        lnM=log(Mb);
        % cell means and sds of log medex by year (for the residuals)
        cmean=nan(1,nY); csd=nan(1,nY);
        for y=1:nY
            ind=ind0 & mobs_in(:,y);
            if sum(ind)>=cellmin
                cmean(y)=mean(lnM(ind,y));
                csd(y)=sqrt(mean((lnM(ind,y)-cmean(y)).^2));
            end
        end
        for y=2:nY % 1996 is the initial condition, never a moment
            ind=ind0 & mobs_in(:,y);
            n=sum(ind);
            MedexCells.Cnt(c,hh,y)=n;
            if n<cellmin
                continue
            end
            MedexCells.Mean(c,hh,y)=mean(Mb(ind,y));
            MedexCells.P90(c,hh,y)=quantile(Mb(ind,y),0.9);
            % bootstrap standard errors
            bm=nan(1,Bboot); bp=nan(1,Bboot);
            for bb=1:Bboot
                d=draws(:,bb); indb=ind0(d) & mobs_in(d,y);
                if sum(indb)>0
                    bm(bb)=mean(Mb(d(indb),y)); bp(bb)=quantile(Mb(d(indb),y),0.9);
                end
            end
            MedexList(end+1,:)={c,hh,years_in(y),'mean',MedexCells.Mean(c,hh,y),n,std(bm,'omitnan')}; %#ok<SAGROW>
            MedexList(end+1,:)={c,hh,years_in(y),'p90',MedexCells.P90(c,hh,y),n,std(bp,'omitnan')}; %#ok<SAGROW>
            for k=1:2 % autocorrelations with year y-k
                if y-k<1 || isnan(cmean(y)) || isnan(cmean(y-k))
                    continue
                end
                indk=ind & mobs_in(:,y-k);
                nk=sum(indk);
                if k==1, MedexCells.CntAC1(c,hh,y)=nk; else, MedexCells.CntAC2(c,hh,y)=nk; end
                if nk<cellmin
                    continue
                end
                l1=lnM(indk,y); l0=lnM(indk,y-k); % the pair population: observed at both dates
                ac=mean(((l1-mean(l1))/std(l1,1)).*((l0-mean(l0))/std(l0,1))); % toolkit definition: the pair population's own means and sds
                acdfj=mean(((l1-cmean(y))/csd(y)).*((l0-cmean(y-k))/csd(y-k))); % DFJ2010 definition: the full cell's means and sds at each date
                bac=nan(1,Bboot); bacdfj=nan(1,Bboot);
                for bb=1:Bboot
                    d=draws(:,bb); indb=ind0(d) & mobs_in(d,y) & mobs_in(d,y-k);
                    if sum(indb)>=2
                        l1=lnM(d(indb),y); l0=lnM(d(indb),y-k);
                        bac(bb)=mean(((l1-mean(l1))/std(l1,1)).*((l0-mean(l0))/std(l0,1)));
                        f1=lnM(d(ind0(d) & mobs_in(d,y)),y); f0=lnM(d(ind0(d) & mobs_in(d,y-k)),y-k); % the resample's full cells at the two dates
                        bacdfj(bb)=mean(((l1-mean(f1))/std(f1,1)).*((l0-mean(f0))/std(f0,1)));
                    end
                end
                if k==1, MedexCells.AC1(c,hh,y)=ac; MedexCells.AC1dfj(c,hh,y)=acdfj; else, MedexCells.AC2(c,hh,y)=ac; MedexCells.AC2dfj(c,hh,y)=acdfj; end
                MedexList(end+1,:)={c,hh,years_in(y),['ac',num2str(k)],ac,nk,std(bac,'omitnan')}; %#ok<SAGROW>
                MedexList(end+1,:)={c,hh,years_in(y),['ac',num2str(k),'dfj'],acdfj,nk,std(bacdfj,'omitnan')}; %#ok<SAGROW>
            end
        end
    end
end
MedexList=cell2table(MedexList,'VariableNames',{'cohort','half','year','stat','value','n','se'});
fprintf('DFJ2010_MedexTargets: %d rows (%d means, %d 90th percentiles, %d first and %d second autocorrelations, each autocorrelation in two definitions; DFJ2010: 249-101=148 moments in total)\n',height(MedexList),sum(strcmp(MedexList.stat,'mean')),sum(strcmp(MedexList.stat,'p90')),sum(strcmp(MedexList.stat,'ac1')),sum(strcmp(MedexList.stat,'ac2')))
fprintf('   cohort 1 (74 in 1996), lower/upper PI half: mean medex 1998..2006 = %s / %s; 90th pct = %s / %s; AC1 = %s / %s\n',mat2str(round(squeeze(MedexCells.Mean(1,1,2:end))'/1000,1)),mat2str(round(squeeze(MedexCells.Mean(1,2,2:end))'/1000,1)),mat2str(round(squeeze(MedexCells.P90(1,1,2:end))'/1000,1)),mat2str(round(squeeze(MedexCells.P90(1,2,2:end))'/1000,1)),mat2str(round(squeeze(MedexCells.AC1(1,1,2:end))',2)),mat2str(round(squeeze(MedexCells.AC1(1,2,2:end))',2)))
