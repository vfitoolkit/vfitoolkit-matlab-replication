% DFJ2010_Targets: median-asset targets, moment matrix and covariance, mirroring the DFJ2010 GAUSS code
%
% This is a script that operates on arrays already in the workspace (so that DFJ2010_DataWork can
% run it both on the rebuilt HRS panel and on DFJ's own wlthmat13b.out and compare):
%   A_in      N x nY  assets by person and moment-year (1996,1998,...,2006), already top-coded
%   obs_in    N x nY  logical: alive and assets observed (GAUSS obsdat)
%   AGE_in    N x nY  age at each moment-year (NaN where not observed)
%   cohort_in N x 1   cohort 1..5 (0 = not in any cohort), from getchrt()
%   PIq_in    N x 1   permanent-income quintile 1..5 (0 = missing), from PIquant2()
%   years_in  1 x nY  calendar years of the columns
% It produces:
%   MedA      5 x 5 x nY  median assets by cohort x PI quintile x year among those alive (GAUSS simqunts, comptype=0)
%   Cnt       5 x 5 x nY  cell counts
%   AvgAge    5 x 5 x nY  mean age in the cell
%   MedA_surv 5 x 5 x nY  same medians but only for people alive in the first AND last year (comptype=1)
%   MomentList table: cohort, PIq, year, avgage, median, n  (one row per moment condition)
%   mmtmtx    N x J     individual-level moment functions (1{a<=median}-0.5)*1{in cell}, GAUSS makemmts_q
%   VCV       J x J     covariance of the data moments, GAUSS compvinv: demeaned mmtmtx'*mmtmtx/N
%   W         J x J     optimal weighting matrix inv(VCV) (GAUSS optwgts=2)
%   Wdiag     J x J     diagonal weighting matrix (optwgts=1)
%   VCV_medians J x J   bootstrap covariance of the cell medians (resampling persons), SE_medians its sqrt diag;
%             MomentList.se holds the same standard errors
%
% GAUSS rules mirrored (singlsav15_wave8.gau + gmmcrit25.g):
%   quantmin=10   fewer observations than this in a cell-year: no median for that year and the
%                 loop over years is abandoned for that cell (simqunts 'break')
%   cellmin=10    a moment is dropped if its cell-year has fewer observations, and so are all later
%                 years for that cell (makemmts_q)
%   keepinit=0    the first year (1996) is never a moment: it is the initial condition
%   The median is taken among survivors in that year (the unbalanced panel), see paper Sec. IV.B.

quantmin=10;
cellmin=10;
N=size(A_in,1);
nY=size(A_in,2);

MedA=nan(5,5,nY);
Cnt=zeros(5,5,nY);
AvgAge=nan(5,5,nY);
MedA_surv=nan(5,5,nY);
mmtmtx=zeros(N,0);
MomentList=zeros(0,6);
survivor=obs_in(:,1) & obs_in(:,end); % comptype=1: alive at first and last moment-year

for c=1:5
    for q=1:5
        ind0=(cohort_in==c & PIq_in==q);
        % medians among those alive in each year (simqunts, comptype=0)
        for y=1:nY
            ind=ind0 & obs_in(:,y);
            n=sum(ind);
            Cnt(c,q,y)=n;
            if n<quantmin
                break
            end
            MedA(c,q,y)=median(A_in(ind,y));
            AvgAge(c,q,y)=mean(AGE_in(ind,y));
        end
        % medians among survivors (comptype=1), for the paper's Figure 2 style comparison
        for y=1:nY
            ind=ind0 & obs_in(:,y) & survivor;
            if sum(ind)<quantmin
                break
            end
            MedA_surv(c,q,y)=median(A_in(ind,y));
        end
        % moment conditions (makemmts_q)
        skiprest=false;
        for y=1:nY
            ind=ind0 & obs_in(:,y);
            n=sum(ind);
            if n<cellmin
                skiprest=true;
            end
            if skiprest || y==1 || isnan(MedA(c,q,y))
                continue
            end
            m=((A_in(:,y)<=MedA(c,q,y))-0.5).*ind;
            m(~ind)=0;
            mmtmtx=[mmtmtx,m]; %#ok<AGROW>
            MomentList=[MomentList; c,q,years_in(y),AvgAge(c,q,y),MedA(c,q,y),n]; %#ok<AGROW>
        end
    end
end
MomentList=array2table(MomentList,'VariableNames',{'cohort','PIq','year','avgage','median','n'});
J=size(mmtmtx,2);

% covariance of the data moments (compvinv)
mm=mmtmtx-mean(mmtmtx,1);
VCV=(mm'*mm)/N;
W=inv(VCV);
Wdiag=diag(1./diag(VCV));

% Bootstrap covariance of the cell MEDIANS themselves (resampling persons), for use with the
% toolkit's estimator, whose moments are the median levels rather than DFJ's indicator moments
% (DFJ's VCV above is of 1{a<=median}-0.5; converting it to medians needs the density at the
% median, GAUSS getpdf/Koning kernel; the bootstrap avoids that).
Bboot=500;
rng(20100118,'twister'); % fixed seed so the weighting matrix is reproducible
MedBoot=nan(Bboot,J);
cellc=MomentList.cohort; cellq=MomentList.PIq; celly=zeros(J,1);
for jj=1:J
    celly(jj)=find(years_in==MomentList.year(jj));
end
for bb=1:Bboot
    draw=randi(N,N,1);
    Ab=A_in(draw,:); obsb=obs_in(draw,:); cohb=cohort_in(draw); PIb=PIq_in(draw);
    for jj=1:J
        ind=(cohb==cellc(jj) & PIb==cellq(jj) & obsb(:,celly(jj)));
        if sum(ind)>0
            MedBoot(bb,jj)=median(Ab(ind,celly(jj)));
        end
    end
end
VCV_medians=cov(MedBoot,'omitrows');
SE_medians=sqrt(diag(VCV_medians));
MomentList.se=SE_medians;

fprintf('DFJ2010_Targets: %d persons, %d moment conditions (DFJ: 101); cells with a 1996 median: %d of 25\n',N,J,sum(~isnan(MedA(:,:,1)),'all'))
fprintf('   bootstrap (B=%d) std errors of the cell medians: median %.0f, range [%.0f, %.0f]\n',Bboot,median(SE_medians),min(SE_medians),max(SE_medians))
