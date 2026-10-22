%% HVY2011 Table 1: estimation of human capital shocks from the wages of older male workers
% HVY2011 Section II.B. Late in the working life human capital production is negligible, so log wage differences are due to
% rental rate differences, the human capital shocks z~N(mu,sigma^2), and measurement error with std dev sigma_u (sigma_epsilon
% in the paper). The variances and covariances of 1, 2 and 3 year log wage differences identify (sigma, sigma_u).
%
% This follows the authors' code (data/hvy_Table: masterfileHVY_EST.m, FirstStageRegressionDataPRP.m,
% FirstStageRegression1_EXPN.m, estimate_HVY_ALLT.m) line by line, using the PSID extracts that come with it
% (data/hvy_Table/out69.txt,...,out100.txt; each panel year follows people for four years; after 1993 the PSID is biennial, so
% the later panel years have two or three of the four years). The wage is real labor earnings divided by hours (hourly wage).
% Selection (as in the authors' code): male household heads in the age range in the first year, not SEO, earnings in each
% year at least emin (1968 dollars; $1,000 if age 30 or below) and at most $1.8M, at least 520 hours at the federal minimum
% wage, hours between 520 (260 if age 30 or below) and 5820, and growth of the hourly wage between years between 1/20 and 20.
% First stage: log hourly wages are regressed on a constant, education, and a cubic in age (pooled), and the residuals used.
% Estimation: two-step GMM (identity weighting matrix, then the inverse of the variance of the moments), with fminsearch.
% Standard errors: sqrt(diag(inv(D'*inv(S)*D)/Tl)) with D the Jacobian of the moments and Tl the smallest number of
% observations in a panel year (the authors' stderrv=1), with D by central differences (the authors use fsolve's Jacobian,
% which is the same here). The point estimates reproduce the authors' exactly, but the standard errors are larger than theirs
% (by a factor of about 1.1 to 1.5, more so the more post-1993 panel years are used), for a reason we have not found.

if ~exist('./SavedOutput','dir'); mkdir('./SavedOutput'); end
datadirT1='./data/hvy_Table/';

% CPI (equal to 1 in 1969) and the federal minimum wage (1969 onwards)
cpi_defT1=readmatrix([datadirT1,'CPI_BLS.xlw'],'FileType','spreadsheet','Range','D5:D41'); % 1969-2005
minwageT1=readmatrix([datadirT1,'MinimumWage.xlsx'],'Range','B1:B35'); % 1969-2003

% Read the PSID extracts once
% Column order (authors' readme): famid68(1) persnr sex(3) educ(4) educftr seq1(6) seq2(7) seq3 seq4 age1(10) ... hrs1(18) hrs2 hrs3
% hrs4(21) inc1(22)..inc4(25) ... wge1(34) wge2(35) wge3(36) wge4(37) ...
allpanelyears=69:1:100;
PSIDT1=cell(1,length(allpanelyears));
for pp=1:length(allpanelyears)
    PSIDT1{pp}=readmatrix([datadirT1,'out',num2str(allpanelyears(pp)),'.txt'],'FileType','text');
end

% The rows of Table 1: [minage, maxage, first panel year, last panel year, emin]
% (the periods are of panel years, i.e. the first of the four years a panel follows people: 1969-1981 and 1982-2004; this
% reproduces the median number of observations in HVY2011 Table 1)
Table1Specs=[55,65,69,100,2000;
             50,60,69,100,2000;
             55,65,69,100,3000;
             55,65,69,100,2000;
             55,65,69,100,1500;
             55,65,69,100,1000;
             55,65,69,81,2000;
             50,60,69,81,2000;
             55,65,82,100,2000;
             50,60,82,100,2000];
nspecs=size(Table1Specs,1);
Table1=struct();
Table1.specs=Table1Specs;
Table1.N=zeros(nspecs,1);
Table1.sigma=zeros(nspecs,1);
Table1.se_sigma=zeros(nspecs,1);
Table1.sigma_u=zeros(nspecs,1);
Table1.se_sigma_u=zeros(nspecs,1);

earnup=1800000; % maximum earnings (1968 dollars)
earnlbY=1000;   % minimum earnings, age 30 and below
upgr=20;        % maximum growth of the hourly wage between years
dngr=1/upgr;
minhours=520;
minhoursY=260;
maxhours=5820;

for ss=1:nspecs
    minage=Table1Specs(ss,1);
    maxage=Table1Specs(ss,2);
    panelyears=Table1Specs(ss,3):1:Table1Specs(ss,4);
    earnlb=Table1Specs(ss,5);
    npy=length(panelyears);

    % Selected samples: log real hourly wages (NaN for years without data), and the regressors of the first stage
    logwage=cell(1,npy);
    educ=cell(1,npy);
    age=cell(1,npy);
    Nselected=zeros(1,npy);
    for pp=1:npy
        jyear=panelyears(pp);
        temp=PSIDT1{jyear-68};
        jy=jyear-68; % index of the calendar year for the CPI and the minimum wage
        earningup=earnup*cpi_defT1(jy:jy+3)';
        earninglb=earnlb*cpi_defT1(jy:jy+3)';
        earninglbY=earnlbY*cpi_defT1(jy:jy+3)';
        agep=temp(:,10);
        % the years of the four that have data: before 1994 all four; 1994: 1,2,3; 1995: 1,2,4; 1996, 1998, 2000: 1,3; 1997, 1999: 2,4
        if jyear<=93
            yrs=[1,2,3,4];
        elseif jyear==94
            yrs=[1,2,3];
        elseif jyear==95
            yrs=[1,2,4];
        elseif jyear==97 || jyear==99
            yrs=[2,4];
        else
            yrs=[1,3];
        end
        select=(agep>=minage & agep<=maxage) & (temp(:,1)>=1 & temp(:,1)<=2930) & (temp(:,4)<=4) & (temp(:,3)==1);
        if jyear==96 || jyear==98 || jyear==100
            select=select & (temp(:,7)==1); % head of household (the sequence number variable moves for these years)
        else
            select=select & (temp(:,6)==1);
        end
        for yy=yrs
            wge=temp(:,33+yy);
            hrs=temp(:,17+yy);
            select=select & (wge~=999999) & (wge~=99999) & (wge~=9999999) & (wge>0) ...
                & (wge>minhours*minwageT1(jy+yy-1)) & (wge<=earningup(yy)) ...
                & ((agep<=30).*(wge>=earninglbY(yy))+(agep>30).*(wge>=earninglb(yy))) ...
                & (hrs<=maxhours & ((agep>30).*(hrs>=minhours)+(agep<=30).*(hrs>=minhoursY)));
        end
        for yy=2:length(yrs) % growth of the hourly wage between consecutive available years
            growth=(temp(:,33+yrs(yy))./temp(:,17+yrs(yy)))./(temp(:,33+yrs(yy-1))./temp(:,17+yrs(yy-1))+.001);
            select=select & (growth<upgr) & (growth>=dngr);
        end
        Nselected(pp)=sum(select);
        lw=nan(Nselected(pp),4);
        for yy=yrs
            lw(:,yy)=log((temp(select,33+yy)/cpi_defT1(jy+yy-1))./temp(select,17+yy)); % log real hourly wage
        end
        logwage{pp}=lw;
        educ{pp}=temp(select,4);
        age{pp}=agep(select);
    end

    % First stage: pooled regression of log hourly wages on a constant, education and a cubic in age (age in that year)
    Yfs=[];
    Xfs=[];
    for pp=1:npy
        for yy=1:4
            if ~isnan(logwage{pp}(1,yy))
                agey=age{pp}+yy-1;
                Yfs=[Yfs; logwage{pp}(:,yy)];
                Xfs=[Xfs; ones(Nselected(pp),1), educ{pp}, agey, agey.^2/100, agey.^3/1000];
            end
        end
    end
    bfs=(Xfs'*Xfs)\(Xfs'*Yfs);
    U=cell(1,npy);
    for pp=1:npy
        U{pp}=nan(Nselected(pp),4);
        for yy=1:4
            if ~isnan(logwage{pp}(1,yy))
                agey=age{pp}+yy-1;
                U{pp}(:,yy)=logwage{pp}(:,yy)-[ones(Nselected(pp),1), educ{pp}, agey, agey.^2/100, agey.^3/1000]*bfs;
            end
        end
    end
    Tl=min(Nselected);

    % Two-step GMM
    [~,gn0]=HVY2011_Table1_Moments([0.1;0.1],U,panelyears,1,Tl);
    W1=eye(length(gn0));
    opts=optimset('MaxIter',2900,'MaxFunEvals',2900,'TolFun',1e-10,'TolX',1e-10);
    A1=fminsearch(@(A) HVY2011_Table1_Moments(A,U,panelyears,W1,Tl),[0.1;0.1],opts);
    [~,~,VarMatrix1]=HVY2011_Table1_Moments(A1,U,panelyears,W1,Tl);
    W2=inv(VarMatrix1);
    opts=optimset('MaxIter',1900,'MaxFunEvals',1900,'TolFun',1e-10,'TolX',1e-10);
    A2=fminsearch(@(A) HVY2011_Table1_Moments(A,U,panelyears,W2,Tl),A1,opts);
    [~,gn2,Sw]=HVY2011_Table1_Moments(A2,U,panelyears,W2,Tl);
    % Jacobian of the moments, by central differences
    Dg=zeros(length(gn2),2);
    for kk=1:2
        hstep=zeros(2,1); hstep(kk)=1e-5;
        [~,gnplus]=HVY2011_Table1_Moments(A2+hstep,U,panelyears,W2,Tl);
        [~,gnminus]=HVY2011_Table1_Moments(A2-hstep,U,panelyears,W2,Tl);
        Dg(:,kk)=(gnplus-gnminus)/(2e-5);
    end
    VarCov=sqrt(diag(inv(Dg'*inv(Sw)*Dg)/Tl));

    Table1.N(ss)=median(Nselected);
    Table1.sigma(ss)=abs(A2(1));
    Table1.sigma_u(ss)=abs(A2(2));
    Table1.se_sigma(ss)=VarCov(1);
    Table1.se_sigma_u(ss)=VarCov(2);
    fprintf('Table 1: ages %i-%i, panel years 19%i-%i, emin %i: N %5.1f, sigma %6.4f (%6.4f), sigma_u %6.4f (%6.4f) \n',minage,maxage,Table1Specs(ss,3),Table1Specs(ss,4)+1900,earnlb,Table1.N(ss),Table1.sigma(ss),Table1.se_sigma(ss),Table1.sigma_u(ss),Table1.se_sigma_u(ss))
end
fprintf('Authors'' output for the first row (data/hvy_Table/outestimationHVY.txt): N 103, sigma 0.1112 (0.0070), sigma_u 0.1371 (0.0047) \n')

save ./SavedOutput/HVY2011_Table1.mat Table1
