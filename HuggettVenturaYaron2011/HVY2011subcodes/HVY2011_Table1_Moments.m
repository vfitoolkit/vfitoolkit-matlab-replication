function [J,gn,VarMatrix]=HVY2011_Table1_Moments(A,U,panelyears,W,Tl)
% GMM objective, moments and variance of moments for HVY2011 Table 1 (the estimation of human capital shocks).
% Follows the authors' code estimate_HVY_ALLT.m (data/hvy_Table).
%
% A=[sigma; sigma_u]: std dev of the human capital shock, and of measurement error in log wages
% U{pp}: [N,4] residual log hourly wages of panel year panelyears(pp), for the four years the panel follows people
%        (NaN for the years without data: after 1993 the PSID is biennial)
% W: weighting matrix, Tl: scale of the objective (the smallest number of observations in a panel year)
%
% With log wage differences y_n over n years, the model implies (HVY2011 Section II.B)
%   var(y_n)=n*sigma^2+2*sigma_u^2,  cov(y_n,y_m)=m*sigma^2+sigma_u^2 for m<n.
% For panel years to 1993 the moments are var(y_1), cov(y_2,y_1), cov(y_3,y_1), var(y_2), var(y_3); later panel years use
% whichever of these the available years allow. The moments are ordered as in the authors' code (all panel years to 1993
% moment by moment, then the later panel years), and the variance of the moments is block diagonal by panel year (with the
% off-diagonal terms within a panel year, except for 1995, as in the authors' code).

sigma2=A(1)^2;
sigmau2=A(2)^2;

early=find(panelyears<=93);
late=find(panelyears>93);
npanels=length(early);

% Number of moments
num_mom=5*npanels;
for pp=late
    if panelyears(pp)==94 || panelyears(pp)==95
        num_mom=num_mom+2;
    else
        num_mom=num_mom+1;
    end
end
gn=zeros(num_mom,1);
VarMatrix=zeros(num_mom,num_mom);

% Panel years to 1993: all four years
for jj=1:npanels
    u=U{early(jj)};
    d1=u(:,2)-u(:,1); d1=d1-mean(d1);
    d2=u(:,3)-u(:,1); d2=d2-mean(d2);
    d3=u(:,4)-u(:,1); d3=d3-mean(d3);
    F=[d1.^2-(sigma2+2*sigmau2), d2.*d1-(sigma2+sigmau2), d3.*d1-(sigma2+sigmau2), d2.^2-(2*sigma2+2*sigmau2), d3.^2-(3*sigma2+2*sigmau2)];
    momindex=jj+npanels*(0:4);
    gn(momindex)=mean(F,1)';
    Fd=F-mean(F,1);
    VarMatrix(momindex,momindex)=(Fd'*Fd)/size(F,1);
end

% Later panel years
kk=5*npanels;
for pp=late
    u=U{pp};
    if panelyears(pp)==94 % years 1, 2, 3
        d1=u(:,2)-u(:,1); d1=d1-mean(d1);
        d2=u(:,3)-u(:,1); d2=d2-mean(d2);
        F=[d1.^2-(sigma2+2*sigmau2), d2.*d1-(sigma2+sigmau2)];
        momindex=kk+(1:2);
        gn(momindex)=mean(F,1)';
        Fd=F-mean(F,1);
        VarMatrix(momindex,momindex)=(Fd'*Fd)/size(F,1);
        kk=kk+2;
    elseif panelyears(pp)==95 % years 1, 2, 4
        d1=u(:,2)-u(:,1); d1=d1-mean(d1);
        d3=u(:,4)-u(:,1); d3=d3-mean(d3);
        F=[d1.^2-(sigma2+2*sigmau2), d3.*d1-(sigma2+sigmau2)];
        momindex=kk+(1:2);
        gn(momindex)=mean(F,1)';
        Fd=F-mean(F,1);
        VarMatrix(momindex,momindex)=diag(mean(Fd.^2,1)); % no off-diagonal terms for 1995, as in the authors' code
        kk=kk+2;
    elseif panelyears(pp)==97 || panelyears(pp)==99 % years 2 and 4
        d=u(:,4)-u(:,2); d=d-mean(d);
        F=d.^2-(2*sigma2+2*sigmau2);
        gn(kk+1)=mean(F);
        VarMatrix(kk+1,kk+1)=mean((F-mean(F)).^2);
        kk=kk+1;
    else % 1996, 1998, 2000: years 1 and 3
        d=u(:,3)-u(:,1); d=d-mean(d);
        F=d.^2-(2*sigma2+2*sigmau2);
        gn(kk+1)=mean(F);
        VarMatrix(kk+1,kk+1)=mean((F-mean(F)).^2);
        kk=kk+1;
    end
end

J=Tl*gn'*W*gn;

end
