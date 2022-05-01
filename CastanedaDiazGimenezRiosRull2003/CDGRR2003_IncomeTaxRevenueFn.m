function [IncomeTaxRevenue]=CDGRR2003_IncomeTaxRevenueFn(h,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3)

w=(1-theta)*(((r+delta)/(theta))^(theta/(theta-1)));

% Income
% If retired
e=0;
y=k*r+e*h*w+omega;
if s<=J %If working
    e=e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4);
    y=k*r+e*h*w;
end

%Income tax revenue
IncomeTaxRevenue=a0*(y-(y^(-a1)+a2)^(-1/a1))+a3*y;

%c=0;
%c=y+a_val-d1_val;

end