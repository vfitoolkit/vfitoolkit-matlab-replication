function F=HuggettVenturaYaron2011_ReturnFn(l,aprime,a,h,r,gamma,agej,Jr,J,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)
% l is the fraction of time working (1-l is the time learning, HVY2011 call it s)
% a is assets (HVY2011 call it k), h is human capital
%
% The rental rate of human capital grows at rate g; we follow a cohort, so the rental rate at age j
% is Rj=(1+g)^(j-1) (normalized to one at age 23, as in the authors' code). Earnings are Rj*h*l.
% Tax system (HVY2011 Table 2, and the authors' code hcapshocks-ge.f):
%   social security: payroll tax tau_ss on earnings before retirement; a transfer ssbenefit in retirement
%   income tax: average tax rate (in percent) taxw1*x-taxw2*x^2 on income (earnings plus interest) before
%   retirement, and taxr1*x-taxr2*x^2 on income (transfer plus interest) in retirement, where x is income
%   in tens of thousands of dollars (taxfactor converts model units to dollars, and dividing by Rj
%   indexes the tax schedule to the rental rate). No income tax on non-positive income.

F=-Inf;

Rj=(1+g)^(agej-1); % rental rate of human capital at age j

if agej<Jr % working age
    earnings=Rj*h*l;
    income=earnings+r*a;
    incometax=0;
    if income>0
        x=income*taxfactor/(Rj*10000);
        incometax=((taxw1*x-taxw2*x^2)/100)*income;
    end
    c=earnings*(1-tau_ss)+(1+r)*a-incometax-aprime;
else % retired (Jr is the first period of retirement)
    income=r*a+ssbenefit;
    incometax=0;
    if income>0
        x=income*taxfactor/(Rj*10000);
        incometax=((taxr1*x-taxr2*x^2)/100)*income;
    end
    c=ssbenefit+(1+r)*a-incometax-aprime;
end

if c>0
    F=(c^(1-gamma))/(1-gamma);
end

% Indifferent about l when retired
% Not a problem, but gives crazy policy for l during retirement
% Clean it up by making l=1 when retired (no time learning, so human capital just gets the shocks)
if agej>=Jr
    if l<1
        F=-Inf;
    end
end

% No debt at death: assets must be non-negative at the end of the last period (as in the authors' code). Together
% with the -Inf return to negative consumption, this imposes the natural borrowing limit at every age.
if agej==J
    if aprime<0
        F=-Inf;
    end
end

end
