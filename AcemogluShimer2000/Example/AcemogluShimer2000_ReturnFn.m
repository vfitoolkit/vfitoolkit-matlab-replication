function F=AcemogluShimer2000_ReturnFn(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,r,Gamma,hbar,eta,theta)
% Acemoglu & Shimer (2000), eqs (3)-(7), normalized by productivity Gamma^t.
% semiz: 1=good job, 2=bad job, 3=unemployed insured, 4=unemployed uninsured
% h: search hours (only used when unemployed)
% s: only enters the SemiExoStateFn (job type sought / quit), not the return
% uincfloor: tiny income of the uninsured, and of the insured when uib=0 (1e-7), see AcemogluShimer2000_Example.m
%
% Budget: a'=((1+r)/Gamma)*(a+y-c), so c=a+y-a'*Gamma/(1+r)
%
% Employed hours do not affect transitions, so instead of putting them on a grid we
% use the static first-order condition, eq (6): eta*c/(hbar-h)=(1-tau)*w. Substituting
% c=a+(1-tau)*w*h-a'*Gamma/(1+r) gives a closed form for h given (a,a').
% (The h decision is then irrelevant for employed workers; it is a harmless tie.)

F=-Inf;

w=0;
if semiz==1
    w=wg;
elseif semiz==2
    w=wb;
end

if semiz<=2
    at=(1-tau)*w; % after-tax wage
    hours=(hbar-eta*(a-aprime*Gamma/(1+r))/at)/(1+eta);
    hours=min(max(hours,0),hbar);
    c=a+at*hours-aprime*Gamma/(1+r);
else
    hours=h;
    if semiz==3
        c=a+max(uib,uincfloor)-aprime*Gamma/(1+r); % floor only bites when uib=0 (No UI)
    else
        c=a+uincfloor-aprime*Gamma/(1+r);
    end
end

leisure=hbar-hours;

if c>0 && leisure>0
    if theta==1
        F=log(c)+eta*log(leisure);
    else
        F=((c*leisure^eta)^(1-theta))/(1-theta);
    end
end

end
