function F=AS2000_ReturnFn(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,R,Gamma,hbar,eta,theta,ojs)
% Acemoglu & Shimer (2000), eqs (3)-(7), normalized by productivity Gamma^t.
% semiz: 1=good job, 2=bad job, 3=unemployed insured, 4=unemployed uninsured
% h: search hours (used when unemployed, and in bad jobs when there is on-the-job search, ojs=1)
% s: only enters the SemiExoStateFn (job type sought / quit), not the return
% uincfloor: tiny income of the uninsured, and of the insured when uib=0 (1e-7), see the note in AcemogluShimer2000.m
% ojs: =1 on-the-job search in bad jobs (Section 5.4, Table 6)
%
% Budget: a'=(R/Gamma)*(a+y-c), so c=a+y-a'*Gamma/R
%
% Employed work hours do not affect transitions, so instead of putting them on a grid we
% use the static first-order condition, eq (6): eta*c/(hbar-hs-h)=(1-tau)*w, where hs is search
% hours on the job (zero unless ojs=1 and in a bad job). Substituting c=a+(1-tau)*w*h-a'*Gamma/R
% gives a closed form for h given (a,a',hs).
% (In good jobs, and in bad jobs without ojs, the h decision is then irrelevant; it is a harmless tie.)

F=-Inf;

w=0;
if semiz==1
    w=wg;
elseif semiz==2
    w=wb;
end

hs=0; % search hours on the job
if semiz==2 && ojs==1
    hs=h;
end

if semiz<=2
    at=(1-tau)*w; % after-tax wage
    hours=(hbar-hs-eta*(a-aprime*Gamma/R)/at)/(1+eta);
    hours=min(max(hours,0),hbar-hs);
    c=a+at*hours-aprime*Gamma/R;
else
    hours=h;
    if semiz==3
        c=a+max(uib,uincfloor)-aprime*Gamma/R; % floor only bites when uib=0 (No UI)
    else
        c=a+uincfloor-aprime*Gamma/R;
    end
end

leisure=hbar-hours-hs;

if c>0 && leisure>0
    if theta==1
        F=log(c)+eta*log(leisure);
    else
        F=((c*leisure^eta)^(1-theta))/(1-theta);
    end
end

end
