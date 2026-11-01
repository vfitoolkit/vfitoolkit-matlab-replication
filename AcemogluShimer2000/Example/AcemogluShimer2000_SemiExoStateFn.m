function prob=AcemogluShimer2000_SemiExoStateFn(semiz,semizprime,h,s,deltag,deltab,mg,mb,phi)
% Transition probabilities of the employment state, which depend on the decisions (h,s).
% semiz: 1=good job, 2=bad job, 3=unemployed insured, 4=unemployed uninsured
%
% Employed (semiz=1,2):  s=1 quit (next period uninsured unemployed), s=0 stay
%   stay: job ends with prob delta, worker becomes insured unemployed
% Unemployed (semiz=3,4): s=1 seek a good job, s=0 seek a bad job
%   find a job with prob m*h (m=mg or mb)
%   insured who do not find a job lose eligibility with prob phi

prob=0;
delta=0;
m=0;
jobstate=0;

if semiz<=2
    if semiz==1
        delta=deltag;
    else
        delta=deltab;
    end
    if s==1 % quit
        if semizprime==4
            prob=1;
        end
    else
        if semizprime==semiz
            prob=1-delta;
        elseif semizprime==3
            prob=delta;
        end
    end
else
    if s==1
        m=mg;
        jobstate=1;
    else
        m=mb;
        jobstate=2;
    end
    p=m*h;
    if semizprime==jobstate
        prob=p;
    elseif semiz==3
        if semizprime==3
            prob=(1-p)*(1-phi);
        elseif semizprime==4
            prob=(1-p)*phi;
        end
    else % semiz==4
        if semizprime==4
            prob=1-p;
        end
    end
end

end
