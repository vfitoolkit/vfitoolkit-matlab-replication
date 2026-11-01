function hours=AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)
% Actual hours (work hours when employed, search hours when unemployed).
% Employed hours come from the static FOC, exactly as in AS2000_ReturnFn.
% With on-the-job search (ojs=1), h is search hours in bad jobs, and this returns the work hours.
% Used by the FnsToEvaluate.

if semiz<=2
    w=wb;
    if semiz==1
        w=wg;
    end
    hs=0;
    if semiz==2 && ojs==1
        hs=h;
    end
    at=(1-tau)*w;
    hours=(hbar-hs-eta*(a-aprime*Gamma/R)/at)/(1+eta);
    hours=min(max(hours,0),hbar-hs);
else
    hours=h;
end

end
