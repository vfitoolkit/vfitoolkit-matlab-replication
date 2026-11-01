function hours=AcemogluShimer2000_Hours(h,aprime,a,semiz,wg,wb,tau,r,Gamma,hbar,eta)
% Actual hours (work hours when employed, search hours when unemployed).
% Employed hours come from the static FOC, exactly as in AcemogluShimer2000_ReturnFn.
% Used by the FnsToEvaluate.

if semiz<=2
    w=wb;
    if semiz==1
        w=wg;
    end
    at=(1-tau)*w;
    hours=(hbar-eta*(a-aprime*Gamma/(1+r))/at)/(1+eta);
    hours=min(max(hours,0),hbar);
else
    hours=h;
end

end
