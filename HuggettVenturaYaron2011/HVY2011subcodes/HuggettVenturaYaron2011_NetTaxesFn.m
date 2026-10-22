function T=HuggettVenturaYaron2011_NetTaxesFn(l,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)
% Net taxes paid (payroll tax plus income tax, less the retirement transfer), divided by the
% rental rate Rj=(1+g)^(j-1) so that summing across ages gives the cross-section at a point in
% time (HVY2011 Appendix A, step 3). HVY2011 require total net taxes to be positive.
% (tax system as in HuggettVenturaYaron2011_ReturnFn)

Rj=(1+g)^(agej-1); % rental rate of human capital at age j

if agej<Jr % working age
    earnings=Rj*h*l;
    income=earnings+r*a;
    incometax=0;
    if income>0
        x=income*taxfactor/(Rj*10000);
        incometax=((taxw1*x-taxw2*x^2)/100)*income;
    end
    T=(tau_ss*earnings+incometax)/Rj;
else % retired
    income=r*a+ssbenefit;
    incometax=0;
    if income>0
        x=income*taxfactor/(Rj*10000);
        incometax=((taxr1*x-taxr2*x^2)/100)*income;
    end
    T=(incometax-ssbenefit)/Rj;
end

end
