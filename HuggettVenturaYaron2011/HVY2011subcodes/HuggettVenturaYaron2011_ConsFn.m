function c=HuggettVenturaYaron2011_ConsFn(l,aprime,a,h,r,agej,Jr,g,tau_ss,ssbenefit,taxfactor,taxw1,taxw2,taxr1,taxr2)
% Consumption, from the budget constraint in HuggettVenturaYaron2011_ReturnFn (see there for the tax system)

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
else % retired
    income=r*a+ssbenefit;
    incometax=0;
    if income>0
        x=income*taxfactor/(Rj*10000);
        incometax=((taxr1*x-taxr2*x^2)/100)*income;
    end
    c=ssbenefit+(1+r)*a-incometax-aprime;
end

end
