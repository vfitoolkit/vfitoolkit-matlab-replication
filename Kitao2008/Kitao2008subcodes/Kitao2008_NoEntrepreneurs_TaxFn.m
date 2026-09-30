function TaxRevenue=Kitao2008_NoEntrepreneurs_TaxFn(aprime,a,eta,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k)
% Revenue of: Income tax plus consumption tax

TaxRevenue=0;
TaxI=0;

I=w*eta+r*a;
if taxincome==1 % This is baseline model where income is taxed
    if I>0
        TaxI=tau_a0*(I-(tau_a2+I^(-tau_a1))^(-1/tau_a1))+tau_I*I;
    end
elseif taxincome==2 % Seperate taxation of labor income and capital income
    Itilde=I-r*a; % Itilde=w*eta
    if Itilde>0
        TaxI=tau_k*r*a+tau_a0*(Itilde-(tau_a2+Itilde^(-tau_a1))^(-1/tau_a1))+tau_I*Itilde;
    end
% Note: taxincome=3 is not relevant to 'no entrepreneur' economy
end
c=(w*eta+(1+r)*a-TaxI-aprime)/(1+tau_c);

TaxRevenue=tau_c*c+TaxI;


end