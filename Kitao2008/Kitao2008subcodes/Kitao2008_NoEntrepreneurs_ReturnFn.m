function F=Kitao2008_NoEntrepreneurs_ReturnFn(aprime,a,eta,sigma,r,w,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k)
% Note: Essentially just copy-paste of main return function, then delete everything which is only relevant to entrepreneurs.

F=-Inf;

I=w*eta+r*a;
TaxI=0;
if taxincome==1 % This is baseline model where income is taxed
    if I>0
        TaxI=tau_a0*(I-(tau_a2+I^(-tau_a1))^(-1/tau_a1))+tau_I*I;
    end
elseif taxincome==2 % Seperate taxation of labor income and capital income
    Itilde=I-r*a;
    if Itilde>0
        TaxI=tau_k*r*a+tau_a0*(Itilde-(tau_a2+Itilde^(-tau_a1))^(-1/tau_a1))+tau_I*Itilde;
    end
% Note: taxincome=3 is not relevant to 'no entrepreneur' economy
end
c=(w*eta+(1+r)*a-TaxI-aprime)/(1+tau_c);

if c>0
    F=(c^(1-sigma))/(1-sigma); % CES utility fn
end

if aprime<0
    F=-Inf; % Borrowing constraint
end

end