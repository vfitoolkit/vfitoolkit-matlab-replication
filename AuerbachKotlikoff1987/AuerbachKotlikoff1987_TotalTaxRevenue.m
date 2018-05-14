function TotalRevenue=AuerbachKotlikoff1987_TotalTaxRevenue(l,aprime,a,z,wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w,tau_k)
% z is unused as model is deterministic

%c=(1+rt)*a+wt*ej*(1-l)-aprime; % in absence of taxes

% Revenue from earnings tax:
Revenue_w=tau_w*wt*ej*(1-l);

% Revenue from capital income tax
Revenue_capital=tau_k*(1/(1-tau_k)*rt)*a; %rt is after-capital-tax (but before income tax). So 1/(1-tau_k)*rt give the before-capital-tax rate of return on capital

TaxableIncome=rt*a+(1-tau_w)*wt*ej*(1-l);
% Revenue from income tax:
Revenue_income=(tau_phi+tau_pi*TaxableIncome/2)*TaxableIncome; % ((tau_phi+tau_pi*TaxableIncome/2) is average tax rate on taxable income

c=(1/(1+tau_c))*((1-(tau_phi+tau_pi*TaxableIncome/2))*TaxableIncome+a-aprime);
% Revenue from consumption tax:
Revenue_c=tau_c*c;

TotalRevenue=Revenue_w+Revenue_capital+Revenue_income+Revenue_c;

end