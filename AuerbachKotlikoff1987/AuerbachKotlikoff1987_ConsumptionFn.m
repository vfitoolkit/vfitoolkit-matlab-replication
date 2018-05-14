function c=AuerbachKotlikoff1987_ConsumptionFn(l,aprime,a,z,wt,rt,ej,tau_phi,tau_pi,tau_c,tau_w)
% z is unused as model is deterministic
% F=-Inf;

%c=(1+rt)*a+wt*ej*(1-l)-aprime; % in absence of taxes
TaxableIncome=rt*a+(1-tau_w)*wt*ej*(1-l);
c=(1/(1+tau_c))*((1-(tau_phi+tau_pi*TaxableIncome/2))*TaxableIncome+a-aprime);
% ((tau_phi+tau_pi*TaxableIncome/2) is average tax rate on taxable income
% 
% if c>0
%     Finner=(c^(1-1/rho)+alpha*l^(1-1/rho))^(1/(1-1/rho));
%     F=(Finner^(1-1/gamma))/(1-1/gamma);
% end
% 
% if jj==55 && aprime<0
%     F=-Inf; % Impose a_J>=0
% end

end