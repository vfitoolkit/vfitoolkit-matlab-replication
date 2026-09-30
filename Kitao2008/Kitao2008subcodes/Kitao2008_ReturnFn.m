function F=Kitao2008_ReturnFn(aprime,eprime,a,e,eta,theta,sigma,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)

F=-Inf;

% Taxable income and the income tax on it. For entrepreneurs this internally solves the static
% production problem (Kitao2008_StaticEntrepreneurProblem), for workers it is just w*eta+r*a.
[TaxI,TaxI_nonlinearpart_unused,I]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,tau_E1);

% Once we have I and TaxI, the budget constraint is the same for workers and entrepreneurs:
%   workers:       (1+tau_c)c+aprime = w*eta+(1+r)a-TaxI = I+a-TaxI   [eqns (3) and (4) of Kitao (2008)]
%   entrepreneurs: (1+tau_c)c+aprime = profit            = I+a-TaxI   [eqns (6), (7) and (8) of Kitao (2008)]
% For the entrepreneur: eqn (7) has +(1-delta)k-(1+rbar)(k-a); we already have -delta*k-rbar*(k-a)
% in I (eqn 8), so we are left with +k-(k-a), which becomes +a
c=(I+a-TaxI-aprime)/(1+tau_c);

if c>0
    F=(c^(1-sigma))/(1-sigma); % CES utility fn
end

if aprime<0
    F=-Inf; % Borrowing constraint
end

end
