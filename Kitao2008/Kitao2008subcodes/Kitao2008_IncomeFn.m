function I=Kitao2008_IncomeFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)
% Taxable income of the household.
% Note: the tax liability itself is discarded here, only the income I is returned. But the tax
% parameters still have to be passed through rather than stubbed out, because under taxincome=2 the
% entrepreneur's capital choice k is made after tax (see (iii) at the top of
% Kitao2008_StaticEntrepreneurProblem), so income does depend on the tax system through k.
% tau_E1 is passed as 0 because it enters only the tax liability, never the income.

[TaxI_unused,TaxI_nonlinearpart_unused,I]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,0);

end
