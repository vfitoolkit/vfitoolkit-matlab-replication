function TaxI=Kitao2008_IncomeTaxFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1)
% Income tax revenue raised from this household (the whole income tax, so both the non-linear
% tau_a0*() part and the flat tau_I, tau_k and tau_E1 parts). Does not include the consumption tax.
% See Kitao2008_IncomeTaxFn_NonLinearPartOnly for just the non-linear part.

[TaxI,TaxI_nonlinearpart_unused,I_unused]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,tau_E1);

end
