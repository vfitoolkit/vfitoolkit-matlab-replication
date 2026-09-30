function TaxRevenue=Kitao2008_TaxFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, tau_c, taxincome, tau_k, tau_E1)
% Revenue of: Income tax plus consumption tax

[TaxI,TaxI_nonlinearpart_unused,I]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,tau_E1);

c=(I+a-TaxI-aprime)/(1+tau_c); % Same budget constraint as in Kitao2008_ReturnFn

TaxRevenue=tau_c*c+TaxI;

end
