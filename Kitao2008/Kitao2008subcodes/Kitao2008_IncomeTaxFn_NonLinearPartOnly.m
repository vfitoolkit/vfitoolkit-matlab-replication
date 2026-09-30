function TaxI_nonlinearpart=Kitao2008_IncomeTaxFn_NonLinearPartOnly(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0, tau_a1, tau_a2, tau_I, taxincome, tau_k, tau_E1)
% Revenue raised from this household by the NON-LINEAR part of the income tax only, that is by the
%    tau_a0*( base - (tau_a2+base^(-tau_a1))^(-1/tau_a1) )
% part of the tax function, and not by the flat tau_I, tau_k or tau_E1 parts.
%
% This is the object that Kitao (2008), Section 3.2, calibrates tau_a2 to: "The parameter a2 is
% pinned down in equilibrium so that the share of the government expenditures raised by the
% non-linear part of the function equals 65%, the fraction of government tax revenues raised by
% the income tax in data (OECD, 2003)." So the 65% entry in Table 2 is this divided by total tax
% revenue, and not Kitao2008_IncomeTaxFn divided by total tax revenue.

[TaxI_unused,TaxI_nonlinearpart,I_unused]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,tau_E1);

end
