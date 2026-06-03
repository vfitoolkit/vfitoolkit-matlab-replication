function d=GourioMiao2010_DividendFn(kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)
% Optimal dividend: given (kprime,k,z), with d and s pinned down by the firm
% budget constraint together with d>=0 and s>=0 (no share repurchases), the
% optimum is d = max(A, 0), where A is the dividend the firm would pay if
% it issued no equity.
%
% Uses profit = (1-alpha_l)*y with l substituted out via the static labor FOC.

y=(z*k^alpha_k)^(1/(1-alpha_l)) * (alpha_l/w)^(alpha_l/(1-alpha_l));
profit=(1-alpha_l)*y;
invest=kprime-(1-delta)*k;
capitaladjcost=(capadjconstant/2)*(invest^2)/k;
T=profit-delta*k-phi*capitaladjcost;
A=profit-tau_corp*T-invest-capitaladjcost;

d=max(A,0);

end
