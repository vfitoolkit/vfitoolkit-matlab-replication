function earnings=GourioMiao2010_EarningsFn(kprime,k,z,w,alpha_k,alpha_l)
% Before-tax earnings = profit. Static labor FOC w*l = alpha_l*y gives
% profit = (1-alpha_l)*y with l substituted out. This holds for any w
% (in or out of GE) since it follows from the firm's static optimization,
% not from labor-market clearing.

y=(z*k^alpha_k)^(1/(1-alpha_l)) * (alpha_l/w)^(alpha_l/(1-alpha_l));
earnings=(1-alpha_l)*y;

end
