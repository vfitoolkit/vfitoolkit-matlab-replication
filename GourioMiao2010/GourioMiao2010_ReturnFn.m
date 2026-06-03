function F=GourioMiao2010_ReturnFn(kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg)
% Decision variable: kprime only. Dividend d and new equity s are not
% separate decisions. Given (kprime,k,z) the firm budget constraint together
% with the constraints d>=0 and s>=0 (no share repurchases) pin them down.
% Let A be the dividend the firm would pay if it issued no equity. If A>=0
% the firm pays d=A and sets s=0; otherwise it pays d=0 and issues s=-A.
% The (1-tau_d)/(1-tau_cg) tax wedge favouring retained funds together with
% d>=0 force this choice.
%
% Static labor FOC implies w*l = alpha_l*y, so profit = (1-alpha_l)*y with l
% substituted out. This holds for any w (in or out of GE) since it follows
% from the firm's static optimization, not from labor-market clearing.
%
% Note: r is not needed here, it enters the firm via the discount factor.

% Output and profit (l substituted out using the static labor FOC)
y=(z*k^alpha_k)^(1/(1-alpha_l)) * (alpha_l/w)^(alpha_l/(1-alpha_l));
profit=(1-alpha_l)*y;

% Investment and capital-adjustment costs
invest=kprime-(1-delta)*k;
capitaladjcost=(capadjconstant/2)*(invest^2)/k;

% Taxable corporate income (-delta*k is investment expensing; phi=0 in GM2010 baseline)
T=profit-delta*k-phi*capitaladjcost;

% A = dividend if s=0, equivalently -(equity issuance) if d=0
A=profit-tau_corp*T-invest-capitaladjcost;

% F = ((1-tau_d)/(1-tau_cg))*d - s, with optimal (d,s) = (max(A,0), max(-A,0))
F=((1-tau_d)/(1-tau_cg))*max(A,0) + min(A,0);

end
