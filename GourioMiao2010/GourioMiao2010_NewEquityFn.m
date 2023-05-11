function s=GourioMiao2010_NewEquityFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output
y=z*(k^alpha_k)*(l^alpha_l);

% Profit
profit=y-w*l;

% Investment
invest=kprime-(1-delta)*k;

% Capital-adjustment costs
capitaladjcost=(capadjconstant/2)*(invest^2)/k; 

% Taxable corporate income
T=profit-delta*k-phi*capitaladjcost;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes (=0 in GM2010 baseline)

% Firms financing constraint gives the new equity issuance
s=dividend+invest+capitaladjcost-(profit-tau_corp*T);

end
