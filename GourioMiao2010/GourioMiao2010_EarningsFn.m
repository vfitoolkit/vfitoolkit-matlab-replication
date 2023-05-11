function earnings=GourioMiao2010_EarningsFn(dividend,kprime,k,z,w,alpha_k,alpha_l)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output
y=z*(k^alpha_k)*(l^alpha_l);

% Profit
profit=y-w*l;

earnings=profit;

% Note: earnings is "before tax earnings" and is defined on line 67 of GM2010 code divGEstatdis.m.
% earnings = exp(z)*(Anorm*k^alpha)^(1/(1-nu)) *(nu/wage)^(nu/(1-nu)) *(1-nu)
% nu is their alpha_l, alpha is there alpha_k.
% This is the same as what I call profit (their formula actually assumes
% you are in general eqm so profit=(1-alpha_l)*y. My formulation does not assume this.

end
