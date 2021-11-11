function F=ConesaKrueger1999_Ptype_ReturnFn(l,aprime,a,u,zeta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma)
% We don't need w and tau as these can be imputed from the others.

% Note that epsilon_j and I_j are both dependent on j

% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivN=((r+delta)/(theta*alpha))^(1/(alpha-1));
% We know w=MPL (MPL is marginal product of labour)
w=theta*(1-alpha)*(KdivN^alpha); % wage rate (per effective labour unit)

% Difference in Ptype is exp(u+zeta) instead of eta directly (note that eta=exp(u+zeta) and without Ptypes zeta would just be zero)
eta=exp(u+zeta);

% Budget constraint
c=I_j*(1-tau)*w*epsilon_j*eta*l+(1+r)*(a+Tr)+(1-I_j)*SS-aprime;

F=-Inf;

% both aprime>=0 and 0<=l<=1 are built into the grids
if c>0
    F=(((c^gamma)*((1-l)^(1-gamma))).^(1-sigma))/(1-sigma);
end

end