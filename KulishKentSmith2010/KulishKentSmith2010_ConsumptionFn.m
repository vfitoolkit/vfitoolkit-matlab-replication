function c=KulishKentSmith2010_ConsumptionFn(l,aprime,a,ej,r,delta,A,alpha,g_e)

% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
% We know w=MPL (MPL is marginal product of labour)
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

c=(1+r)*a+(1-l)*ej*w-(1+g_e)*aprime;

end