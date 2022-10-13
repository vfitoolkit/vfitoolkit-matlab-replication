function F=Huggett1996_ReturnFn(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw)
% We don't need w and tau as these can be imputed from the others.

% Note that ybarj and b are both dependent on j

% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
% From K/Y, substituting the production function gives
KdivY=(KdivL^(1-alpha))/A;
% We know w=MPL (MPL is marginal product of labour)
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)
% Huggett (1996) calibrates tau to the following (see pg 478 for explanation)
tau=0.195/(1-delta*KdivY);

% Pensions:

% The borrowing constraint 
% (will be -w if bc_equalsminusw==1, 0 if bc_equalsminusw==0)
bc=-w*bc_equalsminusw;

% What Hugget calls e(z,t)
e_jz=exp(z+ybarj);

% Budget constraint
c=(1+r*(1-tau))*a+(1-tau-theta)*w*e_jz+b*bvec+T-aprime;

F=-Inf; %(l by aprime)

if c>0 && aprime>=bc
    F=(c.^(1-sigma)-1)/(1-sigma);
end

end
