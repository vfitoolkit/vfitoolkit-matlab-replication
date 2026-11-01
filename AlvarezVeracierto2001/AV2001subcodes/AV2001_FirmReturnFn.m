function F=AV2001_FirmReturnFn(n,e,s,wstar,r,theta,gamma,epsilon,rho)
% Establishment problem, eq (3.1) of AV2001.
% n: current employment (the 'aprime' choice), e: last period employment (the
% endogenous state), s: idiosyncratic productivity (s=0 is a dead establishment).
% Capital is a static choice given (s,n): k=(theta*s*n^gamma/r)^(1/(1-theta)),
% so it is substituted out here rather than being a decision variable.
% Firing cost, eq (3.2): epsilon*wstar per worker dismissed below (1-rho)*e.

k=(theta*s*n^gamma/r)^(1/(1-theta)); % zero when s=0 or n=0
y=s*k^theta*n^gamma;

F=y-wstar*n-r*k-epsilon*wstar*max((1-rho)*e-n,0);

end
