function F=RestucciaRogerson2008_ReturnFn(aprime_val, a_val, s_val, tau_val,w,r,alpha,gamma,taurate, subsidyrate,cf)
% Note that neither aprime_val nor a_val is actually used for anything in
% this model. But VFI toolkit is not set up to handle that they do not exist.

F=-Inf;

tau=taurate*(tau_val>=0)-subsidyrate*(tau_val<0); % convert rate into the actual tax/subsidy which depends on firm type (as captured by tau_val)

% RR2008, pg 711 explains how to derive the optimal decisions in closed form (conditional on remaining operational)

% Physical capital:
kbar=(alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(s_val*(1-tau))^(1/(1-alpha-gamma));

% Labour
nbar=(((1-tau)*s_val*gamma)/w)^(1/(1-gamma)) *kbar^(alpha/(1-gamma));

% Substitute kbar, nbar in pi and compute W(s,kbar(s,theta))
pibar=(1-tau)*s_val*(kbar^alpha)*(nbar^gamma)-w*nbar-r*kbar-cf;

F=pibar;

end