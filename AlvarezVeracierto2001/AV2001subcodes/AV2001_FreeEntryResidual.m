function [FEresidual,V_firm,Policy_firm]=AV2001_FreeEntryResidual(wstar,Params,Grids,vfoptions)
% Free entry condition, eq (3.7): m = (1/(1+i)) sum_s V(0,s) psi(s).
% Returns m minus the expected discounted value of an entrant; increasing in wstar.

Params.wstar=wstar;

FirmReturnFn=@(n,e,s,wstar,r,theta,gamma,epsilon,rho) AV2001_FirmReturnFn(n,e,s,wstar,r,theta,gamma,epsilon,rho);

[V_firm,Policy_firm]=ValueFnIter_InfHorz(0,Grids.n_e,Grids.n_s,[],Grids.e_grid,Grids.s_grid,Grids.Q,FirmReturnFn,Params,{'firmdiscount'},[],vfoptions);

% Entrants have e=0 (first grid point) and draw s in {1,sbar} (grid points 2 and 3) from psi
EV0=Params.psi(1)*gather(V_firm(1,2))+Params.psi(2)*gather(V_firm(1,3));

FEresidual=Params.m-EV0/(1+Params.i);

end
