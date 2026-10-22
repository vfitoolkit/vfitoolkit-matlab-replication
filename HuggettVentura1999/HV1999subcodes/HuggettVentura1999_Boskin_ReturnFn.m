function F=HuggettVentura1999_Boskin_ReturnFn(l,a1prime,a1,a2,logprod,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,C)
% logprod is LOG productivity, and it occupies the slot that is a per-type
% PARAMETER in the no-idiosyncratic spec (logabilityi, ability fixed at birth)
% and an exogenous STATE in the idiosyncratic spec (z, the AR(1) deviation from
% the age-mean log endowment). Same position, same meaning, so one file serves
% both; only the driver's anonymous-function wrapper differs.
% HV1999 Boskin proposal: budget eqn (8), utility eqn (2), benefit eqn (11).
% Decision: l. Endogenous states: a1 (private assets), a2 (accumulated taxes).
% theta_j = theta for j<Jr, 0 for j>=Jr (paid only while working).

F=-Inf;

if l<0 || l>1
    return
end

if agej<Jr
    earnings=exp(logprod)*kappa_j*l*w;
    theta_j=theta;
    benefit=0;
else
    earnings=0;
    theta_j=0;
    b_annuity=C*a2/((1+g)^(agej-Jr));
    benefit=max(b_lower,b_annuity);
end

c=a1*(1+r*(1-tau))+(1-tau-theta_j)*earnings+benefit-(1+g)*a1prime;

if c>0 && l<1
    F=((c^nu*(1-l)^(1-nu))^(1-sigma))/(1-sigma);
end

end
