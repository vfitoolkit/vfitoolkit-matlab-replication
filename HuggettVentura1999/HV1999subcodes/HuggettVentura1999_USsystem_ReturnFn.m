function F=HuggettVentura1999_USsystem_ReturnFn(l,aprime,a,ebar,logprod,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage)
% logprod is LOG productivity, and it occupies the slot that is a per-type
% PARAMETER in the no-idiosyncratic spec (logabilityi, ability fixed at birth)
% and an exogenous STATE in the idiosyncratic spec (z, the AR(1) deviation from
% the age-mean log endowment). Same position, same meaning, so one file serves
% both; only the driver's anonymous-function wrapper differs.
% HV1999 U.S. system: budget eqn (5), utility eqn (2), PIA eqn (15).
% Decision: l (labor). Endogenous states: a, ebar (experience asset).
% Shock z is the productivity log-deviation; e(z,j)=exp(logprod)*kappa_j.
% ebar is in "ratio to mean working-age earnings" units: ebar := AIME/(w*ybar_mean_unwage),
% so bend points are {0.20, 1.24} and the AIME cap is 2.47 (HV1999 Section 3.3).

F=-Inf;

if l<0 || l>1
    return
end

if agej<Jr
    earnings=exp(logprod)*kappa_j*l*w;
    benefit=0;
else
    earnings=0;
    if ebar<=0.20
        PIA=0.90*ebar;
    elseif ebar<=1.24
        PIA=0.90*0.20+0.32*(ebar-0.20);
    else
        PIA=0.90*0.20+0.32*(1.24-0.20)+0.15*(ebar-1.24);
    end
    benefit=b_common+w*ybar_mean_unwage*PIA/((1+g)^(agej-Jr));
end

c=a*(1+r*(1-tau))+(1-tau-theta)*earnings+benefit-(1+g)*aprime;

if c>0 && l<1
    F=((c^nu*(1-l)^(1-nu))^(1-sigma))/(1-sigma);
end

end
