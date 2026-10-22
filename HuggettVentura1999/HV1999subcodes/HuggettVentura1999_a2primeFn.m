function a2prime=HuggettVentura1999_a2primeFn(l,a2,logprod,agej,Jr,g,r,w,theta,kappa_j)
% logprod is LOG productivity, and it occupies the slot that is a per-type
% PARAMETER in the no-idiosyncratic spec (logabilityi, ability fixed at birth)
% and an exogenous STATE in the idiosyncratic spec (z, the AR(1) deviation from
% the age-mean log endowment). Same position, same meaning, so one file serves
% both; only the driver's anonymous-function wrapper differs.
% HV1999 eqn (10): (1+g)*a2' = (a2 + theta*e(z,j)*w*l)*(1+r)   for j<Jr,
%                  a2'       = a2                              for j>=Jr.

a2prime=a2;

if agej<Jr
    a2prime=((a2+theta*exp(logprod)*kappa_j*l*w)*(1+r))/(1+g);
end

end
