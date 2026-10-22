function F=HuggettVentura1999_BoskinMed_ReturnFn(l,a1prime,a1,a2,logprod,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,b_common,C)
% logprod is LOG productivity, and it occupies the slot that is a per-type
% PARAMETER in the no-idiosyncratic spec (logabilityi, ability fixed at birth)
% and an exogenous STATE in the idiosyncratic spec (z, the AR(1) deviation from
% the age-mean log endowment). Same position, same meaning, so one file serves
% both; only the driver's anonymous-function wrapper differs.
%
% The Boskin proposal WITH the hospital/medical transfer retained (HV1999 Fig 4
% and Table VI, specified in footnote 29 p.519). Relative to the plain Boskin
% ReturnFn this differs in exactly one place: the common transfer b_common sits
% inside the max, so the floor applies to TOTAL retirement income from the
% system rather than to the annuity alone:
%
%   b(x,j) = max(b_lower, b_common + C*a2/(1+g)^(j-R))     j >= R
%
% Footnote 29's own wording is "Benefits at retirement are then equal to
% max{b_bar, b + b(a2,j)}", where b_bar is the floor and b the common transfer.
% The alternative reading b_common + max(b_lower, annuity) would make the floor
% apply only to the annuity; that is NOT what is implemented here.
%
% The payroll tax in the agent's budget is the FULL theta: footnote 29 splits
% theta = theta1 + theta2 only to say what each part finances (theta1 the common
% transfer, theta2 the annuity accumulation). Workers pay theta either way. The
% split enters through the aprimeFn, which accumulates a2 at theta2, not theta --
% the driver passes theta2 into that slot.
%
% Decision: l. Endogenous states: a1 (private assets), a2 (accumulated taxes).

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
    benefit=max(b_lower,b_common+b_annuity);
end

c=a1*(1+r*(1-tau))+(1-tau-theta_j)*earnings+benefit-(1+g)*a1prime;

if c>0 && l<1
    F=((c^nu*(1-l)^(1-nu))^(1-sigma))/(1-sigma);
end

end
